#!/usr/bin/env python3
"""Build donor/sample-level PD axis scores from GEO series matrices.

Only rows that map to NeuroFate axis genes are retained. Probe identifiers require an
explicit platform mapping; ambiguous genome-wide conversion is intentionally avoided.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", errors="replace", newline="")


def clean(value: str | None) -> str:
    text = str(value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1]
    return text.strip()


def norm(value: str | None) -> str:
    return clean(value).casefold().replace("_", " ").replace("-", " ")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def axis_gene_map(axes: list[dict[str, str]]) -> dict[str, list[str]]:
    return {axis["axis_id"]: parse_genes(axis.get("gene_members", "")) for axis in axes}


def all_axis_genes(axes: list[dict[str, str]]) -> set[str]:
    genes: set[str] = set()
    for axis in axes:
        genes.update(parse_genes(axis.get("gene_members", "")))
    return genes


def alias_to_symbol(path: Path) -> dict[str, str]:
    aliases: dict[str, str] = {}
    if not path.exists():
        return aliases
    for row in read_tsv(path):
        symbol = clean(row.get("gene_symbol", "")).upper()
        if not symbol:
            continue
        for key in ["ensembl_gene_id", "alias"]:
            value = clean(row.get(key, ""))
            if value:
                aliases[value.upper().split(".", 1)[0]] = symbol
    return aliases


def load_platform_mapping(path: Path | None) -> dict[str, list[str]]:
    if not path:
        return {}
    mapping: dict[str, list[str]] = {}
    rows = read_tsv(path)
    for row in rows:
        probe = clean(row.get("probe_id", "") or row.get("ID", "") or row.get("id", ""))
        symbol = clean(row.get("gene_symbol", "") or row.get("Gene Symbol", "") or row.get("symbol", ""))
        if probe and symbol:
            mapping.setdefault(probe, []).append(symbol.upper())
    return mapping


def parse_series_metadata(path: Path) -> dict[str, dict[str, str]]:
    sample_rows: dict[str, list[list[str]]] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = clean(row[0])
            if key.startswith("!Sample_"):
                sample_rows.setdefault(key, []).append([clean(value) for value in row[1:]])
            if key == "!series_matrix_table_begin":
                break
    titles = sample_rows.get("!Sample_title", [[]])[0]
    accessions = sample_rows.get("!Sample_geo_accession", [[]])[0]
    sources = sample_rows.get("!Sample_source_name_ch1", [[]])[0]
    count = max([len(titles), len(accessions), len(sources), *[len(row) for row in sample_rows.get("!Sample_characteristics_ch1", [])]] or [0])
    samples: dict[str, dict[str, str]] = {}
    for index in range(count):
        geo = accessions[index] if index < len(accessions) else ""
        title = titles[index] if index < len(titles) else geo
        sample = {
            "sample_id": geo or title,
            "geo_accession": geo,
            "sample_title": title,
            "source_name": sources[index] if index < len(sources) else "",
            "disease_state": "",
        }
        characteristic_values = []
        for values in sample_rows.get("!Sample_characteristics_ch1", []):
            if index < len(values):
                characteristic_values.append(values[index])
                if ":" in values[index]:
                    key, value = values[index].split(":", 1)
                    if any(token in norm(key) for token in ["disease", "diagnosis", "condition", "status"]):
                        sample["disease_state"] = clean(value)
        sample["all_metadata_text"] = " ".join([sample["sample_title"], sample["source_name"], sample["disease_state"], *characteristic_values])
        for key in {geo, title}:
            if key:
                samples[key] = sample
    return samples


def infer_pd_label(sample: dict[str, str]) -> str:
    text = norm(sample.get("all_metadata_text", ""))
    if "control" in text or "unaffected" in text or "normal" in text:
        return "0"
    if "parkinson" in text or " pd" in f" {text} " or "pdd" in text:
        return "1"
    return ""


def map_feature(feature_id: str, axis_genes: set[str], aliases: dict[str, str], platform_mapping: dict[str, list[str]]) -> list[str]:
    feature = clean(feature_id)
    upper = feature.upper()
    stripped = upper.split(".", 1)[0]
    if upper in axis_genes:
        return [upper]
    if stripped in aliases and aliases[stripped] in axis_genes:
        return [aliases[stripped]]
    mapped = [gene for gene in platform_mapping.get(feature, []) + platform_mapping.get(upper, []) if gene in axis_genes]
    return sorted(set(mapped))


def to_float(value: str) -> float:
    try:
        return float(clean(value))
    except ValueError:
        return math.nan


def read_axis_expression(
    path: Path,
    axis_genes: set[str],
    aliases: dict[str, str],
    platform_mapping: dict[str, list[str]],
) -> tuple[list[str], dict[str, dict[str, list[float]]]]:
    in_table = False
    sample_columns: list[str] = []
    values: dict[str, dict[str, list[float]]] = {}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            key = clean(row[0])
            if key == "!series_matrix_table_begin":
                in_table = True
                continue
            if key == "!series_matrix_table_end":
                break
            if not in_table:
                continue
            if key in {"ID_REF", "ID", "Gene", "gene_symbol", "symbol"}:
                sample_columns = [clean(value) for value in row[1:]]
                continue
            if not sample_columns:
                continue
            genes = map_feature(row[0], axis_genes, aliases, platform_mapping)
            if not genes:
                continue
            for sample, raw_value in zip(sample_columns, row[1:], strict=False):
                value = to_float(raw_value)
                if math.isnan(value):
                    continue
                for gene in genes:
                    values.setdefault(sample, {}).setdefault(gene, []).append(value)
    if not sample_columns:
        raise SystemExit("GEO series matrix expression table was not found; use processed supplementary files or platform mapping.")
    return sample_columns, values


def mean(values: list[float]) -> float:
    observed = [value for value in values if not math.isnan(value)]
    return sum(observed) / len(observed) if observed else math.nan


def standardize(values: list[float]) -> list[float]:
    observed = [value for value in values if not math.isnan(value)]
    if not observed:
        return [math.nan for _ in values]
    mu = sum(observed) / len(observed)
    sd = math.sqrt(sum((value - mu) ** 2 for value in observed) / max(1, len(observed) - 1)) or 1.0
    return [(value - mu) / sd if not math.isnan(value) else math.nan for value in values]


def build_axis_scores(
    sample_columns: list[str],
    gene_values: dict[str, dict[str, list[float]]],
    metadata: dict[str, dict[str, str]],
    axes: list[dict[str, str]],
    cohort_id: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    axis_map = axis_gene_map(axes)
    sample_gene_means = {
        sample: {gene: mean(vals) for gene, vals in genes.items()}
        for sample, genes in gene_values.items()
    }
    usable_samples = [sample for sample in sample_columns if sample in sample_gene_means and sample in metadata and infer_pd_label(metadata[sample]) in {"0", "1"}]
    raw_scores: dict[str, list[float]] = {axis["axis_id"]: [] for axis in axes}
    for sample in usable_samples:
        for axis_id, genes in axis_map.items():
            raw_scores[axis_id].append(mean([sample_gene_means[sample].get(gene, math.nan) for gene in genes]))
    standardized = {axis_id: standardize(vals) for axis_id, vals in raw_scores.items()}
    rows: list[dict[str, str]] = []
    label_counts: dict[str, int] = {}
    for index, sample in enumerate(usable_samples):
        meta = metadata[sample]
        label = infer_pd_label(meta)
        label_counts[label] = label_counts.get(label, 0) + 1
        row = {
            "cohort_id": cohort_id,
            "sample_id": sample,
            "geo_accession": meta.get("geo_accession", ""),
            "sample_title": meta.get("sample_title", ""),
            "disease_state": meta.get("disease_state", ""),
            "label__pd_vs_control": label,
        }
        for axis in axes:
            value = standardized[axis["axis_id"]][index]
            row[f"axis__{axis['axis_id']}"] = "" if math.isnan(value) else f"{value:.8g}"
        rows.append(row)
    coverage: list[dict[str, str]] = []
    found_universe = {gene for genes in sample_gene_means.values() for gene in genes}
    for axis in axes:
        genes = axis_map[axis["axis_id"]]
        found = sorted(set(genes) & found_universe)
        missing = sorted(set(genes) - set(found))
        coverage.append(
            {
                "cohort_id": cohort_id,
                "axis_id": axis["axis_id"],
                "genes_requested": str(len(genes)),
                "genes_found": str(len(found)),
                "genes_missing": str(len(missing)),
                "found_gene_members": ";".join(found),
                "missing_gene_members": ";".join(missing),
                "status": "ok" if found else "insufficient_coverage",
            }
        )
    labels = [{"cohort_id": cohort_id, "label__pd_vs_control": key, "count": str(value)} for key, value in sorted(label_counts.items())]
    return rows, coverage, labels


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build PD axis scores from GEO series matrix expression tables.")
    parser.add_argument("--series-matrix", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--alias-table", type=Path, default=Path("metadata/neurofate_axis_gene_aliases.tsv"))
    parser.add_argument("--platform-annotation", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--coverage-output", type=Path)
    parser.add_argument("--label-summary-output", type=Path)
    parser.add_argument("--log-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/163_{args.cohort_id}_axis_scores.log")
    configure_logging(log_file)
    axes = read_tsv(args.axis_registry)
    axis_genes = all_axis_genes(axes)
    sample_columns, expression = read_axis_expression(args.series_matrix, axis_genes, alias_to_symbol(args.alias_table), load_platform_mapping(args.platform_annotation))
    if not any(expression.values()):
        raise SystemExit("No NeuroFate axis genes were found. Provide a reviewed platform annotation mapping before proceeding.")
    rows, coverage, labels = build_axis_scores(sample_columns, expression, parse_series_metadata(args.series_matrix), axes, args.cohort_id)
    output = args.output or Path(f"results/tables/phase33_{args.cohort_id}_axis_scores.tsv")
    coverage_output = args.coverage_output or Path(f"results/tables/phase33_{args.cohort_id}_axis_feature_coverage.tsv")
    label_output = args.label_summary_output or Path(f"results/tables/phase33_{args.cohort_id}_label_summary.tsv")
    columns = ["cohort_id", "sample_id", "geo_accession", "sample_title", "disease_state", "label__pd_vs_control"] + [f"axis__{axis['axis_id']}" for axis in axes]
    write_tsv(output, rows, columns)
    write_tsv(coverage_output, coverage, ["cohort_id", "axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    write_tsv(label_output, labels, ["cohort_id", "label__pd_vs_control", "count"])
    logging.info("Built Phase 33 PD axis scores cohort=%s samples=%d", args.cohort_id, len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
