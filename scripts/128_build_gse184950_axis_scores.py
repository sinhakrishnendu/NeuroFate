#!/usr/bin/env python3
"""Build sample-level GSE184950 NeuroFate-Axis scores from sample-level gene summaries."""

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


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "") else math.nan
    except ValueError:
        return math.nan


def axis_gene_map(axis_registry: Path) -> dict[str, list[str]]:
    axes = {}
    for row in read_tsv(axis_registry):
        axes[row["axis_id"]] = [gene.strip().upper() for gene in row.get("gene_members", "").replace(",", ";").split(";") if gene.strip()]
    return axes


def standardize(values: list[float]) -> list[float]:
    observed = [value for value in values if not math.isnan(value)]
    if not observed:
        return [math.nan for _ in values]
    mean = sum(observed) / len(observed)
    variance = sum((value - mean) ** 2 for value in observed) / max(1, len(observed) - 1)
    sd = math.sqrt(variance) or 1.0
    return [(value - mean) / sd if not math.isnan(value) else math.nan for value in values]


def metadata_label(row: dict[str, str]) -> str:
    if row.get("label__pd_pdd_vs_control") in {"0", "1"}:
        return row["label__pd_pdd_vs_control"]
    disease = (row.get("disease_state") or "").lower().replace("'", "")
    if "control" in disease:
        return "0"
    if "parkinson" in disease:
        return "1"
    return ""


def expression_value(row: dict[str, str]) -> float:
    for column in ("mean_expression", "gene_sum"):
        value = to_float(row.get(column))
        if not math.isnan(value):
            return value
    return math.nan


def build_scores(expression_rows: list[dict[str, str]], metadata_rows: list[dict[str, str]], axes: dict[str, list[str]]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    metadata = {row.get("sample_name", "") or row.get("sample_id", ""): row for row in metadata_rows if row.get("sample_name") or row.get("sample_id")}
    valid_samples = {sample for sample, row in metadata.items() if row.get("disease_state") and metadata_label(row) in {"0", "1"}}
    sample_gene_values: dict[str, dict[str, float]] = {}
    for row in expression_rows:
        sample = row.get("sample_id", "")
        if sample not in valid_samples:
            continue
        gene = row.get("gene_symbol", "").upper()
        value = expression_value(row)
        if sample and gene and not math.isnan(value):
            sample_gene_values.setdefault(sample, {})[gene] = value
    samples = sorted(sample for sample in valid_samples if sample in sample_gene_values)
    raw_scores: dict[str, list[float]] = {}
    coverage: list[dict[str, str]] = []
    for axis_id, genes in axes.items():
        values = []
        found = sorted({gene for sample_values in sample_gene_values.values() for gene in genes if gene in sample_values})
        missing = sorted(set(genes) - set(found))
        for sample in samples:
            observed = [sample_gene_values[sample].get(gene, math.nan) for gene in genes]
            observed = [value for value in observed if not math.isnan(value)]
            values.append(sum(observed) / len(observed) if observed else math.nan)
        raw_scores[axis_id] = standardize(values)
        coverage.append(
            {
                "axis_id": axis_id,
                "genes_requested": str(len(genes)),
                "genes_found": str(len(found)),
                "genes_missing": str(len(missing)),
                "found_gene_members": ";".join(found),
                "missing_gene_members": ";".join(missing),
                "status": "ok" if found else "insufficient_coverage",
            }
        )
    rows: list[dict[str, str]] = []
    label_counts: dict[str, int] = {}
    for index, sample in enumerate(samples):
        meta = metadata.get(sample, {})
        label = metadata_label(meta)
        if label not in {"0", "1"}:
            continue
        label_counts[label] = label_counts.get(label, 0) + 1
        out = {
            "sample_id": sample,
            "disease_state": meta.get("disease_state", ""),
            "label__pd_pdd_vs_control": label,
        }
        for axis_id in axes:
            value = raw_scores[axis_id][index]
            out[f"axis__{axis_id}"] = "" if math.isnan(value) else f"{value:.8g}"
        rows.append(out)
    label_summary = [{"label__pd_pdd_vs_control": label, "count": str(count)} for label, count in sorted(label_counts.items())]
    return rows, coverage, label_summary


def preferred_metadata(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if rows:
        if len(rows) <= 2:
            logging.warning("GSE184950 metadata has only %d rows; use Phase 25 series metadata for replication.", len(rows))
        return rows
    fallback = Path("results/tables/phase24_gse184950_sample_metadata.tsv")
    rows = read_tsv(fallback)
    if rows and len(rows) <= 2:
        logging.warning("Falling back to incomplete Phase 24 workbook metadata with %d rows.", len(rows))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build GSE184950 sample-level axis scores.")
    parser.add_argument("--expression", type=Path, default=Path("data/interim/external/gse184950_pd_sn/gse184950_axis_gene_sample_summary.tsv.gz"))
    parser.add_argument("--sample-metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase27_gse184950_axis_scores_clean.tsv"))
    parser.add_argument("--coverage-output", type=Path, default=Path("results/tables/phase27_gse184950_axis_feature_coverage_clean.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/phase27_gse184950_axis_label_summary_clean.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/128_build_gse184950_axis_scores.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows, coverage, label_summary = build_scores(read_gzip_tsv(args.expression), preferred_metadata(args.sample_metadata), axis_gene_map(args.axis_registry))
    axis_columns = [column for column in (rows[0] if rows else {}) if column.startswith("axis__")]
    write_tsv(args.output, rows, ["sample_id", "disease_state", "label__pd_pdd_vs_control", *axis_columns])
    write_tsv(args.coverage_output, coverage, ["axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    write_tsv(args.label_summary_output, label_summary, ["label__pd_pdd_vs_control", "count"])
    logging.info("Built GSE184950 Phase 25 axis scores samples=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
