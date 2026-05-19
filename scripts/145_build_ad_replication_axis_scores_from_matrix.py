#!/usr/bin/env python3
"""Build AD replication NeuroFate-Axis scores from sample-level CSV/TSV matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
from pathlib import Path


METADATA_HINTS = ("label", "diagnosis", "disease", "condition", "status", "pathology", "sample", "donor", "age", "sex", "pmi")


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", newline="")


def delimiter_for(path: Path) -> str:
    return "\t" if any(path.name.endswith(suffix) for suffix in [".tsv", ".tsv.gz", ".txt", ".txt.gz"]) else ","


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
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


def norm(value: str | None) -> str:
    return str(value or "").strip().casefold()


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def axis_gene_map(axes: list[dict[str, str]]) -> dict[str, list[str]]:
    return {axis["axis_id"]: parse_genes(axis.get("gene_members", "")) for axis in axes}


def all_axis_genes(axes: list[dict[str, str]]) -> set[str]:
    genes: set[str] = set()
    for axis in axes:
        genes.update(parse_genes(axis.get("gene_members", "")))
    return genes


def is_metadata_column(column: str) -> bool:
    lowered = column.lower()
    return any(hint in lowered for hint in METADATA_HINTS)


def read_metadata(path: Path, sample_col: str, label_col: str, positive: str, negative: str) -> dict[str, dict[str, str]]:
    delim = delimiter_for(path)
    metadata: dict[str, dict[str, str]] = {}
    with opener(path) as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        for row in reader:
            sample = row.get(sample_col, "")
            if not sample:
                continue
            label = row.get(label_col, "")
            if norm(label) not in {norm(positive), norm(negative)}:
                continue
            metadata[sample] = {"sample_id": sample, label_col: label, "label__ad_replication_binary": "1" if norm(label) == norm(positive) else "0"}
    return metadata


def detect_orientation(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    with opener(path) as handle:
        header = next(csv.reader(handle, delimiter=delimiter_for(path)), [])
    first = header[0].strip().lower() if header else ""
    return "genes_rows" if first in {"gene", "genes", "gene_symbol", "symbol", "feature"} else "samples_rows"


def read_genes_rows(path: Path, genes: set[str]) -> dict[str, dict[str, float]]:
    sample_values: dict[str, dict[str, float]] = {}
    with opener(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter_for(path))
        header = next(reader)
        samples = header[1:]
        for row in reader:
            if not row:
                continue
            gene = row[0].strip().upper()
            if gene not in genes:
                continue
            for sample, value in zip(samples, row[1:], strict=False):
                sample_values.setdefault(sample, {})[gene] = to_float(value)
    return sample_values


def read_samples_rows(path: Path, genes: set[str], sample_id_column: str) -> dict[str, dict[str, float]]:
    sample_values: dict[str, dict[str, float]] = {}
    with opener(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(path))
        fieldnames = reader.fieldnames or []
        gene_columns = [column for column in fieldnames if column.upper() in genes and not is_metadata_column(column)]
        for row in reader:
            sample = row.get(sample_id_column, "")
            if not sample:
                continue
            sample_values[sample] = {gene.upper(): to_float(row.get(gene)) for gene in gene_columns}
    return sample_values


def mean(values: list[float]) -> float:
    observed = [value for value in values if not math.isnan(value)]
    return sum(observed) / len(observed) if observed else math.nan


def standardize(values: list[float]) -> list[float]:
    observed = [value for value in values if not math.isnan(value)]
    if not observed:
        return [math.nan for _ in values]
    mu = sum(observed) / len(observed)
    variance = sum((value - mu) ** 2 for value in observed) / max(1, len(observed) - 1)
    sd = math.sqrt(variance) or 1.0
    return [(value - mu) / sd if not math.isnan(value) else math.nan for value in values]


def build_scores(sample_values: dict[str, dict[str, float]], metadata: dict[str, dict[str, str]], axes: list[dict[str, str]], label_col: str) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    axis_map = axis_gene_map(axes)
    samples = sorted(sample for sample in sample_values if sample in metadata)
    raw_scores = {axis["axis_id"]: [] for axis in axes}
    for sample in samples:
        values = sample_values[sample]
        for axis_id, genes in axis_map.items():
            raw_scores[axis_id].append(mean([values.get(gene, math.nan) for gene in genes]))
    standardized = {axis_id: standardize(values) for axis_id, values in raw_scores.items()}
    rows: list[dict[str, str]] = []
    label_counts: dict[str, int] = {}
    for index, sample in enumerate(samples):
        meta = metadata[sample]
        label = meta["label__ad_replication_binary"]
        label_counts[label] = label_counts.get(label, 0) + 1
        row = {"sample_id": sample, label_col: meta[label_col], "label__ad_replication_binary": label}
        for axis in axes:
            value = standardized[axis["axis_id"]][index]
            row[f"axis__{axis['axis_id']}"] = "" if math.isnan(value) else f"{value:.8g}"
        rows.append(row)
    coverage: list[dict[str, str]] = []
    for axis in axes:
        genes = axis_map[axis["axis_id"]]
        found = sorted({gene for values in sample_values.values() for gene in genes if gene in values})
        missing = sorted(set(genes) - set(found))
        coverage.append(
            {
                "axis_id": axis["axis_id"],
                "genes_requested": str(len(genes)),
                "genes_found": str(len(found)),
                "genes_missing": str(len(missing)),
                "found_gene_members": ";".join(found),
                "missing_gene_members": ";".join(missing),
                "status": "ok" if found else "insufficient_coverage",
            }
        )
    labels = [{"label__ad_replication_binary": key, "count": str(value)} for key, value in sorted(label_counts.items())]
    return rows, coverage, labels


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build AD replication axis scores from sample-level expression matrices.")
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--sample-id-column", required=True)
    parser.add_argument("--label-column", required=True)
    parser.add_argument("--positive-class", required=True)
    parser.add_argument("--negative-class", required=True)
    parser.add_argument("--orientation", choices=["auto", "genes_rows", "samples_rows"], default="auto")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--coverage-output", type=Path)
    parser.add_argument("--label-summary-output", type=Path)
    parser.add_argument("--log-file", type=Path)
    parser.add_argument("--max-file-size-mb", type=float, default=250.0)
    parser.add_argument("--allow-large-matrix", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/145_{args.cohort_id}_axis_scores.log")
    configure_logging(log_file)
    if args.matrix.stat().st_size > args.max_file_size_mb * 1024 * 1024 and not args.allow_large_matrix:
        raise SystemExit("Matrix is above the lightweight limit. Use a reviewed manual command with --allow-large-matrix only if appropriate.")
    axes = read_tsv(args.axis_registry)
    genes = all_axis_genes(axes)
    orientation = detect_orientation(args.matrix, args.orientation)
    sample_values = read_genes_rows(args.matrix, genes) if orientation == "genes_rows" else read_samples_rows(args.matrix, genes, args.sample_id_column)
    rows, coverage, labels = build_scores(sample_values, read_metadata(args.metadata, args.sample_id_column, args.label_column, args.positive_class, args.negative_class), axes, args.label_column)
    output = args.output or Path(f"results/tables/phase28_{args.cohort_id}_axis_scores.tsv")
    coverage_output = args.coverage_output or Path(f"results/tables/phase28_{args.cohort_id}_axis_feature_coverage.tsv")
    label_output = args.label_summary_output or Path(f"results/tables/phase28_{args.cohort_id}_axis_label_summary.tsv")
    columns = ["sample_id", args.label_column, "label__ad_replication_binary"] + [f"axis__{axis['axis_id']}" for axis in axes]
    write_tsv(output, rows, columns)
    write_tsv(coverage_output, coverage, ["axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    write_tsv(label_output, labels, ["label__ad_replication_binary", "count"])
    logging.info("Built AD replication axis scores cohort=%s samples=%d orientation=%s", args.cohort_id, len(rows), orientation)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
