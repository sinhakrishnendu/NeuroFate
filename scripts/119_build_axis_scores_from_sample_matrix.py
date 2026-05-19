#!/usr/bin/env python3
"""Build NeuroFate-Axis scores from a donor/sample-level CSV/TSV matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
from pathlib import Path


LABEL_HINTS = ("label", "diagnosis", "disease", "status", "pathology", "sample_id", "donor_id")


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


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def axis_gene_map(axes: list[dict[str, str]]) -> dict[str, list[str]]:
    return {axis["axis_id"]: parse_genes(axis.get("gene_members", "")) for axis in axes}


def all_axis_genes(axes: list[dict[str, str]]) -> set[str]:
    genes: set[str] = set()
    for axis in axes:
        genes.update(parse_genes(axis.get("gene_members", "")))
    return genes


def is_label_or_metadata_column(column: str) -> bool:
    lowered = column.lower()
    return any(hint in lowered for hint in LABEL_HINTS)


def read_metadata(path: Path, sample_col: str, label_col: str) -> dict[str, dict[str, str]]:
    delim = delimiter_for(path)
    with opener(path) as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        return {row.get(sample_col, ""): {sample_col: row.get(sample_col, ""), label_col: row.get(label_col, "")} for row in reader if row.get(sample_col)}


def detect_orientation(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    delim = delimiter_for(path)
    with opener(path) as handle:
        reader = csv.reader(handle, delimiter=delim)
        header = next(reader, [])
    first = header[0].lower() if header else ""
    return "genes_rows" if first in {"gene", "genes", "gene_symbol", "symbol", "feature"} else "samples_rows"


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


def read_genes_rows(path: Path, genes: set[str]) -> dict[str, dict[str, float]]:
    delim = delimiter_for(path)
    sample_values: dict[str, dict[str, float]] = {}
    with opener(path) as handle:
        reader = csv.reader(handle, delimiter=delim)
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
    delim = delimiter_for(path)
    sample_values: dict[str, dict[str, float]] = {}
    with opener(path) as handle:
        reader = csv.DictReader(handle, delimiter=delim)
        if reader.fieldnames is None:
            return sample_values
        gene_columns = [column for column in reader.fieldnames if column.upper() in genes and not is_label_or_metadata_column(column)]
        for row in reader:
            sample = row.get(sample_id_column, "")
            if not sample:
                continue
            sample_values[sample] = {gene.upper(): to_float(row.get(gene)) for gene in gene_columns}
    return sample_values


def build_scores(sample_values: dict[str, dict[str, float]], metadata: dict[str, dict[str, str]], axes: list[dict[str, str]], label_col: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    axis_genes = axis_gene_map(axes)
    raw_scores: dict[str, list[float]] = {axis["axis_id"]: [] for axis in axes}
    samples = sorted(sample_values)
    for sample in samples:
        values = sample_values[sample]
        for axis_id, genes in axis_genes.items():
            raw_scores[axis_id].append(mean([values.get(gene, math.nan) for gene in genes]))
    standardized = {axis_id: standardize(values) for axis_id, values in raw_scores.items()}
    rows: list[dict[str, str]] = []
    for index, sample in enumerate(samples):
        out = {"sample_id": sample, label_col: metadata.get(sample, {}).get(label_col, "")}
        for axis in axes:
            value = standardized[axis["axis_id"]][index]
            out[f"axis__{axis['axis_id']}"] = "" if math.isnan(value) else f"{value:.8g}"
        rows.append(out)
    coverage = []
    for axis in axes:
        genes = axis_genes[axis["axis_id"]]
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
    return rows, coverage


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build axis scores from a sample-level expression matrix.")
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--sample-id-column", required=True)
    parser.add_argument("--label-column", required=True)
    parser.add_argument("--orientation", choices=["auto", "genes_rows", "samples_rows"], default="auto")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--coverage-output", type=Path)
    parser.add_argument("--log-file", type=Path)
    parser.add_argument("--max-file-size-mb", type=float, default=250.0)
    parser.add_argument("--allow-large-matrix", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/119_{args.cohort_id}_axis_scores.log")
    configure_logging(log_file)
    if args.matrix.stat().st_size > args.max_file_size_mb * 1024 * 1024 and not args.allow_large_matrix:
        raise SystemExit("Matrix is above the lightweight size limit. Use a manual reviewed command with --allow-large-matrix only when appropriate.")
    axes = read_tsv(args.axis_registry)
    genes = all_axis_genes(axes)
    orientation = detect_orientation(args.matrix, args.orientation)
    sample_values = read_genes_rows(args.matrix, genes) if orientation == "genes_rows" else read_samples_rows(args.matrix, genes, args.sample_id_column)
    rows, coverage = build_scores(sample_values, read_metadata(args.metadata, args.sample_id_column, args.label_column), axes, args.label_column)
    output = args.output or Path(f"results/tables/phase23_{args.cohort_id}_axis_scores.tsv")
    coverage_output = args.coverage_output or Path(f"results/tables/phase23_{args.cohort_id}_axis_feature_coverage.tsv")
    columns = ["sample_id", args.label_column] + [f"axis__{axis['axis_id']}" for axis in axes]
    write_tsv(output, rows, columns)
    write_tsv(coverage_output, coverage, ["axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    logging.info("Built axis scores for cohort=%s samples=%d orientation=%s", args.cohort_id, len(rows), orientation)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
