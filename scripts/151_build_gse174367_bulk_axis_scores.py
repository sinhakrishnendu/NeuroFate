#!/usr/bin/env python3
"""Build sample-level GSE174367 bulk RNA NeuroFate-Axis scores."""

from __future__ import annotations

import argparse
import csv
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


def read_matrix(path: Path) -> tuple[list[str], dict[str, dict[str, float]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        samples = [field for field in (reader.fieldnames or []) if field != "gene_symbol"]
        values = {sample: {} for sample in samples}
        for row in reader:
            gene = row.get("gene_symbol", "").upper()
            if not gene:
                continue
            for sample in samples:
                values[sample][gene] = to_float(row.get(sample))
    return samples, values


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


def build_scores(
    matrix_values: dict[str, dict[str, float]],
    sample_map: list[dict[str, str]],
    axes: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    valid = [row for row in sample_map if row.get("match_status") == "matched" and row.get("label__ad_vs_control") in {"0", "1"} and row.get("expression_sample_id") in matrix_values]
    raw: dict[str, list[float]] = {axis["axis_id"]: [] for axis in axes}
    for row in valid:
        values = matrix_values[row["expression_sample_id"]]
        for axis in axes:
            genes = parse_genes(axis.get("gene_members", ""))
            raw[axis["axis_id"]].append(mean([values.get(gene, math.nan) for gene in genes]))
    standardized = {axis_id: standardize(values) for axis_id, values in raw.items()}
    out_rows: list[dict[str, str]] = []
    label_counts: dict[str, int] = {}
    for index, row in enumerate(valid):
        label = row["label__ad_vs_control"]
        label_counts[label] = label_counts.get(label, 0) + 1
        out = {
            "sample_id": row["expression_sample_id"],
            "sample_title": row.get("sample_title", ""),
            "geo_accession": row.get("geo_accession", ""),
            "inferred_ad_endpoint": row.get("inferred_ad_endpoint", ""),
            "label__ad_vs_control": label,
        }
        for axis in axes:
            value = standardized[axis["axis_id"]][index]
            out[f"axis__{axis['axis_id']}"] = "" if math.isnan(value) else f"{value:.8g}"
        out_rows.append(out)
    coverage: list[dict[str, str]] = []
    available = {gene for values in matrix_values.values() for gene in values}
    for axis in axes:
        genes = parse_genes(axis.get("gene_members", ""))
        found = sorted(set(genes) & available)
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
    labels = [{"label__ad_vs_control": key, "count": str(value)} for key, value in sorted(label_counts.items())]
    return out_rows, coverage, labels


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build GSE174367 bulk RNA axis scores.")
    parser.add_argument("--axis-gene-matrix", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_matrix.tsv"))
    parser.add_argument("--sample-map", type=Path, default=Path("results/tables/phase31_gse174367_bulk_sample_map.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_scores.tsv"))
    parser.add_argument("--coverage-output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_feature_coverage.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_label_summary.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/151_build_gse174367_bulk_axis_scores_phase31.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    _samples, matrix_values = read_matrix(args.axis_gene_matrix)
    axes = read_tsv(args.axis_registry)
    rows, coverage, labels = build_scores(matrix_values, read_tsv(args.sample_map), axes)
    columns = ["sample_id", "sample_title", "geo_accession", "inferred_ad_endpoint", "label__ad_vs_control"] + [f"axis__{axis['axis_id']}" for axis in axes]
    write_tsv(args.output, rows, columns)
    write_tsv(args.coverage_output, coverage, ["axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    write_tsv(args.label_summary_output, labels, ["label__ad_vs_control", "count"])
    logging.info("Built GSE174367 bulk axis scores samples=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
