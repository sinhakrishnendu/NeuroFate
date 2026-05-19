#!/usr/bin/env python3
"""Build donor/sample-level NeuroFate biological axis scores from existing feature tables."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


LABEL_PREFIXES = ("label__",)
METADATA_COLUMNS = {
    "dataset_id",
    "dataset_unit_id",
    "unit_type",
    "sample_id",
    "donor_id",
    "diagnosis",
    "cognitive_status",
    "cognitive status",
    "apoe_genotype",
    "sex",
    "age",
    "pmi",
    "rin",
    "braak",
    "cerad",
    "lewy_body_midbrain",
    "lewy_body_limbic",
    "lewy_body_neocortical",
}
FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "cell_fraction__",
    "celltype_gene_mean__",
    "celltype_gene_detection__",
    "index__",
)


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


def is_label_or_metadata_column(column: str) -> bool:
    lowered = column.lower()
    return lowered.startswith(LABEL_PREFIXES) or lowered in METADATA_COLUMNS


def feature_columns(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    columns = list(rows[0])
    return [
        column
        for column in columns
        if column.startswith(FEATURE_PREFIXES)
        and not is_label_or_metadata_column(column)
        and any(not math.isnan(to_float(row.get(column))) for row in rows)
    ]


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def field_matches_gene(field: str, gene: str) -> bool:
    upper_field = field.upper()
    token = gene.upper()
    return (
        upper_field.endswith(f"__{token}")
        or f"__{token}__" in upper_field
        or upper_field.endswith(f"_{token}")
        or upper_field == token
    )


def axis_feature_map(features: list[str], genes: list[str]) -> tuple[list[str], list[str], list[str]]:
    used = [feature for feature in features if any(field_matches_gene(feature, gene) for gene in genes)]
    found_genes = sorted({gene for gene in genes if any(field_matches_gene(feature, gene) for feature in used)})
    missing_genes = sorted(set(genes) - set(found_genes))
    return used, found_genes, missing_genes


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


def metadata_columns(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    return [
        column
        for column in rows[0]
        if column.startswith("label__") or column.lower() in METADATA_COLUMNS
    ]


def build_axis_scores(
    rows: list[dict[str, str]],
    cohort: str,
    axes: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    features = feature_columns(rows)
    meta_cols = metadata_columns(rows)
    score_rows: list[dict[str, str]] = []
    coverage_rows: list[dict[str, str]] = []
    raw_axis_values: dict[str, list[float]] = {}
    axis_feature_counts: dict[str, tuple[list[str], list[str], list[str]]] = {}
    for axis in axes:
        genes = parse_genes(axis.get("gene_members", ""))
        used, found_genes, missing_genes = axis_feature_map(features, genes)
        raw_values = [mean([to_float(row.get(feature)) for feature in used]) for row in rows]
        raw_axis_values[axis["axis_id"]] = standardize(raw_values)
        axis_feature_counts[axis["axis_id"]] = (used, found_genes, missing_genes)
        coverage_rows.append(
            {
                "cohort": cohort,
                "axis_id": axis["axis_id"],
                "axis_name": axis.get("axis_name", ""),
                "genes_requested": str(len(genes)),
                "genes_found": str(len(found_genes)),
                "genes_missing": str(len(missing_genes)),
                "features_used": str(len(used)),
                "found_gene_members": ";".join(found_genes),
                "missing_gene_members": ";".join(missing_genes),
                "status": "ok" if used else "insufficient_coverage",
            }
        )
    for row_index, row in enumerate(rows):
        out = {
            "cohort": cohort,
            "unit_id": row.get("dataset_unit_id") or row.get("donor_id") or row.get("sample_id") or str(row_index + 1),
        }
        for column in meta_cols:
            out[column] = row.get(column, "")
        for axis in axes:
            value = raw_axis_values[axis["axis_id"]][row_index]
            out[f"axis__{axis['axis_id']}"] = "" if math.isnan(value) else f"{value:.8g}"
        score_rows.append(out)
    return score_rows, coverage_rows


def output_columns(rows: list[dict[str, str]], axes: list[dict[str, str]]) -> list[str]:
    base = ["cohort", "unit_id"]
    meta = sorted({column for row in rows for column in row if column not in base and not column.startswith("axis__")})
    axis_cols = [f"axis__{axis['axis_id']}" for axis in axes]
    return base + meta + axis_cols


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build NeuroFate donor/sample-level biological axis scores.")
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--pd-features", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/106_build_neurofate_axis_scores.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    axes = read_tsv(args.axis_registry)
    sea_rows, sea_coverage = build_axis_scores(read_tsv(args.sea_ad_features), "sea_ad", axes)
    pd_rows, pd_coverage = build_axis_scores(read_tsv(args.pd_features), "gse243639_pd_snpc", axes)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.outdir / "phase21_sea_ad_axis_scores.tsv", sea_rows, output_columns(sea_rows, axes))
    write_tsv(args.outdir / "phase21_gse243639_axis_scores.tsv", pd_rows, output_columns(pd_rows, axes))
    write_tsv(
        args.outdir / "phase21_axis_feature_coverage.tsv",
        [*sea_coverage, *pd_coverage],
        [
            "cohort",
            "axis_id",
            "axis_name",
            "genes_requested",
            "genes_found",
            "genes_missing",
            "features_used",
            "found_gene_members",
            "missing_gene_members",
            "status",
        ],
    )
    logging.info("Axis scoring completed for SEA-AD rows=%d and PD rows=%d", len(sea_rows), len(pd_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
