#!/usr/bin/env python3
"""Build Phase 5 donor-level feature table from sparse target-gene expression.

This script uses only the sparse target-gene expression TSV and decoded metadata
TSV. It streams expression rows and writes compact donor-level features for
interpretable modeling.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


DEFAULT_MAX_ROW_CHUNK = 1_000_000
MAX_ALLOWED_ROW_CHUNK = 5_000_000

DONOR_FIELD = "Donor ID"
CELLTYPE_FIELD = "Subclass"
LABEL_FIELDS = [
    "Cognitive Status",
    "Overall AD neuropathological Change",
    "Braak",
    "CERAD score",
    "APOE Genotype",
    "Highest Lewy Body Disease",
    "LATE",
    "Overall CAA Score",
    "Neurotypical reference",
]
METADATA_FIELDS = [DONOR_FIELD, CELLTYPE_FIELD, *LABEL_FIELDS]

MICROGLIAL_GENES = ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1"]
ASTROCYTE_GENES = ["GFAP"]
NEURONAL_GENES = ["SLC17A7", "SST", "PVALB", "LAMP5"]
NEURODEGENERATION_GENES = ["PINK1", "PRKN", "SNCA", "MAPT", "APOE"]
INFLAMMATORY_GENES = ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1", "IL1B", "TNF", "NFKB1", "B2M"]
MITOCHONDRIAL_GENES = ["PINK1", "PRKN", "LRRK2"]

DONOR_INDICES = {
    "MAI": MICROGLIAL_GENES,
    "ASI": ASTROCYTE_GENES,
    "NVI": NEURODEGENERATION_GENES,
    "inflammatory_index": INFLAMMATORY_GENES,
    "mitochondrial_index": MITOCHONDRIAL_GENES,
    "neuronal_index": NEURONAL_GENES,
}


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
        ],
    )


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def clean_label(value: str | None) -> str:
    if value is None:
        return "missing"
    value = value.strip()
    return value if value else "missing"


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return cleaned or "missing"


def mode_label(counter: Counter[str]) -> str:
    if not counter:
        return "missing"
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0][0]


def build_gene_to_index_lookup() -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = defaultdict(list)
    for index_name, genes in DONOR_INDICES.items():
        for gene in genes:
            lookup[gene].append(index_name)
    return lookup


def load_metadata(
    metadata_path: Path,
) -> tuple[
    list[str],
    list[str],
    dict[str, int],
    dict[str, Counter[str]],
    dict[str, Counter[str]],
    set[str],
]:
    row_donors: list[str] = []
    row_celltypes: list[str] = []
    donor_cell_counts: dict[str, int] = defaultdict(int)
    donor_label_counts: dict[str, Counter[str]] = defaultdict(Counter)
    donor_celltype_counts: dict[str, Counter[str]] = defaultdict(Counter)
    celltypes: set[str] = set()

    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            donor = clean_label(row.get(DONOR_FIELD) or f"missing_donor_{row_number}")
            celltype = clean_label(row.get(CELLTYPE_FIELD))
            row_donors.append(donor)
            row_celltypes.append(celltype)
            donor_cell_counts[donor] += 1
            donor_celltype_counts[donor][celltype] += 1
            celltypes.add(celltype)
            for field in LABEL_FIELDS:
                donor_label_counts[donor][f"label__{safe_name(field)}::{clean_label(row.get(field))}"] += 1

    logging.info("Loaded decoded metadata rows: %d", len(row_donors))
    logging.info("Detected donors: %d", len(donor_cell_counts))
    logging.info("Detected cell subclasses: %d", len(celltypes))
    return (
        row_donors,
        row_celltypes,
        donor_cell_counts,
        donor_label_counts,
        donor_celltype_counts,
        celltypes,
    )


def extract_label_modes(donor_label_counts: dict[str, Counter[str]]) -> dict[str, dict[str, str]]:
    donor_labels: dict[str, dict[str, str]] = defaultdict(dict)
    for donor, counter in donor_label_counts.items():
        grouped: dict[str, Counter[str]] = defaultdict(Counter)
        for combined, count in counter.items():
            field, label = combined.split("::", 1)
            grouped[field][label] += count
        for field, field_counter in grouped.items():
            donor_labels[donor][field] = mode_label(field_counter)
    return donor_labels


def stream_expression_features(
    expression_path: Path,
    row_donors: list[str],
    row_celltypes: list[str],
    max_row_chunk: int,
) -> tuple[
    set[str],
    dict[tuple[str, str], float],
    dict[tuple[str, str], int],
    dict[tuple[str, str], float],
    dict[tuple[str, str, str], float],
]:
    if max_row_chunk > MAX_ALLOWED_ROW_CHUNK:
        raise ValueError(f"max_row_chunk may not exceed {MAX_ALLOWED_ROW_CHUNK}")

    gene_to_indices = build_gene_to_index_lookup()
    genes_seen: set[str] = set()
    donor_gene_sum: dict[tuple[str, str], float] = defaultdict(float)
    donor_gene_nonzero: dict[tuple[str, str], int] = defaultdict(int)
    donor_index_sum: dict[tuple[str, str], float] = defaultdict(float)
    donor_celltype_index_sum: dict[tuple[str, str, str], float] = defaultdict(float)
    processed = 0

    with open_text(expression_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            processed += 1
            if processed % max_row_chunk == 0:
                logging.info("Processed sparse expression rows: %d", processed)

            row_index = int(row["row_index"])
            if row_index >= len(row_donors):
                logging.warning("Skipping out-of-range row_index=%d", row_index)
                continue
            donor = row_donors[row_index]
            celltype = row_celltypes[row_index]
            gene = row.get("gene_symbol", "")
            value = float(row["expression_value"])
            genes_seen.add(gene)
            donor_gene_sum[(donor, gene)] += value
            donor_gene_nonzero[(donor, gene)] += 1
            for index_name in gene_to_indices.get(gene, []):
                donor_index_sum[(donor, index_name)] += value
                donor_celltype_index_sum[(donor, celltype, index_name)] += value

    logging.info("Total sparse expression rows processed: %d", processed)
    logging.info("Genes observed in sparse expression table: %d", len(genes_seen))
    return genes_seen, donor_gene_sum, donor_gene_nonzero, donor_index_sum, donor_celltype_index_sum


def present_gene_count(index_name: str, genes_seen: set[str]) -> int:
    genes = DONOR_INDICES[index_name]
    return max(1, len([gene for gene in genes if gene in genes_seen]))


def build_feature_rows(
    row_donors: list[str],
    row_celltypes: list[str],
    donor_cell_counts: dict[str, int],
    donor_label_counts: dict[str, Counter[str]],
    donor_celltype_counts: dict[str, Counter[str]],
    celltypes: set[str],
    genes_seen: set[str],
    donor_gene_sum: dict[tuple[str, str], float],
    donor_gene_nonzero: dict[tuple[str, str], int],
    donor_index_sum: dict[tuple[str, str], float],
    donor_celltype_index_sum: dict[tuple[str, str, str], float],
) -> tuple[list[dict[str, str]], list[str]]:
    donor_labels = extract_label_modes(donor_label_counts)
    donors = sorted(donor_cell_counts)
    genes = sorted(genes_seen)
    sorted_celltypes = sorted(celltypes)

    label_columns = [f"label__{safe_name(field)}" for field in LABEL_FIELDS]
    gene_columns = [f"gene_mean__{safe_name(gene)}" for gene in genes] + [
        f"gene_detection__{safe_name(gene)}" for gene in genes
    ]
    index_columns = [f"index__{safe_name(index_name)}" for index_name in sorted(DONOR_INDICES)]
    cell_fraction_columns = [f"cell_fraction__{safe_name(celltype)}" for celltype in sorted_celltypes]
    celltype_index_columns = [
        f"celltype_index__{safe_name(index_name)}__{safe_name(celltype)}"
        for index_name in sorted(DONOR_INDICES)
        for celltype in sorted_celltypes
    ]
    fieldnames = [
        "donor_id",
        "n_cells",
        *label_columns,
        *gene_columns,
        *index_columns,
        *cell_fraction_columns,
        *celltype_index_columns,
    ]

    rows: list[dict[str, str]] = []
    for donor in donors:
        cell_count = donor_cell_counts[donor]
        row: dict[str, str] = {"donor_id": donor, "n_cells": str(cell_count)}
        for column in label_columns:
            row[column] = donor_labels.get(donor, {}).get(column, "missing")

        for gene in genes:
            row[f"gene_mean__{safe_name(gene)}"] = f"{donor_gene_sum.get((donor, gene), 0.0) / cell_count:.8g}"
            row[f"gene_detection__{safe_name(gene)}"] = f"{donor_gene_nonzero.get((donor, gene), 0) / cell_count:.8g}"

        for index_name in sorted(DONOR_INDICES):
            denominator = cell_count * present_gene_count(index_name, genes_seen)
            row[f"index__{safe_name(index_name)}"] = f"{donor_index_sum.get((donor, index_name), 0.0) / denominator:.8g}"

        for celltype in sorted_celltypes:
            celltype_count = donor_celltype_counts[donor].get(celltype, 0)
            row[f"cell_fraction__{safe_name(celltype)}"] = f"{celltype_count / cell_count:.8g}"
            for index_name in sorted(DONOR_INDICES):
                denominator = celltype_count * present_gene_count(index_name, genes_seen)
                value = (
                    donor_celltype_index_sum.get((donor, celltype, index_name), 0.0) / denominator
                    if denominator
                    else 0.0
                )
                row[f"celltype_index__{safe_name(index_name)}__{safe_name(celltype)}"] = f"{value:.8g}"
        rows.append(row)

    logging.info("Prepared donor feature rows: %d", len(rows))
    logging.info("Feature columns prepared: %d", len(fieldnames) - 2 - len(label_columns))
    return rows, fieldnames


def write_feature_table(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote donor feature table: %s", path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Phase 5 donor-level feature table.")
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/phase5_donor_feature_table.tsv"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/22_build_donor_feature_table.log"),
    )
    parser.add_argument("--max-row-chunk", type=int, default=DEFAULT_MAX_ROW_CHUNK)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting Phase 5 donor-level feature table build.")
    logging.info("Expression input: %s", args.expression)
    logging.info("Metadata input: %s", args.metadata)

    (
        row_donors,
        row_celltypes,
        donor_cell_counts,
        donor_label_counts,
        donor_celltype_counts,
        celltypes,
    ) = load_metadata(args.metadata)
    (
        genes_seen,
        donor_gene_sum,
        donor_gene_nonzero,
        donor_index_sum,
        donor_celltype_index_sum,
    ) = stream_expression_features(args.expression, row_donors, row_celltypes, args.max_row_chunk)

    rows, fieldnames = build_feature_rows(
        row_donors,
        row_celltypes,
        donor_cell_counts,
        donor_label_counts,
        donor_celltype_counts,
        celltypes,
        genes_seen,
        donor_gene_sum,
        donor_gene_nonzero,
        donor_index_sum,
        donor_celltype_index_sum,
    )
    write_feature_table(args.output, rows, fieldnames)
    logging.info("No H5AD file, full matrix, or pipeline-engine analysis was used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
