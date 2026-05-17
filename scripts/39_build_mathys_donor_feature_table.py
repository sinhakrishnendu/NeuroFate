#!/usr/bin/env python3
"""Build future Mathys 2019 donor-level feature table from sparse target genes.

This script consumes only lightweight metadata TSVs and sparse target-gene TSVs.
It aligns output columns with the Phase 5 donor feature schema where possible.
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


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "index__",
    "cell_fraction__",
    "celltype_index__",
)

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

LABEL_FIELDS = [
    "Cognitive Status",
    "Overall AD neuropathological Change",
    "Braak",
    "CERAD score",
    "APOE Genotype",
    "Highest Lewy Body Disease",
    "LATE",
    "Overall CAA Score",
    "diagnosis",
    "disease_status",
]


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def clean_label(value: str | None) -> str:
    value = (value or "").strip()
    return value if value else "missing"


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return cleaned or "missing"


def mode_label(counter: Counter[str]) -> str:
    if not counter:
        return "missing"
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0][0]


def gene_to_index_lookup() -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = defaultdict(list)
    for index_name, genes in DONOR_INDICES.items():
        for gene in genes:
            lookup[gene].append(index_name)
    return lookup


def load_metadata(
    metadata_path: Path,
    donor_field: str,
    celltype_field: str,
) -> tuple[list[str], list[str], dict[str, int], dict[str, Counter[str]], dict[str, Counter[str]], set[str]]:
    row_donors: list[str] = []
    row_celltypes: list[str] = []
    donor_cell_counts: dict[str, int] = defaultdict(int)
    donor_label_counts: dict[str, Counter[str]] = defaultdict(Counter)
    donor_celltype_counts: dict[str, Counter[str]] = defaultdict(Counter)
    celltypes: set[str] = set()
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            donor = clean_label(row.get(donor_field) or row.get("donor") or row.get("individual") or f"mathys_missing_donor_{row_number}")
            celltype = clean_label(row.get(celltype_field) or row.get("cell_type") or row.get("celltype") or row.get("cluster"))
            row_donors.append(donor)
            row_celltypes.append(celltype)
            donor_cell_counts[donor] += 1
            donor_celltype_counts[donor][celltype] += 1
            celltypes.add(celltype)
            for field in LABEL_FIELDS:
                if field in row:
                    donor_label_counts[donor][f"label__{safe_name(field)}::{clean_label(row.get(field))}"] += 1
    return row_donors, row_celltypes, donor_cell_counts, donor_label_counts, donor_celltype_counts, celltypes


def stream_sparse_features(
    expression_path: Path,
    row_donors: list[str],
    row_celltypes: list[str],
) -> tuple[set[str], dict[tuple[str, str], float], dict[tuple[str, str], int], dict[tuple[str, str], float], dict[tuple[str, str, str], float]]:
    lookup = gene_to_index_lookup()
    genes_seen: set[str] = set()
    donor_gene_sum: dict[tuple[str, str], float] = defaultdict(float)
    donor_gene_nonzero: dict[tuple[str, str], int] = defaultdict(int)
    donor_index_sum: dict[tuple[str, str], float] = defaultdict(float)
    donor_celltype_index_sum: dict[tuple[str, str, str], float] = defaultdict(float)
    with open_text(expression_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_index = int(row["row_index"])
            if row_index >= len(row_donors):
                continue
            donor = row_donors[row_index]
            celltype = row_celltypes[row_index]
            gene = row.get("gene_symbol", "")
            value = float(row.get("expression_value", "0"))
            genes_seen.add(gene)
            donor_gene_sum[(donor, gene)] += value
            donor_gene_nonzero[(donor, gene)] += 1
            for index_name in lookup.get(gene, []):
                donor_index_sum[(donor, index_name)] += value
                donor_celltype_index_sum[(donor, celltype, index_name)] += value
    return genes_seen, donor_gene_sum, donor_gene_nonzero, donor_index_sum, donor_celltype_index_sum


def present_gene_count(index_name: str, genes_seen: set[str]) -> int:
    return max(1, len([gene for gene in DONOR_INDICES[index_name] if gene in genes_seen]))


def build_rows(
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
    donors = sorted(donor_cell_counts)
    genes = sorted(genes_seen)
    sorted_celltypes = sorted(celltypes)
    label_columns = sorted({combined.split("::", 1)[0] for counts in donor_label_counts.values() for combined in counts})
    fieldnames = [
        "donor_id",
        "n_cells",
        *label_columns,
        *[f"gene_mean__{safe_name(gene)}" for gene in genes],
        *[f"gene_detection__{safe_name(gene)}" for gene in genes],
        *[f"index__{safe_name(index_name)}" for index_name in sorted(DONOR_INDICES)],
        *[f"cell_fraction__{safe_name(celltype)}" for celltype in sorted_celltypes],
        *[
            f"celltype_index__{safe_name(index_name)}__{safe_name(celltype)}"
            for index_name in sorted(DONOR_INDICES)
            for celltype in sorted_celltypes
        ],
    ]
    rows: list[dict[str, str]] = []
    for donor in donors:
        cell_count = donor_cell_counts[donor]
        row = {"donor_id": donor, "n_cells": str(cell_count)}
        for column in label_columns:
            choices = Counter()
            for combined, count in donor_label_counts[donor].items():
                field, label = combined.split("::", 1)
                if field == column:
                    choices[label] += count
            row[column] = mode_label(choices)
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
                value = donor_celltype_index_sum.get((donor, celltype, index_name), 0.0) / denominator if denominator else 0.0
                row[f"celltype_index__{safe_name(index_name)}__{safe_name(celltype)}"] = f"{value:.8g}"
        rows.append(row)
    return rows, fieldnames


def schema_alignment_rows(mathys_fieldnames: list[str], phase5_schema_path: Path) -> list[dict[str, str]]:
    if not phase5_schema_path.exists():
        return [
            {
                "feature": "phase5_schema_missing",
                "status": "not_checked",
                "notes": f"SEA-AD Phase 5 table not found: {phase5_schema_path}",
            }
        ]
    phase5_rows = read_tsv(phase5_schema_path)
    phase5_fields = set(phase5_rows[0]) if phase5_rows else set()
    mathys_fields = set(mathys_fieldnames)
    rows: list[dict[str, str]] = []
    for feature in sorted(phase5_fields | mathys_fields):
        rows.append(
            {
                "feature": feature,
                "status": "shared" if feature in phase5_fields and feature in mathys_fields else "schema_specific",
                "notes": "align before cross-cohort validation",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Mathys donor-level NeuroFate feature table.")
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/interim/external/mathys_2019/mathys_obs_metadata_decoded.tsv"),
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/interim/external/mathys_2019/sparse_gene_panel_expression.tsv.gz"),
    )
    parser.add_argument("--donor-field", default="donor_id")
    parser.add_argument("--celltype-field", default="cell_type")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"),
    )
    parser.add_argument(
        "--schema-output",
        type=Path,
        default=Path("results/tables/mathys_2019_feature_schema_alignment.tsv"),
    )
    parser.add_argument(
        "--phase5-schema",
        type=Path,
        default=Path("results/tables/phase5_donor_feature_table.tsv"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/39_build_mathys_donor_feature_table.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    metadata = load_metadata(args.metadata, args.donor_field, args.celltype_field)
    row_donors, row_celltypes, donor_cell_counts, donor_label_counts, donor_celltype_counts, celltypes = metadata
    sparse = stream_sparse_features(args.expression, row_donors, row_celltypes)
    genes_seen, donor_gene_sum, donor_gene_nonzero, donor_index_sum, donor_celltype_index_sum = sparse
    rows, fieldnames = build_rows(
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
    write_tsv(args.output, rows, fieldnames)
    write_tsv(args.schema_output, schema_alignment_rows(fieldnames, args.phase5_schema), ["feature", "status", "notes"])
    logging.info("Mathys donor feature table prepared from sparse target-gene TSV only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
