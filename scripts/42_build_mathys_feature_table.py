#!/usr/bin/env python3
"""Build Mathys 2019 donor/sample-level NeuroFate feature table from CSV-derived sparse data."""

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

MATHYS_CELLTYPE_FIELD = "oupSample.cellType"
MATHYS_BATCH_COND_FIELD = "oupSample.batchCond"
MATHYS_SUBCLUST_ID_FIELD = "oupSample.subclustID"
MATHYS_SUBCLUST_COND_FIELD = "oupSample.subclustCond"

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


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")
    return cleaned or "missing"


def clean_label(value: str | None) -> str:
    value = (value or "").strip()
    return value if value else "missing"


def mode_label(counter: Counter[str]) -> str:
    if not counter:
        return "missing"
    return sorted(counter.items(), key=lambda item: (-item[1], item[0]))[0][0]


def infer_sample_id(cell_id: str) -> str:
    match = re.search(r"(AD\d+[_-]AD\d+|AD\d+|CT\d+|Control\d+|Ctrl\d+)", cell_id, flags=re.IGNORECASE)
    if match:
        return match.group(1).replace("-", "_")
    pieces = re.split(r"[_:\-.]", cell_id)
    return pieces[-1] if len(pieces) > 1 and pieces[-1] else "mathys_pseudo_donor"


def diagnosis_from_covariates(row: dict[str, str]) -> str:
    for field in [MATHYS_BATCH_COND_FIELD, MATHYS_SUBCLUST_COND_FIELD]:
        value = clean_label(row.get(field))
        if value != "missing":
            lowered = value.lower()
            if "ad" in lowered:
                return "AD"
            if "ct" in lowered or "control" in lowered:
                return "Control"
            return value
    return "missing"


def gene_to_index_lookup() -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = defaultdict(list)
    for index_name, genes in DONOR_INDICES.items():
        for gene in genes:
            lookup[gene].append(index_name)
    return lookup


def load_covariates(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, Counter[str]], dict[str, Counter[str]], dict[str, Counter[str]], set[str]]:
    cell_to_metadata: dict[str, dict[str, str]] = {}
    sample_cell_counts: dict[str, Counter[str]] = defaultdict(Counter)
    sample_label_counts: dict[str, Counter[str]] = defaultdict(Counter)
    sample_celltype_counts: dict[str, Counter[str]] = defaultdict(Counter)
    celltypes: set[str] = set()
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Mathys covariates file has no header: {path}")
        cell_id_field = reader.fieldnames[0]
        for row in reader:
            cell_id = row.get(cell_id_field, "")
            sample_id = infer_sample_id(cell_id)
            diagnosis = diagnosis_from_covariates(row)
            cell_type = clean_label(row.get(MATHYS_CELLTYPE_FIELD))
            subcluster = clean_label(row.get(MATHYS_SUBCLUST_ID_FIELD))
            cell_to_metadata[cell_id] = {
                "sample_id": sample_id,
                "diagnosis": diagnosis,
                "cell_type": cell_type,
                "subcluster": subcluster,
            }
            sample_cell_counts[sample_id]["n_cells"] += 1
            sample_label_counts[sample_id][f"label__diagnosis::{diagnosis}"] += 1
            sample_label_counts[sample_id][f"label__subcluster::{subcluster}"] += 1
            sample_celltype_counts[sample_id][cell_type] += 1
            celltypes.add(cell_type)
    return cell_to_metadata, sample_cell_counts, sample_label_counts, sample_celltype_counts, celltypes


def stream_expression(
    expression_path: Path,
    cell_to_metadata: dict[str, dict[str, str]],
) -> tuple[set[str], dict[tuple[str, str], float], dict[tuple[str, str], int], dict[tuple[str, str], float], dict[tuple[str, str, str], float], int]:
    lookup = gene_to_index_lookup()
    genes_seen: set[str] = set()
    sample_gene_sum: dict[tuple[str, str], float] = defaultdict(float)
    sample_gene_nonzero: dict[tuple[str, str], int] = defaultdict(int)
    sample_index_sum: dict[tuple[str, str], float] = defaultdict(float)
    sample_celltype_index_sum: dict[tuple[str, str, str], float] = defaultdict(float)
    missing_cells = 0
    with open_text(expression_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cell_id = row["cell_id"]
            metadata = cell_to_metadata.get(cell_id)
            if metadata is None:
                missing_cells += 1
                continue
            sample_id = metadata["sample_id"]
            cell_type = metadata["cell_type"]
            gene = row["gene_symbol"]
            value = float(row["expression_value"])
            genes_seen.add(gene)
            sample_gene_sum[(sample_id, gene)] += value
            sample_gene_nonzero[(sample_id, gene)] += 1
            for index_name in lookup.get(gene, []):
                sample_index_sum[(sample_id, index_name)] += value
                sample_celltype_index_sum[(sample_id, cell_type, index_name)] += value
    return genes_seen, sample_gene_sum, sample_gene_nonzero, sample_index_sum, sample_celltype_index_sum, missing_cells


def present_gene_count(index_name: str, genes_seen: set[str]) -> int:
    return max(1, len([gene for gene in DONOR_INDICES[index_name] if gene in genes_seen]))


def build_feature_rows(
    sample_cell_counts: dict[str, Counter[str]],
    sample_label_counts: dict[str, Counter[str]],
    sample_celltype_counts: dict[str, Counter[str]],
    celltypes: set[str],
    genes_seen: set[str],
    sample_gene_sum: dict[tuple[str, str], float],
    sample_gene_nonzero: dict[tuple[str, str], int],
    sample_index_sum: dict[tuple[str, str], float],
    sample_celltype_index_sum: dict[tuple[str, str, str], float],
) -> tuple[list[dict[str, str]], list[str]]:
    genes = sorted(genes_seen)
    sorted_celltypes = sorted(celltypes)
    fieldnames = [
        "donor_id",
        "n_cells",
        "label__diagnosis",
        "label__Overall_AD_neuropathological_Change",
        "label__subcluster",
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
    for sample_id in sorted(sample_cell_counts):
        n_cells = sample_cell_counts[sample_id]["n_cells"]
        diagnosis = mode_label(Counter({key.split("::", 1)[1]: count for key, count in sample_label_counts[sample_id].items() if key.startswith("label__diagnosis::")}))
        row = {
            "donor_id": sample_id,
            "n_cells": str(n_cells),
            "label__diagnosis": diagnosis,
            "label__Overall_AD_neuropathological_Change": "high" if diagnosis == "AD" else "none_or_low",
            "label__subcluster": mode_label(Counter({key.split("::", 1)[1]: count for key, count in sample_label_counts[sample_id].items() if key.startswith("label__subcluster::")})),
        }
        for gene in genes:
            row[f"gene_mean__{safe_name(gene)}"] = f"{sample_gene_sum.get((sample_id, gene), 0.0) / n_cells:.8g}"
            row[f"gene_detection__{safe_name(gene)}"] = f"{sample_gene_nonzero.get((sample_id, gene), 0) / n_cells:.8g}"
        for index_name in sorted(DONOR_INDICES):
            denominator = n_cells * present_gene_count(index_name, genes_seen)
            row[f"index__{safe_name(index_name)}"] = f"{sample_index_sum.get((sample_id, index_name), 0.0) / denominator:.8g}"
        for celltype in sorted_celltypes:
            celltype_count = sample_celltype_counts[sample_id].get(celltype, 0)
            row[f"cell_fraction__{safe_name(celltype)}"] = f"{celltype_count / n_cells:.8g}"
            for index_name in sorted(DONOR_INDICES):
                denominator = celltype_count * present_gene_count(index_name, genes_seen)
                value = sample_celltype_index_sum.get((sample_id, celltype, index_name), 0.0) / denominator if denominator else 0.0
                row[f"celltype_index__{safe_name(index_name)}__{safe_name(celltype)}"] = f"{value:.8g}"
        rows.append(row)
    return rows, fieldnames


def schema_alignment(rows: list[dict[str, str]], phase5_schema_path: Path) -> list[dict[str, str]]:
    mathys_fields = set(rows[0]) if rows else set()
    phase5_rows = read_tsv(phase5_schema_path) if phase5_schema_path.exists() else []
    phase5_fields = set(phase5_rows[0]) if phase5_rows else set()
    all_fields = sorted(mathys_fields | phase5_fields)
    return [
        {
            "feature": field,
            "status": "shared" if field in mathys_fields and field in phase5_fields else "schema_specific",
            "notes": "Mathys CSV feature alignment to SEA-AD Phase 5 schema",
        }
        for field in all_fields
    ]


def label_summary(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary: Counter[tuple[str, str]] = Counter()
    for row in rows:
        for field in ["label__diagnosis", "label__Overall_AD_neuropathological_Change", "label__subcluster"]:
            summary[(field, row.get(field, "missing"))] += 1
    return [
        {"label_field": field, "label": label, "sample_count": str(count)}
        for (field, label), count in sorted(summary.items())
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Mathys CSV donor/sample-level feature table.")
    parser.add_argument("--expression", type=Path, default=Path("data/interim/external/mathys_2019/mathys_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--covariates", type=Path, default=Path("data/raw/external/mathys_2019/GSE138852_covariates.csv.gz"))
    parser.add_argument("--phase5-schema", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"))
    parser.add_argument("--schema-output", type=Path, default=Path("results/tables/mathys_2019_feature_schema_alignment.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/mathys_2019_label_summary.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/42_build_mathys_feature_table.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    cell_to_metadata, sample_cell_counts, sample_label_counts, sample_celltype_counts, celltypes = load_covariates(args.covariates)
    sparse = stream_expression(args.expression, cell_to_metadata)
    genes_seen, sample_gene_sum, sample_gene_nonzero, sample_index_sum, sample_celltype_index_sum, missing_cells = sparse
    rows, fieldnames = build_feature_rows(
        sample_cell_counts,
        sample_label_counts,
        sample_celltype_counts,
        celltypes,
        genes_seen,
        sample_gene_sum,
        sample_gene_nonzero,
        sample_index_sum,
        sample_celltype_index_sum,
    )
    if len(rows) == 1 and rows[0]["donor_id"] == "mathys_pseudo_donor":
        logging.warning("Could not infer sample IDs; using a cohort-level pseudo-donor.")
    if missing_cells:
        logging.warning("Sparse expression rows with no covariate match: %d", missing_cells)
    write_tsv(args.output, rows, fieldnames)
    write_tsv(args.schema_output, schema_alignment(rows, args.phase5_schema), ["feature", "status", "notes"])
    write_tsv(args.label_summary_output, label_summary(rows), ["label_field", "label", "sample_count"])
    logging.info("Mathys CSV feature table complete. No AnnData/H5AD object was created.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
