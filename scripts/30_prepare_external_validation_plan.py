#!/usr/bin/env python3
"""Prepare a lightweight external validation plan for NeuroFate.

This script reads only the external validation registry and target gene panel.
It does not open expression files or process external datasets.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


REGISTRY_COLUMNS = [
    "dataset_id",
    "disease",
    "modality",
    "species",
    "brain_region",
    "cells_or_donors",
    "accession",
    "download_status",
    "processing_status",
    "notes",
]

REQUIRED_METADATA_FIELDS = {
    "donor_id": "donor; individual; subject; specimen_id",
    "cohort_id": "cohort; dataset_id; study",
    "diagnosis": "diagnosis; disease_status; clinical_diagnosis",
    "cognitive_status": "cognitive_status; dementia; clinical_dementia_rating",
    "ad_pathology": "braak; cerad; overall_ad_neuropathological_change",
    "cell_type": "cell_type; subclass; cluster; annotation",
    "brain_region": "brain_region; tissue; dissection",
    "age": "age; age_at_death",
    "sex": "sex; gender",
    "apoe_genotype": "apoe; apoe_genotype",
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_columns(rows: list[dict[str, str]], required: list[str], label: str) -> None:
    if not rows:
        raise RuntimeError(f"{label} has no rows.")
    missing = [column for column in required if column not in rows[0]]
    if missing:
        raise RuntimeError(f"{label} missing columns: {', '.join(missing)}")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def build_gene_overlap_rows(
    registry_rows: list[dict[str, str]],
    panel_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset in registry_rows:
        for gene in panel_rows:
            rows.append(
                {
                    "dataset_id": dataset["dataset_id"],
                    "gene_symbol": gene["gene_symbol"],
                    "priority_tier": gene.get("priority_tier", ""),
                    "panel_status": "required_target_panel_gene",
                    "external_gene_status": "unknown_until_var_metadata_available",
                    "harmonization_action": "compare external gene symbols and gene IDs before extraction",
                    "notes": "planning only; no expression values loaded",
                }
            )
    return rows


def metadata_status(dataset: dict[str, str], field: str) -> str:
    modality = dataset["modality"].lower()
    disease = dataset["disease"].lower()
    if field == "cell_type" and "bulk" in modality and "single" not in modality:
        return "not_required_for_bulk_donor_level_validation"
    if field == "ad_pathology" and "alzheimer" in disease:
        return "required"
    if field == "apoe_genotype" and "alzheimer" in disease:
        return "preferred"
    if field == "diagnosis":
        return "required"
    if field in {"donor_id", "cohort_id", "brain_region"}:
        return "required"
    return "preferred"


def build_metadata_overlap_rows(registry_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for dataset in registry_rows:
        for field, aliases in REQUIRED_METADATA_FIELDS.items():
            rows.append(
                {
                    "dataset_id": dataset["dataset_id"],
                    "required_field": field,
                    "expected_aliases": aliases,
                    "compatibility_status": metadata_status(dataset, field),
                    "harmonization_action": "map source column to canonical NeuroFate field before modeling",
                    "notes": "planning only; validate once metadata TSV exists",
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare external validation planning tables.")
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("metadata/external_validation_registry.tsv"),
    )
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/30_prepare_external_validation_plan.log"),
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting lightweight Phase 7 external validation planning.")
    registry_rows = read_tsv(args.registry)
    panel_rows = read_tsv(args.panel)
    validate_columns(registry_rows, REGISTRY_COLUMNS, "external validation registry")
    validate_columns(panel_rows, ["gene_symbol"], "target gene panel")

    gene_rows = build_gene_overlap_rows(registry_rows, panel_rows)
    metadata_rows = build_metadata_overlap_rows(registry_rows)
    logging.info("External cohorts planned: %d", len(registry_rows))
    logging.info("Target genes planned for overlap checks: %d", len(panel_rows))

    if args.dry_run:
        logging.info("Dry run requested; planning tables were validated but not written.")
        return 0

    write_tsv(
        args.tables_dir / "external_validation_gene_overlap.tsv",
        gene_rows,
        [
            "dataset_id",
            "gene_symbol",
            "priority_tier",
            "panel_status",
            "external_gene_status",
            "harmonization_action",
            "notes",
        ],
    )
    write_tsv(
        args.tables_dir / "external_validation_metadata_overlap.tsv",
        metadata_rows,
        [
            "dataset_id",
            "required_field",
            "expected_aliases",
            "compatibility_status",
            "harmonization_action",
            "notes",
        ],
    )
    logging.info("No external expression file was opened or processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
