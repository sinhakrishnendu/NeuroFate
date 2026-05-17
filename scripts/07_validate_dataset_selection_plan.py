#!/usr/bin/env python3
"""LIGHTWEIGHT dataset selection plan validation for NeuroFate Phase 1D.

This script validates planning TSV files only. It never downloads datasets,
accesses remote URLs, opens h5ad/HDF5 files, runs Scanpy, or processes data.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


SELECTION_COLUMNS = [
    "claim_id",
    "manuscript_claim",
    "disease_area",
    "required_dataset_id",
    "data_modality",
    "minimum_required_fields",
    "preferred_source",
    "access_type",
    "estimated_size_category",
    "manual_download_priority",
    "phase",
    "blocking_risk",
    "notes",
]

MVDP_COLUMNS = [
    "phase",
    "dataset_id",
    "why_needed",
    "minimum_file_needed",
    "can_start_without_it",
    "first_analysis_enabled",
    "manual_action_needed",
    "notes",
]

EXPECTED_CLAIMS = [
    "claim_ad_single_cell_state_reconstruction",
    "claim_pd_single_cell_state_reconstruction",
    "claim_gut_brain_metabolite_network_layer",
    "claim_ppi_network_biology_layer",
    "claim_evolutionary_conservation_layer",
    "claim_positive_selection_layer",
    "claim_multimodal_neurofate_score",
    "claim_interpretability_reporting_layer",
]


def configure_logging(log_file: Path | None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, mode="w", encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def missing_columns(fieldnames: list[str], required: list[str]) -> list[str]:
    return [column for column in required if column not in fieldnames]


def split_dataset_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def validate_unique_ids(rows: list[dict[str, str]], column: str, label: str) -> list[str]:
    errors: list[str] = []
    seen: set[str] = set()
    for index, row in enumerate(rows, start=2):
        value = row.get(column, "").strip()
        if not value:
            errors.append(f"{label} line {index}: missing {column}")
            continue
        if value in seen:
            errors.append(f"{label} line {index}: duplicate {column}: {value}")
        seen.add(value)
    return errors


def priority_key(value: str) -> tuple[int, str]:
    normalized = value.strip().upper()
    if normalized.startswith("P") and normalized[1:].isdigit():
        return int(normalized[1:]), normalized
    return 999, normalized


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT NeuroFate dataset selection validator.")
    parser.add_argument(
        "--dataset-selection-plan",
        type=Path,
        default=Path("metadata/dataset_selection_plan.tsv"),
    )
    parser.add_argument(
        "--minimum-viable-plan",
        type=Path,
        default=Path("metadata/minimum_viable_dataset_plan.tsv"),
    )
    parser.add_argument("--dataset-registry", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/07_validate_dataset_selection_plan.log"),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate planning tables only. No dataset or remote access occurs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting LIGHTWEIGHT Phase 1D dataset selection validation.")
    logging.info("Dry run: %s", args.dry_run)
    logging.info("Dataset selection plan: %s", args.dataset_selection_plan)
    logging.info("Minimum viable dataset plan: %s", args.minimum_viable_plan)
    logging.info("Dataset registry: %s", args.dataset_registry)

    errors: list[str] = []
    try:
        selection_columns, selection_rows = read_tsv(args.dataset_selection_plan)
        mvd_columns, mvd_rows = read_tsv(args.minimum_viable_plan)
        _, registry_rows = read_tsv(args.dataset_registry)
    except FileNotFoundError as exc:
        logging.error("Required planning file missing: %s", exc)
        return 2

    for label, columns, required in [
        ("dataset_selection_plan", selection_columns, SELECTION_COLUMNS),
        ("minimum_viable_dataset_plan", mvd_columns, MVDP_COLUMNS),
    ]:
        missing = missing_columns(columns, required)
        if missing:
            errors.append(f"{label} missing required columns: {', '.join(missing)}")

    errors.extend(validate_unique_ids(selection_rows, "claim_id", "dataset_selection_plan"))

    registry_ids = {row.get("dataset_id", "").strip() for row in registry_rows}
    for row in selection_rows:
        claim_id = row.get("claim_id", "")
        for dataset_id in split_dataset_ids(row.get("required_dataset_id", "")):
            if dataset_id not in registry_ids:
                errors.append(
                    f"dataset_selection_plan claim {claim_id} references unknown dataset_id: {dataset_id}"
                )

    for row in mvd_rows:
        dataset_id = row.get("dataset_id", "").strip()
        if dataset_id not in registry_ids:
            errors.append(f"minimum_viable_dataset_plan references unknown dataset_id: {dataset_id}")
        if row.get("can_start_without_it", "").strip().lower() not in {"true", "false"}:
            errors.append(
                "minimum_viable_dataset_plan "
                f"{dataset_id} has non-boolean can_start_without_it"
            )

    present_claims = {row.get("claim_id", "") for row in selection_rows}
    missing_claims = [claim for claim in EXPECTED_CLAIMS if claim not in present_claims]
    if missing_claims:
        errors.append("dataset_selection_plan missing expected claims: " + ", ".join(missing_claims))

    logging.info("Selection plan rows: %d", len(selection_rows))
    logging.info("Minimum viable dataset rows: %d", len(mvd_rows))
    logging.info("Priority order:")
    for row in sorted(selection_rows, key=lambda item: priority_key(item["manual_download_priority"])):
        logging.info(
            "  %s | %s | dataset=%s | phase=%s | risk=%s",
            row["manual_download_priority"],
            row["claim_id"],
            row["required_dataset_id"],
            row["phase"],
            row["blocking_risk"],
        )

    logging.info("Minimum viable dataset order:")
    for row in sorted(mvd_rows, key=lambda item: priority_key(item["phase"])):
        logging.info(
            "  %s | %s | can_start_without_it=%s | first_analysis=%s",
            row["phase"],
            row["dataset_id"],
            row["can_start_without_it"],
            row["first_analysis_enabled"],
        )

    logging.info("Errors: %d", len(errors))
    if errors:
        for error in errors:
            logging.error(error)
        logging.error("Dataset selection plan validation failed.")
        logging.info("No downloads, remote access, dataset opening, or processing occurred.")
        return 2

    logging.info("Dataset selection plan validation passed.")
    logging.info("No downloads, remote access, dataset opening, or processing occurred.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
