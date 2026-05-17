#!/usr/bin/env python3
"""LIGHTWEIGHT registry validation for NeuroFate Phase 1B.

This script validates metadata TSV schemas and placeholder path rules only.
It never downloads data, opens h5ad/HDF5 files, runs Scanpy, or processes datasets.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


DATASET_REQUIRED_COLUMNS = [
    "dataset_id",
    "disease_area",
    "modality",
    "source",
    "expected_file_type",
    "local_path",
    "status",
    "heavy_to_download",
    "notes",
]

FEATURE_REQUIRED_COLUMNS = [
    "feature_id",
    "feature_group",
    "biological_layer",
    "expected_input",
    "expected_output",
    "status",
    "notes",
]

MODULE_MAP_REQUIRED_COLUMNS = [
    "module_id",
    "manuscript_section",
    "scientific_layer",
    "expected_data_source",
    "expected_pipeline_script",
    "expected_result_table",
    "expected_result_figure",
    "status",
    "notes",
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
        fieldnames = reader.fieldnames or []
        rows = [{key: value for key, value in row.items()} for row in reader]
    return fieldnames, rows


def missing_columns(fieldnames: list[str], required_columns: list[str]) -> list[str]:
    return [column for column in required_columns if column not in fieldnames]


def validate_unique_ids(rows: list[dict[str, str]], id_column: str, label: str) -> list[str]:
    errors: list[str] = []
    seen: set[str] = set()
    for index, row in enumerate(rows, start=2):
        value = row.get(id_column, "").strip()
        if not value:
            errors.append(f"{label}: missing {id_column} at TSV line {index}")
            continue
        if value in seen:
            errors.append(f"{label}: duplicate {id_column}: {value}")
        seen.add(value)
    return errors


def validate_dataset_paths(rows: list[dict[str, str]]) -> list[str]:
    errors: list[str] = []
    for row in rows:
        dataset_id = row.get("dataset_id", "").strip()
        status = row.get("status", "").strip().lower()
        local_path = row.get("local_path", "").strip()
        if not local_path:
            errors.append(f"dataset_registry: {dataset_id} has an empty local_path")
            continue
        if status == "placeholder":
            continue
        path = Path(local_path)
        if not path.exists():
            errors.append(
                "dataset_registry: "
                f"{dataset_id} has status={status!r} but local_path is missing: {local_path}"
            )
    return errors


def validate_registry(
    path: Path,
    required_columns: list[str],
    id_column: str,
    label: str,
) -> tuple[list[dict[str, str]], list[str]]:
    errors: list[str] = []
    fieldnames, rows = read_tsv(path)
    missing = missing_columns(fieldnames, required_columns)
    if missing:
        errors.append(f"{label}: missing required columns: {', '.join(missing)}")
    errors.extend(validate_unique_ids(rows, id_column, label))
    logging.info("%s rows: %d", label, len(rows))
    logging.info("%s columns: %s", label, ", ".join(fieldnames))
    return rows, errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT NeuroFate registry validator.")
    parser.add_argument("--dataset-registry", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument("--feature-registry", type=Path, default=Path("metadata/feature_registry.tsv"))
    parser.add_argument(
        "--module-map",
        type=Path,
        default=Path("metadata/manuscript_module_map.tsv"),
    )
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/04_validate_registries.log"))
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate metadata only. This script never processes datasets in any mode.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting LIGHTWEIGHT Phase 1B registry validation.")
    logging.info("Dry run: %s", args.dry_run)
    logging.info("Dataset registry: %s", args.dataset_registry)
    logging.info("Feature registry: %s", args.feature_registry)
    logging.info("Manuscript module map: %s", args.module_map)

    errors: list[str] = []
    try:
        dataset_rows, dataset_errors = validate_registry(
            args.dataset_registry,
            DATASET_REQUIRED_COLUMNS,
            "dataset_id",
            "dataset_registry",
        )
        errors.extend(dataset_errors)
        errors.extend(validate_dataset_paths(dataset_rows))
    except FileNotFoundError as exc:
        errors.append(f"dataset_registry: file not found: {exc}")
        dataset_rows = []

    try:
        feature_rows, feature_errors = validate_registry(
            args.feature_registry,
            FEATURE_REQUIRED_COLUMNS,
            "feature_id",
            "feature_registry",
        )
        errors.extend(feature_errors)
    except FileNotFoundError as exc:
        errors.append(f"feature_registry: file not found: {exc}")
        feature_rows = []

    try:
        module_rows, module_errors = validate_registry(
            args.module_map,
            MODULE_MAP_REQUIRED_COLUMNS,
            "module_id",
            "manuscript_module_map",
        )
        errors.extend(module_errors)
    except FileNotFoundError as exc:
        errors.append(f"manuscript_module_map: file not found: {exc}")
        module_rows = []

    logging.info("Validation summary:")
    logging.info("  dataset rows: %d", len(dataset_rows))
    logging.info("  feature rows: %d", len(feature_rows))
    logging.info("  manuscript module rows: %d", len(module_rows))
    logging.info("  errors: %d", len(errors))

    if errors:
        for error in errors:
            logging.error(error)
        logging.error("Registry validation failed.")
        logging.info("No datasets were opened, downloaded, or processed.")
        return 2

    logging.info("Registry validation passed.")
    logging.info("No datasets were opened, downloaded, or processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
