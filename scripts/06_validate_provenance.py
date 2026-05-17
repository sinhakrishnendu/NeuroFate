#!/usr/bin/env python3
"""LIGHTWEIGHT provenance template validation for NeuroFate Phase 1C.

This script validates provenance metadata only. It never accesses remote URLs,
opens large files, computes checksums over real files, or processes datasets.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


PROVENANCE_COLUMNS = [
    "provenance_id",
    "dataset_id",
    "source_name",
    "source_url_or_accession",
    "source_database",
    "download_status",
    "download_command_manual_only",
    "date_accessed",
    "license_or_terms",
    "original_file_name",
    "local_expected_path",
    "checksum_algorithm",
    "checksum_value",
    "file_size_expected",
    "file_size_observed",
    "verified",
    "notes",
]

COMPLETE_STATUSES = {"complete", "completed", "downloaded", "verified"}
MISSING_CHECKSUM_VALUES = {"", "pending", "not_computed", "na", "n/a"}


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT NeuroFate provenance validator.")
    parser.add_argument(
        "--provenance-template",
        type=Path,
        default=Path("metadata/provenance_template.tsv"),
    )
    parser.add_argument("--dataset-registry", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/06_validate_provenance.log"))
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate metadata only. This script never computes checksums in this phase.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting LIGHTWEIGHT Phase 1C provenance validation.")
    logging.info("Dry run: %s", args.dry_run)
    logging.info("Provenance template: %s", args.provenance_template)
    logging.info("Dataset registry: %s", args.dataset_registry)

    errors: list[str] = []
    try:
        provenance_columns, provenance_rows = read_tsv(args.provenance_template)
    except FileNotFoundError:
        logging.error("Provenance template not found: %s", args.provenance_template)
        return 2

    try:
        _, dataset_rows = read_tsv(args.dataset_registry)
    except FileNotFoundError:
        logging.error("Dataset registry not found: %s", args.dataset_registry)
        return 2

    missing_columns = [column for column in PROVENANCE_COLUMNS if column not in provenance_columns]
    if missing_columns:
        errors.append("provenance_template missing columns: " + ", ".join(missing_columns))

    dataset_ids = {row.get("dataset_id", "") for row in dataset_rows}
    provenance_dataset_ids = {row.get("dataset_id", "") for row in provenance_rows}
    missing_dataset_ids = sorted(dataset_ids - provenance_dataset_ids)
    if missing_dataset_ids:
        errors.append(
            "provenance_template missing dataset_id rows: " + ", ".join(missing_dataset_ids)
        )

    seen_provenance_ids: set[str] = set()
    for index, row in enumerate(provenance_rows, start=2):
        provenance_id = row.get("provenance_id", "").strip()
        dataset_id = row.get("dataset_id", "").strip()
        download_status = row.get("download_status", "").strip().lower()
        checksum_algorithm = row.get("checksum_algorithm", "").strip().lower()
        checksum_value = row.get("checksum_value", "").strip().lower()
        file_size_observed = row.get("file_size_observed", "").strip().lower()
        verified = row.get("verified", "").strip().lower()
        command = row.get("download_command_manual_only", "").strip()

        if not provenance_id:
            errors.append(f"provenance_template line {index}: empty provenance_id")
        elif provenance_id in seen_provenance_ids:
            errors.append(f"provenance_template line {index}: duplicate provenance_id {provenance_id}")
        seen_provenance_ids.add(provenance_id)

        if dataset_id not in dataset_ids:
            errors.append(
                f"provenance_template line {index}: dataset_id not in registry: {dataset_id}"
            )

        if not command.startswith("MANUAL_ONLY:"):
            errors.append(
                f"provenance_template line {index}: download command must start with MANUAL_ONLY:"
            )

        if download_status in COMPLETE_STATUSES:
            if checksum_algorithm not in {"sha256", "sha512", "md5"}:
                errors.append(
                    f"provenance_template line {index}: complete row lacks checksum algorithm"
                )
            if checksum_value in MISSING_CHECKSUM_VALUES:
                errors.append(
                    f"provenance_template line {index}: complete row lacks checksum value"
                )
            if file_size_observed in {"", "pending", "na", "n/a"}:
                errors.append(
                    f"provenance_template line {index}: complete row lacks observed file size"
                )
            if verified != "true":
                errors.append(
                    f"provenance_template line {index}: complete row must have verified=true"
                )

    logging.info("Provenance rows: %d", len(provenance_rows))
    logging.info("Dataset registry rows: %d", len(dataset_rows))
    logging.info("Provenance columns: %s", ", ".join(provenance_columns))
    logging.info("Errors: %d", len(errors))

    if errors:
        for error in errors:
            logging.error(error)
        logging.error("Provenance validation failed.")
        logging.info("No remote URLs were accessed and no files were opened.")
        return 2

    logging.info("Provenance validation passed.")
    logging.info("No remote URLs were accessed.")
    logging.info("No files were opened for checksum or content validation.")
    logging.info("No datasets were downloaded or processed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
