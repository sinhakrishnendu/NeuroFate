#!/usr/bin/env python3
"""LIGHTWEIGHT real dataset source validation for NeuroFate Phase 2A.

This script validates local planning metadata and manual download templates only.
It never accesses remote URLs, downloads data, opens dataset files, reads h5ad/HDF5,
or executes manual download scripts.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


SOURCE_COLUMNS = [
    "dataset_id",
    "real_source_name",
    "official_resource",
    "accession_or_resource_id",
    "access_type",
    "manual_download_required",
    "controlled_access",
    "recommended_priority",
    "expected_format",
    "expected_size_category",
    "notes",
]

EXPECTED_MANUAL_SCRIPTS = {
    "sea_ad_single_nucleus": Path("scripts/manual_downloads/download_sea_ad_manual.sh"),
    "string_ppi_placeholder": Path("scripts/manual_downloads/download_string_manual.sh"),
    "mathys_2019_ad_single_nucleus": Path(
        "scripts/manual_downloads/download_mathys_synapse_manual.sh"
    ),
    "rosmap_ad_transcriptomics": Path(
        "scripts/manual_downloads/download_rosmap_synapse_manual.sh"
    ),
}

CONTROLLED_DATASETS = {
    "mathys_2019_ad_single_nucleus",
    "rosmap_ad_transcriptomics",
}


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


def validate_script_template(path: Path) -> list[str]:
    errors: list[str] = []
    if not path.exists():
        return [f"manual download script missing: {path}"]

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    if not lines or lines[0].strip() != "set -euo pipefail":
        errors.append(f"{path}: first line must be 'set -euo pipefail'")
    for required_text in [
        "DO NOT RUN FROM CODEX",
        "MANUAL_HEAVY",
        "RUN_MANUAL_DOWNLOAD",
        'RUN_MANUAL_DOWNLOAD}" != "YES"',
    ]:
        if required_text not in text:
            errors.append(f"{path}: missing guard text {required_text!r}")
    if "data/raw/" not in text and "data/external/" not in text:
        errors.append(f"{path}: must write only under data/raw/ or data/external/")
    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT NeuroFate real source validator.")
    parser.add_argument(
        "--real-dataset-sources",
        type=Path,
        default=Path("metadata/real_dataset_sources.tsv"),
    )
    parser.add_argument("--dataset-registry", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument(
        "--manual-download-dir",
        type=Path,
        default=Path("scripts/manual_downloads"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/08_validate_real_dataset_sources.log"),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate local metadata and script templates only.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting LIGHTWEIGHT Phase 2A real dataset source validation.")
    logging.info("Dry run: %s", args.dry_run)
    logging.info("Real dataset sources: %s", args.real_dataset_sources)
    logging.info("Dataset registry: %s", args.dataset_registry)
    logging.info("Manual download directory: %s", args.manual_download_dir)

    errors: list[str] = []
    try:
        source_columns, source_rows = read_tsv(args.real_dataset_sources)
        _, registry_rows = read_tsv(args.dataset_registry)
    except FileNotFoundError as exc:
        logging.error("Required local metadata file missing: %s", exc)
        return 2

    missing = missing_columns(source_columns, SOURCE_COLUMNS)
    if missing:
        errors.append("real_dataset_sources missing required columns: " + ", ".join(missing))

    registry_ids = {row.get("dataset_id", "").strip() for row in registry_rows}
    source_ids = {row.get("dataset_id", "").strip() for row in source_rows}
    for dataset_id in sorted(source_ids):
        if dataset_id not in registry_ids:
            errors.append(f"real_dataset_sources references unknown dataset_id: {dataset_id}")

    for dataset_id in CONTROLLED_DATASETS:
        row = next((item for item in source_rows if item.get("dataset_id") == dataset_id), None)
        if row is None:
            errors.append(f"controlled-access dataset missing from source table: {dataset_id}")
            continue
        if row.get("controlled_access", "").strip().lower() != "true":
            errors.append(f"{dataset_id} must have controlled_access=true")

    for dataset_id, script_path in EXPECTED_MANUAL_SCRIPTS.items():
        if dataset_id not in source_ids:
            errors.append(f"manual script expected for dataset missing from source table: {dataset_id}")
        errors.extend(validate_script_template(script_path))

    logging.info("Real dataset source rows: %d", len(source_rows))
    logging.info("Manual download scripts expected: %d", len(EXPECTED_MANUAL_SCRIPTS))
    logging.info("Priority order:")
    for row in sorted(source_rows, key=lambda item: item.get("recommended_priority", "")):
        logging.info(
            "  %s | %s | %s | controlled_access=%s",
            row.get("recommended_priority", ""),
            row.get("dataset_id", ""),
            row.get("real_source_name", ""),
            row.get("controlled_access", ""),
        )

    logging.info("Errors: %d", len(errors))
    if errors:
        for error in errors:
            logging.error(error)
        logging.error("Real dataset source validation failed.")
        logging.info("No remote URLs were accessed and no downloads were attempted.")
        return 2

    logging.info("Real dataset source validation passed.")
    logging.info("No remote URLs were accessed.")
    logging.info("No downloads were attempted.")
    logging.info("No dataset files, h5ad files, or HDF5 files were opened.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
