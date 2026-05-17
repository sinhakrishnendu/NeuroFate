#!/usr/bin/env python3
"""LIGHTWEIGHT dataset registry preparation.

Reads dataset metadata from a YAML file and previews or writes a small TSV registry.
This script does not download datasets or inspect data file contents.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any


FIELDS = [
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


def load_yaml(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise SystemExit("PyYAML is required for this script. Install the environment first.") from exc

    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected mapping in {path}")
    return data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT dataset registry preparation.")
    parser.add_argument("--datasets", type=Path, default=Path("configs/datasets.yaml"))
    parser.add_argument("--output", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write the registry TSV. Without this flag the script is a dry run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview only. This is the default when --write is not supplied.",
    )
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/01_prepare_dataset_registry.log"))
    return parser.parse_args()


def stringify_field(value: object) -> str:
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    effective_dry_run = args.dry_run or not args.write

    logging.info("Starting LIGHTWEIGHT dataset registry preparation.")
    logging.info("Dry run: %s", effective_dry_run)
    logging.info("Input dataset config: %s", args.datasets)
    logging.info("Output registry: %s", args.output)

    if not args.datasets.exists():
        logging.error("Dataset config not found: %s", args.datasets)
        return 1

    config = load_yaml(args.datasets)
    records = config.get("datasets", [])
    if not isinstance(records, list):
        logging.error("Expected 'datasets' to be a list.")
        return 1

    normalized: list[dict[str, str]] = []
    for record in records:
        if not isinstance(record, dict):
            logging.warning("Skipping non-mapping record: %r", record)
            continue
        normalized.append({field: stringify_field(record.get(field, "")) for field in FIELDS})

    logging.info("Prepared %d registry rows.", len(normalized))
    for row in normalized:
        logging.info("Registry row preview: %s | %s | %s", row["dataset_id"], row["modality"], row["status"])

    if effective_dry_run:
        logging.info("Dry run complete. No files were written.")
        logging.info("No dataset downloads or data inspection were executed.")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(normalized)

    logging.info("Wrote registry: %s", args.output)
    logging.info("No dataset downloads or data inspection were executed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
