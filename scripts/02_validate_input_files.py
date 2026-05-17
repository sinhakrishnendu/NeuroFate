#!/usr/bin/env python3
"""LIGHTWEIGHT input path validation for registered datasets.

This script checks path existence and basic file sizes only. It never opens HDF5,
AnnData, zarr, MTX, parquet, or other large biological data files.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT registered input validation.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/dataset_registry.tsv"))
    parser.add_argument("--max-lightweight-file-mb", type=float, default=50.0)
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/02_validate_input_files.log"))
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report checks only. No files are opened beyond metadata stat calls.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting LIGHTWEIGHT input file validation.")
    logging.info("Dry run: %s", args.dry_run)
    logging.info("Registry: %s", args.registry)

    if not args.registry.exists():
        logging.error("Registry not found: %s", args.registry)
        return 1

    missing_paths = 0
    placeholder_paths = 0
    oversized_paths = 0

    with args.registry.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            dataset_id = row.get("dataset_id", "")
            local_path = row.get("local_path", "")
            if "PLACEHOLDER" in local_path:
                placeholder_paths += 1
                logging.info("Placeholder path for %s: %s", dataset_id, local_path)
                continue

            path = Path(local_path)
            if not path.exists():
                missing_paths += 1
                logging.warning("Missing path for %s: %s", dataset_id, path)
                continue

            if path.is_file():
                size_mb = path.stat().st_size / (1024 * 1024)
                logging.info("Found file for %s: %s (%.2f MB)", dataset_id, path, size_mb)
                if size_mb > args.max_lightweight_file_mb:
                    oversized_paths += 1
                    logging.warning(
                        "File exceeds lightweight threshold and was not opened: %s (%.2f MB)",
                        path,
                        size_mb,
                    )
            elif path.is_dir():
                logging.info("Found directory for %s: %s", dataset_id, path)
            else:
                logging.warning("Path exists but is neither file nor directory: %s", path)

    logging.info("Placeholder paths: %d", placeholder_paths)
    logging.info("Missing paths: %d", missing_paths)
    logging.info("Oversized files not opened: %d", oversized_paths)
    logging.info("No dataset contents were read or processed.")
    return 0 if missing_paths == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
