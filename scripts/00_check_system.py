#!/usr/bin/env python3
"""LIGHTWEIGHT system check for the NeuroFate skeleton.

This script performs import availability checks only. It does not download data,
open large files, run Scanpy-based analysis, initialize PyTorch models, or train anything.
"""

from __future__ import annotations

import argparse
import importlib.util
import logging
import platform
import sys
from pathlib import Path


REQUIRED_MODULES = [
    "anndata",
    "h5py",
    "matplotlib",
    "networkx",
    "numpy",
    "pandas",
    "polars",
    "pyarrow",
    "yaml",
    "scanpy",
    "sklearn",
    "scipy",
    "zarr",
]

OPTIONAL_MODULES = ["torch"]


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


def module_available(module_name: str) -> bool:
    return importlib.util.find_spec(module_name) is not None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT NeuroFate system check.")
    parser.add_argument("--config", type=Path, default=Path("configs/project_config.yaml"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/00_check_system.log"))
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report checks only. This is the intended mode for project initiation.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting LIGHTWEIGHT system check.")
    logging.info("Dry run: %s", args.dry_run)
    logging.info("Python executable: %s", sys.executable)
    logging.info("Python version: %s", sys.version.replace("\n", " "))
    logging.info("Platform: %s", platform.platform())
    logging.info("Machine: %s", platform.machine())

    if not (sys.version_info >= (3, 11) and sys.version_info < (3, 13)):
        logging.warning("Python 3.11 or 3.12 is recommended.")

    if not args.config.exists():
        logging.warning("Config file not found: %s", args.config)
    else:
        logging.info("Config file found: %s", args.config)

    missing_required: list[str] = []
    for module_name in REQUIRED_MODULES:
        if module_available(module_name):
            logging.info("Required module available: %s", module_name)
        else:
            missing_required.append(module_name)
            logging.warning("Required module missing: %s", module_name)

    for module_name in OPTIONAL_MODULES:
        if module_available(module_name):
            logging.info("Optional module available: %s", module_name)
        else:
            logging.info("Optional module not installed: %s", module_name)

    logging.info("No dataset processing, downloads, or training were executed.")
    if missing_required:
        logging.warning("Missing required modules: %s", ", ".join(missing_required))
        logging.warning("Create/activate the environment before analysis steps.")
        return 2

    logging.info("System check completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
