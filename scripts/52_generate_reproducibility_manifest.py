#!/usr/bin/env python3
"""Generate a NeuroFate reproducibility manifest."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from pathlib import Path


PACKAGE_NAMES = ["numpy", "scipy", "pandas", "polars", "pyarrow", "h5py", "zarr", "scikit-learn", "matplotlib", "networkx", "pyyaml", "torch"]
INPUT_PATTERNS = ["data/raw/**/*", "data/interim/**/*", "metadata/*.tsv", "configs/*.yaml"]
OUTPUT_PATTERNS = ["results/tables/*", "results/figures/*", "results/models/*", "results/reports/*"]


def package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for package in PACKAGE_NAMES:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    return versions


def git_commit() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "not_available"
    return result.stdout.strip()


def sha256(path: Path, size_limit_mb: int) -> str:
    if path.stat().st_size > size_limit_mb * 1024 * 1024:
        return f"skipped_file_larger_than_{size_limit_mb}_mb"
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inventory(patterns: list[str], checksum_limit_mb: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for pattern in patterns:
        for path in sorted(Path(".").glob(pattern)):
            if not path.is_file():
                continue
            rows.append(
                {
                    "path": str(path),
                    "size_bytes": str(path.stat().st_size),
                    "sha256": sha256(path, checksum_limit_mb),
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate NeuroFate reproducibility manifest.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/reproducibility_manifest.json"))
    parser.add_argument("--checksum-size-limit-mb", type=int, default=512)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = {
        "python_version": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "package_versions": package_versions(),
        "git_commit": git_commit(),
        "input_file_inventory": inventory(INPUT_PATTERNS, args.checksum_size_limit_mb),
        "output_inventory": inventory(OUTPUT_PATTERNS, args.checksum_size_limit_mb),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
