#!/usr/bin/env python3
"""Create a clean NeuroFate source-code release ZIP."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


INCLUDE_PATHS = [
    "neurofate",
    "scripts",
    "tests",
    "configs",
    "metadata",
    "docs",
    "examples/tiny_demo",
    ".github/workflows/ci.yml",
    "README.md",
    "LICENSE",
    "PYPI_RELEASE_CHECKLIST.md",
    "environment.yml",
    "CITATION.cff",
    "codemeta.json",
    "CHANGELOG.md",
    "RELEASE_CHECKLIST.md",
    "CONTRIBUTING.md",
    "CODE_OF_CONDUCT.md",
    "MANIFEST.in",
    "pyproject.toml",
    "RESULTS_INTERPRETATION.md",
]

EXCLUDED_PREFIXES = [
    "data/",
    "results/models/",
    "results/figures/",
    "results/tables/",
    "results/logs/",
    "dist/",
    "release_artifacts/",
]
EXCLUDED_PARTS = {"__pycache__", ".DS_Store"}
EXCLUDED_SUFFIXES = [".h5ad", ".pt", ".pth", ".pyc"]


def should_include(path: Path) -> tuple[bool, str]:
    text = path.as_posix()
    if any(part in EXCLUDED_PARTS for part in path.parts):
        return False, "cache_or_system_file"
    if any(text.startswith(prefix) for prefix in EXCLUDED_PREFIXES):
        return False, "excluded_release_artifact_or_data"
    if any(text.endswith(suffix) for suffix in EXCLUDED_SUFFIXES):
        return False, "excluded_binary_or_raw_suffix"
    if text.endswith(".gz") and text != "examples/tiny_demo/tiny_sparse_expression.tsv.gz":
        return False, "excluded_compressed_file"
    return True, "included_source_release"


def candidate_files() -> list[Path]:
    files: list[Path] = []
    for path_text in INCLUDE_PATHS:
        path = Path(path_text)
        if path.is_dir():
            files.extend(child for child in path.rglob("*") if child.is_file())
        elif path.is_file():
            files.append(path)
    return sorted(set(files))


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "size_bytes", "included", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build NeuroFate source release package.")
    parser.add_argument("--output-dir", type=Path, default=Path("release_artifacts"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    zip_path = args.output_dir / f"neurofate_source_release_{timestamp}.zip"
    manifest_path = args.output_dir / "source_release_manifest.tsv"
    rows: list[dict[str, str]] = []
    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as archive:
        for path in candidate_files():
            include, reason = should_include(path)
            rows.append(
                {
                    "path": path.as_posix(),
                    "size_bytes": str(path.stat().st_size),
                    "included": str(include).lower(),
                    "reason": reason,
                }
            )
            if include:
                archive.write(path, path.as_posix())
    write_manifest(manifest_path, rows)
    print(f"Wrote {zip_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
