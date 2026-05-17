#!/usr/bin/env python3
"""Inventory NeuroFate output artifacts without recomputing analyses."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path


SEARCH_ROOTS = [Path("results/tables"), Path("results/figures"), Path("results/models"), Path("results/reports"), Path("results/demo")]


def artifact_type(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".png", ".jpg", ".jpeg", ".svg", ".pdf"}:
        return "figure"
    if suffix in {".tsv", ".csv", ".json"}:
        return "table_or_manifest"
    if suffix in {".md", ".html", ".txt"}:
        return "report"
    if suffix in {".pt", ".pkl", ".joblib"}:
        return "model"
    return "other"


def infer_phase(path: Path) -> str:
    text = str(path).lower()
    for phase in [f"phase{index}" for index in range(1, 12)]:
        if phase in text:
            return phase
    if "figure" in text:
        return "figure"
    if "demo" in text:
        return "demo"
    if "manifest" in text or "audit" in text or "validation" in text:
        return "platform"
    return "unassigned"


def is_user_facing(path: Path) -> str:
    return str(path.suffix.lower() in {".png", ".html", ".md", ".txt", ".tsv", ".json"}).lower()


def inventory_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for root in SEARCH_ROOTS:
        if not root.exists():
            continue
        for path in sorted(root.rglob("*")):
            if not path.is_file() or path.name == ".gitkeep":
                continue
            stat = path.stat()
            rows.append(
                {
                    "path": str(path),
                    "size_bytes": str(stat.st_size),
                    "modified_time": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                    "artifact_type": artifact_type(path),
                    "phase": infer_phase(path),
                    "user_facing": is_user_facing(path),
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a NeuroFate output inventory.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/output_inventory.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["path", "size_bytes", "modified_time", "artifact_type", "phase", "user_facing"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(inventory_rows())
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
