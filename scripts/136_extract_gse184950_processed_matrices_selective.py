#!/usr/bin/env python3
"""Guarded extraction of selected processed 10x matrix files from GSE184950 nested archives."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import tarfile
from pathlib import Path


ALLOWED_NAMES = {"matrix.mtx.gz", "features.tsv.gz", "genes.tsv.gz", "barcodes.tsv.gz"}
AUDIT_COLUMNS = ["sample_id", "outer_member_path", "nested_member_path", "output_path", "status"]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def safe_basename(path: str) -> str:
    return Path(path).name


def selected_inventory_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected = []
    for row in rows:
        name = safe_basename(row.get("nested_member_path", ""))
        if name in ALLOWED_NAMES:
            selected.append(row)
    return selected


def guarded() -> bool:
    return os.environ.get("RUN_MANUAL_GSE184950_EXTRACTION") == "YES"


def extract_selected(raw_tar: Path, inventory: list[dict[str, str]], output_root: Path, execute: bool) -> list[dict[str, str]]:
    selected = selected_inventory_rows(inventory)
    selected_by_outer: dict[str, list[dict[str, str]]] = {}
    for row in selected:
        selected_by_outer.setdefault(row.get("outer_member_path", ""), []).append(row)
    audit: list[dict[str, str]] = []
    if not execute:
        for row in selected:
            sample = row.get("sample_id", "unknown")
            audit.append(
                {
                    "sample_id": sample,
                    "outer_member_path": row.get("outer_member_path", ""),
                    "nested_member_path": row.get("nested_member_path", ""),
                    "output_path": str(output_root / sample / safe_basename(row.get("nested_member_path", ""))),
                    "status": "planned_not_executed",
                }
            )
        return audit
    if not guarded():
        raise SystemExit("Set RUN_MANUAL_GSE184950_EXTRACTION=YES to extract reviewed processed matrix files.")
    output_root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw_tar, "r:*") as outer:
        outer_members = {member.name: member for member in outer.getmembers()}
        for outer_path, wanted_rows in selected_by_outer.items():
            outer_member = outer_members.get(outer_path)
            if outer_member is None:
                continue
            nested_handle = outer.extractfile(outer_member)
            if nested_handle is None:
                continue
            with tarfile.open(fileobj=nested_handle, mode="r:*") as nested:
                nested_members = {member.name: member for member in nested.getmembers()}
                for row in wanted_rows:
                    nested_path = row.get("nested_member_path", "")
                    member = nested_members.get(nested_path)
                    if member is None:
                        continue
                    name = safe_basename(nested_path)
                    if name not in ALLOWED_NAMES:
                        continue
                    sample = row.get("sample_id", "unknown")
                    destination = output_root / sample / name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    source = nested.extractfile(member)
                    if source is None:
                        continue
                    with source, destination.open("wb") as handle:
                        while True:
                            chunk = source.read(1024 * 1024)
                            if not chunk:
                                break
                            handle.write(chunk)
                    audit.append(
                        {
                            "sample_id": sample,
                            "outer_member_path": outer_path,
                            "nested_member_path": nested_path,
                            "output_path": str(destination),
                            "status": "extracted_processed_matrix_file_only",
                        }
                    )
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Selectively extract GSE184950 processed matrix files.")
    parser.add_argument("--raw-tar", type=Path, default=Path("data/raw/external/gse184950_pd_sn/GSE184950_RAW.tar"))
    parser.add_argument("--nested-inventory", type=Path, default=Path("results/tables/phase26_gse184950_nested_archive_inventory.tsv"))
    parser.add_argument("--output-root", type=Path, default=Path("data/interim/external/gse184950_pd_sn/processed_matrices"))
    parser.add_argument("--audit-output", type=Path, default=Path("results/tables/phase26_gse184950_selected_extraction_audit.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/136_extract_gse184950_processed_matrices_selective.log"))
    parser.add_argument("--execute", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    audit = extract_selected(args.raw_tar, read_tsv(args.nested_inventory), args.output_root, args.execute)
    write_tsv(args.audit_output, audit)
    logging.info("GSE184950 selective extraction audit rows=%d execute=%s", len(audit), args.execute)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
