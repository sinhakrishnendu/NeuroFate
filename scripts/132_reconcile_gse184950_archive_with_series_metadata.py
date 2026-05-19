#!/usr/bin/env python3
"""Reconcile listed GSE184950 RAW archive members against series-matrix sample files."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


OUTPUT_COLUMNS = [
    "sample_name",
    "expected_archive_member",
    "processed_tar_name",
    "found_in_archive",
    "archive_member_path",
    "contains_matrix_mtx",
    "contains_barcodes",
    "contains_features_or_genes",
    "contains_fastq",
    "processed_matrix_availability",
]


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


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def truth(value: bool) -> str:
    return str(value).lower()


def basename(path: str) -> str:
    return Path(path).name


def reconcile(inventory: list[dict[str, str]], manifest: list[dict[str, str]]) -> list[dict[str, str]]:
    member_paths = [row.get("member_path", "") for row in inventory]
    member_basenames = {basename(path): path for path in member_paths}
    lower_paths = [path.lower() for path in member_paths]
    rows: list[dict[str, str]] = []
    for item in manifest:
        expected = item.get("expected_archive_member") or item.get("processed_tar_name", "")
        processed_tar = item.get("processed_tar_name", expected)
        found_path = member_basenames.get(expected, "")
        if not found_path:
            found_path = next((path for path in member_paths if expected and expected in path), "")
        sample_prefix = processed_tar.replace(".tar.gz", "").replace(".tgz", "")
        related = [path for path in lower_paths if sample_prefix.lower() and sample_prefix.lower() in path]
        has_matrix = any("matrix.mtx" in path for path in related)
        has_barcodes = any("barcodes.tsv" in path for path in related)
        has_features = any("features.tsv" in path or "genes.tsv" in path for path in related)
        has_fastq = any("fastq" in path for path in related)
        availability = "processed_archive_present" if found_path else "missing_expected_archive"
        if has_matrix and has_barcodes and has_features:
            availability = "processed_10x_members_visible"
        rows.append(
            {
                "sample_name": item.get("sample_name", sample_prefix),
                "expected_archive_member": expected,
                "processed_tar_name": processed_tar,
                "found_in_archive": truth(bool(found_path)),
                "archive_member_path": found_path,
                "contains_matrix_mtx": truth(has_matrix),
                "contains_barcodes": truth(has_barcodes),
                "contains_features_or_genes": truth(has_features),
                "contains_fastq": truth(has_fastq),
                "processed_matrix_availability": availability,
            }
        )
    return rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    found = sum(row["found_in_archive"] == "true" for row in rows)
    visible = sum(row["processed_matrix_availability"] == "processed_10x_members_visible" for row in rows)
    fastq = sum(row["contains_fastq"] == "true" for row in rows)
    lines = [
        "# Phase 25 GSE184950 Archive-Series Reconciliation",
        "",
        "Archive inventory was compared against series-matrix supplementary tar names only. No archive members were extracted.",
        "",
        f"- Expected per-sample processed archives: {len(rows)}",
        f"- Found expected archives: {found}",
        f"- Samples with visible processed 10x members: {visible}",
        f"- Samples with FASTQ-like members detected: {fastq}",
        "",
        "FASTQ/SRA processing remains out of scope for this route.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reconcile GSE184950 RAW archive listing with series metadata.")
    parser.add_argument("--archive-inventory", type=Path, default=Path("results/tables/phase24_gse184950_raw_archive_inventory.tsv"))
    parser.add_argument("--series-manifest", type=Path, default=Path("results/tables/phase25_gse184950_series_processed_file_manifest.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase25_gse184950_archive_series_reconciliation.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/reports/phase25_gse184950_archive_series_reconciliation.md"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/132_reconcile_gse184950_archive_with_series_metadata.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = reconcile(read_tsv(args.archive_inventory), read_tsv(args.series_manifest))
    write_tsv(args.output, rows, OUTPUT_COLUMNS)
    write_summary(args.summary_output, rows)
    logging.info("Reconciled GSE184950 expected archives=%d without extraction", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
