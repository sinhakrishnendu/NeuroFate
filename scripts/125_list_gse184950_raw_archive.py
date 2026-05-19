#!/usr/bin/env python3
"""List GSE184950_RAW.tar contents safely without extracting files."""

from __future__ import annotations

import argparse
import csv
import logging
import tarfile
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def extension_for(name: str) -> str:
    lower = name.lower()
    for suffix in [".matrix.mtx.gz", ".barcodes.tsv.gz", ".features.tsv.gz", ".genes.tsv.gz", ".fastq.gz", ".tar.gz", ".csv.gz", ".tsv.gz"]:
        if lower.endswith(suffix):
            return suffix
    return Path(lower).suffix


def likely_role(name: str) -> str:
    lower = name.lower()
    if lower.endswith(".tar.gz"):
        return "per_sample_processed_archive"
    if "matrix.mtx" in lower:
        return "tenx_matrix"
    if "barcodes.tsv" in lower:
        return "tenx_barcodes"
    if "features.tsv" in lower or "genes.tsv" in lower:
        return "tenx_features"
    if "fastq" in lower:
        return "raw_fastq_do_not_process_here"
    if lower.endswith((".xlsx", ".xls")):
        return "metadata_workbook"
    return "unknown"


def sample_prefix(name: str) -> str:
    base = Path(name).name
    for token in [".tar.gz", "_matrix", "_barcodes", "_features", "_genes", "_R1", "_R2"]:
        if token in base:
            return base.split(token)[0]
    return base.split(".")[0]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["member_path", "size_bytes", "extension", "likely_role", "sample_prefix"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["likely_role"]] = counts.get(row["likely_role"], 0) + 1
    lines = ["# Phase 24 GSE184950 RAW Archive Summary", "", "Archive was listed only; no members were extracted.", ""]
    for role, count in sorted(counts.items()):
        lines.append(f"- {role}: {count}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="List GSE184950 RAW archive contents without extraction.")
    parser.add_argument("--tar", type=Path, default=Path("data/raw/external/gse184950_pd_sn/GSE184950_RAW.tar"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase24_gse184950_raw_archive_inventory.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/reports/phase24_gse184950_raw_archive_summary.md"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/125_list_gse184950_raw_archive.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows: list[dict[str, str]] = []
    with tarfile.open(args.tar, "r:*") as archive:
        for member in archive.getmembers():
            rows.append(
                {
                    "member_path": member.name,
                    "size_bytes": str(member.size),
                    "extension": extension_for(member.name),
                    "likely_role": likely_role(member.name),
                    "sample_prefix": sample_prefix(member.name),
                }
            )
    write_tsv(args.output, rows)
    write_summary(args.summary_output, rows)
    logging.info("Listed GSE184950 archive members=%d without extraction", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
