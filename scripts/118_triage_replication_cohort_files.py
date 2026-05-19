#!/usr/bin/env python3
"""Triage local replication cohort directories without opening large biological files."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FORMAT_SUFFIXES = {
    ".h5ad": "h5ad_container",
    ".h5": "hdf5_container",
    ".mtx": "matrix_market",
    ".csv": "csv_table",
    ".tsv": "tsv_table",
    ".txt": "text_table",
    ".xlsx": "excel_metadata",
    ".xls": "excel_metadata",
    ".tar": "tar_archive",
    ".gz": "gzip_compressed",
    ".rds": "r_object",
    ".loom": "loom_container",
    ".fastq": "fastq",
    ".sra": "sra_archive",
}


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def detect_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".csv.gz"):
        return "csv_gzip_table"
    if name.endswith(".tsv.gz") or name.endswith(".txt.gz"):
        return "text_gzip_table"
    for suffix, label in FORMAT_SUFFIXES.items():
        if name.endswith(suffix):
            return label
    return "unknown"


def recommend_action(format_counts: dict[str, int], dir_exists: bool) -> str:
    if not dir_exists:
        return "manual_download_required"
    if any(key in format_counts for key in ["csv_table", "csv_gzip_table", "tsv_table", "text_gzip_table", "excel_metadata"]):
        return "inspect_metadata_and_matrix_headers"
    if any(key in format_counts for key in ["tar_archive", "sra_archive", "fastq"]):
        return "manual_preprocessing_plan_required_do_not_run_sra_from_codex"
    if any(key in format_counts for key in ["h5ad_container", "hdf5_container", "loom_container"]):
        return "safe_metadata_inspection_only_no_matrix_loading"
    return "manual_file_review_required"


def triage_row(cohort: dict[str, str]) -> tuple[list[dict[str, str]], dict[str, str]]:
    raw_dir = Path(cohort["local_raw_dir"])
    files = sorted(path for path in raw_dir.glob("*") if path.is_file()) if raw_dir.exists() else []
    inventory: list[dict[str, str]] = []
    format_counts: dict[str, int] = {}
    for path in files:
        fmt = detect_format(path)
        format_counts[fmt] = format_counts.get(fmt, 0) + 1
        inventory.append(
            {
                "cohort_id": cohort["cohort_id"],
                "geo_accession": cohort["geo_accession"],
                "local_raw_dir": str(raw_dir),
                "dir_exists": str(raw_dir.exists()).lower(),
                "file_name": path.name,
                "file_size_bytes": str(path.stat().st_size),
                "detected_format": fmt,
                "large_file_flag": str(path.stat().st_size > 500_000_000).lower(),
                "recommended_next_action": "",
            }
        )
    recommendation = recommend_action(format_counts, raw_dir.exists())
    if not inventory:
        inventory.append(
            {
                "cohort_id": cohort["cohort_id"],
                "geo_accession": cohort["geo_accession"],
                "local_raw_dir": str(raw_dir),
                "dir_exists": str(raw_dir.exists()).lower(),
                "file_name": "",
                "file_size_bytes": "",
                "detected_format": "none",
                "large_file_flag": "false",
                "recommended_next_action": recommendation,
            }
        )
    for row in inventory:
        row["recommended_next_action"] = recommendation
    summary = {
        "cohort_id": cohort["cohort_id"],
        "geo_accession": cohort["geo_accession"],
        "file_count": str(len(files)),
        "formats": ";".join(f"{key}:{value}" for key, value in sorted(format_counts.items())) or "none",
        "recommended_next_action": recommendation,
    }
    return inventory, summary


def write_markdown(path: Path, summaries: list[dict[str, str]]) -> None:
    lines = ["# Phase 23 Replication Next Actions", ""]
    for row in summaries:
        lines.append(f"- `{row['cohort_id']}` ({row['geo_accession']}): {row['recommended_next_action']} ({row['file_count']} local files; formats={row['formats']}).")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Triage local files for Phase 23 replication cohorts.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase23_replication_cohort_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/118_triage_replication_cohort_files.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    inventory_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, str]] = []
    for cohort in read_tsv(args.registry):
        inventory, summary = triage_row(cohort)
        inventory_rows.extend(inventory)
        summary_rows.append(summary)
    write_tsv(
        args.outdir / "phase23_replication_file_inventory.tsv",
        inventory_rows,
        ["cohort_id", "geo_accession", "local_raw_dir", "dir_exists", "file_name", "file_size_bytes", "detected_format", "large_file_flag", "recommended_next_action"],
    )
    write_markdown(args.outdir / "phase23_replication_next_actions.md", summary_rows)
    logging.info("Triaged replication cohorts=%d", len(summary_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
