#!/usr/bin/env python3
"""Triage local files for Phase 33 PD replication cohorts without opening large matrices."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


FORMAT_SUFFIXES = {
    "series_matrix": ("series_matrix.txt", "series_matrix.txt.gz"),
    "soft": (".soft", ".soft.gz"),
    "supplementary_archive": (".tar", ".tar.gz", ".tgz", ".zip"),
    "expression_table": (".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz"),
    "platform_annotation": ("GPL", ".annot", ".annot.gz"),
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


def classify_file(path: Path) -> str:
    name = path.name.lower()
    if "series_matrix" in name:
        return "series_matrix"
    if name.endswith((".soft", ".soft.gz")):
        return "soft"
    if name.endswith((".tar", ".tar.gz", ".tgz", ".zip")):
        return "supplementary_archive"
    if name.startswith("gpl") or "platform" in name or "annot" in name:
        return "platform_annotation"
    if name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "expression_or_metadata_table"
    return "unknown"


def likely_next_action(files: list[dict[str, str]]) -> str:
    roles = {row["likely_role"] for row in files}
    if "series_matrix" in roles:
        return "Run scripts/163_build_axis_scores_from_geo_series_matrix.py on the series matrix; add platform annotation only if probe IDs cannot be mapped."
    if "expression_or_metadata_table" in roles:
        return "Inspect table headers manually, then use a sample-level matrix builder or platform mapping if needed."
    if "supplementary_archive" in roles:
        return "List archive contents safely before extracting; prefer processed sample-level expression or matrix files."
    return "Manual acquisition required: download series matrix and processed supplementary files only."


def inspect_cohort(row: dict[str, str]) -> tuple[list[dict[str, str]], dict[str, str]]:
    cohort_id = row["cohort_id"]
    raw_dir = Path(row["local_raw_dir"])
    files: list[dict[str, str]] = []
    if raw_dir.exists():
        for path in sorted(raw_dir.rglob("*")):
            if path.is_file():
                files.append(
                    {
                        "cohort_id": cohort_id,
                        "geo_accession": row.get("geo_accession", ""),
                        "local_raw_dir": str(raw_dir),
                        "file_path": str(path),
                        "file_name": path.name,
                        "size_bytes": str(path.stat().st_size),
                        "likely_role": classify_file(path),
                    }
                )
    status = "local_files_present" if files else "missing_local_files"
    summary = {
        "cohort_id": cohort_id,
        "geo_accession": row.get("geo_accession", ""),
        "priority": row.get("priority", ""),
        "status": status,
        "file_count": str(len(files)),
        "next_action": likely_next_action(files),
    }
    if not files:
        files.append(
            {
                "cohort_id": cohort_id,
                "geo_accession": row.get("geo_accession", ""),
                "local_raw_dir": str(raw_dir),
                "file_path": "",
                "file_name": "",
                "size_bytes": "0",
                "likely_role": "missing",
            }
        )
    return files, summary


def write_next_actions(path: Path, summaries: list[dict[str, str]]) -> None:
    lines = ["# Phase 33 PD Replication Next Actions", ""]
    for row in summaries:
        lines.append(f"## {row['cohort_id']}")
        lines.append(f"- GEO accession: {row['geo_accession']}")
        lines.append(f"- Status: {row['status']}")
        lines.append(f"- Local file count: {row['file_count']}")
        lines.append(f"- Next action: {row['next_action']}")
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Triage local files for Phase 33 PD replication candidates.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase33_pd_replication_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/162_triage_pd_replication_geo_files.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    inventory: list[dict[str, str]] = []
    summaries: list[dict[str, str]] = []
    for row in read_tsv(args.registry):
        files, summary = inspect_cohort(row)
        inventory.extend(files)
        summaries.append(summary)
    write_tsv(
        args.outdir / "phase33_pd_replication_file_inventory.tsv",
        inventory,
        ["cohort_id", "geo_accession", "local_raw_dir", "file_path", "file_name", "size_bytes", "likely_role"],
    )
    write_next_actions(args.outdir / "phase33_pd_replication_next_actions.md", summaries)
    logging.info("Triaged %d PD replication cohorts", len(summaries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
