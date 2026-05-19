#!/usr/bin/env python3
"""Triage local GEO files for Phase 34 PD microarray/bulk replication."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


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


def classify(path: Path) -> str:
    name = path.name.lower()
    if "series_matrix" in name:
        return "series_matrix"
    if name.endswith((".soft", ".soft.gz")):
        return "soft"
    if name.endswith((".xml", ".xml.gz", ".xml.tgz")) or "miniml" in name:
        return "miniml"
    if name.startswith("gpl") or "platform" in name or "annot" in name:
        return "platform_annotation"
    if name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "supplementary_expression_or_metadata_table"
    if name.endswith((".tar", ".tar.gz", ".tgz", ".zip")):
        return "supplementary_archive"
    return "unknown"


def readiness(files: list[dict[str, str]]) -> str:
    roles = {row["likely_role"] for row in files}
    if "series_matrix" in roles and "platform_annotation" in roles:
        return "ready_for_metadata_and_probe_mapping"
    if "series_matrix" in roles:
        return "metadata_ready_platform_mapping_may_be_needed"
    if "supplementary_expression_or_metadata_table" in roles and "platform_annotation" in roles:
        return "processed_table_ready_for_review"
    if files:
        return "local_files_need_manual_review"
    return "manual_acquisition_required"


def next_action(status: str) -> str:
    if status == "ready_for_metadata_and_probe_mapping":
        return "Parse series metadata, build platform probe map, then build axis scores from processed expression."
    if status == "metadata_ready_platform_mapping_may_be_needed":
        return "Parse series metadata and acquire platform annotation if expression rows are probe IDs."
    if status == "processed_table_ready_for_review":
        return "Inspect table header and run Phase 34 axis-score builder with the probe map."
    if status == "local_files_need_manual_review":
        return "Inspect filenames and acquire missing series matrix or platform annotation."
    return "Run the guarded manual download template for the priority cohort."


def main() -> int:
    parser = argparse.ArgumentParser(description="Triage Phase 34 PD GEO files.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase34_pd_replication_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/162_triage_phase34_pd_geo_files.log"))
    args = parser.parse_args()
    configure_logging(args.log_file)

    inventory: list[dict[str, str]] = []
    summaries: list[dict[str, str]] = []
    for row in read_tsv(args.registry):
        raw_dir = Path(row["local_raw_dir"])
        files: list[dict[str, str]] = []
        if raw_dir.exists():
            for path in sorted(raw_dir.rglob("*")):
                if path.is_file():
                    item = {
                        "cohort_id": row["cohort_id"],
                        "geo_accession": row["geo_accession"],
                        "file_path": str(path),
                        "file_name": path.name,
                        "size_bytes": str(path.stat().st_size),
                        "likely_role": classify(path),
                    }
                    files.append(item)
                    inventory.append(item)
        if not files:
            inventory.append(
                {
                    "cohort_id": row["cohort_id"],
                    "geo_accession": row["geo_accession"],
                    "file_path": "",
                    "file_name": "",
                    "size_bytes": "0",
                    "likely_role": "missing",
                }
            )
        status = readiness(files)
        summaries.append({"cohort_id": row["cohort_id"], "geo_accession": row["geo_accession"], "status": status, "file_count": str(len(files)), "next_action": next_action(status)})

    write_tsv(args.outdir / "phase34_pd_geo_file_inventory.tsv", inventory, ["cohort_id", "geo_accession", "file_path", "file_name", "size_bytes", "likely_role"])
    lines = ["# Phase 34 PD GEO Next Actions", ""]
    for row in summaries:
        lines.extend([f"## {row['cohort_id']}", f"- GEO accession: {row['geo_accession']}", f"- Status: {row['status']}", f"- File count: {row['file_count']}", f"- Next action: {row['next_action']}", ""])
    (args.outdir / "phase34_pd_geo_next_actions.md").write_text("\n".join(lines), encoding="utf-8")
    logging.info("Triaged Phase 34 cohorts=%d", len(summaries))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
