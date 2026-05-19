#!/usr/bin/env python3
"""Triage local files for Phase 28 AD replication cohorts without opening large matrices."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


INVENTORY_COLUMNS = ["cohort_id", "path", "size_bytes", "detected_role", "recommended_next_action"]


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


def role_for(path: Path) -> str:
    name = path.name.lower()
    if "series_matrix" in name and name.endswith((".txt.gz", ".txt")):
        return "geo_series_matrix"
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        return "raw_or_processed_archive_manual_inspection_needed"
    if name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "sample_or_count_matrix_candidate"
    if name.endswith((".h5", ".h5ad")):
        return "h5_container_metadata_or_manual_sparse_plan_only"
    if name.endswith((".xlsx", ".xls")):
        return "supplementary_metadata_workbook"
    if "matrix.mtx" in name or name.endswith(".mtx.gz"):
        return "processed_mtx_matrix_candidate"
    if "features.tsv" in name or "genes.tsv" in name or "barcodes.tsv" in name:
        return "processed_10x_component_candidate"
    return "manual_review_required"


def action_for(role: str) -> str:
    if role == "geo_series_matrix":
        return "parse_with_143_parse_geo_series_matrix_generic"
    if role == "sample_or_count_matrix_candidate":
        return "review_header_then_run_145_if_sample_level_or_bulk_matrix"
    if role == "raw_or_processed_archive_manual_inspection_needed":
        return "list_archive_only_then_plan_processed_matrix_route"
    if role == "h5_container_metadata_or_manual_sparse_plan_only":
        return "metadata_inspection_only_no_full_matrix"
    if role == "supplementary_metadata_workbook":
        return "inspect_metadata_schema_only"
    if "10x" in role or "mtx" in role:
        return "plan_processed_10x_axis_extraction_with_146"
    return "manual_review_required"


def triage(registry: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for cohort in registry:
        cohort_id = cohort["cohort_id"]
        raw_dir = Path(cohort["local_raw_dir"])
        files = sorted(path for path in raw_dir.rglob("*") if path.is_file()) if raw_dir.exists() else []
        if not files:
            rows.append(
                {
                    "cohort_id": cohort_id,
                    "path": str(raw_dir),
                    "size_bytes": "",
                    "detected_role": "local_raw_dir_missing_or_empty",
                    "recommended_next_action": "run_guarded_manual_acquisition_template",
                }
            )
            continue
        for path in files:
            role = role_for(path)
            rows.append(
                {
                    "cohort_id": cohort_id,
                    "path": str(path),
                    "size_bytes": str(path.stat().st_size),
                    "detected_role": role,
                    "recommended_next_action": action_for(role),
                }
            )
    return rows


def write_actions(path: Path, rows: list[dict[str, str]]) -> None:
    lines = ["# Phase 28 AD Replication Next Actions", "", "No files were downloaded or opened as expression matrices.", ""]
    for row in rows:
        lines.append(f"- `{row['cohort_id']}`: {row['recommended_next_action']} ({row['detected_role']})")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Triage local AD replication cohort files.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase28_ad_replication_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/144_triage_ad_replication_files.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = triage(read_tsv(args.registry))
    write_tsv(args.outdir / "phase28_ad_replication_file_inventory.tsv", rows, INVENTORY_COLUMNS)
    write_actions(args.outdir / "phase28_ad_replication_next_actions.md", rows)
    logging.info("Triaged Phase 28 AD replication files rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
