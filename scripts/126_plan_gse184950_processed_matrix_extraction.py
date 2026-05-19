#!/usr/bin/env python3
"""Plan selective GSE184950 processed 10x matrix extraction without extracting expression."""

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
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_preferred_metadata(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if rows:
        if len(rows) <= 2 and "phase24" in str(path):
            logging.warning("Only %d workbook metadata rows are available; prefer the Phase 25 series matrix metadata.", len(rows))
        return rows
    phase25 = Path("results/tables/phase25_gse184950_series_sample_metadata.tsv")
    if path != phase25 and phase25.exists():
        logging.info("Using Phase 25 series metadata: %s", phase25)
        return read_tsv(phase25)
    phase24 = Path("results/tables/phase24_gse184950_sample_metadata.tsv")
    if path != phase24 and phase24.exists():
        rows = read_tsv(phase24)
        if len(rows) <= 2:
            logging.warning("Falling back to incomplete Phase 24 workbook metadata with %d rows.", len(rows))
        return rows
    return []


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sample_from_processed(file_name: str) -> str:
    base = Path(file_name).name
    return base.replace(".tar.gz", "").replace(".tgz", "")


def build_plan(inventory: list[dict[str, str]], metadata: list[dict[str, str]], axis_registry: Path) -> list[dict[str, str]]:
    read_tsv(axis_registry)
    declared = {}
    for row in metadata:
        for column in ("processed_data_file", "processed_tar_name", "expected_archive_member"):
            if row.get(column):
                declared[row[column]] = row
    processed = [row for row in inventory if row.get("likely_role") == "per_sample_processed_archive"]
    matrix_members = [row for row in inventory if row.get("likely_role") == "tenx_matrix"]
    rows: list[dict[str, str]] = []
    if processed:
        for row in processed:
            file_name = Path(row["member_path"]).name
            sample = sample_from_processed(file_name)
            rows.append(
                {
                    "sample_name": sample,
                    "processed_archive": row["member_path"],
                    "metadata_declared": str(file_name in declared or row["member_path"] in declared).lower(),
                    "matrix_status": "archive_member_requires_manual_selective_extraction",
                    "recommended_route": "manual_extract_per_sample_tar_then_run_127_on_processed_matrices",
                    "fastq_processing_status": "not_used",
                }
            )
    elif matrix_members:
        for row in matrix_members:
            rows.append(
                {
                    "sample_name": row.get("sample_prefix", ""),
                    "processed_archive": "",
                    "metadata_declared": "unknown",
                    "matrix_status": "tenx_matrix_member_visible",
                    "recommended_route": "manual_selective_extract_matrix_barcodes_features_only",
                    "fastq_processing_status": "not_used",
                }
            )
    else:
        rows.append(
            {
                "sample_name": "",
                "processed_archive": "",
                "metadata_declared": "false",
                "matrix_status": "processed_10x_matrix_not_detected",
                "recommended_route": "manual_review_required_fastq_processing_not_performed_by_neurofate",
                "fastq_processing_status": "blocked",
            }
        )
    return rows


def write_manual_template(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: GSE184950 selective processed-matrix extraction template."
echo "Do not run from Codex. Manual user execution only."
RUN_MANUAL_EXTRACTION="${RUN_MANUAL_EXTRACTION:-NO}"
if [[ "${RUN_MANUAL_EXTRACTION}" != "YES" ]]; then
  echo "Set RUN_MANUAL_EXTRACTION=YES only after reviewing phase24_gse184950_axis_extraction_plan.tsv."
  exit 1
fi

RAW_ARCHIVE="data/raw/external/gse184950_pd_sn/GSE184950_RAW.tar"
OUT_ROOT="data/interim/external/gse184950_pd_sn/processed_matrices"
mkdir -p "${OUT_ROOT}"

# MANUAL_HEAVY examples only. Prefer processed 10x matrices; do not process FASTQ here.
# tar -tf "${RAW_ARCHIVE}" | grep -E 'matrix.mtx.gz|features.tsv.gz|genes.tsv.gz|barcodes.tsv.gz|\\.tar\\.gz$'
# tar -xf "${RAW_ARCHIVE}" -C "${OUT_ROOT}" <reviewed_sample_archive_or_matrix_members>
# python scripts/127_extract_gse184950_axis_genes_from_10x.py --matrix-dir-root "${OUT_ROOT}"
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan GSE184950 processed matrix extraction.")
    parser.add_argument("--archive-inventory", type=Path, default=Path("results/tables/phase24_gse184950_raw_archive_inventory.tsv"))
    parser.add_argument("--sample-metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output-plan", type=Path, default=Path("results/tables/phase24_gse184950_axis_extraction_plan.tsv"))
    parser.add_argument("--manual-script-output", type=Path, default=Path("results/logs/manual_phase24_gse184950_axis_extraction_template.sh"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/126_plan_gse184950_processed_matrix_extraction.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = build_plan(read_tsv(args.archive_inventory), read_preferred_metadata(args.sample_metadata), args.axis_registry)
    write_tsv(args.output_plan, rows, ["sample_name", "processed_archive", "metadata_declared", "matrix_status", "recommended_route", "fastq_processing_status"])
    write_manual_template(args.manual_script_output)
    logging.info("Wrote GSE184950 matrix extraction plan rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
