#!/usr/bin/env python3
"""Plan AD snRNA-seq replication extraction routes without extracting or processing data."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


PLAN_COLUMNS = ["cohort_id", "file_name", "file_size_bytes", "detected_format", "suggested_route", "heavy_execution_status"]


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
    if name.endswith((".h5ad", ".h5")):
        return "h5_container"
    if "matrix.mtx" in name or name.endswith(".mtx.gz"):
        return "mtx_matrix"
    if "features.tsv" in name or "genes.tsv" in name or "barcodes.tsv" in name:
        return "tenx_component"
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        return "archive_manual_listing_needed"
    if name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "csv_tsv_candidate"
    if name.endswith((".xlsx", ".xls")):
        return "metadata_workbook"
    return "manual_review"


def route_for(fmt: str) -> str:
    routes = {
        "h5_container": "metadata_inspection_then_manual_sparse_axis_extraction_no_full_matrix",
        "mtx_matrix": "manual_processed_10x_axis_gene_extraction_plan",
        "tenx_component": "collect_matching_matrix_features_barcodes_before_manual_axis_extraction",
        "archive_manual_listing_needed": "list_archive_only_then_select_processed_matrix_files",
        "csv_tsv_candidate": "review_header_then_use_145_for_bulk_or_sample_matrix",
        "metadata_workbook": "metadata_schema_review_only",
    }
    return routes.get(fmt, "manual_review_required")


def write_manual_template(path: Path, cohort_id: str, input_dir: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: Phase 28 {cohort_id} AD snRNA extraction planning template."
echo "Do not run from Codex. Manual user execution only."
RUN_MANUAL_EXTRACTION="${{RUN_MANUAL_EXTRACTION:-NO}}"
if [[ "${{RUN_MANUAL_EXTRACTION}}" != "YES" ]]; then
  echo "Set RUN_MANUAL_EXTRACTION=YES only after reviewing the Phase 28 plan."
  exit 1
fi

INPUT_DIR="{input_dir}"
echo "Review files under ${{INPUT_DIR}}. Prefer processed 10x matrices or sample-level matrices."
# MANUAL_HEAVY templates only:
# python scripts/145_build_ad_replication_axis_scores_from_matrix.py --matrix <sample_matrix.tsv.gz> --metadata <metadata.tsv> --cohort-id {cohort_id} --sample-id-column <sample_id> --label-column <diagnosis> --positive-class <AD> --negative-class <Control>
# python scripts/73_prepare_external_sparse_extraction_plan.py --dataset-id {cohort_id} --format <format> --input-matrix <matrix> --metadata-file <metadata> --feature-file <features> --panel metadata/target_gene_panel_v1.tsv --output-plan results/tables/phase28_{cohort_id}_sparse_plan.tsv --manual-script-output results/logs/manual_phase28_{cohort_id}_sparse_template.sh
""",
        encoding="utf-8",
    )


def build_plan(cohort_id: str, input_dir: Path, axis_registry: Path) -> list[dict[str, str]]:
    read_tsv(axis_registry)
    files = sorted(path for path in input_dir.rglob("*") if path.is_file()) if input_dir.exists() else []
    if not files:
        return [
            {
                "cohort_id": cohort_id,
                "file_name": "",
                "file_size_bytes": "",
                "detected_format": "no_local_files",
                "suggested_route": "manual_acquisition_required",
                "heavy_execution_status": "manual_template_only",
            }
        ]
    rows: list[dict[str, str]] = []
    for path in files:
        fmt = detect_format(path)
        rows.append(
            {
                "cohort_id": cohort_id,
                "file_name": str(path),
                "file_size_bytes": str(path.stat().st_size),
                "detected_format": fmt,
                "suggested_route": route_for(fmt),
                "heavy_execution_status": "manual_template_only",
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan AD snRNA replication extraction.")
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output-plan", type=Path)
    parser.add_argument("--manual-script-output", type=Path)
    parser.add_argument("--log-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/146_{args.cohort_id}_snrna_plan.log")
    configure_logging(log_file)
    output = args.output_plan or Path(f"results/tables/phase28_{args.cohort_id}_snrna_extraction_plan.tsv")
    manual = args.manual_script_output or Path(f"results/logs/manual_phase28_{args.cohort_id}_snrna_extraction.sh")
    rows = build_plan(args.cohort_id, args.input_dir, args.axis_registry)
    write_tsv(output, rows, PLAN_COLUMNS)
    write_manual_template(manual, args.cohort_id, args.input_dir)
    logging.info("Wrote Phase 28 AD snRNA plan cohort=%s rows=%d", args.cohort_id, len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
