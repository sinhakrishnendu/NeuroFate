#!/usr/bin/env python3
"""Prepare manual snRNA-seq axis extraction plans for replication cohorts."""

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


def detect_route(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".h5ad", ".h5", ".loom")):
        return "safe_metadata_inspection_then_manual_sparse_axis_extraction_no_full_matrix"
    if name.endswith((".mtx", ".mtx.gz")):
        return "manual_mtx_target_gene_axis_extraction_plan"
    if name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")):
        return "manual_csv_tsv_axis_extraction_or_sample_matrix_axis_scoring"
    if name.endswith((".tar", ".sra", ".fastq", ".fastq.gz")):
        return "manual_preprocessing_required_do_not_run_sra_from_codex"
    if name.endswith((".xlsx", ".xls")):
        return "metadata_or_annotation_review"
    return "manual_review_required"


def write_manual_template(path: Path, cohort_id: str, input_dir: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: Phase 23 {cohort_id} axis extraction template."
echo "Do not run from Codex. Manual user execution only."
RUN_MANUAL_EXTRACTION="${{RUN_MANUAL_EXTRACTION:-NO}}"
if [[ "${{RUN_MANUAL_EXTRACTION}}" != "YES" ]]; then
  echo "Set RUN_MANUAL_EXTRACTION=YES after reviewing the plan."
  exit 1
fi

INPUT_DIR="{input_dir}"
echo "Review files under ${{INPUT_DIR}} and choose a safe extraction route."
# MANUAL_HEAVY template examples only:
# python scripts/119_build_axis_scores_from_sample_matrix.py --matrix <sample_matrix.csv.gz> --metadata <metadata.tsv> --cohort-id {cohort_id} --sample-id-column <sample_id> --label-column <label>
# python scripts/73_prepare_external_sparse_extraction_plan.py --dataset-id {cohort_id} --format <format> --input-matrix <matrix> --metadata-file <metadata> --feature-file <features> --panel metadata/target_gene_panel_v1.tsv --output-plan results/tables/phase23_{cohort_id}_sparse_plan.tsv --manual-script-output results/logs/manual_phase23_{cohort_id}_sparse_template.sh
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a manual snRNA axis extraction plan for a replication cohort.")
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output-plan", type=Path)
    parser.add_argument("--manual-script-output", type=Path)
    parser.add_argument("--log-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/120_{args.cohort_id}_snrna_axis_plan.log")
    configure_logging(log_file)
    read_tsv(args.axis_registry)
    files = sorted(path for path in args.input_dir.glob("*") if path.is_file()) if args.input_dir.exists() else []
    rows = [
        {
            "cohort_id": args.cohort_id,
            "file_name": path.name,
            "file_size_bytes": str(path.stat().st_size),
            "suggested_route": detect_route(path),
            "heavy_execution_status": "manual_template_only",
        }
        for path in files
    ]
    if not rows:
        rows.append(
            {
                "cohort_id": args.cohort_id,
                "file_name": "",
                "file_size_bytes": "",
                "suggested_route": "manual_acquisition_or_file_inventory_required",
                "heavy_execution_status": "manual_template_only",
            }
        )
    output_plan = args.output_plan or Path(f"results/tables/phase23_{args.cohort_id}_snrna_axis_plan.tsv")
    manual_script = args.manual_script_output or Path(f"results/logs/manual_phase23_{args.cohort_id}_axis_extraction_template.sh")
    write_tsv(output_plan, rows, ["cohort_id", "file_name", "file_size_bytes", "suggested_route", "heavy_execution_status"])
    write_manual_template(manual_script, args.cohort_id, args.input_dir)
    logging.info("Wrote snRNA replication plan for cohort=%s files=%d", args.cohort_id, len(files))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
