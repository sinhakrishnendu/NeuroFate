#!/usr/bin/env python3
"""Create manual sparse extraction plans for external cohorts without executing extraction."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


SUPPORTED_FORMATS = ["h5ad_csr", "csv_genes_as_rows", "csv_cells_as_rows", "mtx_features_barcodes", "bulk_matrix"]


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=path, level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def plan_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    return [
        {
            "dataset_id": args.dataset_id,
            "format": args.format,
            "input_matrix": str(args.input_matrix),
            "metadata_file": str(args.metadata_file),
            "feature_file": str(args.feature_file),
            "panel": str(args.panel),
            "execution_mode": "manual_template_only",
            "safety": "chunked_target_gene_only_no_dense_conversion",
        }
    ]


def write_plan(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset_id", "format", "input_matrix", "metadata_file", "feature_file", "panel", "execution_mode", "safety"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def manual_template(args: argparse.Namespace) -> str:
    output = f"data/interim/external/{args.dataset_id}/{args.dataset_id}_sparse_gene_panel_expression.tsv.gz"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'echo "Do not run from Codex. Manual user execution only."',
        'if [[ "${RUN_MANUAL_EXTRACTION:-NO}" != "YES" ]]; then',
        '  echo "Set RUN_MANUAL_EXTRACTION=YES after reviewing memory limits and official files."',
        "  exit 1",
        "fi",
        "",
        f"# Dataset: {args.dataset_id}",
        f"# Format: {args.format}",
        "# MANUAL_HEAVY extraction templates below are intentionally commented.",
    ]
    if args.format == "h5ad_csr":
        lines.append(f"# python scripts/15_sparse_gene_extraction_safe.py --input {args.input_matrix} --gene-panel {args.panel} --out {output} --chunk-size 5000 --max-genes 64")
    elif args.format == "csv_genes_as_rows":
        lines.append(f"# python scripts/41_extract_mathys_target_gene_panel.py --counts {args.input_matrix} --gene-panel {args.panel} --orientation genes_as_rows --out {output}")
    elif args.format == "csv_cells_as_rows":
        lines.append(f"# python scripts/41_extract_mathys_target_gene_panel.py --counts {args.input_matrix} --gene-panel {args.panel} --orientation cells_as_rows --out {output}")
    elif args.format == "mtx_features_barcodes":
        lines.append(f"# python scripts/FUTURE_extract_mtx_target_gene_panel.py --matrix {args.input_matrix} --features {args.feature_file} --metadata {args.metadata_file} --panel {args.panel} --out {output}")
    elif args.format == "bulk_matrix":
        lines.append(f"# python scripts/FUTURE_extract_bulk_target_gene_panel.py --matrix {args.input_matrix} --metadata {args.metadata_file} --features {args.feature_file} --panel {args.panel} --out {output}")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare external sparse extraction manual plan.")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--format", choices=SUPPORTED_FORMATS, required=True)
    parser.add_argument("--input-matrix", type=Path, required=True)
    parser.add_argument("--metadata-file", type=Path, required=True)
    parser.add_argument("--feature-file", type=Path, required=True)
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--output-plan", type=Path, required=True)
    parser.add_argument("--manual-script-output", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log_file)
    write_plan(args.output_plan, plan_rows(args))
    args.manual_script_output.parent.mkdir(parents=True, exist_ok=True)
    args.manual_script_output.write_text(manual_template(args), encoding="utf-8")
    logging.info("Wrote manual extraction plan for %s", args.dataset_id)
    print(f"Wrote {args.output_plan}")
    print(f"Wrote {args.manual_script_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
