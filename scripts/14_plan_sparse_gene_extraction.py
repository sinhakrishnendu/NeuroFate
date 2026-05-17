#!/usr/bin/env python3
"""Plan safe sparse extraction for a small SEA-AD gene panel.

This script reads only lightweight TSV files: the target gene panel and the
metadata-only var gene table. It does not open the H5AD file and never reads
expression values.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


PRESENCE_OUTPUT = "target_gene_panel_presence.tsv"
MANUAL_TEMPLATE_OUTPUT = "manual_sparse_gene_extraction_template.sh"


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
        ],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_var_index(var_path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(var_path)
    index: dict[str, dict[str, str]] = {}
    for row_number, row in enumerate(rows):
        gene_symbol = row.get("gene_symbol", "")
        gene_id = row.get("gene_id", "")
        for key in {gene_symbol, gene_id}:
            if key and key not in index:
                index[key] = {
                    "gene_symbol": gene_symbol,
                    "gene_id": gene_id,
                    "var_index": str(row_number),
                }
    return index


def write_presence_table(
    panel_rows: list[dict[str, str]],
    var_index: dict[str, dict[str, str]],
    output_path: Path,
) -> tuple[list[str], list[str]]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    found: list[str] = []
    missing: list[str] = []
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "gene_symbol",
            "present_in_var",
            "matched_gene_symbol",
            "matched_gene_id",
            "var_index",
            "priority_tier",
            "biological_role",
            "manuscript_relevance",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in panel_rows:
            gene = row["gene_symbol"]
            match = var_index.get(gene)
            present = match is not None
            if present:
                found.append(gene)
            else:
                missing.append(gene)
            writer.writerow(
                {
                    "gene_symbol": gene,
                    "present_in_var": str(present).lower(),
                    "matched_gene_symbol": match["gene_symbol"] if match else "",
                    "matched_gene_id": match["gene_id"] if match else "",
                    "var_index": match["var_index"] if match else "",
                    "priority_tier": row.get("priority_tier", ""),
                    "biological_role": row.get("biological_role", ""),
                    "manuscript_relevance": row.get("manuscript_relevance", ""),
                }
            )
    return found, missing


def write_manual_template(template_path: Path, h5ad_path: Path, panel_path: Path) -> None:
    template_path.parent.mkdir(parents=True, exist_ok=True)
    content = f"""#!/usr/bin/env bash
set -euo pipefail

# DO NOT RUN FROM CODEX.
# MANUAL_SPARSE_EXPRESSION_EXTRACTION: targeted gene panel only.
# This template runs a dry run first. It does not request dense output.

python scripts/15_sparse_gene_extraction_safe.py \\
  --input {h5ad_path} \\
  --var data/interim/sea_ad/sea_ad_var_genes.tsv \\
  --panel {panel_path} \\
  --output data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \\
  --chunk-size 5000 \\
  --max-genes 64 \\
  --dry-run

# After reviewing the dry-run output, run manually only if intended:
# RUN_MANUAL_EXTRACTION=YES python scripts/15_sparse_gene_extraction_safe.py \\
#   --input {h5ad_path} \\
#   --var data/interim/sea_ad/sea_ad_var_genes.tsv \\
#   --panel {panel_path} \\
#   --output data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz \\
#   --max-genes 64 \\
#   --m5-max-profile \\
#   --execute
"""
    template_path.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan sparse SEA-AD target gene extraction.")
    parser.add_argument("--var", type=Path, default=Path("data/interim/sea_ad/sea_ad_var_genes.tsv"))
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/14_sparse_gene_extraction_plan.log"),
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=Path("data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting sparse gene extraction planning.")
    logging.info("Var genes TSV: %s", args.var)
    logging.info("Target panel: %s", args.panel)

    panel_rows = read_tsv(args.panel)
    var_index = load_var_index(args.var)
    presence_path = args.tables_dir / PRESENCE_OUTPUT
    found, missing = write_presence_table(panel_rows, var_index, presence_path)
    template_path = args.log_file.parent / MANUAL_TEMPLATE_OUTPUT
    write_manual_template(template_path, args.h5ad, args.panel)

    logging.info("Panel genes: %d", len(panel_rows))
    logging.info("Genes found in var: %d", len(found))
    logging.info("Genes missing from var: %d", len(missing))
    logging.info("Found genes: %s", ", ".join(found) if found else "none")
    logging.info("Missing genes: %s", ", ".join(missing) if missing else "none")
    logging.info("Wrote %s", presence_path)
    logging.info("Wrote %s", template_path)
    logging.info("No H5AD file was opened and no expression values were read.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
