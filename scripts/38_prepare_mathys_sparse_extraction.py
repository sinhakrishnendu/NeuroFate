#!/usr/bin/env python3
"""Prepare a guarded sparse extraction plan for Mathys 2019 target genes.

This planning script reads only gene-overlap TSVs and writes manual extraction
commands. It never opens expression files.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


MAX_GENES = 64
DEFAULT_CHUNK_SIZE = 5000
MAX_CHUNK_SIZE = 50000
DEFAULT_MEMORY_LIMIT_MB = 4096


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


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def plan_rows(overlap_rows: list[dict[str, str]], chunk_size: int, memory_limit_mb: int) -> list[dict[str, str]]:
    present = [row for row in overlap_rows if row.get("mathys_status") == "present"]
    if len(present) > MAX_GENES:
        raise ValueError(f"Mathys target genes exceed hard limit {MAX_GENES}.")
    if chunk_size > MAX_CHUNK_SIZE:
        raise ValueError(f"chunk_size may not exceed {MAX_CHUNK_SIZE}.")
    rows: list[dict[str, str]] = []
    for row in present:
        rows.append(
            {
                "dataset_id": "mathys_2019_ad",
                "gene_symbol": row["gene_symbol"],
                "mathys_gene_symbol": row.get("mathys_gene_symbol", ""),
                "mathys_gene_id": row.get("mathys_gene_id", ""),
                "mathys_var_index": row.get("mathys_var_index", ""),
                "extraction_status": "planned",
                "chunk_size": str(chunk_size),
                "memory_limit_mb": str(memory_limit_mb),
                "safety_notes": "manual guarded sparse extraction only; no dense conversion",
            }
        )
    return rows


def write_manual_template(
    path: Path,
    input_path: Path,
    var_path: Path,
    output_path: Path,
    chunk_size: int,
    memory_limit_mb: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"""#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: Mathys 2019 sparse target-gene extraction."
echo "DO NOT RUN FROM CODEX. Review paths and memory settings before executing."

: "${{RUN_MANUAL_EXTERNAL_EXTRACTION:=NO}}"
if [[ "${{RUN_MANUAL_EXTERNAL_EXTRACTION}}" != "YES" ]]; then
  echo "RUN_MANUAL_EXTERNAL_EXTRACTION is not YES. Exiting without extraction."
  exit 0
fi

python scripts/31_sparse_external_gene_extraction.py \\
  --dataset-id mathys_2019_ad \\
  --input {input_path} \\
  --var {var_path} \\
  --panel metadata/target_gene_panel_v1.tsv \\
  --output {output_path} \\
  --log-file results/logs/31_sparse_external_gene_extraction_mathys.log \\
  --chunk-size {chunk_size} \\
  --memory-limit-mb {memory_limit_mb} \\
  --execute
"""
    path.write_text(text, encoding="utf-8")
    logging.info("Wrote manual extraction template: %s", path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Mathys sparse extraction plan.")
    parser.add_argument("--overlap", type=Path, default=Path("results/tables/mathys_gene_overlap.tsv"))
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/raw/external/mathys_2019/mathys2019_external.h5ad"),
    )
    parser.add_argument(
        "--var",
        type=Path,
        default=Path("data/interim/external/mathys_2019/mathys_var_genes.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/interim/external/mathys_2019/sparse_gene_panel_expression.tsv.gz"),
    )
    parser.add_argument(
        "--plan-output",
        type=Path,
        default=Path("results/tables/mathys_sparse_extraction_plan.tsv"),
    )
    parser.add_argument(
        "--manual-script-output",
        type=Path,
        default=Path("results/logs/manual_mathys_sparse_extraction_template.sh"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/38_prepare_mathys_sparse_extraction.log"),
    )
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--memory-limit-mb", type=int, default=DEFAULT_MEMORY_LIMIT_MB)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    overlap_rows = read_tsv(args.overlap)
    rows = plan_rows(overlap_rows, args.chunk_size, args.memory_limit_mb)
    write_tsv(
        args.plan_output,
        rows,
        [
            "dataset_id",
            "gene_symbol",
            "mathys_gene_symbol",
            "mathys_gene_id",
            "mathys_var_index",
            "extraction_status",
            "chunk_size",
            "memory_limit_mb",
            "safety_notes",
        ],
    )
    write_manual_template(
        args.manual_script_output,
        args.input,
        args.var,
        args.output,
        args.chunk_size,
        args.memory_limit_mb,
    )
    logging.info("Sparse extraction plan complete. No expression file was opened.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
