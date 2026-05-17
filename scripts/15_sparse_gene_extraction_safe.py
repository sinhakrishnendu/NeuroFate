#!/usr/bin/env python3
"""Safely extract a small targeted gene panel from SEA-AD CSR expression.

This script is for future manual use only. It reads selected CSR rows in chunks
and writes nonzero values for a small gene panel. It never creates dense matrices
and does not run single-cell analysis workflows.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
from pathlib import Path
from typing import Any

import h5py


MAX_GENES_DEFAULT = 64
CHUNK_SIZE_DEFAULT = 5000
MAX_CHUNK_SIZE = 50000
DEFAULT_MEMORY_LIMIT_MB = 4096
M5_MAX_HIGH_MEMORY_LIMIT_MB = 32768
OUTPUT_COLUMNS = ["obs_index", "row_index", "gene_symbol", "gene_id", "var_index", "expression_value"]


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


def decode_scalar(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item"):
        try:
            value = value.item()
        except ValueError:
            pass
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return "" if value is None else str(value)


def read_dataset_slice(dataset: h5py.Dataset, start: int, stop: int) -> list[str]:
    return [decode_scalar(value) for value in dataset[start:stop]]


def obs_index_name(obs_group: h5py.Group) -> str:
    raw = obs_group.attrs.get("_index", "_index")
    return decode_scalar(raw)


def load_selected_genes(var_path: Path, panel_path: Path, max_genes: int) -> list[dict[str, str]]:
    if max_genes > MAX_GENES_DEFAULT:
        raise ValueError(f"max_genes may not exceed {MAX_GENES_DEFAULT}")
    panel_genes = [row["gene_symbol"] for row in read_tsv(panel_path)]
    if len(panel_genes) > max_genes:
        raise ValueError(f"Panel has {len(panel_genes)} genes but max_genes={max_genes}")

    selected: list[dict[str, str]] = []
    wanted = set(panel_genes)
    for index, row in enumerate(read_tsv(var_path)):
        gene_symbol = row.get("gene_symbol", "")
        gene_id = row.get("gene_id", "")
        if gene_symbol in wanted or gene_id in wanted:
            selected.append(
                {
                    "gene_symbol": gene_symbol,
                    "gene_id": gene_id,
                    "var_index": str(index),
                }
            )
    return selected


def estimate_chunk_bytes(nnz: int, row_count: int) -> int:
    return nnz * (8 + 4) + (row_count + 1) * 8


def validate_runtime_limits(selected_genes: list[dict[str, str]], chunk_size: int, memory_limit_mb: int) -> None:
    if not selected_genes:
        raise ValueError("No selected genes found in var table.")
    if len(selected_genes) > MAX_GENES_DEFAULT:
        raise ValueError(f"Selected genes exceed hard limit of {MAX_GENES_DEFAULT}.")
    if chunk_size > MAX_CHUNK_SIZE:
        raise ValueError(f"chunk_size may not exceed {MAX_CHUNK_SIZE}.")
    if memory_limit_mb <= 0:
        raise ValueError("memory_limit_mb must be positive.")


def write_sparse_panel(
    h5ad_path: Path,
    selected_genes: list[dict[str, str]],
    output_path: Path,
    chunk_size: int,
    memory_limit_mb: int,
) -> int:
    selected_by_index = {int(row["var_index"]): row for row in selected_genes}
    written_rows = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5ad_path, "r") as handle, gzip.open(output_path, "wt", encoding="utf-8", newline="") as out:
        matrix = handle["X"]
        indptr = matrix["indptr"]
        indices = matrix["indices"]
        data = matrix["data"]
        obs_group = handle["obs"]
        obs_key = obs_index_name(obs_group)
        obs_index_dataset = obs_group[obs_key]
        total_rows = int(indptr.shape[0]) - 1

        writer = csv.DictWriter(out, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()

        for start in range(0, total_rows, chunk_size):
            stop = min(start + chunk_size, total_rows)
            chunk_indptr = indptr[start : stop + 1]
            data_start = int(chunk_indptr[0])
            data_stop = int(chunk_indptr[-1])
            nnz = data_stop - data_start
            estimated_mb = estimate_chunk_bytes(nnz, stop - start) / (1024 * 1024)
            if estimated_mb > memory_limit_mb:
                raise MemoryError(
                    f"Chunk {start}-{stop} requires estimated {estimated_mb:.2f} MB, "
                    f"above limit {memory_limit_mb} MB."
                )

            chunk_indices = indices[data_start:data_stop]
            chunk_values = data[data_start:data_stop]
            obs_ids = read_dataset_slice(obs_index_dataset, start, stop)

            for row_offset in range(stop - start):
                row_start = int(chunk_indptr[row_offset]) - data_start
                row_stop = int(chunk_indptr[row_offset + 1]) - data_start
                row_index = start + row_offset
                obs_id = obs_ids[row_offset]
                for local_position in range(row_start, row_stop):
                    var_index = int(chunk_indices[local_position])
                    gene = selected_by_index.get(var_index)
                    if gene is None:
                        continue
                    writer.writerow(
                        {
                            "obs_index": obs_id,
                            "row_index": row_index,
                            "gene_symbol": gene["gene_symbol"],
                            "gene_id": gene["gene_id"],
                            "var_index": var_index,
                            "expression_value": float(chunk_values[local_position]),
                        }
                    )
                    written_rows += 1

            logging.info("Processed sparse rows %d-%d of %d", start + 1, stop, total_rows)
    return written_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Safely extract a small sparse SEA-AD gene panel.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--var", type=Path, default=Path("data/interim/sea_ad/sea_ad_var_genes.tsv"))
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz"),
    )
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/15_sparse_gene_extraction_safe.log"))
    parser.add_argument("--max-genes", type=int, default=MAX_GENES_DEFAULT)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE_DEFAULT)
    parser.add_argument("--memory-limit-mb", type=int, default=DEFAULT_MEMORY_LIMIT_MB)
    parser.add_argument(
        "--m5-max-profile",
        action="store_true",
        help="Use a manual high-memory profile: chunk-size 50000 and memory-limit-mb 32768.",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--execute", action="store_true", help="Required for writing expression output.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    if args.m5_max_profile:
        args.chunk_size = MAX_CHUNK_SIZE
        args.memory_limit_mb = M5_MAX_HIGH_MEMORY_LIMIT_MB
        logging.info("Using M5 Max high-memory profile.")
    selected_genes = load_selected_genes(args.var, args.panel, args.max_genes)
    validate_runtime_limits(selected_genes, args.chunk_size, args.memory_limit_mb)

    logging.info("Starting safe sparse gene extraction.")
    logging.info("Input H5AD: %s", args.input)
    logging.info("Selected genes found in var: %d", len(selected_genes))
    logging.info("Chunk size: %d", args.chunk_size)
    logging.info("Memory limit MB: %d", args.memory_limit_mb)
    logging.info("Output: %s", args.output)

    if args.dry_run or not args.execute:
        logging.info("Dry run only. No expression values were read and no output was written.")
        logging.info("Selected gene symbols: %s", ", ".join(row["gene_symbol"] for row in selected_genes))
        return 0

    if os.environ.get("RUN_MANUAL_EXTRACTION") != "YES":
        logging.error("Refusing to execute. Set RUN_MANUAL_EXTRACTION=YES for manual extraction.")
        return 2

    written_rows = write_sparse_panel(
        args.input,
        selected_genes,
        args.output,
        args.chunk_size,
        args.memory_limit_mb,
    )
    logging.info("Wrote nonzero sparse expression rows: %d", written_rows)
    logging.info("No dense matrix was created.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
