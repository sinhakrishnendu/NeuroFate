#!/usr/bin/env python3
"""Manual sparse external-cohort target-gene extraction template.

This script is guarded for future manual use. Dry runs validate target genes and
runtime limits without opening external expression files. Execution requires both
`--execute` and `RUN_MANUAL_EXTERNAL_EXTRACTION=YES`.
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


MAX_GENES = 64
DEFAULT_CHUNK_SIZE = 5000
MAX_CHUNK_SIZE = 50000
DEFAULT_MEMORY_LIMIT_MB = 4096
OUTPUT_COLUMNS = [
    "dataset_id",
    "obs_index",
    "row_index",
    "gene_symbol",
    "gene_id",
    "var_index",
    "expression_value",
]


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
    return "" if value is None else str(value)


def obs_index_name(obs_group: h5py.Group) -> str:
    raw = obs_group.attrs.get("_index", "_index")
    return decode_scalar(raw)


def read_obs_ids(dataset: h5py.Dataset, start: int, stop: int) -> list[str]:
    return [decode_scalar(value) for value in dataset[start:stop]]


def load_selected_genes(var_path: Path, panel_path: Path, max_genes: int) -> list[dict[str, str]]:
    if max_genes > MAX_GENES:
        raise ValueError(f"max_genes may not exceed {MAX_GENES}")
    panel_genes = [row["gene_symbol"] for row in read_tsv(panel_path)]
    if len(panel_genes) > max_genes:
        raise ValueError(f"Target panel has {len(panel_genes)} genes but max_genes={max_genes}")
    wanted = set(panel_genes)
    selected: list[dict[str, str]] = []
    for index, row in enumerate(read_tsv(var_path)):
        gene_symbol = row.get("gene_symbol", "")
        gene_id = row.get("gene_id", "")
        if gene_symbol in wanted or gene_id in wanted:
            selected.append({"gene_symbol": gene_symbol, "gene_id": gene_id, "var_index": str(index)})
    return selected


def validate_limits(selected_genes: list[dict[str, str]], chunk_size: int, memory_limit_mb: int) -> None:
    if not selected_genes:
        raise ValueError("No target genes were found in the external var metadata TSV.")
    if len(selected_genes) > MAX_GENES:
        raise ValueError(f"Selected gene count exceeds {MAX_GENES}.")
    if chunk_size > MAX_CHUNK_SIZE:
        raise ValueError(f"chunk_size may not exceed {MAX_CHUNK_SIZE}.")
    if memory_limit_mb <= 0:
        raise ValueError("memory_limit_mb must be positive.")


def estimate_chunk_mb(nnz: int, row_count: int) -> float:
    return (nnz * (8 + 4) + (row_count + 1) * 8) / (1024 * 1024)


def write_sparse_panel(
    dataset_id: str,
    input_path: Path,
    selected_genes: list[dict[str, str]],
    output_path: Path,
    chunk_size: int,
    memory_limit_mb: int,
) -> int:
    selected_by_index = {int(row["var_index"]): row for row in selected_genes}
    written_rows = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(input_path, "r") as handle, gzip.open(output_path, "wt", encoding="utf-8", newline="") as out:
        matrix = handle["X"]
        indptr = matrix["indptr"]
        indices = matrix["indices"]
        data = matrix["data"]
        obs_group = handle["obs"]
        obs_dataset = obs_group[obs_index_name(obs_group)]
        row_count = int(indptr.shape[0]) - 1
        writer = csv.DictWriter(out, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()

        for start in range(0, row_count, chunk_size):
            stop = min(start + chunk_size, row_count)
            chunk_indptr = indptr[start : stop + 1]
            data_start = int(chunk_indptr[0])
            data_stop = int(chunk_indptr[-1])
            nnz = data_stop - data_start
            chunk_mb = estimate_chunk_mb(nnz, stop - start)
            if chunk_mb > memory_limit_mb:
                raise MemoryError(f"Chunk estimate {chunk_mb:.2f} MB exceeds limit {memory_limit_mb} MB.")
            chunk_indices = indices[data_start:data_stop]
            chunk_values = data[data_start:data_stop]
            obs_ids = read_obs_ids(obs_dataset, start, stop)

            for row_offset in range(stop - start):
                row_start = int(chunk_indptr[row_offset]) - data_start
                row_stop = int(chunk_indptr[row_offset + 1]) - data_start
                row_index = start + row_offset
                for local_position in range(row_start, row_stop):
                    var_index = int(chunk_indices[local_position])
                    gene = selected_by_index.get(var_index)
                    if gene is None:
                        continue
                    writer.writerow(
                        {
                            "dataset_id": dataset_id,
                            "obs_index": obs_ids[row_offset],
                            "row_index": row_index,
                            "gene_symbol": gene["gene_symbol"],
                            "gene_id": gene["gene_id"],
                            "var_index": var_index,
                            "expression_value": float(chunk_values[local_position]),
                        }
                    )
                    written_rows += 1
            logging.info("Processed sparse rows %d-%d of %d", start + 1, stop, row_count)
    return written_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Manual sparse external target-gene extraction.")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--var", type=Path, required=True, help="External var metadata TSV with gene_symbol/gene_id.")
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/31_sparse_external_gene_extraction.log"),
    )
    parser.add_argument("--max-genes", type=int, default=MAX_GENES)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--memory-limit-mb", type=int, default=DEFAULT_MEMORY_LIMIT_MB)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--execute", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    selected_genes = load_selected_genes(args.var, args.panel, args.max_genes)
    validate_limits(selected_genes, args.chunk_size, args.memory_limit_mb)
    logging.info("Dataset ID: %s", args.dataset_id)
    logging.info("Selected genes found in external var metadata: %d", len(selected_genes))
    logging.info("Chunk size: %d", args.chunk_size)
    logging.info("Memory limit MB: %d", args.memory_limit_mb)

    if args.dry_run or not args.execute:
        logging.info("Dry run only. No external expression file was opened.")
        logging.info("Selected genes: %s", ", ".join(row["gene_symbol"] for row in selected_genes))
        return 0
    if os.environ.get("RUN_MANUAL_EXTERNAL_EXTRACTION") != "YES":
        raise RuntimeError("Set RUN_MANUAL_EXTERNAL_EXTRACTION=YES and pass --execute for manual extraction.")

    written = write_sparse_panel(
        args.dataset_id,
        args.input,
        selected_genes,
        args.output,
        args.chunk_size,
        args.memory_limit_mb,
    )
    logging.info("Wrote sparse target-gene rows: %d", written)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
