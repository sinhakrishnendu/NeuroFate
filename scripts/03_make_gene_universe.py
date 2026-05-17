#!/usr/bin/env python3
"""LIGHTWEIGHT gene universe builder from small text files.

This script accepts plain TSV/CSV/text inputs only. It refuses large files by default
and never reads HDF5, h5ad, zarr, loom, MTX, or parquet single-cell data.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


BLOCKED_SUFFIXES = {".h5", ".hdf5", ".h5ad", ".zarr", ".loom", ".mtx", ".parquet"}


def configure_logging(log_file: Path | None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, mode="w", encoding="utf-8"))
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LIGHTWEIGHT gene universe builder.")
    parser.add_argument(
        "--input",
        type=Path,
        action="append",
        required=True,
        help="Small TSV/CSV/text file. May be provided multiple times.",
    )
    parser.add_argument("--output", type=Path, default=Path("data/interim/gene_universe.tsv"))
    parser.add_argument("--gene-column", default="gene_symbol")
    parser.add_argument("--max-lightweight-file-mb", type=float, default=50.0)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write output. Without this flag the script is a dry run.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview only. This is the default when --write is not supplied.",
    )
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/03_make_gene_universe.log"))
    return parser.parse_args()


def delimiter_for(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def read_gene_symbols(path: Path, gene_column: str, max_mb: float) -> set[str]:
    suffix = path.suffix.lower()
    if suffix in BLOCKED_SUFFIXES:
        raise ValueError(f"Refusing heavy or binary input: {path}")
    if not path.exists():
        raise FileNotFoundError(path)
    size_mb = path.stat().st_size / (1024 * 1024)
    if size_mb > max_mb:
        raise ValueError(f"Refusing file above lightweight threshold ({size_mb:.2f} MB): {path}")

    genes: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.readline()
        handle.seek(0)
        if gene_column in sample:
            reader = csv.DictReader(handle, delimiter=delimiter_for(path))
            for row in reader:
                value = (row.get(gene_column) or "").strip()
                if value:
                    genes.add(value)
        else:
            for line in handle:
                value = line.strip().split(delimiter_for(path))[0]
                if value and not value.startswith("#"):
                    genes.add(value)
    return genes


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    effective_dry_run = args.dry_run or not args.write

    logging.info("Starting LIGHTWEIGHT gene universe preparation.")
    logging.info("Dry run: %s", effective_dry_run)
    logging.info("Output: %s", args.output)

    all_genes: set[str] = set()
    for input_path in args.input:
        logging.info("Inspecting lightweight text input: %s", input_path)
        try:
            genes = read_gene_symbols(input_path, args.gene_column, args.max_lightweight_file_mb)
        except (FileNotFoundError, ValueError) as exc:
            logging.warning("%s", exc)
            continue
        logging.info("Parsed %d unique candidate genes from %s", len(genes), input_path)
        all_genes.update(genes)

    logging.info("Total unique candidate genes: %d", len(all_genes))
    if effective_dry_run:
        logging.info("Dry run complete. No files were written.")
        logging.info("No single-cell matrices or large files were opened.")
        return 0

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_symbol"])
        for gene in sorted(all_genes):
            writer.writerow([gene])

    logging.info("Wrote gene universe: %s", args.output)
    logging.info("No single-cell matrices or large files were opened.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
