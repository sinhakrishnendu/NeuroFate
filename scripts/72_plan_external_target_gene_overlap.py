#!/usr/bin/env python3
"""Plan target-gene overlap for external cohorts without loading count matrices."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


def setup_logging(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(filename=path, level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")


def read_panel(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [row["gene_symbol"].strip() for row in csv.DictReader(handle, delimiter="\t") if row.get("gene_symbol")]


def read_feature_file(path: Path) -> set[str]:
    opener = open
    mode = "r"
    if path.suffix == ".gz":
        import gzip

        opener = gzip.open
        mode = "rt"
    genes: set[str] = set()
    with opener(path, mode, encoding="utf-8", newline="") as handle:
        sample = handle.readline()
        delimiter = "\t" if "\t" in sample else ","
        handle.seek(0)
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, [])
        lower_header = [col.lower() for col in header]
        gene_index = 0
        for candidate in ["gene_symbol", "gene", "symbol", "features", "feature_name", "gene_name"]:
            if candidate in lower_header:
                gene_index = lower_header.index(candidate)
                break
        for row in reader:
            if row and len(row) > gene_index and row[gene_index]:
                genes.add(row[gene_index].strip())
    return genes


def write_outputs(dataset_id: str, panel: list[str], genes: set[str], output: Path, missing_output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    missing_output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset_id", "gene_symbol", "status"], delimiter="\t")
        writer.writeheader()
        for gene in panel:
            writer.writerow({"dataset_id": dataset_id, "gene_symbol": gene, "status": "found" if gene in genes else "missing"})
    with missing_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset_id", "gene_symbol"], delimiter="\t")
        writer.writeheader()
        for gene in panel:
            if gene not in genes:
                writer.writerow({"dataset_id": dataset_id, "gene_symbol": gene})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan external target-gene overlap.")
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--feature-file", type=Path, required=True)
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--missing-output", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_logging(args.log_file)
    panel = read_panel(args.panel)
    genes = read_feature_file(args.feature_file)
    write_outputs(args.dataset_id, panel, genes, args.output, args.missing_output)
    logging.info("Dataset %s gene overlap: found=%s missing=%s", args.dataset_id, len(set(panel) & genes), len(set(panel) - genes))
    print(f"Wrote {args.output}")
    print(f"Wrote {args.missing_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
