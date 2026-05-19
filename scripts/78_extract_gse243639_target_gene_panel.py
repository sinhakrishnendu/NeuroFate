#!/usr/bin/env python3
"""Stream GSE243639 target-gene counts into a sparse-like long table."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path
from typing import TextIO


MAX_TARGET_GENES = 64
OUTPUT_COLUMNS = ["cell_id", "sample_id", "gene_symbol", "expression_value"]


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


def open_text(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def sample_id_from_cell_id(cell_id: str) -> str:
    return cell_id.split("_", 1)[0].strip()


def read_panel_genes(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        genes = [row["gene_symbol"].strip() for row in csv.DictReader(handle, delimiter="\t") if row.get("gene_symbol")]
    if len(genes) > MAX_TARGET_GENES:
        raise ValueError(f"Target panel has {len(genes)} genes; maximum allowed is {MAX_TARGET_GENES}.")
    return genes


def is_nonzero(value: str) -> bool:
    try:
        return float(value) != 0.0
    except ValueError:
        return False


def write_cell_map(path: Path, cell_ids: list[str]) -> dict[str, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    sample_counts: dict[str, int] = {}
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_id", "sample_id"], delimiter="\t")
        writer.writeheader()
        for cell_id in cell_ids:
            sample_id = sample_id_from_cell_id(cell_id)
            sample_counts[sample_id] = sample_counts.get(sample_id, 0) + 1
            writer.writerow({"cell_id": cell_id, "sample_id": sample_id})
    return sample_counts


def write_audit(
    path: Path,
    requested: list[str],
    found: set[str],
    cell_count: int,
    sparse_rows: int,
    sample_counts: dict[str, int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    requested_upper = {gene.upper() for gene in requested}
    found_upper = {gene.upper() for gene in found}
    missing = sorted(requested_upper - found_upper)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset_id",
                "requested_target_genes",
                "requested_gene_symbols",
                "extracted_target_genes",
                "extracted_gene_symbols",
                "missing_target_genes",
                "missing_gene_symbols",
                "cell_columns",
                "sample_units",
                "sparse_expression_rows",
                "status",
                "notes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "dataset_id": "gse243639_pd_snpc",
                "requested_target_genes": str(len(requested)),
                "requested_gene_symbols": ",".join(requested),
                "extracted_target_genes": str(len(found_upper)),
                "extracted_gene_symbols": ",".join(sorted(found)),
                "missing_target_genes": str(len(missing)),
                "missing_gene_symbols": ",".join(missing),
                "cell_columns": str(cell_count),
                "sample_units": str(len(sample_counts)),
                "sparse_expression_rows": str(sparse_rows),
                "status": "ok" if found else "no_target_genes_extracted",
                "notes": "Streaming genes-as-rows CSV extraction; nonzero target-gene values only.",
            }
        )


def extract_target_genes(
    counts_path: Path,
    panel_path: Path,
    output_path: Path,
    audit_output: Path,
    cell_map_output: Path,
    dry_run: bool = False,
) -> tuple[int, int, set[str]]:
    requested = read_panel_genes(panel_path)
    requested_upper = {gene.upper() for gene in requested}
    sparse_rows = 0
    found: set[str] = set()
    with open_text(counts_path) as source:
        reader = csv.reader(source)
        header = next(reader)
        if len(header) < 2:
            raise RuntimeError("GSE243639 count table header does not contain cell columns.")
        cell_ids = header[1:]
        sample_counts = write_cell_map(cell_map_output, cell_ids)
        if dry_run:
            write_audit(audit_output, requested, found, len(cell_ids), 0, sample_counts)
            return len(cell_ids), 0, found
        with open_output(output_path) as target:
            writer = csv.DictWriter(target, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
            writer.writeheader()
            for row in reader:
                if not row:
                    continue
                gene_symbol = row[0].strip()
                if gene_symbol.upper() not in requested_upper:
                    continue
                found.add(gene_symbol)
                for cell_id, value in zip(cell_ids, row[1:], strict=False):
                    if not is_nonzero(value):
                        continue
                    writer.writerow(
                        {
                            "cell_id": cell_id,
                            "sample_id": sample_id_from_cell_id(cell_id),
                            "gene_symbol": gene_symbol,
                            "expression_value": value,
                        }
                    )
                    sparse_rows += 1
    write_audit(audit_output, requested, found, len(cell_ids), sparse_rows, sample_counts)
    return len(cell_ids), sparse_rows, found


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract GSE243639 target genes using streaming CSV logic.")
    parser.add_argument("--counts", type=Path, default=Path("data/raw/external/gse243639_pd_snpc/GSE243639_Filtered_count_table.csv.gz"))
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--output", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--audit-output", type=Path, default=Path("results/tables/phase16_gse243639_gene_extraction_audit.tsv"))
    parser.add_argument("--cell-map-output", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_cell_sample_map.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/78_extract_gse243639_target_gene_panel.log"))
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    cell_count, sparse_rows, found = extract_target_genes(
        args.counts,
        args.panel,
        args.output,
        args.audit_output,
        args.cell_map_output,
        args.dry_run,
    )
    logging.info("GSE243639 cell columns: %s", cell_count)
    logging.info("GSE243639 extracted target genes: %s", len(found))
    logging.info("GSE243639 sparse nonzero rows: %s", sparse_rows)
    print(f"Wrote {args.audit_output}")
    if not args.dry_run:
        print(f"Wrote {args.output}")
    print(f"Wrote {args.cell_map_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
