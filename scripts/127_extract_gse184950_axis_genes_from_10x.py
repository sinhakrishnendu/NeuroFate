#!/usr/bin/env python3
"""Manual-only sample-level extraction of NeuroFate axis genes from processed sparse 10x directories."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from pathlib import Path


EXPRESSION_COLUMNS = ["sample_id", "gene_symbol", "mean_expression", "detection_rate", "n_cells", "gene_sum", "detected_cells"]
AUDIT_COLUMNS = ["sample_id", "genes_requested", "genes_found", "found_gene_symbols", "missing_gene_symbols", "n_cells", "matrix_path", "features_path", "barcodes_path", "status"]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def axis_genes(axis_registry: Path) -> set[str]:
    genes: set[str] = set()
    for row in read_tsv(axis_registry):
        genes.update(gene.strip().upper() for gene in row.get("gene_members", "").replace(",", ";").split(";") if gene.strip())
    return genes


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.name.endswith(".gz") else path.open("r", encoding="utf-8")


def find_file(root: Path, names: tuple[str, ...]) -> Path | None:
    preferred_roots = [root / "filtered_feature_bc_matrix", root]
    for base in preferred_roots:
        for name in names:
            candidate = base / name
            if candidate.exists():
                return candidate
    for name in names:
        matches = list(root.rglob(name))
        if matches:
            return matches[0]
    return None


def count_barcodes(path: Path | None) -> int:
    if path is None:
        return 0
    with open_text(path) as handle:
        return sum(1 for line in handle if line.strip())


def read_features(path: Path, wanted: set[str]) -> dict[int, str]:
    found: dict[int, str] = {}
    with open_text(path) as handle:
        for index, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            symbols = [part.upper() for part in parts[:2]]
            match = next((symbol for symbol in symbols if symbol in wanted), None)
            if match:
                found[index] = match
    return found


def read_matrix_for_genes(matrix: Path, gene_rows: dict[int, str], barcode_count: int) -> tuple[dict[str, float], dict[str, int], int]:
    sums = {gene: 0.0 for gene in gene_rows.values()}
    detected = {gene: 0 for gene in gene_rows.values()}
    dimensions_seen = False
    n_cells = barcode_count
    with open_text(matrix) as handle:
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            if not dimensions_seen:
                dimensions_seen = True
                if not n_cells:
                    n_cells = int(parts[1])
                continue
            row_index = int(parts[0])
            if row_index not in gene_rows:
                continue
            value = float(parts[2])
            gene = gene_rows[row_index]
            sums[gene] += value
            if value > 0:
                detected[gene] += 1
    return sums, detected, n_cells


def candidate_sample_dirs(root: Path, metadata_by_sample: set[str]) -> list[Path]:
    if not root.exists():
        return []
    dirs = []
    for sample_id in sorted(metadata_by_sample):
        sample_dir = root / sample_id
        if sample_dir.exists() and sample_dir.is_dir():
            dirs.append(sample_dir)
    return dirs


def extract_sample(sample_dir: Path, wanted: set[str], metadata_samples: set[str]) -> tuple[list[dict[str, str]], dict[str, str]]:
    sample_id = sample_dir.name
    if sample_id not in metadata_samples:
        return [], {
            "sample_id": sample_id,
            "genes_requested": str(len(wanted)),
            "genes_found": "0",
            "found_gene_symbols": "",
            "missing_gene_symbols": ";".join(sorted(wanted)),
            "n_cells": "0",
            "matrix_path": "",
            "features_path": "",
            "barcodes_path": "",
            "status": "invalid_non_metadata_sample_skipped",
        }
    if not sample_dir.exists():
        return [], {
            "sample_id": sample_id,
            "genes_requested": str(len(wanted)),
            "genes_found": "0",
            "found_gene_symbols": "",
            "missing_gene_symbols": ";".join(sorted(wanted)),
            "n_cells": "0",
            "matrix_path": "",
            "features_path": "",
            "barcodes_path": "",
            "status": "missing_sample_directory",
        }
    matrix = find_file(sample_dir, ("matrix.mtx.gz", "matrix.mtx"))
    features = find_file(sample_dir, ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
    barcodes = find_file(sample_dir, ("barcodes.tsv.gz", "barcodes.tsv"))
    if matrix is None or features is None:
        return [], {
            "sample_id": sample_id,
            "genes_requested": str(len(wanted)),
            "genes_found": "0",
            "found_gene_symbols": "",
            "missing_gene_symbols": ";".join(sorted(wanted)),
            "n_cells": "0",
            "matrix_path": str(matrix or ""),
            "features_path": str(features or ""),
            "barcodes_path": str(barcodes or ""),
            "status": "missing_processed_10x_files",
        }
    gene_rows = read_features(features, wanted)
    found_genes = sorted(set(gene_rows.values()))
    missing_genes = sorted(wanted - set(found_genes))
    barcode_count = count_barcodes(barcodes)
    sums, detected, n_cells = read_matrix_for_genes(matrix, gene_rows, barcode_count)
    rows = []
    for gene in sorted(sums):
        denom = n_cells or 1
        rows.append(
            {
                "sample_id": sample_id,
                "gene_symbol": gene,
                "mean_expression": f"{sums[gene] / denom:.8g}",
                "detection_rate": f"{detected[gene] / denom:.8g}",
                "n_cells": str(n_cells),
                "gene_sum": f"{sums[gene]:.8g}",
                "detected_cells": str(detected[gene]),
            }
        )
    audit = {
        "sample_id": sample_id,
        "genes_requested": str(len(wanted)),
        "genes_found": str(len(sums)),
        "found_gene_symbols": ";".join(found_genes),
        "missing_gene_symbols": ";".join(missing_genes),
        "n_cells": str(n_cells),
        "matrix_path": str(matrix),
        "features_path": str(features),
        "barcodes_path": str(barcodes or ""),
        "status": "extracted_axis_genes_sample_level_only" if sums else "no_axis_genes_found",
    }
    return rows, audit


def write_gzip_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPRESSION_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract GSE184950 axis genes from manually prepared processed 10x directories.")
    parser.add_argument("--matrix-dir-root", type=Path, default=Path("data/interim/external/gse184950_pd_sn/processed_matrices"))
    parser.add_argument("--sample-metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output", type=Path, default=Path("data/interim/external/gse184950_pd_sn/gse184950_axis_gene_sample_summary.tsv.gz"))
    parser.add_argument("--audit-output", type=Path, default=Path("results/tables/phase27_gse184950_axis_gene_extraction_audit_clean.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/127_extract_gse184950_axis_genes_from_10x.log"))
    parser.add_argument("--run-manual-extraction", choices=["YES", "NO"], default="NO")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    if args.run_manual_extraction != "YES":
        raise SystemExit("Manual extraction guard active. Re-run with --run-manual-extraction YES only after selective processed matrices are prepared.")
    metadata = read_tsv(args.sample_metadata)
    if len(metadata) <= 2:
        logging.warning("GSE184950 metadata has only %d rows; prefer Phase 25 series metadata for the full cohort.", len(metadata))
    metadata_samples = {row.get("sample_name", "") for row in metadata if row.get("sample_name")}
    wanted = axis_genes(args.axis_registry)
    all_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    for sample in sorted(metadata_samples):
        sample_dir = args.matrix_dir_root / sample
        rows, audit = extract_sample(sample_dir, wanted, metadata_samples)
        all_rows.extend(rows)
        audit_rows.append(audit)
    write_gzip_tsv(args.output, all_rows)
    write_tsv(args.audit_output, audit_rows)
    logging.info("Extracted sample-level GSE184950 axis genes for samples=%d rows=%d", len(audit_rows), len(all_rows))
    logging.info("No raw-sequence, single-cell toolkit, embedding, clustering, or full dense matrix workflow was used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
