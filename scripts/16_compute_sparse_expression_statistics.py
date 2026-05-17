#!/usr/bin/env python3
"""Compute Phase 3 sparse expression summaries from extracted target genes.

Inputs are the sparse long-format target-gene expression TSV and decoded metadata
TSV only. This script never opens H5AD, never creates dense matrices, and never
runs single-cell analysis pipelines.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from typing import TextIO


DEFAULT_MAX_ROW_CHUNK = 1_000_000
MAX_ALLOWED_ROW_CHUNK = 5_000_000

CELLTYPE_FIELD = "Subclass"
PATHOLOGY_FIELDS = ["Braak", "CERAD score", "Overall AD neuropathological Change"]
COGNITIVE_FIELD = "Cognitive Status"
METADATA_FIELDS = [
    "Class",
    "Subclass",
    "Supertype",
    "Braak",
    "CERAD score",
    "Overall AD neuropathological Change",
    "Cognitive Status",
]

MICROGLIAL_GENES = ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1"]
ASTROCYTE_GENES = ["GFAP"]
NEURONAL_GENES = ["SLC17A7", "SST", "PVALB", "LAMP5"]
NEURODEGENERATION_GENES = ["PINK1", "PRKN", "SNCA", "MAPT", "APOE"]

SIGNATURES = {
    "microglial_activation": MICROGLIAL_GENES,
    "astrocyte_stress": ASTROCYTE_GENES,
    "neuronal_signature": NEURONAL_GENES,
    "neurodegeneration_signature": NEURODEGENERATION_GENES,
}

SIGNATURE_OUTPUTS = {
    "microglial_activation": "microglial_activation_signature.tsv",
    "astrocyte_stress": "astrocyte_stress_signature.tsv",
    "neuronal_signature": "neuronal_signature_summary.tsv",
    "neurodegeneration_signature": "neurodegeneration_signature_summary.tsv",
}


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
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_panel_genes(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [row["gene_symbol"] for row in csv.DictReader(handle, delimiter="\t")]


def clean_label(value: str | None) -> str:
    if value is None:
        return "missing"
    value = value.strip()
    return value if value else "missing"


def load_metadata_vectors(metadata_path: Path) -> tuple[dict[str, list[str]], dict[tuple[str, str], int]]:
    vectors: dict[str, list[str]] = {field: [] for field in METADATA_FIELDS}
    vectors["cell_type"] = []
    group_counts: dict[tuple[str, str], int] = defaultdict(int)

    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for field in METADATA_FIELDS:
                vectors[field].append(clean_label(row.get(field)))
            cell_type = clean_label(row.get(CELLTYPE_FIELD) or row.get("Class") or row.get("Supertype"))
            vectors["cell_type"].append(cell_type)
            group_counts[("cell_type", cell_type)] += 1
            for field in PATHOLOGY_FIELDS:
                group_counts[(field, clean_label(row.get(field)))] += 1
            group_counts[(COGNITIVE_FIELD, clean_label(row.get(COGNITIVE_FIELD)))] += 1

    return vectors, group_counts


def empty_stat() -> dict[str, float]:
    return {"nonzero_count": 0.0, "total_expression": 0.0}


def update_stat(stats: dict[tuple[str, ...], dict[str, float]], key: tuple[str, ...], value: float) -> None:
    stats[key]["nonzero_count"] += 1
    stats[key]["total_expression"] += value


def stream_expression_statistics(
    expression_path: Path,
    target_genes: list[str],
    metadata_vectors: dict[str, list[str]],
    max_row_chunk: int,
) -> tuple[
    dict[tuple[str, str], dict[str, float]],
    dict[tuple[str, str, str], dict[str, float]],
    dict[tuple[str, str], dict[str, float]],
    dict[tuple[str, str, str], dict[str, float]],
    set[str],
]:
    if max_row_chunk > MAX_ALLOWED_ROW_CHUNK:
        raise ValueError(f"max_row_chunk may not exceed {MAX_ALLOWED_ROW_CHUNK}")

    target_set = set(target_genes)
    signature_lookup: dict[str, list[str]] = defaultdict(list)
    for signature_name, genes in SIGNATURES.items():
        for gene in genes:
            signature_lookup[gene].append(signature_name)

    gene_celltype_stats: dict[tuple[str, str], dict[str, float]] = defaultdict(empty_stat)
    gene_pathology_stats: dict[tuple[str, str, str], dict[str, float]] = defaultdict(empty_stat)
    gene_cognitive_stats: dict[tuple[str, str], dict[str, float]] = defaultdict(empty_stat)
    signature_stats: dict[tuple[str, str, str], dict[str, float]] = defaultdict(empty_stat)
    seen_genes: set[str] = set()
    processed = 0

    with open_text(expression_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            processed += 1
            if processed % max_row_chunk == 0:
                logging.info("Processed sparse expression rows: %d", processed)

            gene = row.get("gene_symbol", "")
            if gene not in target_set:
                continue
            row_index = int(row["row_index"])
            if row_index >= len(metadata_vectors["cell_type"]):
                logging.warning("Skipping expression row with out-of-range row_index=%d", row_index)
                continue
            value = float(row["expression_value"])
            seen_genes.add(gene)

            cell_type = metadata_vectors["cell_type"][row_index]
            update_stat(gene_celltype_stats, (gene, cell_type), value)

            for field in PATHOLOGY_FIELDS:
                label = metadata_vectors[field][row_index]
                update_stat(gene_pathology_stats, (gene, field, label), value)
                for signature_name in signature_lookup.get(gene, []):
                    update_stat(signature_stats, (signature_name, field, label), value)

            cognitive_label = metadata_vectors[COGNITIVE_FIELD][row_index]
            update_stat(gene_cognitive_stats, (gene, cognitive_label), value)
            for signature_name in signature_lookup.get(gene, []):
                update_stat(signature_stats, (signature_name, "cell_type", cell_type), value)
                update_stat(signature_stats, (signature_name, COGNITIVE_FIELD, cognitive_label), value)

    logging.info("Total sparse expression rows processed: %d", processed)
    return gene_celltype_stats, gene_pathology_stats, gene_cognitive_stats, signature_stats, seen_genes


def write_gene_by_celltype(
    path: Path,
    target_genes: list[str],
    group_counts: dict[tuple[str, str], int],
    stats: dict[tuple[str, str], dict[str, float]],
) -> None:
    cell_types = sorted(label for variable, label in group_counts if variable == "cell_type")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_symbol",
                "cell_type",
                "cell_count",
                "expressing_cell_count",
                "detection_rate",
                "total_expression",
                "mean_expression_all_cells",
                "mean_expression_expressing_cells",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for gene in target_genes:
            for cell_type in cell_types:
                denominator = group_counts[("cell_type", cell_type)]
                stat = stats.get((gene, cell_type), empty_stat())
                write_gene_stat(writer, gene, "cell_type", cell_type, denominator, stat)


def write_gene_by_pathology(
    path: Path,
    target_genes: list[str],
    group_counts: dict[tuple[str, str], int],
    stats: dict[tuple[str, str, str], dict[str, float]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_symbol",
                "pathology_variable",
                "pathology_label",
                "cell_count",
                "expressing_cell_count",
                "detection_rate",
                "total_expression",
                "mean_expression_all_cells",
                "mean_expression_expressing_cells",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for gene in target_genes:
            for variable in PATHOLOGY_FIELDS:
                labels = sorted(label for field, label in group_counts if field == variable)
                for label in labels:
                    denominator = group_counts[(variable, label)]
                    stat = stats.get((gene, variable, label), empty_stat())
                    write_gene_stat(writer, gene, variable, label, denominator, stat)


def write_gene_by_cognitive_status(
    path: Path,
    target_genes: list[str],
    group_counts: dict[tuple[str, str], int],
    stats: dict[tuple[str, str], dict[str, float]],
) -> None:
    labels = sorted(label for field, label in group_counts if field == COGNITIVE_FIELD)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "gene_symbol",
                "cognitive_status",
                "cell_count",
                "expressing_cell_count",
                "detection_rate",
                "total_expression",
                "mean_expression_all_cells",
                "mean_expression_expressing_cells",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for gene in target_genes:
            for label in labels:
                denominator = group_counts[(COGNITIVE_FIELD, label)]
                stat = stats.get((gene, label), empty_stat())
                write_gene_stat(writer, gene, COGNITIVE_FIELD, label, denominator, stat)


def write_gene_stat(
    writer: csv.DictWriter,
    gene: str,
    variable: str,
    label: str,
    denominator: int,
    stat: dict[str, float],
) -> None:
    nonzero = int(stat["nonzero_count"])
    total = float(stat["total_expression"])
    row = {
        "gene_symbol": gene,
        "cell_count": denominator,
        "expressing_cell_count": nonzero,
        "detection_rate": f"{nonzero / denominator:.8f}" if denominator else "0.00000000",
        "total_expression": f"{total:.6f}",
        "mean_expression_all_cells": f"{total / denominator:.8f}" if denominator else "0.00000000",
        "mean_expression_expressing_cells": f"{total / nonzero:.8f}" if nonzero else "0.00000000",
    }
    if variable == "cell_type":
        row["cell_type"] = label
    elif variable == COGNITIVE_FIELD:
        row["cognitive_status"] = label
    else:
        row["pathology_variable"] = variable
        row["pathology_label"] = label
    writer.writerow(row)


def write_signature_table(
    path: Path,
    signature_name: str,
    signature_genes: list[str],
    group_counts: dict[tuple[str, str], int],
    stats: dict[tuple[str, str, str], dict[str, float]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "signature_name",
                "grouping_variable",
                "group_label",
                "cell_count",
                "gene_count",
                "genes",
                "nonzero_gene_cell_count",
                "detection_fraction",
                "total_expression",
                "mean_expression_per_gene_per_cell",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for variable in ["cell_type", *PATHOLOGY_FIELDS, COGNITIVE_FIELD]:
            labels = sorted(label for field, label in group_counts if field == variable)
            for label in labels:
                denominator = group_counts[(variable, label)]
                stat = stats.get((signature_name, variable, label), empty_stat())
                gene_count = len(signature_genes)
                denominator_gene_cells = denominator * gene_count
                nonzero = int(stat["nonzero_count"])
                total = float(stat["total_expression"])
                writer.writerow(
                    {
                        "signature_name": signature_name,
                        "grouping_variable": variable,
                        "group_label": label,
                        "cell_count": denominator,
                        "gene_count": gene_count,
                        "genes": ",".join(signature_genes),
                        "nonzero_gene_cell_count": nonzero,
                        "detection_fraction": (
                            f"{nonzero / denominator_gene_cells:.8f}"
                            if denominator_gene_cells
                            else "0.00000000"
                        ),
                        "total_expression": f"{total:.6f}",
                        "mean_expression_per_gene_per_cell": (
                            f"{total / denominator_gene_cells:.8f}"
                            if denominator_gene_cells
                            else "0.00000000"
                        ),
                    }
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute sparse target-gene Phase 3 statistics.")
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv"),
    )
    parser.add_argument("--panel", type=Path, default=Path("metadata/target_gene_panel_v1.tsv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/16_compute_sparse_expression_statistics.log"),
    )
    parser.add_argument("--max-row-chunk", type=int, default=DEFAULT_MAX_ROW_CHUNK)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting Phase 3 sparse expression statistics.")
    logging.info("Expression TSV: %s", args.expression)
    logging.info("Metadata TSV: %s", args.metadata)
    logging.info("Target panel: %s", args.panel)

    args.tables_dir.mkdir(parents=True, exist_ok=True)
    target_genes = read_panel_genes(args.panel)
    metadata_vectors, group_counts = load_metadata_vectors(args.metadata)
    (
        gene_celltype_stats,
        gene_pathology_stats,
        gene_cognitive_stats,
        signature_stats,
        seen_genes,
    ) = stream_expression_statistics(args.expression, target_genes, metadata_vectors, args.max_row_chunk)

    write_gene_by_celltype(
        args.tables_dir / "gene_by_celltype_summary.tsv",
        target_genes,
        group_counts,
        gene_celltype_stats,
    )
    write_gene_by_pathology(
        args.tables_dir / "gene_by_ad_pathology.tsv",
        target_genes,
        group_counts,
        gene_pathology_stats,
    )
    write_gene_by_cognitive_status(
        args.tables_dir / "gene_by_cognitive_status.tsv",
        target_genes,
        group_counts,
        gene_cognitive_stats,
    )
    for signature_name, genes in SIGNATURES.items():
        write_signature_table(
            args.tables_dir / SIGNATURE_OUTPUTS[signature_name],
            signature_name,
            genes,
            group_counts,
            signature_stats,
        )

    logging.info("Target panel genes: %d", len(target_genes))
    logging.info("Genes observed in sparse expression file: %d", len(seen_genes))
    logging.info("Observed genes: %s", ", ".join(sorted(seen_genes)) if seen_genes else "none")
    logging.info("No dense matrices, Scanpy pipelines, clustering, or model training were used.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
