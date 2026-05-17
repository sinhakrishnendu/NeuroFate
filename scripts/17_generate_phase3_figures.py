#!/usr/bin/env python3
"""Generate Phase 3 manuscript figures from sparse summary TSVs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def save_bar(labels: list[str], values: list[float], title: str, ylabel: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 5.5))
    positions = range(len(labels))
    plt.bar(positions, values, color="#476A6F")
    plt.xticks(positions, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def figure1_celltype_composition(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "gene_by_celltype_summary.tsv")
    cell_counts: dict[str, int] = {}
    for row in rows:
        cell_type = row["cell_type"]
        cell_counts.setdefault(cell_type, int(float(row["cell_count"])))
    top = sorted(cell_counts.items(), key=lambda item: item[1], reverse=True)[:20]
    save_bar(
        [label for label, _ in top],
        [value for _, value in top],
        "SEA-AD Cell-Type Composition",
        "Cell/nucleus count",
        figures_dir / "figure1_celltype_composition.png",
    )


def figure2_microglial_activation(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "microglial_activation_signature.tsv")
        if row["grouping_variable"] == "Overall AD neuropathological Change"
    ]
    rows = sorted(rows, key=lambda row: row["group_label"])
    save_bar(
        [row["group_label"] for row in rows],
        [to_float(row["mean_expression_per_gene_per_cell"]) for row in rows],
        "Microglial Activation Signature Across AD Pathology",
        "Mean expression per gene per cell",
        figures_dir / "figure2_microglial_activation.png",
    )


def figure3_neurodegeneration_signatures(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "neurodegeneration_signature_summary.tsv")
        if row["grouping_variable"] == "Overall AD neuropathological Change"
    ]
    rows = sorted(rows, key=lambda row: row["group_label"])
    save_bar(
        [row["group_label"] for row in rows],
        [to_float(row["mean_expression_per_gene_per_cell"]) for row in rows],
        "Neurodegeneration Signature Across AD Pathology",
        "Mean expression per gene per cell",
        figures_dir / "figure3_neurodegeneration_signatures.png",
    )


def figure4_ad_pathology_gene_trends(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "gene_by_ad_pathology.tsv")
        if row["pathology_variable"] == "Overall AD neuropathological Change"
    ]
    by_gene: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_gene.setdefault(row["gene_symbol"], []).append(row)
    ranked: list[tuple[str, float]] = []
    for gene, gene_rows in by_gene.items():
        values = [to_float(row["mean_expression_all_cells"]) for row in gene_rows]
        if values:
            ranked.append((gene, max(values) - min(values)))
    top_genes = [gene for gene, _ in sorted(ranked, key=lambda item: item[1], reverse=True)[:8]]
    labels: list[str] = []
    values: list[float] = []
    for gene in top_genes:
        gene_rows = sorted(by_gene[gene], key=lambda row: row["pathology_label"])
        if not gene_rows:
            continue
        labels.append(gene)
        values.append(max(to_float(row["mean_expression_all_cells"]) for row in gene_rows))
    save_bar(
        labels,
        values,
        "Top AD Pathology-Associated Target Genes",
        "Maximum mean expression across AD pathology groups",
        figures_dir / "figure4_ad_pathology_gene_trends.png",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 3 sparse-expression figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure1_celltype_composition(args.tables_dir, args.figures_dir)
    figure2_microglial_activation(args.tables_dir, args.figures_dir)
    figure3_neurodegeneration_signatures(args.tables_dir, args.figures_dir)
    figure4_ad_pathology_gene_trends(args.tables_dir, args.figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
