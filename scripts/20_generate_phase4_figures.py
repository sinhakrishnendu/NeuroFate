#!/usr/bin/env python3
"""Generate Phase 4 manuscript figures from Phase 4 statistical TSVs."""

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
    except (TypeError, ValueError):
        return 0.0


def save_bar(
    labels: list[str],
    values: list[float],
    title: str,
    ylabel: str,
    out_path: Path,
    color: str = "#476A6F",
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(10, 5.8))
    positions = range(len(labels))
    plt.bar(positions, values, color=color)
    plt.xticks(positions, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def figure5_braak_associations(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "phase4_gene_statistics.tsv")
        if row["test_variable"] == "Braak"
    ]
    rows.sort(
        key=lambda row: (
            to_float(row["fdr_p_value"]) if row["fdr_p_value"] != "nan" else 1.0,
            -abs(to_float(row["effect_size_max_minus_min"])),
        )
    )
    top = rows[:12]
    save_bar(
        [row["gene_symbol"] for row in top],
        [to_float(row["effect_size_max_minus_min"]) for row in top],
        "Braak-Associated Target Genes",
        "Max-minus-min donor mean expression",
        figures_dir / "figure5_braak_associations.png",
        "#6E5E4F",
    )


def figure6_apoe_microglia(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "phase4_apoe_analysis.tsv")
        if row["index_id"] == "MAI"
    ]
    rows.sort(key=lambda row: row["group_label"])
    save_bar(
        [row["group_label"] for row in rows],
        [to_float(row["donor_mean_index"]) for row in rows],
        "APOE Genotype and Microglial Activation Index",
        "Donor mean MAI",
        figures_dir / "figure6_apoe_microglia.png",
        "#4D6A8A",
    )


def figure7_celltype_vulnerability_heatmap(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase4_celltype_vulnerability.tsv")
    signatures = [
        "neurodegeneration_signature",
        "inflammatory_signature",
        "mitochondrial_dysfunction_signature",
    ]
    top_celltypes = [
        row["cell_subclass"]
        for row in sorted(
            [item for item in rows if item["signature_id"] == "neurodegeneration_signature"],
            key=lambda item: int(item["rank"]),
        )[:18]
    ]
    value_lookup = {
        (row["signature_id"], row["cell_subclass"]): to_float(row["donor_mean_index"])
        for row in rows
    }
    grid = [
        [value_lookup.get((signature, celltype), 0.0) for celltype in top_celltypes]
        for signature in signatures
    ]

    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(12, 4.8))
    image = plt.imshow(grid, aspect="auto", cmap="viridis")
    plt.colorbar(image, label="Donor mean index")
    plt.xticks(range(len(top_celltypes)), top_celltypes, rotation=45, ha="right", fontsize=8)
    plt.yticks(range(len(signatures)), signatures, fontsize=8)
    plt.title("Cell-Type Vulnerability Signatures")
    plt.tight_layout()
    plt.savefig(figures_dir / "figure7_celltype_vulnerability_heatmap.png", dpi=300)
    plt.close()


def figure8_composite_indices(tables_dir: Path, figures_dir: Path) -> None:
    rows = [
        row for row in read_tsv(tables_dir / "phase4_composite_indices.tsv")
        if row["grouping_variable"] == "Overall AD neuropathological Change"
    ]
    labels: list[str] = []
    values: list[float] = []
    for index_id in ["MAI", "ASI", "NVI"]:
        index_rows = sorted(
            [row for row in rows if row["index_id"] == index_id],
            key=lambda row: row["group_label"],
        )
        if not index_rows:
            continue
        labels.append(index_id)
        values.append(max(to_float(row["donor_mean_index"]) for row in index_rows))
    save_bar(
        labels,
        values,
        "Composite Indices Across AD Neuropathology",
        "Maximum donor mean index",
        figures_dir / "figure8_composite_indices.png",
        "#8A5A44",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 4 statistical biology figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure5_braak_associations(args.tables_dir, args.figures_dir)
    figure6_apoe_microglia(args.tables_dir, args.figures_dir)
    figure7_celltype_vulnerability_heatmap(args.tables_dir, args.figures_dir)
    figure8_composite_indices(args.tables_dir, args.figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
