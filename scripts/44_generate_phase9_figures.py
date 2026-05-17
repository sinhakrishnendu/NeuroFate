#!/usr/bin/env python3
"""Generate Phase 9 Mathys CSV external-validation figures."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def save_bar(labels: list[str], values: list[float], title: str, ylabel: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(9.5, 5.5))
    positions = range(len(labels))
    plt.bar(positions, values, color="#4D6A8A")
    plt.xticks(positions, labels, rotation=45, ha="right", fontsize=8)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def figure20_gene_overlap(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "mathys_gene_overlap.tsv")
    counts = Counter(row.get("mathys_status", "missing") for row in rows)
    save_bar(
        ["present", "missing"],
        [counts.get("present", 0), counts.get("missing", 0)],
        "Mathys Target-Gene Panel Overlap",
        "Gene count",
        figures_dir / "figure20_mathys_gene_overlap.png",
    )


def figure21_external_validation(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "phase9_mathys_external_validation_metrics.tsv")
    save_bar(
        [row["validation_mode"] for row in rows],
        [to_float(row["auroc"]) for row in rows],
        "Mathys External Validation",
        "AUROC",
        figures_dir / "figure21_mathys_external_validation.png",
    )


def figure22_celltype_composition(tables_dir: Path, figures_dir: Path) -> None:
    rows = read_tsv(tables_dir / "mathys_2019_phase5_donor_feature_table.tsv")
    if not rows:
        save_bar([], [], "Mathys Cell-Type Composition", "Mean fraction", figures_dir / "figure22_mathys_celltype_composition.png")
        return
    celltype_columns = [field for field in rows[0] if field.startswith("cell_fraction__")]
    means = [
        (column.removeprefix("cell_fraction__"), sum(to_float(row.get(column, "0")) for row in rows) / len(rows))
        for column in celltype_columns
    ]
    means.sort(key=lambda item: item[1], reverse=True)
    top = means[:20]
    save_bar(
        [label for label, _ in top],
        [value for _, value in top],
        "Mathys Cell-Type Composition",
        "Mean sample fraction",
        figures_dir / "figure22_mathys_celltype_composition.png",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 9 Mathys figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure20_gene_overlap(args.tables_dir, args.figures_dir)
    figure21_external_validation(args.tables_dir, args.figures_dir)
    figure22_celltype_composition(args.tables_dir, args.figures_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
