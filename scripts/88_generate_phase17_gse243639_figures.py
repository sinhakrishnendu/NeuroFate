#!/usr/bin/env python3
"""Generate Phase 17 GSE243639 cell-type-aware PD figures using matplotlib."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

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


def save_bar(path: Path, labels: list[str], values: list[float], title: str, ylabel: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels or ["unavailable"], values or [0.0], color="#4f7f6f")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=30, labelsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def celltype_composition(summary: Path, output: Path) -> None:
    counts: dict[str, float] = {}
    for row in read_tsv(summary):
        cell_type = row.get("cell_type", "unannotated")
        counts[cell_type] = counts.get(cell_type, 0.0) + to_float(row.get("cell_count", "0"))
    top = sorted(counts.items(), key=lambda item: item[1], reverse=True)[:20]
    save_bar(output, [item[0] for item in top], [item[1] for item in top], "GSE243639 cell-type composition", "cell count")


def feature_importance(importance: Path, output: Path) -> None:
    rows = read_tsv(importance)[:20]
    save_bar(output, [row.get("feature", "") for row in rows], [to_float(row.get("importance", "0")) for row in rows], "Phase 17 feature importance", "importance")


def phase_comparison(comparison: Path, output: Path) -> None:
    rows = read_tsv(comparison)
    save_bar(output, [row.get("phase", "") for row in rows], [to_float(row.get("auroc", "0")) for row in rows], "Phase 16 vs Phase 17 PD signal", "AUROC")


def reliability(metrics: Path, output: Path) -> None:
    counts: dict[str, float] = {}
    for row in read_tsv(metrics):
        flag = row.get("reliability_flag", "unavailable")
        counts[flag] = counts.get(flag, 0.0) + 1.0
    save_bar(output, list(counts), list(counts.values()), "Phase 17 reliability categories", "row count")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 17 GSE243639 figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    celltype_composition(args.tables_dir / "phase17_gse243639_cell_annotation_summary.tsv", args.figures_dir / "figure39_gse243639_celltype_composition.png")
    feature_importance(args.tables_dir / "phase17_gse243639_celltype_feature_importance.tsv", args.figures_dir / "figure40_gse243639_celltype_feature_importance.png")
    phase_comparison(args.tables_dir / "phase17_pd_validation_comparison.tsv", args.figures_dir / "figure41_gse243639_phase16_vs_phase17.png")
    reliability(args.tables_dir / "phase17_gse243639_celltype_validation_metrics.tsv", args.figures_dir / "figure42_gse243639_pd_signal_reliability.png")
    print(f"Wrote Phase 17 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
