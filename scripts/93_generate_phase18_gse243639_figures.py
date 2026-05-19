#!/usr/bin/env python3
"""Generate Phase 18 GSE243639 repair figures with matplotlib only."""

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


def to_float(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "", "nan", "unavailable") else 0.0
    except ValueError:
        return 0.0


def save_bar(path: Path, labels: list[str], values: list[float], title: str, ylabel: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels or ["unavailable"], values or [0.0], color="#5f7f9f")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=30, labelsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def cell_id_audit(audit: Path, output: Path) -> None:
    rows = read_tsv(audit)
    keep = [
        "direct_overlap_count",
        "overlap_after_removing_sample_prefix",
        "overlap_after_adding_sample_prefix",
        "overlap_after_removing_trailing_dot_suffix",
        "overlap_after_normalizing_punctuation",
    ]
    values = {row.get("metric", ""): to_float(row.get("value")) for row in rows}
    save_bar(output, keep, [values.get(metric, 0.0) for metric in keep], "GSE243639 cell-ID matching audit", "overlap count")


def feature_groups(counts: Path, output: Path) -> None:
    rows = read_tsv(counts)
    save_bar(output, [row.get("feature_group", "") for row in rows], [to_float(row.get("feature_count")) for row in rows], "Repaired cell-type feature groups", "feature count")


def phase_comparison(comparison: Path, output: Path) -> None:
    rows = read_tsv(comparison)
    save_bar(output, [row.get("phase", "") for row in rows], [to_float(row.get("auroc")) for row in rows], "Phase 16/17/18 PD validation", "AUROC")


def repaired_validation(metrics: Path, output: Path) -> None:
    rows = read_tsv(metrics)
    save_bar(output, [row.get("validation_mode", "") for row in rows], [to_float(row.get("balanced_accuracy")) for row in rows], "Repaired PD validation", "balanced accuracy")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 18 GSE243639 figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cell_id_audit(args.tables_dir / "phase18_gse243639_cell_id_matching_audit.tsv", args.figures_dir / "figure43_gse243639_cell_id_matching_audit.png")
    feature_groups(args.tables_dir / "phase18_gse243639_feature_group_counts.tsv", args.figures_dir / "figure44_gse243639_repaired_celltype_features.png")
    phase_comparison(args.tables_dir / "phase18_pd_validation_comparison.tsv", args.figures_dir / "figure45_gse243639_phase16_17_18_comparison.png")
    repaired_validation(args.tables_dir / "phase18_gse243639_celltype_validation_metrics.tsv", args.figures_dir / "figure46_gse243639_repaired_pd_validation.png")
    print(f"Wrote Phase 18 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
