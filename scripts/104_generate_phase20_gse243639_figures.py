#!/usr/bin/env python3
"""Generate Phase 20 GSE243639 safe-map repair figures with matplotlib."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "") else 0.0
    except ValueError:
        return 0.0


def save_bar(labels: list[str], values: list[float], title: str, ylabel: str, path: Path) -> None:
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, values, color="#4C78A8")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def figure_schema_audit(rows: list[dict[str, str]], path: Path) -> None:
    overlap_rows = [row for row in rows if row.get("audit_item") == "candidate_join_column_overlap"]
    labels = [row.get("candidate_join_column", "") for row in overlap_rows]
    values = [to_float(row.get("overlap_rate")) for row in overlap_rows]
    save_bar(labels or ["unavailable"], values or [0.0], "Phase 20 Annotation Map Join Overlap", "Overlap rate", path)


def figure_feature_groups(rows: list[dict[str, str]], path: Path) -> None:
    labels = [row.get("feature_group", "") for row in rows]
    values = [to_float(row.get("feature_count")) for row in rows]
    save_bar(labels or ["unavailable"], values or [0.0], "Phase 20 Feature Groups", "Feature count", path)


def figure_comparison(rows: list[dict[str, str]], path: Path) -> None:
    labels = [row.get("phase", "").replace("_", "\n") for row in rows]
    values = [to_float(row.get("auroc")) for row in rows]
    save_bar(labels or ["unavailable"], values or [0.0], "Phase 16/17/18/20 PD Validation", "AUROC", path)


def figure_validation(rows: list[dict[str, str]], path: Path) -> None:
    labels = [row.get("model", "") + "\n" + row.get("validation_mode", "") for row in rows]
    values = [to_float(row.get("balanced_accuracy")) for row in rows]
    save_bar(labels or ["unavailable"], values or [0.0], "Phase 20 PD Validation", "Balanced accuracy", path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 20 GSE243639 figures.")
    parser.add_argument("--schema-audit", type=Path, default=Path("results/tables/phase20_safe_annotation_map_schema_audit.tsv"))
    parser.add_argument("--feature-groups", type=Path, default=Path("results/tables/phase20_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--comparison", type=Path, default=Path("results/tables/phase20_pd_validation_comparison.tsv"))
    parser.add_argument("--metrics", type=Path, default=Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure_schema_audit(read_tsv(args.schema_audit), args.figures_dir / "figure47_phase20_annotation_map_audit.png")
    figure_feature_groups(read_tsv(args.feature_groups), args.figures_dir / "figure48_phase20_celltype_feature_groups.png")
    figure_comparison(read_tsv(args.comparison), args.figures_dir / "figure49_phase16_17_18_20_comparison.png")
    figure_validation(read_tsv(args.metrics), args.figures_dir / "figure50_phase20_pd_validation.png")
    print(f"Wrote Phase 20 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
