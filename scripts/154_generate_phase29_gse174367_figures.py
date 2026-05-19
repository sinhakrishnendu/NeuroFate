#!/usr/bin/env python3
"""Generate Phase 29 GSE174367 bulk AD replication figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIGURES = {
    "sample_summary": Path("results/figures/figure67_gse174367_bulk_sample_summary.png"),
    "axis_coverage": Path("results/figures/figure68_gse174367_bulk_axis_coverage.png"),
    "effects": Path("results/figures/figure69_gse174367_ad_axis_replication_effects.png"),
    "readiness": Path("results/figures/figure70_updated_pnas_ad_replication_status.png"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def write_skipped(path: Path, reason: str) -> None:
    marker = path.with_suffix(path.suffix + ".skipped.txt")
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(reason + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 29 GSE174367 figures.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase29_gse174367_bulk_axis_scores.tsv"))
    parser.add_argument("--coverage", type=Path, default=Path("results/tables/phase29_gse174367_bulk_axis_feature_coverage.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        for path in FIGURES.values():
            write_skipped(path, "matplotlib unavailable; Phase 29 figure skipped.")
        return 0

    labels: dict[str, int] = {}
    for row in read_tsv(args.axis_scores):
        key = row.get("label__ad_vs_control", "missing") or "missing"
        labels[key] = labels.get(key, 0) + 1
    FIGURES["sample_summary"].parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(5, 3))
    plt.bar(list(labels), list(labels.values()), color=["#6b7280", "#b91c1c", "#9ca3af"][: len(labels)])
    plt.title("GSE174367 AD endpoint")
    plt.ylabel("Samples")
    plt.tight_layout()
    plt.savefig(FIGURES["sample_summary"], dpi=200)
    plt.close()

    coverage = read_tsv(args.coverage)
    axes = [row.get("axis_id", "") for row in coverage]
    found = [to_float(row.get("genes_found")) for row in coverage]
    plt.figure(figsize=(8, 4))
    plt.bar(range(len(axes)), found, color="#2563eb")
    plt.xticks(range(len(axes)), axes, rotation=75, ha="right", fontsize=7)
    plt.ylabel("Genes found")
    plt.title("Bulk RNA axis-gene coverage")
    plt.tight_layout()
    plt.savefig(FIGURES["axis_coverage"], dpi=200)
    plt.close()

    replication = read_tsv(args.replication)
    axes = [row.get("axis_id", "") for row in replication]
    effects = [to_float(row.get("effect_size")) for row in replication]
    colors = ["#15803d" if row.get("evidence_label") == "statistically_supported_ad_replication" else "#f59e0b" for row in replication]
    plt.figure(figsize=(8, 4))
    plt.bar(range(len(axes)), effects, color=colors)
    plt.axhline(0, color="#111827", linewidth=0.8)
    plt.xticks(range(len(axes)), axes, rotation=75, ha="right", fontsize=7)
    plt.ylabel("Rank-biserial effect")
    plt.title("GSE174367 endpoint-locked AD replication")
    plt.tight_layout()
    plt.savefig(FIGURES["effects"], dpi=200)
    plt.close()

    readiness = read_tsv(args.readiness)
    statuses: dict[str, int] = {}
    for row in readiness:
        status = row.get("status", "missing")
        statuses[status] = statuses.get(status, 0) + 1
    plt.figure(figsize=(7, 3.5))
    plt.bar(range(len(statuses)), list(statuses.values()), color="#4b5563")
    plt.xticks(range(len(statuses)), list(statuses), rotation=35, ha="right", fontsize=7)
    plt.ylabel("Criteria")
    plt.title("PNAS readiness status")
    plt.tight_layout()
    plt.savefig(FIGURES["readiness"], dpi=200)
    plt.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
