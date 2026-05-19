#!/usr/bin/env python3
"""Generate Phase 32 cross-cohort evidence figures with matplotlib only."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIGURES = [
    Path("results/figures/figure71_crosscohort_axis_evidence_heatmap.png"),
    Path("results/figures/figure72_neuronal_vulnerability_axis_replication.png"),
    Path("results/figures/figure73_ad_vs_pd_replication_status.png"),
    Path("results/figures/figure74_pnas_readiness_summary.png"),
]


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


def write_skipped(reason: str) -> None:
    for path in FIGURES:
        marker = path.with_suffix(path.suffix + ".skipped.txt")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(reason + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 32 cross-cohort figures.")
    parser.add_argument("--summary", type=Path, default=Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        write_skipped("matplotlib unavailable; Phase 32 figure skipped.")
        return 0
    rows = read_tsv(args.summary)
    readiness = read_tsv(args.readiness)
    for path in FIGURES:
        path.parent.mkdir(parents=True, exist_ok=True)
    axes = [row.get("axis_id", "") for row in rows]
    matrix = [
        [to_float(row.get("sea_ad_effect")), to_float(row.get("gse174367_effect")), to_float(row.get("gse184950_effect"))]
        for row in rows
    ]
    plt.figure(figsize=(7, max(3, len(rows) * 0.35)))
    plt.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-0.6, vmax=0.6)
    plt.yticks(range(len(axes)), axes, fontsize=7)
    plt.xticks(range(3), ["SEA-AD", "GSE174367", "GSE184950"], rotation=20)
    plt.colorbar(label="Effect")
    plt.title("Cross-cohort axis evidence")
    plt.tight_layout()
    plt.savefig(FIGURES[0], dpi=200)
    plt.close()

    neuronal = next((row for row in rows if row.get("axis_id") == "neuronal_vulnerability_axis"), {})
    plt.figure(figsize=(5, 3))
    plt.bar(["SEA-AD", "GSE174367", "GSE184950"], [to_float(neuronal.get("sea_ad_effect")), to_float(neuronal.get("gse174367_effect")), to_float(neuronal.get("gse184950_effect"))], color=["#1d4ed8", "#15803d", "#f59e0b"])
    plt.axhline(0, color="#111827", linewidth=0.8)
    plt.ylabel("Effect")
    plt.title("Neuronal vulnerability axis")
    plt.tight_layout()
    plt.savefig(FIGURES[1], dpi=200)
    plt.close()

    class_counts: dict[str, int] = {}
    for row in rows:
        key = row.get("crosscohort_evidence_class", "missing")
        class_counts[key] = class_counts.get(key, 0) + 1
    plt.figure(figsize=(7, 3.5))
    plt.bar(range(len(class_counts)), list(class_counts.values()), color="#4b5563")
    plt.xticks(range(len(class_counts)), list(class_counts), rotation=35, ha="right", fontsize=7)
    plt.ylabel("Axes")
    plt.title("AD vs PD replication status")
    plt.tight_layout()
    plt.savefig(FIGURES[2], dpi=200)
    plt.close()

    statuses: dict[str, int] = {}
    for row in readiness:
        status = row.get("status", "missing")
        statuses[status] = statuses.get(status, 0) + 1
    plt.figure(figsize=(7, 3.5))
    plt.bar(range(len(statuses)), list(statuses.values()), color="#0f766e")
    plt.xticks(range(len(statuses)), list(statuses), rotation=35, ha="right", fontsize=7)
    plt.ylabel("Criteria")
    plt.title("PNAS readiness summary")
    plt.tight_layout()
    plt.savefig(FIGURES[3], dpi=200)
    plt.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
