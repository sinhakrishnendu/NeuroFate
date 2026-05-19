#!/usr/bin/env python3
"""Generate Phase 38 GSE7621 direction-aware figures."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


FIGURES = [
    "figure75_gse7621_axis_effects.png",
    "figure76_gse7621_synuclein_mitochondrial_axis.png",
    "figure77_gse7621_neuronal_vulnerability_axis.png",
    "figure78_updated_pd_replication_status.png",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def skipped(outdir: Path, reason: str) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for name in FIGURES:
        (outdir / f"{name}.skipped.txt").write_text(reason + "\n", encoding="utf-8")


def axis_values(axis_scores: list[dict[str, str]], axis_id: str) -> tuple[list[float], list[float]]:
    col = f"axis__{axis_id}"
    controls = [to_float(row.get(col)) for row in axis_scores if row.get("label__pd_vs_control") == "0"]
    pds = [to_float(row.get(col)) for row in axis_scores if row.get("label__pd_vs_control") == "1"]
    return [v for v in controls if not math.isnan(v)], [v for v in pds if not math.isnan(v)]


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Phase 38 GSE7621 figures.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_axis_scores.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_pd_axis_replication_statistics.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/figures"))
    args = parser.parse_args()
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        skipped(args.outdir, f"matplotlib unavailable: {exc}")
        return 0

    args.outdir.mkdir(parents=True, exist_ok=True)
    stats = read_tsv(args.replication)
    scores = read_tsv(args.axis_scores)
    readiness = read_tsv(args.readiness)

    axes = [row.get("axis_id", "") for row in stats]
    effects = [to_float(row.get("effect_size"), 0.0) for row in stats]
    colors = ["#7b3294" if row.get("evidence_label") == "opposite_direction" and (to_float(row.get("pvalue"), 1.0) < 0.05 or to_float(row.get("fdr"), 1.0) < 0.1) else "#4575b4" for row in stats]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(range(len(axes)), effects, color=colors)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(range(len(axes)))
    ax.set_xticklabels([axis.replace("_axis", "").replace("_", "\n") for axis in axes], fontsize=8)
    ax.set_ylabel("Rank-biserial effect")
    ax.set_title("GSE7621 PD vs Control axis effects")
    fig.tight_layout()
    fig.savefig(args.outdir / "figure75_gse7621_axis_effects.png", dpi=200)
    plt.close(fig)

    for axis_id, filename, title in [
        ("synuclein_mitochondrial_axis", "figure76_gse7621_synuclein_mitochondrial_axis.png", "Synuclein-mitochondrial axis"),
        ("neuronal_vulnerability_axis", "figure77_gse7621_neuronal_vulnerability_axis.png", "Neuronal vulnerability axis"),
    ]:
        control, pd = axis_values(scores, axis_id)
        fig, ax = plt.subplots(figsize=(4.8, 4.2))
        ax.boxplot([control, pd], labels=["Control", "PD"], patch_artist=True, boxprops={"facecolor": "#d9f0d3"}, medianprops={"color": "black"})
        ax.scatter([1] * len(control), control, color="#4575b4", s=22, alpha=0.8)
        ax.scatter([2] * len(pd), pd, color="#d73027", s=22, alpha=0.8)
        ax.set_ylabel("Standardized axis score")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(args.outdir / filename, dpi=200)
        plt.close(fig)

    status_rows = [row for row in readiness if row.get("criterion") in {"independent_ad_replication", "independent_pd_replication", "pd_divergent_axis_candidate", "shared_ad_pd_axis_claim", "pnas_biological_claim"}]
    fig, ax = plt.subplots(figsize=(8, 4))
    labels = [row["criterion"].replace("_", "\n") for row in status_rows]
    status_to_value = {
        "statistically_supported": 3,
        "nominally_supported": 2,
        "present": 2,
        "mixed_pd_evidence": 1.5,
        "available_but_preliminary": 1,
        "promising_but_requires_pd_resolution": 1,
        "promising_but_not_ready": 1,
        "not_ready": 0,
        "missing_or_pending": 0,
    }
    values = [status_to_value.get(row.get("status", ""), 0) for row in status_rows]
    ax.bar(range(len(labels)), values, color=["#74add1", "#fdae61", "#7b3294", "#d73027", "#f46d43"][: len(labels)])
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 3.2)
    ax.set_ylabel("Readiness tier")
    ax.set_title("Updated replication/readiness status")
    fig.tight_layout()
    fig.savefig(args.outdir / "figure78_updated_pd_replication_status.png", dpi=200)
    plt.close(fig)
    print(f"Wrote Phase 38 figures to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
