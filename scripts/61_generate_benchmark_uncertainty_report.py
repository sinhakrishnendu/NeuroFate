#!/usr/bin/env python3
"""Generate a user-facing Phase 12 benchmark uncertainty report."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def table_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def task_stability(summary: pd.DataFrame) -> tuple[list[str], list[str], list[str]]:
    stable: list[str] = []
    unstable: list[str] = []
    failed: list[str] = []
    if summary.empty:
        return stable, unstable, ["benchmark_summary_missing"]
    for _, row in summary.iterrows():
        task = str(row.get("task", "unknown"))
        auroc = float(row.get("auroc_mean", 0) or 0)
        sd = float(row.get("auroc_sd", 999) or 999)
        if auroc >= 0.70 and sd <= 0.08:
            stable.append(task)
        elif auroc >= 0.60:
            unstable.append(task)
        else:
            failed.append(task)
    return sorted(set(stable)), sorted(set(unstable)), sorted(set(failed))


def best_models(summary: pd.DataFrame) -> list[str]:
    if summary.empty or "auroc_mean" not in summary:
        return ["Not available."]
    rows = []
    for task, group in summary.groupby("task"):
        best = group.sort_values("auroc_mean", ascending=False).iloc[0]
        rows.append(f"- {task}: {best['model']} (mean AUROC {best['auroc_mean']:.3f})")
    return rows


def write_report(
    summary_path: Path,
    pvalues_path: Path,
    ablation_path: Path,
    evidence_path: Path,
    output_path: Path,
) -> None:
    summary = table_or_empty(summary_path)
    pvalues = table_or_empty(pvalues_path)
    ablation = table_or_empty(ablation_path)
    evidence = table_or_empty(evidence_path)
    stable, unstable, failed = task_stability(summary)
    lines = [
        "# Phase 12 Benchmark Uncertainty Report",
        "",
        "This report summarizes donor-level benchmark robustness from existing Phase 12 tables.",
        "It does not rerun models or access single-cell matrices.",
        "",
        "## Stable Tasks",
        *(f"- {task}" for task in stable),
        "" if stable else "- None identified.",
        "",
        "## Unstable Tasks",
        *(f"- {task}" for task in unstable),
        "" if unstable else "- None identified.",
        "",
        "## Failed Or Missing Tasks",
        *(f"- {task}" for task in failed),
        "" if failed else "- None identified.",
        "",
        "## Best Models",
        *best_models(summary),
        "",
        "## Null-Model Comparison",
    ]
    if pvalues.empty:
        lines.append("- Permutation p-values are not available.")
    else:
        for _, row in pvalues.iterrows():
            lines.append(
                f"- {row.get('task')}: empirical p={float(row.get('empirical_pvalue', 1)):.3f}"
            )
    lines.extend(["", "## Feature Groups That Matter"])
    if ablation.empty:
        lines.append("- Feature ablation results are not available.")
    else:
        top = ablation.sort_values("delta_auroc_when_removed", ascending=False).head(10)
        for _, row in top.iterrows():
            lines.append(
                f"- {row.get('task')} / {row.get('feature_group')}: delta AUROC {float(row.get('delta_auroc_when_removed', 0)):.3f}"
            )
    lines.extend(["", "## Recommended Claim Strength"])
    if evidence.empty:
        lines.append("- Evidence strength matrix is not available.")
    else:
        for _, row in evidence.iterrows():
            lines.append(f"- {row.get('task')}: {row.get('evidence_category')} - {row.get('rationale')}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 12 uncertainty report.")
    parser.add_argument("--summary", type=Path, default=Path("results/tables/phase12_repeated_benchmark_summary.tsv"))
    parser.add_argument("--pvalues", type=Path, default=Path("results/tables/phase12_empirical_pvalues.tsv"))
    parser.add_argument("--ablation", type=Path, default=Path("results/tables/phase12_feature_group_importance.tsv"))
    parser.add_argument("--evidence", type=Path, default=Path("results/reports/evidence_strength_matrix.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase12_benchmark_uncertainty_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    write_report(args.summary, args.pvalues, args.ablation, args.evidence, args.output)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
