#!/usr/bin/env python3
"""Classify NeuroFate task evidence strength from Phase 12 benchmark outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EVIDENCE_CATEGORIES = [
    "strong_internal",
    "moderate_internal",
    "preliminary_external",
    "insufficient",
    "failed_or_unstable",
]


def table_or_empty(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def high_overclaiming_flags(path: Path) -> int:
    audit = table_or_empty(path)
    if audit.empty or "severity" not in audit:
        return 0
    return int((audit["severity"] == "high").sum())


def external_available(path: Path) -> bool:
    if not path.exists():
        return False
    frame = pd.read_csv(path, sep="\t")
    return len(frame) >= 6


def ablation_consistent(ablation: pd.DataFrame, task: str) -> bool:
    if ablation.empty or "delta_auroc_when_removed" not in ablation:
        return False
    task_rows = ablation[ablation["task"] == task]
    if task_rows.empty:
        return False
    deltas = pd.to_numeric(task_rows["delta_auroc_when_removed"], errors="coerce").dropna()
    return bool((deltas > 0.01).any())


def classify_task(
    row: pd.Series,
    pvalues: pd.DataFrame,
    ablation: pd.DataFrame,
    external_ok: bool,
    overclaiming_flags: int,
) -> dict[str, object]:
    task = str(row.get("task"))
    sample_size = int(row.get("n_samples_min", 0) or 0)
    auroc_mean = float(row.get("auroc_mean", 0) or 0)
    auroc_sd = float(row.get("auroc_sd", 999) or 999)
    p_row = pvalues[pvalues["task"] == task] if not pvalues.empty and "task" in pvalues else pd.DataFrame()
    empirical_pvalue = float(p_row["empirical_pvalue"].iloc[0]) if not p_row.empty else 1.0
    ablation_ok = ablation_consistent(ablation, task)
    if auroc_mean < 0.60 or auroc_sd > 0.15:
        category = "failed_or_unstable"
        rationale = "AUROC is weak or seed instability is high."
    elif sample_size >= 40 and auroc_mean >= 0.75 and auroc_sd <= 0.05 and empirical_pvalue <= 0.05 and ablation_ok and overclaiming_flags == 0:
        category = "strong_internal"
        rationale = "Large enough internal sample, stable AUROC, permutation support, and consistent ablation."
    elif sample_size >= 20 and auroc_mean >= 0.65 and auroc_sd <= 0.10 and empirical_pvalue <= 0.10:
        category = "moderate_internal"
        rationale = "Internal signal is present but not strong enough for high-confidence claims."
    elif external_ok:
        category = "preliminary_external"
        rationale = "External feasibility data exist, but evidence is not strong enough for definitive validation claims."
    else:
        category = "insufficient"
        rationale = "Benchmark evidence is incomplete or below threshold."
    return {
        "task": task,
        "evidence_category": category,
        "sample_size": sample_size,
        "auroc_mean": auroc_mean,
        "auroc_sd": auroc_sd,
        "empirical_pvalue": empirical_pvalue,
        "ablation_consistency": str(ablation_ok).lower(),
        "external_validation_available": str(external_ok).lower(),
        "no_overclaiming_high_flags": overclaiming_flags,
        "rationale": rationale,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify NeuroFate evidence strength.")
    parser.add_argument("--summary", type=Path, default=Path("results/tables/phase12_repeated_benchmark_summary.tsv"))
    parser.add_argument("--pvalues", type=Path, default=Path("results/tables/phase12_empirical_pvalues.tsv"))
    parser.add_argument("--ablation", type=Path, default=Path("results/tables/phase12_feature_group_importance.tsv"))
    parser.add_argument("--external", type=Path, default=Path("results/tables/mathys_2019_phase5_donor_feature_table.tsv"))
    parser.add_argument("--overclaiming", type=Path, default=Path("results/reports/no_overclaiming_audit.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/evidence_strength_matrix.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = table_or_empty(args.summary)
    pvalues = table_or_empty(args.pvalues)
    ablation = table_or_empty(args.ablation)
    external_ok = external_available(args.external)
    flags = high_overclaiming_flags(args.overclaiming)
    rows = []
    if not summary.empty:
        for _, row in summary.iterrows():
            rows.append(classify_task(row, pvalues, ablation, external_ok, flags))
    else:
        rows.append(
            {
                "task": "all_tasks",
                "evidence_category": "insufficient",
                "sample_size": 0,
                "auroc_mean": "",
                "auroc_sd": "",
                "empirical_pvalue": "",
                "ablation_consistency": "false",
                "external_validation_available": str(external_ok).lower(),
                "no_overclaiming_high_flags": flags,
                "rationale": "Repeated benchmark summary is missing.",
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
