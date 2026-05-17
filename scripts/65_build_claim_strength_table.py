#!/usr/bin/env python3
"""Build conservative NeuroFate claim-strength tables from existing metrics."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


CLAIM_COLUMNS = [
    "claim_id",
    "task",
    "model",
    "evidence_layer",
    "internal_auroc",
    "internal_auroc_sd",
    "internal_auprc",
    "balanced_accuracy",
    "brier_score",
    "permutation_empirical_p",
    "feature_ablation_support",
    "external_validation_status",
    "external_sample_units",
    "leakage_status",
    "overclaiming_status",
    "claim_strength",
    "allowed_claim_text",
    "disallowed_claim_text",
    "reviewer_risk",
]

BEST_COLUMNS = [
    "rank",
    "claim",
    "evidence_category",
    "primary_supporting_result",
    "supporting_tables",
    "safe_manuscript_sentence",
    "reviewer_caveat",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def best_repeated_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, dict[str, str]] = {}
    for row in rows:
        task = row.get("task", "unknown")
        current = grouped.get(task)
        if current is None or to_float(row.get("auroc_mean"), -1) > to_float(current.get("auroc_mean"), -1):
            grouped[task] = row
    return list(grouped.values())


def permutation_lookup(rows: list[dict[str, str]]) -> dict[str, float]:
    lookup: dict[str, float] = {}
    for row in rows:
        lookup[row.get("task", "")] = to_float(row.get("empirical_pvalue"), math.nan)
    return lookup


def ablation_support(rows: list[dict[str, str]], task: str) -> str:
    task_rows = [row for row in rows if row.get("task") == task]
    if not task_rows:
        return "unavailable"
    deltas = [to_float(row.get("delta_auroc_when_removed"), 0.0) for row in task_rows]
    return "supported" if any(delta > 0.01 for delta in deltas) else "weak_or_inconsistent"


def leakage_status(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "not_run"
    high = [
        row
        for row in rows
        if row.get("leakage_risk") == "high" and row.get("column_role") not in {"label", "identifier"}
    ]
    if high:
        return "potential_predictor_leakage"
    detected = any(row.get("leakage_risk") in {"high", "medium"} for row in rows)
    return "detected_and_excluded" if detected else "no_leakage_flags"


def overclaiming_status(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "not_run"
    high = [row for row in rows if row.get("severity") == "high" and row.get("allowed") != "true"]
    return "high_flags_present" if high else "no_high_flags"


def external_status(rows: list[dict[str, str]]) -> tuple[str, int]:
    if not rows:
        return "unavailable", 0
    sample_units = 0
    for row in rows:
        sample_units = max(sample_units, int(to_float(row.get("n_test"), 0)))
    if sample_units and sample_units < 20:
        return "preliminary_external_feasibility", sample_units
    if sample_units >= 20:
        return "external_support_requires_review", sample_units
    return "insufficient_external_validation", sample_units


def classify_claim(
    task: str,
    auroc: float,
    auroc_sd: float,
    pvalue: float,
    ablation: str,
    external: str,
    leakage: str,
    overclaiming: str,
) -> tuple[str, str, str, str]:
    no_leakage = leakage in {"detected_and_excluded", "no_leakage_flags", "not_run"}
    no_overclaim = overclaiming in {"no_high_flags", "not_run"}
    if not no_leakage:
        return (
            "failed_or_unstable",
            "NeuroFate predictor claims should be withheld until leakage is resolved.",
            "Leakage-unsafe predictive claims.",
            "high",
        )
    if auroc >= 0.75 and auroc_sd <= 0.05 and pvalue < 0.05 and ablation == "supported" and no_overclaim:
        strength = "strong_internal"
    elif auroc >= 0.68 and auroc_sd <= 0.10 and (math.isnan(pvalue) or pvalue < 0.10):
        strength = "moderate_internal"
    elif auroc >= 0.60 and auroc_sd <= 0.15:
        strength = "exploratory_internal"
    elif external == "preliminary_external_feasibility":
        strength = "preliminary_external_feasibility"
    elif external in {"unavailable", "insufficient_external_validation"}:
        strength = "insufficient_external_validation"
    else:
        strength = "failed_or_unstable"
    if task == "apoe_risk_prediction" and (auroc < 0.65 or (not math.isnan(pvalue) and pvalue >= 0.10)):
        strength = "failed_or_unstable"
    allowed = {
        "strong_internal": "NeuroFate shows strong internal donor-level evidence for this task.",
        "moderate_internal": "NeuroFate shows moderate internal donor-level evidence for this task.",
        "exploratory_internal": "NeuroFate shows exploratory internal signal that needs stronger validation.",
        "preliminary_external_feasibility": "NeuroFate shows preliminary external feasibility only.",
        "insufficient_external_validation": "NeuroFate evidence is insufficient for validation claims.",
        "failed_or_unstable": "NeuroFate results for this task are unstable or unsupported.",
    }[strength]
    disallowed = "Do not claim clinical utility, causality, foundation-model status, or definitive cross-cohort validation."
    risk = "low" if strength == "strong_internal" else "medium" if strength in {"moderate_internal", "exploratory_internal"} else "high"
    return strength, allowed, disallowed, risk


def build_claim_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    repeated = best_repeated_rows(read_tsv(args.repeated_summary))
    pvalues = permutation_lookup(read_tsv(args.pvalues))
    ablation_rows = read_tsv(args.ablation)
    external, external_n = external_status(read_tsv(args.external_metrics))
    leakage = leakage_status(read_tsv(args.leakage_audit))
    overclaiming = overclaiming_status(read_tsv(args.overclaiming_audit))
    claim_rows: list[dict[str, str]] = []
    for index, row in enumerate(repeated, start=1):
        task = row.get("task", "unknown")
        auroc = to_float(row.get("auroc_mean"))
        auroc_sd = to_float(row.get("auroc_sd"))
        auprc = to_float(row.get("auprc_mean"))
        bal = to_float(row.get("balanced_accuracy_mean"))
        brier = to_float(row.get("brier_mean"))
        pvalue = pvalues.get(task, math.nan)
        ablation = ablation_support(ablation_rows, task)
        strength, allowed, disallowed, risk = classify_claim(
            task, auroc, auroc_sd, pvalue, ablation, external, leakage, overclaiming
        )
        claim_rows.append(
            {
                "claim_id": f"claim_{index:03d}",
                "task": task,
                "model": row.get("model", ""),
                "evidence_layer": "donor_level_internal_plus_external_feasibility",
                "internal_auroc": f"{auroc:.6g}" if not math.isnan(auroc) else "",
                "internal_auroc_sd": f"{auroc_sd:.6g}" if not math.isnan(auroc_sd) else "",
                "internal_auprc": f"{auprc:.6g}" if not math.isnan(auprc) else "",
                "balanced_accuracy": f"{bal:.6g}" if not math.isnan(bal) else "",
                "brier_score": f"{brier:.6g}" if not math.isnan(brier) else "",
                "permutation_empirical_p": f"{pvalue:.6g}" if not math.isnan(pvalue) else "unavailable",
                "feature_ablation_support": ablation,
                "external_validation_status": external,
                "external_sample_units": str(external_n),
                "leakage_status": leakage,
                "overclaiming_status": overclaiming,
                "claim_strength": strength,
                "allowed_claim_text": allowed,
                "disallowed_claim_text": disallowed,
                "reviewer_risk": risk,
            }
        )
    if not claim_rows:
        claim_rows.append(
            {column: "" for column in CLAIM_COLUMNS}
            | {
                "claim_id": "claim_001",
                "task": "all_tasks",
                "evidence_layer": "missing_phase12_inputs",
                "external_validation_status": external,
                "external_sample_units": str(external_n),
                "leakage_status": leakage,
                "overclaiming_status": overclaiming,
                "claim_strength": "insufficient_external_validation",
                "allowed_claim_text": "Inputs are missing; no performance claim is supported.",
                "disallowed_claim_text": "Do not claim validation or clinical performance.",
                "reviewer_risk": "high",
            }
        )
    ranked = sorted(
        claim_rows,
        key=lambda row: (
            {"strong_internal": 0, "moderate_internal": 1, "exploratory_internal": 2}.get(
                row["claim_strength"], 3
            ),
            -to_float(row.get("internal_auroc"), 0),
        ),
    )
    best_rows: list[dict[str, str]] = []
    for rank, row in enumerate(ranked[:5], start=1):
        claim = f"{row['task']} supported by {row['model'] or 'available model'}"
        best_rows.append(
            {
                "rank": str(rank),
                "claim": claim,
                "evidence_category": row["claim_strength"],
                "primary_supporting_result": f"AUROC={row['internal_auroc']}, SD={row['internal_auroc_sd']}",
                "supporting_tables": "phase12_repeated_benchmark_summary.tsv;phase12_empirical_pvalues.tsv;phase12_feature_group_importance.tsv",
                "safe_manuscript_sentence": row["allowed_claim_text"],
                "reviewer_caveat": "External evidence remains preliminary when based on six Mathys sample-level units.",
            }
        )
    return claim_rows, best_rows


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build NeuroFate claim-strength tables.")
    parser.add_argument("--phase5-metrics", type=Path, default=Path("results/tables/phase5_model_metrics.tsv"))
    parser.add_argument("--phase6-metrics", type=Path, default=Path("results/tables/phase6_mps_model_metrics.tsv"))
    parser.add_argument("--repeated-summary", type=Path, default=Path("results/tables/phase12_repeated_benchmark_summary.tsv"))
    parser.add_argument("--pvalues", type=Path, default=Path("results/tables/phase12_empirical_pvalues.tsv"))
    parser.add_argument("--ablation", type=Path, default=Path("results/tables/phase12_feature_group_importance.tsv"))
    parser.add_argument("--external-metrics", type=Path, default=Path("results/tables/phase9_mathys_external_validation_metrics.tsv"))
    parser.add_argument("--leakage-audit", type=Path, default=Path("results/reports/feature_leakage_audit.tsv"))
    parser.add_argument("--overclaiming-audit", type=Path, default=Path("results/reports/no_overclaiming_audit.tsv"))
    parser.add_argument("--claim-output", type=Path, default=Path("results/reports/claim_strength_table.tsv"))
    parser.add_argument("--best-output", type=Path, default=Path("results/reports/best_supported_claims.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    claim_rows, best_rows = build_claim_tables(args)
    write_tsv(args.claim_output, claim_rows, CLAIM_COLUMNS)
    write_tsv(args.best_output, best_rows, BEST_COLUMNS)
    print(f"Wrote {args.claim_output}")
    print(f"Wrote {args.best_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
