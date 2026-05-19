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

PHASE17_DELTA_COLUMNS = [
    "task",
    "phase17_status",
    "phase17_sample_units",
    "phase17_auroc",
    "phase17_permutation_pvalue",
    "affected_existing_claim",
    "upgrade_applied",
    "reason",
]

PHASE18_DELTA_COLUMNS = [
    "task",
    "phase18_status",
    "phase18_sample_units",
    "phase18_auroc",
    "phase18_permutation_pvalue",
    "phase18_supersedes_phase17",
    "affected_existing_claim",
    "upgrade_applied",
    "reason",
]

PHASE20_DELTA_COLUMNS = [
    "task",
    "phase20_status",
    "phase20_sample_units",
    "phase20_feature_count",
    "phase20_annotation_match_rate",
    "phase20_auroc",
    "phase20_auprc",
    "phase20_balanced_accuracy",
    "phase20_permutation_pvalue",
    "phase20_supersedes_phase17_phase18",
    "affected_existing_claim",
    "upgrade_applied",
    "reason",
]

PHASE21_DELTA_COLUMNS = [
    "axis_id",
    "axis_classification",
    "claim_strength",
    "affected_existing_claim",
    "upgrade_applied",
    "reason",
]

PHASE22_DELTA_COLUMNS = [
    "axis_id",
    "endpoint_id",
    "axis_claim_class",
    "affected_phase21_claim",
    "upgrade_applied",
    "reason",
]

PHASE27_DELTA_COLUMNS = [
    "cohort_id",
    "statistically_supported_axes",
    "directionally_consistent_preliminary_axes",
    "weak_or_no_replication_axes",
    "claim_upgrade_allowed",
    "reason",
]


def append_phase32_crosscohort_claims(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    rows = read_tsv(args.phase32_crosscohort_evidence)
    if not rows:
        return
    for evidence in rows:
        axis_id = evidence.get("axis_id", "")
        evidence_class = evidence.get("crosscohort_evidence_class", "insufficient_data")
        if not axis_id:
            continue
        if evidence_class == "strong_ad_axis_with_nominal_external_replication":
            strength = "preliminary_external_feasibility"
            risk = "medium"
            allowed = (
                f"{axis_id} has nominal independent AD replication in GSE174367 and directionally consistent SEA-AD support; "
                "describe it as an AD replication candidate, not as a definitive mechanism."
            )
        elif evidence_class == "strong_ad_axis_without_external_replication":
            strength = "axis_level_preliminary_evidence"
            risk = "medium"
            allowed = f"{axis_id} has SEA-AD endpoint-locked support but still requires independent AD replication."
        elif evidence_class == "preliminary_shared_ad_pd_axis_candidate":
            strength = "axis_level_preliminary_evidence"
            risk = "medium"
            allowed = (
                f"{axis_id} may be described as an exploratory cross-disease convergence candidate; "
                "additional PD replication is required."
            )
        else:
            strength = "axis_level_insufficient_validation"
            risk = "high"
            allowed = evidence.get("safe_claim", f"{axis_id} remains insufficiently supported.")
        claim_rows.append(
            {column: "" for column in CLAIM_COLUMNS}
            | {
                "claim_id": f"claim_phase32_{axis_id}",
                "task": f"crosscohort_axis:{axis_id}",
                "model": "endpoint_locked_axis_association",
                "evidence_layer": "phase32_crosscohort_axis_consolidation",
                "permutation_empirical_p": "",
                "feature_ablation_support": "phase22_matched_random_controls_context",
                "external_validation_status": evidence_class,
                "external_sample_units": "90" if evidence.get("gse174367_p") else "",
                "leakage_status": leakage_status(read_tsv(args.leakage_audit)),
                "overclaiming_status": overclaiming_status(read_tsv(args.overclaiming_audit)),
                "claim_strength": strength,
                "allowed_claim_text": allowed,
                "disallowed_claim_text": (
                    "Do not claim clinical utility, diagnostic use, causality, definitive shared AD/PD mechanism, "
                    "or FDR-robust external replication unless future evidence supports it."
                ),
                "reviewer_risk": risk,
            }
        )


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


def phase15_external_status(rows: list[dict[str, str]]) -> tuple[str, int, str]:
    if not rows:
        return "unavailable", 0, "no_phase15_metrics"
    pd_internal = [row for row in rows if row.get("reliability_flag") == "moderate_pd_internal_validation"]
    if pd_internal:
        units = max(int(to_float(row.get("n_samples") or row.get("n_test"), 0)) for row in pd_internal)
        return "moderate_pd_internal_validation", units, "phase16_pd_internal_validation_available"
    reliable = [row for row in rows if row.get("reliability_flag") == "reliable_external_validation"]
    if reliable:
        units = max(int(to_float(row.get("n_test"), 0)) for row in reliable)
        return "reliable_external_validation", units, "phase15_reliable_external_validation_available"
    cross_disease = [row for row in rows if row.get("reliability_flag") == "preliminary_cross_disease_feature_transfer"]
    if cross_disease:
        units = max(int(to_float(row.get("n_samples") or row.get("n_test"), 0)) for row in cross_disease)
        return "preliminary_cross_disease_feature_transfer", units, "phase16_cross_disease_feature_transfer_only"
    preliminary = [row for row in rows if row.get("reliability_flag") == "preliminary_external_feasibility"]
    if preliminary:
        units = max(int(to_float(row.get("n_test"), 0)) for row in preliminary)
        return "preliminary_external_feasibility", units, "phase15_external_feasibility_only"
    units = max([int(to_float(row.get("n_test"), 0)) for row in rows] or [0])
    return "insufficient_external_validation", units, "phase15_external_metrics_not_reliable"


def phase17_pd_status(rows: list[dict[str, str]]) -> tuple[str, int, str, str, str]:
    if not rows:
        return "unavailable", 0, "", "", "no_phase17_metrics"
    preferred = next(
        (
            row
            for row in rows
            if row.get("model") == "logistic_regression"
            and row.get("validation_mode") == "repeated_stratified_split"
        ),
        rows[0],
    )
    status = preferred.get("reliability_flag", "unavailable")
    units = int(to_float(preferred.get("n_samples"), 0))
    auroc = preferred.get("auroc", "")
    pvalue = preferred.get("empirical_permutation_pvalue", "")
    if status == "moderate_pd_internal_validation":
        reason = "phase17_celltype_pd_signal_robust_enough_for_moderate_internal_pd_extension"
    elif status == "preliminary_pd_internal_signal":
        reason = "phase17_celltype_pd_signal_preliminary_only"
    else:
        reason = "phase17_celltype_pd_signal_weak_or_unavailable"
    return status, units, auroc, pvalue, reason


def phase18_pd_status(rows: list[dict[str, str]]) -> tuple[str, int, str, str, str]:
    if not rows:
        return "unavailable", 0, "", "", "no_phase18_metrics"
    preferred = next(
        (
            row
            for row in rows
            if row.get("model") == "logistic_regression"
            and row.get("validation_mode") == "repeated_stratified_split"
        ),
        rows[0],
    )
    status = preferred.get("reliability_flag", "unavailable")
    units = int(to_float(preferred.get("n_samples"), 0))
    auroc = preferred.get("auroc", "")
    pvalue = preferred.get("empirical_permutation_pvalue", "")
    if status == "moderate_pd_internal_validation":
        reason = "phase18_repaired_celltype_pd_signal_robust_enough_for_moderate_internal_pd_extension"
    elif status == "preliminary_pd_internal_signal":
        reason = "phase18_repaired_celltype_pd_signal_preliminary_only"
    elif status == "technical_failure_annotation_join":
        reason = "phase18_annotation_join_technical_failure"
    else:
        reason = "phase18_repaired_celltype_pd_signal_weak_or_unavailable"
    return status, units, auroc, pvalue, reason


def phase20_pd_status(rows: list[dict[str, str]], feature_group_rows: list[dict[str, str]]) -> tuple[str, int, str, str, str, str, str, str, str]:
    if not rows:
        return "unavailable", 0, "", "", "", "", "", "", "no_phase20_metrics"
    preferred = next(
        (
            row
            for row in rows
            if row.get("model") == "logistic_regression"
            and row.get("validation_mode") == "repeated_stratified_split"
        ),
        rows[0],
    )
    status = preferred.get("reliability_flag", "unavailable")
    units = int(to_float(preferred.get("n_samples"), 0))
    feature_count = preferred.get("feature_count", "")
    if not feature_count and feature_group_rows:
        feature_count = str(sum(int(to_float(row.get("feature_count"), 0)) for row in feature_group_rows))
    match_rate = preferred.get("annotation_match_rate", "")
    if not match_rate:
        match_rate = next((row.get("annotation_match_rate", "") for row in feature_group_rows if row.get("annotation_match_rate")), "")
    auroc = preferred.get("auroc", "")
    auprc = preferred.get("auprc", "")
    balanced_accuracy = preferred.get("balanced_accuracy", "")
    pvalue = preferred.get("empirical_permutation_pvalue", "")
    if status == "moderate_pd_internal_validation":
        reason = "phase20_safe_map_celltype_pd_signal_supports_moderate_internal_pd_extension"
    elif status == "preliminary_pd_internal_signal":
        reason = "phase20_safe_map_celltype_pd_signal_preliminary_because_permutation_support_is_not_significant"
    elif status == "technical_failure_annotation_join":
        reason = "phase20_safe_map_annotation_join_technical_failure"
    else:
        reason = "phase20_safe_map_celltype_pd_signal_weak_or_unavailable"
    return status, units, feature_count, match_rate, auroc, auprc, balanced_accuracy, pvalue, reason


def phase27_replication_summary(rows: list[dict[str, str]]) -> tuple[int, int, int, str, str]:
    supported = sum(row.get("evidence_label") == "replicated_statistically_supported" for row in rows)
    directional = sum(row.get("evidence_label") == "directionally_consistent_but_not_significant" for row in rows)
    weak = sum(row.get("evidence_label") in {"weak_or_no_replication", "opposite_direction"} for row in rows)
    upgrade = "true" if supported else "false"
    if supported:
        reason = "Phase 27 contains statistically supported GSE184950 axis replication."
    elif directional:
        reason = "Phase 27 contains directionally consistent GSE184950 signals, but p/FDR support is insufficient for replication claims."
    elif rows:
        reason = "Phase 27 GSE184950 clean replication is weak or unsupported."
    else:
        reason = "Phase 27 GSE184950 clean replication table is unavailable."
    return supported, directional, weak, upgrade, reason


def append_phase20_pd_claim(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> None:
    status, units, feature_count, match_rate, auroc, auprc, balanced_accuracy, pvalue, reason = phase20_pd_status(
        read_tsv(args.phase20_pd_metrics),
        read_tsv(args.phase20_feature_groups),
    )
    if status == "unavailable":
        return
    if any(row.get("task") == "parkinsons_vs_control" and row.get("evidence_layer") == "gse243639_phase20_safe_map_celltype" for row in claim_rows):
        return
    if status == "moderate_pd_internal_validation":
        strength = "moderate_internal"
        allowed = "NeuroFate shows moderate sample-level internal PD evidence in GSE243639 after safe-map cell-type-aware repair."
        risk = "medium"
    elif status == "preliminary_pd_internal_signal":
        strength = "exploratory_internal"
        allowed = (
            "NeuroFate shows a preliminary sample-level cell-type-aware PD internal signal in GSE243639; "
            "the result improves over global features but requires independent replication."
        )
        risk = "medium"
    elif status == "technical_failure_annotation_join":
        strength = "failed_or_unstable"
        allowed = "Phase 20 did not meet annotation-linkage requirements; no cell-type-aware PD claim is supported."
        risk = "high"
    else:
        strength = "failed_or_unstable"
        allowed = "Phase 20 PD evidence is weak or unavailable."
        risk = "high"
    claim_rows.append(
        {
            "claim_id": f"claim_{len(claim_rows) + 1:03d}",
            "task": "parkinsons_vs_control",
            "model": "logistic_regression",
            "evidence_layer": "gse243639_phase20_safe_map_celltype",
            "internal_auroc": auroc,
            "internal_auroc_sd": "",
            "internal_auprc": auprc,
            "balanced_accuracy": balanced_accuracy,
            "brier_score": "",
            "permutation_empirical_p": pvalue,
            "feature_ablation_support": "not_run_for_phase20",
            "external_validation_status": status,
            "external_sample_units": str(units),
            "leakage_status": leakage_status(read_tsv(args.leakage_audit)),
            "overclaiming_status": overclaiming_status(read_tsv(args.overclaiming_audit)),
            "claim_strength": strength,
            "allowed_claim_text": allowed,
            "disallowed_claim_text": "Do not claim clinical PD prediction, diagnostic utility, causality, or cross-disease validation.",
            "reviewer_risk": risk,
        }
    )


def append_phase27_replication_claim(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    rows = read_tsv(args.phase27_gse184950_replication)
    supported, directional, weak, upgrade, reason = phase27_replication_summary(rows)
    if rows:
        claim_rows.append(
            {column: "" for column in CLAIM_COLUMNS}
            | {
                "claim_id": "claim_gse184950_phase27_clean_replication",
                "task": "gse184950_pd_pdd_vs_control_axis_replication",
                "model": "rank_based_axis_association",
                "evidence_layer": "gse184950_phase27_clean_replication",
                "external_validation_status": "statistically_supported_replication" if supported else "directionally_consistent_or_weak_replication",
                "external_sample_units": "34",
                "claim_strength": "preliminary_external_feasibility",
                "allowed_claim_text": "GSE184950 provides a clean independent PD/PDD sample-level replication test; direction-only signals remain preliminary unless p/FDR support is present.",
                "disallowed_claim_text": "Do not claim validated AD/PD shared mechanism, clinical biomarker, diagnostic axis, or causal mechanism from weak FDR support.",
                "reviewer_risk": "medium" if supported else "high",
            }
        )
    return [
        {
            "cohort_id": "gse184950_pd_sn",
            "statistically_supported_axes": str(supported),
            "directionally_consistent_preliminary_axes": str(directional),
            "weak_or_no_replication_axes": str(weak),
            "claim_upgrade_allowed": upgrade,
            "reason": reason,
        }
    ]


def append_phase21_axis_claims(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    if read_tsv(args.phase22_axis_evidence):
        return [
            {
                "axis_id": "all_phase21_axes",
                "axis_classification": "superseded_by_phase22_endpoint_locked_evidence",
                "claim_strength": "exploratory_only",
                "affected_existing_claim": "true",
                "upgrade_applied": "false",
                "reason": "Phase 22 endpoint-locked evidence supersedes Phase 21 largest-effect-across-label comparisons for PNAS-facing claims.",
            }
        ]
    axis_rows = read_tsv(args.phase21_axis_claims)
    delta_rows: list[dict[str, str]] = []
    if not axis_rows:
        return delta_rows
    existing_axis_ids = {
        row.get("task", "").replace("axis:", "")
        for row in claim_rows
        if row.get("evidence_layer") == "phase21_axis_biology"
    }
    for axis_row in axis_rows:
        axis_id = axis_row.get("axis_id", "")
        if not axis_id or axis_id in existing_axis_ids:
            continue
        category = axis_row.get("axis_classification", "inconclusive_axis")
        strength = axis_row.get("claim_strength", "axis_level_insufficient_validation")
        if strength in {"axis_level_preliminary_evidence", "disease_specific_candidate"}:
            allowed = (
                f"NeuroFate supports {axis_id} as a candidate donor/sample-level biological axis; "
                "this is association evidence that requires independent replication."
            )
            risk = "medium"
        else:
            allowed = f"NeuroFate does not yet support a biological claim for {axis_id}."
            risk = "high"
        claim_rows.append(
            {
                "claim_id": f"claim_{len(claim_rows) + 1:03d}",
                "task": f"axis:{axis_id}",
                "model": "axis_score_association",
                "evidence_layer": "phase21_axis_biology",
                "internal_auroc": "",
                "internal_auroc_sd": "",
                "internal_auprc": "",
                "balanced_accuracy": "",
                "brier_score": "",
                "permutation_empirical_p": "",
                "feature_ablation_support": "random_axis_controls_required",
                "external_validation_status": category,
                "external_sample_units": "",
                "leakage_status": leakage_status(read_tsv(args.leakage_audit)),
                "overclaiming_status": overclaiming_status(read_tsv(args.overclaiming_audit)),
                "claim_strength": strength,
                "allowed_claim_text": allowed,
                "disallowed_claim_text": (
                    "Do not claim a causal axis, proven disease mechanism, clinical biomarker, "
                    "definitive shared mechanism, or validation across diseases."
                ),
                "reviewer_risk": risk,
            }
        )
        delta_rows.append(
            {
                "axis_id": axis_id,
                "axis_classification": category,
                "claim_strength": strength,
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": "Phase 21 adds biological axis claim rows without upgrading clinical or predictive claims.",
            }
        )
    return delta_rows


def append_phase22_axis_claims(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    evidence_rows = read_tsv(args.phase22_axis_evidence)
    delta_rows: list[dict[str, str]] = []
    if not evidence_rows:
        return delta_rows
    for evidence in evidence_rows:
        axis_id = evidence.get("axis_id", "")
        endpoint_id = evidence.get("endpoint_id", "")
        if not axis_id or not endpoint_id:
            continue
        axis_claim = evidence.get("axis_claim_class", "exploratory_or_inconclusive_endpoint_locked")
        if "candidate" in axis_claim:
            strength = "axis_level_preliminary_evidence"
            risk = "medium"
        else:
            strength = "axis_level_insufficient_validation"
            risk = "high"
        claim_rows.append(
            {
                "claim_id": f"claim_{len(claim_rows) + 1:03d}",
                "task": f"endpoint_axis:{axis_id}:{endpoint_id}",
                "model": "endpoint_locked_axis_association",
                "evidence_layer": "phase22_endpoint_locked_axis_biology",
                "internal_auroc": "",
                "internal_auroc_sd": "",
                "internal_auprc": "",
                "balanced_accuracy": "",
                "brier_score": "",
                "permutation_empirical_p": evidence.get("empirical_pvalue", ""),
                "feature_ablation_support": "matched_endpoint_random_axis_controls",
                "external_validation_status": axis_claim,
                "external_sample_units": evidence.get("n", ""),
                "leakage_status": leakage_status(read_tsv(args.leakage_audit)),
                "overclaiming_status": overclaiming_status(read_tsv(args.overclaiming_audit)),
                "claim_strength": strength,
                "allowed_claim_text": evidence.get("allowed_claim", ""),
                "disallowed_claim_text": evidence.get("disallowed_claim", ""),
                "reviewer_risk": risk,
            }
        )
        delta_rows.append(
            {
                "axis_id": axis_id,
                "endpoint_id": endpoint_id,
                "axis_claim_class": axis_claim,
                "affected_phase21_claim": "true",
                "upgrade_applied": "false",
                "reason": "Phase 22 replaces exploratory Phase 21 axis claims with endpoint-locked evidence; strong claims are not allowed yet.",
            }
        )
    return delta_rows


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


def apply_phase15_delta(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    phase15_rows = read_tsv(args.phase15_metrics)
    phase16_rows = read_tsv(args.phase16_pd_metrics)
    phase15_status, phase15_units, reason = phase15_external_status([*phase15_rows, *phase16_rows])
    delta_rows: list[dict[str, str]] = []
    for row in claim_rows:
        old_strength = row["claim_strength"]
        new_strength = old_strength
        if phase15_status == "reliable_external_validation" and old_strength in {"strong_internal", "moderate_internal"}:
            new_strength = "strong_internal"
            row["external_validation_status"] = "reliable_external_validation"
            row["external_sample_units"] = str(phase15_units)
            row["claim_strength"] = new_strength
            row["allowed_claim_text"] = "NeuroFate has internal evidence plus reliable external validation for this task."
        elif phase15_status == "preliminary_external_feasibility":
            row["external_validation_status"] = "preliminary_external_feasibility"
            row["external_sample_units"] = str(max(phase15_units, int(to_float(row.get("external_sample_units"), 0))))
        elif phase15_status in {"moderate_pd_internal_validation", "preliminary_cross_disease_feature_transfer"}:
            row["external_validation_status"] = phase15_status
            row["external_sample_units"] = str(max(phase15_units, int(to_float(row.get("external_sample_units"), 0))))
        delta_rows.append(
            {
                "task": row["task"],
                "old_claim_strength": old_strength,
                "new_claim_strength": row["claim_strength"],
                "phase15_external_status": phase15_status,
                "phase15_external_sample_units": str(phase15_units),
                "upgrade_applied": str(old_strength != row["claim_strength"]).lower(),
                "reason": reason,
            }
        )
    return delta_rows


def apply_phase17_delta(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    phase20_rows = read_tsv(args.phase20_pd_metrics)
    if phase20_rows:
        status, units, feature_count, match_rate, auroc, auprc, balanced_accuracy, pvalue, reason = phase20_pd_status(
            phase20_rows,
            read_tsv(args.phase20_feature_groups),
        )
        return [
            {
                "task": "parkinsons_vs_control",
                "phase17_status": "superseded_by_phase20",
                "phase17_sample_units": str(units),
                "phase17_auroc": auroc,
                "phase17_permutation_pvalue": pvalue,
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": f"Phase 17 technical result superseded by Phase 20 safe-map assessment: {reason}; feature_count={feature_count}; match_rate={match_rate}",
            }
        ]
    phase18_rows = read_tsv(args.phase18_pd_metrics)
    if phase18_rows:
        status, units, auroc, pvalue, reason = phase18_pd_status(phase18_rows)
        return [
            {
                "task": "parkinsons_vs_control",
                "phase17_status": "superseded_by_phase18",
                "phase17_sample_units": str(units),
                "phase17_auroc": auroc,
                "phase17_permutation_pvalue": pvalue,
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": f"Phase 17 weak result superseded by Phase 18 repaired assessment: {reason}",
            }
        ]
    status, units, auroc, pvalue, reason = phase17_pd_status(read_tsv(args.phase17_pd_metrics))
    delta_rows: list[dict[str, str]] = []
    affected = False
    for row in claim_rows:
        task = row.get("task", "")
        if task not in {"parkinsons_vs_control", "pd_external_validation"}:
            continue
        affected = True
        old_strength = row["claim_strength"]
        if status == "moderate_pd_internal_validation":
            row["external_validation_status"] = status
            row["external_sample_units"] = str(units)
            row["claim_strength"] = "moderate_internal"
            row["allowed_claim_text"] = "NeuroFate shows moderate sample-level internal PD evidence in GSE243639 with cell-type-aware features."
        elif status == "preliminary_pd_internal_signal":
            row["external_validation_status"] = status
            row["external_sample_units"] = str(units)
            row["claim_strength"] = "exploratory_internal"
            row["allowed_claim_text"] = "NeuroFate shows preliminary PD internal signal in GSE243639 that requires larger validation."
        delta_rows.append(
            {
                "task": task,
                "phase17_status": status,
                "phase17_sample_units": str(units),
                "phase17_auroc": auroc,
                "phase17_permutation_pvalue": pvalue,
                "affected_existing_claim": "true",
                "upgrade_applied": str(old_strength != row["claim_strength"]).lower(),
                "reason": reason,
            }
        )
    if not affected:
        delta_rows.append(
            {
                "task": "parkinsons_vs_control",
                "phase17_status": status,
                "phase17_sample_units": str(units),
                "phase17_auroc": auroc,
                "phase17_permutation_pvalue": pvalue,
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": f"{reason}; no existing AD claim was modified",
            }
        )
    return delta_rows


def apply_phase18_delta(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    phase20_rows = read_tsv(args.phase20_pd_metrics)
    if phase20_rows:
        status, units, feature_count, match_rate, auroc, auprc, balanced_accuracy, pvalue, reason = phase20_pd_status(
            phase20_rows,
            read_tsv(args.phase20_feature_groups),
        )
        return [
            {
                "task": "parkinsons_vs_control",
                "phase18_status": "superseded_by_phase20",
                "phase18_sample_units": str(units),
                "phase18_auroc": auroc,
                "phase18_permutation_pvalue": pvalue,
                "phase18_supersedes_phase17": "false",
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": f"Phase 18 annotation-join failure superseded by Phase 20 safe-map assessment: {reason}; feature_count={feature_count}; match_rate={match_rate}",
            }
        ]
    status, units, auroc, pvalue, reason = phase18_pd_status(read_tsv(args.phase18_pd_metrics))
    delta_rows: list[dict[str, str]] = []
    affected = False
    for row in claim_rows:
        task = row.get("task", "")
        if task not in {"parkinsons_vs_control", "pd_external_validation"}:
            continue
        affected = True
        old_strength = row["claim_strength"]
        if status == "moderate_pd_internal_validation":
            row["external_validation_status"] = status
            row["external_sample_units"] = str(units)
            row["claim_strength"] = "moderate_internal"
            row["allowed_claim_text"] = "NeuroFate shows moderate sample-level internal PD evidence in GSE243639 after repaired annotation matching."
        elif status == "preliminary_pd_internal_signal":
            row["external_validation_status"] = status
            row["external_sample_units"] = str(units)
            row["claim_strength"] = "exploratory_internal"
            row["allowed_claim_text"] = "NeuroFate shows preliminary repaired PD internal signal in GSE243639 that requires larger validation."
        elif status == "technical_failure_annotation_join":
            row["external_validation_status"] = status
            row["claim_strength"] = "failed_or_unstable"
            row["allowed_claim_text"] = "The repaired GSE243639 annotation join failed technical thresholds; no PD biological claim is supported."
        delta_rows.append(
            {
                "task": task,
                "phase18_status": status,
                "phase18_sample_units": str(units),
                "phase18_auroc": auroc,
                "phase18_permutation_pvalue": pvalue,
                "phase18_supersedes_phase17": str(bool(read_tsv(args.phase18_pd_metrics))).lower(),
                "affected_existing_claim": "true",
                "upgrade_applied": str(old_strength != row["claim_strength"]).lower(),
                "reason": reason,
            }
        )
    if not affected:
        delta_rows.append(
            {
                "task": "parkinsons_vs_control",
                "phase18_status": status,
                "phase18_sample_units": str(units),
                "phase18_auroc": auroc,
                "phase18_permutation_pvalue": pvalue,
                "phase18_supersedes_phase17": str(bool(read_tsv(args.phase18_pd_metrics))).lower(),
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": f"{reason}; no existing AD claim was modified",
            }
        )
    return delta_rows


def apply_phase20_delta(claim_rows: list[dict[str, str]], args: argparse.Namespace) -> list[dict[str, str]]:
    status, units, feature_count, match_rate, auroc, auprc, balanced_accuracy, pvalue, reason = phase20_pd_status(
        read_tsv(args.phase20_pd_metrics),
        read_tsv(args.phase20_feature_groups),
    )
    delta_rows: list[dict[str, str]] = []
    affected = False
    for row in claim_rows:
        if row.get("task") != "parkinsons_vs_control":
            continue
        affected = True
        old_strength = row["claim_strength"]
        if status == "moderate_pd_internal_validation":
            row["external_validation_status"] = status
            row["external_sample_units"] = str(units)
            row["claim_strength"] = "moderate_internal"
            row["allowed_claim_text"] = "NeuroFate shows moderate sample-level internal PD evidence in GSE243639 after safe-map cell-type-aware repair."
        elif status == "preliminary_pd_internal_signal":
            row["external_validation_status"] = status
            row["external_sample_units"] = str(units)
            row["claim_strength"] = "exploratory_internal"
            row["allowed_claim_text"] = "NeuroFate shows a preliminary sample-level cell-type-aware PD internal signal in GSE243639 that improves over global features but requires independent replication."
        elif status == "technical_failure_annotation_join":
            row["external_validation_status"] = status
            row["claim_strength"] = "failed_or_unstable"
            row["allowed_claim_text"] = "Phase 20 did not meet annotation-linkage requirements; no cell-type-aware PD claim is supported."
        delta_rows.append(
            {
                "task": "parkinsons_vs_control",
                "phase20_status": status,
                "phase20_sample_units": str(units),
                "phase20_feature_count": feature_count,
                "phase20_annotation_match_rate": match_rate,
                "phase20_auroc": auroc,
                "phase20_auprc": auprc,
                "phase20_balanced_accuracy": balanced_accuracy,
                "phase20_permutation_pvalue": pvalue,
                "phase20_supersedes_phase17_phase18": str(bool(read_tsv(args.phase20_pd_metrics))).lower(),
                "affected_existing_claim": "true",
                "upgrade_applied": str(old_strength != row["claim_strength"]).lower(),
                "reason": reason,
            }
        )
    if not affected:
        delta_rows.append(
            {
                "task": "parkinsons_vs_control",
                "phase20_status": status,
                "phase20_sample_units": str(units),
                "phase20_feature_count": feature_count,
                "phase20_annotation_match_rate": match_rate,
                "phase20_auroc": auroc,
                "phase20_auprc": auprc,
                "phase20_balanced_accuracy": balanced_accuracy,
                "phase20_permutation_pvalue": pvalue,
                "phase20_supersedes_phase17_phase18": str(bool(read_tsv(args.phase20_pd_metrics))).lower(),
                "affected_existing_claim": "false",
                "upgrade_applied": "false",
                "reason": f"{reason}; no existing AD claim was modified",
            }
        )
    return delta_rows


def build_claim_tables(
    args: argparse.Namespace,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
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
    append_phase20_pd_claim(claim_rows, args)
    phase27_delta_rows = append_phase27_replication_claim(claim_rows, args)
    phase22_delta_rows = append_phase22_axis_claims(claim_rows, args)
    phase21_delta_rows = append_phase21_axis_claims(claim_rows, args)
    append_phase32_crosscohort_claims(claim_rows, args)
    delta_rows = apply_phase15_delta(claim_rows, args)
    phase17_delta_rows = apply_phase17_delta(claim_rows, args)
    phase18_delta_rows = apply_phase18_delta(claim_rows, args)
    phase20_delta_rows = apply_phase20_delta(claim_rows, args)
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
                "reviewer_caveat": (
                    "External evidence remains preliminary for Mathys n=6; Phase 20 GSE243639 PD evidence is preliminary internal signal unless independently replicated."
                ),
            }
        )
    return claim_rows, best_rows, delta_rows, phase17_delta_rows, phase18_delta_rows, phase20_delta_rows, phase21_delta_rows, phase22_delta_rows, phase27_delta_rows


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
    parser.add_argument("--phase15-metrics", type=Path, default=Path("results/tables/phase15_multi_external_validation_metrics.tsv"))
    parser.add_argument("--phase16-pd-metrics", type=Path, default=Path("results/tables/phase16_gse243639_external_validation_metrics.tsv"))
    parser.add_argument("--phase17-pd-metrics", type=Path, default=Path("results/tables/phase17_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase18-pd-metrics", type=Path, default=Path("results/tables/phase18_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase20-pd-metrics", type=Path, default=Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase20-feature-groups", type=Path, default=Path("results/tables/phase20_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--phase21-axis-claims", type=Path, default=Path("results/tables/phase21_axis_claim_strength.tsv"))
    parser.add_argument("--phase22-axis-evidence", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv"))
    parser.add_argument("--phase27-gse184950-replication", type=Path, default=Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"))
    parser.add_argument("--phase32-crosscohort-evidence", type=Path, default=Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv"))
    parser.add_argument("--leakage-audit", type=Path, default=Path("results/reports/feature_leakage_audit.tsv"))
    parser.add_argument("--overclaiming-audit", type=Path, default=Path("results/reports/no_overclaiming_audit.tsv"))
    parser.add_argument("--claim-output", type=Path, default=Path("results/reports/claim_strength_table.tsv"))
    parser.add_argument("--best-output", type=Path, default=Path("results/reports/best_supported_claims.tsv"))
    parser.add_argument("--phase15-delta-output", type=Path, default=Path("results/reports/phase15_claim_strength_delta.tsv"))
    parser.add_argument("--phase17-delta-output", type=Path, default=Path("results/reports/phase17_claim_strength_delta.tsv"))
    parser.add_argument("--phase18-delta-output", type=Path, default=Path("results/reports/phase18_claim_strength_delta.tsv"))
    parser.add_argument("--phase20-delta-output", type=Path, default=Path("results/reports/phase20_claim_strength_delta.tsv"))
    parser.add_argument("--phase21-delta-output", type=Path, default=Path("results/reports/phase21_axis_claim_strength_delta.tsv"))
    parser.add_argument("--phase22-delta-output", type=Path, default=Path("results/reports/phase22_endpoint_locked_axis_claim_strength_delta.tsv"))
    parser.add_argument("--phase27-delta-output", type=Path, default=Path("results/reports/phase27_gse184950_claim_strength_delta.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    claim_rows, best_rows, delta_rows, phase17_delta_rows, phase18_delta_rows, phase20_delta_rows, phase21_delta_rows, phase22_delta_rows, phase27_delta_rows = build_claim_tables(args)
    write_tsv(args.claim_output, claim_rows, CLAIM_COLUMNS)
    write_tsv(args.best_output, best_rows, BEST_COLUMNS)
    write_tsv(
        args.phase15_delta_output,
        delta_rows,
        [
            "task",
            "old_claim_strength",
            "new_claim_strength",
            "phase15_external_status",
            "phase15_external_sample_units",
            "upgrade_applied",
            "reason",
        ],
    )
    write_tsv(args.phase17_delta_output, phase17_delta_rows, PHASE17_DELTA_COLUMNS)
    write_tsv(args.phase18_delta_output, phase18_delta_rows, PHASE18_DELTA_COLUMNS)
    write_tsv(args.phase20_delta_output, phase20_delta_rows, PHASE20_DELTA_COLUMNS)
    write_tsv(args.phase21_delta_output, phase21_delta_rows, PHASE21_DELTA_COLUMNS)
    write_tsv(args.phase22_delta_output, phase22_delta_rows, PHASE22_DELTA_COLUMNS)
    write_tsv(args.phase27_delta_output, phase27_delta_rows, PHASE27_DELTA_COLUMNS)
    print(f"Wrote {args.claim_output}")
    print(f"Wrote {args.best_output}")
    print(f"Wrote {args.phase15_delta_output}")
    print(f"Wrote {args.phase17_delta_output}")
    print(f"Wrote {args.phase18_delta_output}")
    print(f"Wrote {args.phase20_delta_output}")
    print(f"Wrote {args.phase21_delta_output}")
    print(f"Wrote {args.phase22_delta_output}")
    print(f"Wrote {args.phase27_delta_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
