#!/usr/bin/env python3
"""Build conservative endpoint-locked NeuroFate-Axis evidence tables and report."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def comparison_lookup(rows: list[dict[str, str]]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in rows:
        if row.get("comparison_type", "").startswith("primary"):
            lookup[row.get("axis_id", "")] = row.get("axis_claim_class", "axis_inconclusive_endpoint_locked")
    return lookup


def empirical_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    return {(row.get("endpoint_id", ""), row.get("axis_id", "")): row for row in rows}


def classify_claim(stat: dict[str, str], axis_class: str, empirical: dict[str, str]) -> tuple[str, str, str, str, str]:
    fdr = to_float(stat.get("fdr"), 1.0)
    empirical_p = to_float(empirical.get("empirical_pvalue"), 1.0)
    if axis_class == "shared_axis_candidate_endpoint_locked":
        base = "candidate_shared_axis_endpoint_locked"
    elif "ad_enriched" in axis_class or "pd_enriched" in axis_class:
        base = "candidate_disease_enriched_axis_endpoint_locked"
    elif axis_class == "insufficient_endpoint_support":
        base = "insufficient_endpoint_support"
    else:
        base = "exploratory_or_inconclusive_endpoint_locked"
    if empirical_p >= 0.05 and base.startswith("candidate"):
        base = "preliminary_" + base
    if fdr > 0.10 and empirical_p > 0.10:
        base = "exploratory_or_inconclusive_endpoint_locked"
    if base.startswith("candidate") or base.startswith("preliminary_candidate"):
        allowed = (
            f"{stat['axis_id']} may be described as a {base.replace('_', ' ')} "
            f"for endpoint {stat['endpoint_id']}, pending replication."
        )
        risk = "medium"
    else:
        allowed = f"{stat['axis_id']} is not sufficiently supported for a PNAS-facing biological claim at endpoint {stat['endpoint_id']}."
        risk = "high"
    disallowed = "Do not claim causal axis, proven mechanism, clinical biomarker, diagnostic utility, validated across diseases, or definitive shared mechanism."
    coverage = "ok" if empirical.get("random_feature_count", "0") not in {"", "0"} else "insufficient_feature_coverage"
    return base, allowed, disallowed, risk, coverage


def build_rows(stats: list[dict[str, str]], comparisons: list[dict[str, str]], empirical_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    axis_classes = comparison_lookup(comparisons)
    empirical = empirical_lookup(empirical_rows)
    rows: list[dict[str, str]] = []
    for stat in stats:
        emp = empirical.get((stat.get("endpoint_id", ""), stat.get("axis_id", "")), {})
        axis_class = axis_classes.get(stat.get("axis_id", ""), "axis_inconclusive_endpoint_locked")
        claim, allowed, disallowed, risk, coverage = classify_claim(stat, axis_class, emp)
        rows.append(
            {
                "axis_id": stat.get("axis_id", ""),
                "endpoint_id": stat.get("endpoint_id", ""),
                "cohort": stat.get("cohort", ""),
                "endpoint_role": stat.get("endpoint_role", ""),
                "effect_size": stat.get("effect_size", ""),
                "standardized_mean_difference": stat.get("standardized_mean_difference", ""),
                "pvalue": stat.get("pvalue", ""),
                "fdr": stat.get("fdr", ""),
                "empirical_pvalue": emp.get("empirical_pvalue", ""),
                "n": stat.get("n", ""),
                "positive_n_or_valid_n": stat.get("positive_n") or stat.get("valid_n", ""),
                "negative_n": stat.get("negative_n", ""),
                "feature_count": emp.get("random_feature_count", ""),
                "coverage_status": coverage,
                "axis_claim_class": claim,
                "allowed_claim": allowed,
                "disallowed_claim": disallowed,
                "reviewer_risk": risk,
            }
        )
    return rows


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    candidates = [row for row in rows if row["axis_claim_class"].startswith(("candidate", "preliminary_candidate"))]
    text = [
        "# Phase 22 Endpoint-Locked Axis Claims",
        "",
        "Phase 22 supersedes the Phase 21 largest-effect-across-label comparison for PNAS-facing biological claims.",
        "Strong claims are not allowed yet. All claims are endpoint-locked and remain candidate or preliminary; no clinical, diagnostic, causal, or validated-across-diseases claim is supported.",
        "",
        "## Candidate Rows",
    ]
    if candidates:
        for row in candidates[:20]:
            text.append(f"- `{row['axis_id']}` at `{row['endpoint_id']}`: {row['axis_claim_class']} (effect={row['effect_size']}, empirical p={row['empirical_pvalue'] or 'NA'}).")
    else:
        text.append("- No endpoint-locked candidate axes are currently supported.")
    text.extend(
        [
            "",
            "## Disallowed Language",
            "",
            "Do not claim causal axes, proven mechanisms, clinical biomarkers, diagnostic utility, definitive shared mechanisms, or validation across diseases.",
        ]
    )
    path.write_text("\n".join(text) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build endpoint-locked NeuroFate axis evidence table.")
    parser.add_argument("--endpoint-stats", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_statistics.tsv"))
    parser.add_argument("--axis-comparison", type=Path, default=Path("results/tables/phase22_endpoint_locked_ad_pd_axis_similarity.tsv"))
    parser.add_argument("--empirical", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_empirical_pvalues.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--report-dir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/115_build_endpoint_locked_axis_evidence_table.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = build_rows(read_tsv(args.endpoint_stats), read_tsv(args.axis_comparison), read_tsv(args.empirical))
    columns = [
        "axis_id",
        "endpoint_id",
        "cohort",
        "endpoint_role",
        "effect_size",
        "standardized_mean_difference",
        "pvalue",
        "fdr",
        "empirical_pvalue",
        "n",
        "positive_n_or_valid_n",
        "negative_n",
        "feature_count",
        "coverage_status",
        "axis_claim_class",
        "allowed_claim",
        "disallowed_claim",
        "reviewer_risk",
    ]
    write_tsv(args.outdir / "phase22_endpoint_locked_axis_evidence_table.tsv", rows, columns)
    write_report(args.report_dir / "phase22_endpoint_locked_axis_claims.md", rows)
    logging.info("Endpoint-locked evidence rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
