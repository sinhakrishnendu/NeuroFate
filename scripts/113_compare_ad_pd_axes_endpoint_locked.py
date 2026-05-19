#!/usr/bin/env python3
"""Compare AD and PD NeuroFate axes only through locked primary/secondary endpoints."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


PRIMARY_AD = "sea_ad_cognitive_dementia"
SECONDARY_AD = "sea_ad_ad_pathology_ordinal"
PRIMARY_PD = "gse243639_pd_diagnosis"


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


def sign(value: float) -> int:
    if math.isnan(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def stat_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    return {(row.get("endpoint_id", ""), row.get("axis_id", "")): row for row in rows}


def classify_axis(ad_effect: float, pd_effect: float, ad_n: int, pd_n: int) -> str:
    if ad_n < 20 or pd_n < 20 or math.isnan(ad_effect) or math.isnan(pd_effect):
        return "insufficient_endpoint_support"
    ad_abs = abs(ad_effect)
    pd_abs = abs(pd_effect)
    if sign(ad_effect) == sign(pd_effect) and sign(ad_effect) != 0 and ad_abs >= 0.20 and pd_abs >= 0.20:
        return "shared_axis_candidate_endpoint_locked"
    if ad_abs >= 0.30 and ad_abs >= 1.5 * max(pd_abs, 0.01):
        return "ad_enriched_axis_candidate_endpoint_locked"
    if pd_abs >= 0.30 and pd_abs >= 1.5 * max(ad_abs, 0.01):
        return "pd_enriched_axis_candidate_endpoint_locked"
    return "axis_inconclusive_endpoint_locked"


def claim_strength(axis_class: str) -> str:
    if axis_class == "shared_axis_candidate_endpoint_locked":
        return "endpoint_locked_candidate_shared_axis"
    if axis_class in {"ad_enriched_axis_candidate_endpoint_locked", "pd_enriched_axis_candidate_endpoint_locked"}:
        return "endpoint_locked_disease_enriched_candidate"
    return "endpoint_locked_insufficient_or_inconclusive"


def comparison_row(
    lookup: dict[tuple[str, str], dict[str, str]],
    axis: str,
    ad_endpoint: str,
    pd_endpoint: str,
    comparison_type: str,
) -> dict[str, str]:
    ad = lookup.get((ad_endpoint, axis), {})
    pd = lookup.get((pd_endpoint, axis), {})
    ad_effect = to_float(ad.get("effect_size"))
    pd_effect = to_float(pd.get("effect_size"))
    ad_n = int(to_float(ad.get("n"), 0))
    pd_n = int(to_float(pd.get("n"), 0))
    category = classify_axis(ad_effect, pd_effect, ad_n, pd_n)
    return {
        "comparison_type": comparison_type,
        "axis_id": axis,
        "ad_endpoint_id": ad_endpoint,
        "pd_endpoint_id": pd_endpoint,
        "ad_effect_size": "" if math.isnan(ad_effect) else f"{ad_effect:.8g}",
        "pd_effect_size": "" if math.isnan(pd_effect) else f"{pd_effect:.8g}",
        "same_direction": str(sign(ad_effect) == sign(pd_effect) and sign(ad_effect) != 0).lower(),
        "ad_n": str(ad_n),
        "pd_n": str(pd_n),
        "ad_pvalue": ad.get("pvalue", ""),
        "pd_pvalue": pd.get("pvalue", ""),
        "ad_fdr": ad.get("fdr", ""),
        "pd_fdr": pd.get("fdr", ""),
        "axis_claim_class": category,
    }


def build_comparisons(stats: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    lookup = stat_lookup(stats)
    axes = sorted({row.get("axis_id", "") for row in stats if row.get("axis_id")})
    similarity_rows: list[dict[str, str]] = []
    for axis in axes:
        similarity_rows.append(comparison_row(lookup, axis, PRIMARY_AD, PRIMARY_PD, "primary_ad_cognition_vs_primary_pd_diagnosis"))
        similarity_rows.append(comparison_row(lookup, axis, SECONDARY_AD, PRIMARY_PD, "secondary_ad_pathology_vs_primary_pd_diagnosis"))
    primary_rows = [row for row in similarity_rows if row["comparison_type"].startswith("primary")]
    shared_rows = [
        {
            "axis_id": row["axis_id"],
            "comparison_type": row["comparison_type"],
            "classification": row["axis_claim_class"],
            "interpretation": (
                "endpoint-locked candidate only; requires replication"
                if "candidate" in row["axis_claim_class"]
                else "insufficient or inconclusive endpoint-locked evidence"
            ),
        }
        for row in similarity_rows
    ]
    claim_rows = [
        {
            "axis_id": row["axis_id"],
            "comparison_type": row["comparison_type"],
            "axis_claim_class": row["axis_claim_class"],
            "claim_strength": claim_strength(row["axis_claim_class"]),
            "allowed_claim_text": (
                f"{row['axis_id']} is an endpoint-locked candidate axis under {row['comparison_type']}."
                if "candidate" in row["axis_claim_class"]
                else f"{row['axis_id']} is inconclusive under endpoint-locked comparison."
            ),
            "disallowed_claim_text": "Do not claim causal axis, proven mechanism, clinical biomarker, diagnostic utility, validated across diseases, or definitive shared mechanism.",
            "reviewer_risk": "medium" if "candidate" in row["axis_claim_class"] else "high",
        }
        for row in primary_rows
    ]
    return similarity_rows, shared_rows, claim_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare AD/PD axes with locked endpoints only.")
    parser.add_argument("--endpoint-stats", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_statistics.tsv"))
    parser.add_argument("--endpoint-registry", type=Path, default=Path("metadata/neurofate_axis_endpoint_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/113_compare_ad_pd_axes_endpoint_locked.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    # Registry is read to make the endpoint dependency explicit and auditable.
    read_tsv(args.endpoint_registry)
    similarity_rows, shared_rows, claim_rows = build_comparisons(read_tsv(args.endpoint_stats))
    similarity_columns = [
        "comparison_type",
        "axis_id",
        "ad_endpoint_id",
        "pd_endpoint_id",
        "ad_effect_size",
        "pd_effect_size",
        "same_direction",
        "ad_n",
        "pd_n",
        "ad_pvalue",
        "pd_pvalue",
        "ad_fdr",
        "pd_fdr",
        "axis_claim_class",
    ]
    write_tsv(args.outdir / "phase22_endpoint_locked_ad_pd_axis_similarity.tsv", similarity_rows, similarity_columns)
    write_tsv(
        args.outdir / "phase22_endpoint_locked_shared_vs_specific_axes.tsv",
        shared_rows,
        ["axis_id", "comparison_type", "classification", "interpretation"],
    )
    write_tsv(
        args.outdir / "phase22_endpoint_locked_axis_claim_strength.tsv",
        claim_rows,
        ["axis_id", "comparison_type", "axis_claim_class", "claim_strength", "allowed_claim_text", "disallowed_claim_text", "reviewer_risk"],
    )
    logging.info("Endpoint-locked AD/PD comparisons rows=%d", len(similarity_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
