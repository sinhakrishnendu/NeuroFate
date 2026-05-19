#!/usr/bin/env python3
"""Compare AD and PD NeuroFate axis association patterns conservatively."""

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


def to_float(value: str | None) -> float:
    try:
        return float(value) if value not in (None, "") else math.nan
    except ValueError:
        return math.nan


def best_axis_effect(rows: list[dict[str, str]], cohort_contains: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        if cohort_contains not in row.get("cohort", ""):
            continue
        axis = row.get("axis_id", "")
        effect = abs(to_float(row.get("effect_size")))
        current = out.get(axis)
        if current is None or effect > abs(to_float(current.get("effect_size"))):
            out[axis] = row
    return out


def sign(value: float) -> int:
    if math.isnan(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def classify_axis(ad_effect: float, pd_effect: float, ad_n: int, pd_n: int) -> str:
    if ad_n < 4 or pd_n < 4 or math.isnan(ad_effect) or math.isnan(pd_effect):
        return "insufficient_coverage"
    ad_abs = abs(ad_effect)
    pd_abs = abs(pd_effect)
    if ad_abs >= 0.20 and pd_abs >= 0.20 and sign(ad_effect) == sign(pd_effect):
        return "shared_ad_pd_candidate"
    if ad_abs >= 0.30 and ad_abs >= 1.5 * max(pd_abs, 0.01):
        return "ad_enriched_axis"
    if pd_abs >= 0.30 and pd_abs >= 1.5 * max(ad_abs, 0.01):
        return "pd_enriched_axis"
    return "inconclusive_axis"


def claim_strength(axis_class: str) -> str:
    if axis_class == "shared_ad_pd_candidate":
        return "axis_level_preliminary_evidence"
    if axis_class in {"ad_enriched_axis", "pd_enriched_axis"}:
        return "disease_specific_candidate"
    return "axis_level_insufficient_validation"


def build_rows(stats: list[dict[str, str]]) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    ad = best_axis_effect(stats, "sea_ad")
    pd = best_axis_effect(stats, "gse243639")
    axes = sorted(set(ad) | set(pd))
    similarity_rows: list[dict[str, str]] = []
    shared_rows: list[dict[str, str]] = []
    claim_rows: list[dict[str, str]] = []
    for axis in axes:
        ad_row = ad.get(axis, {})
        pd_row = pd.get(axis, {})
        ad_effect = to_float(ad_row.get("effect_size"))
        pd_effect = to_float(pd_row.get("effect_size"))
        ad_n = int(to_float(ad_row.get("n")) if not math.isnan(to_float(ad_row.get("n"))) else 0)
        pd_n = int(to_float(pd_row.get("n")) if not math.isnan(to_float(pd_row.get("n"))) else 0)
        category = classify_axis(ad_effect, pd_effect, ad_n, pd_n)
        similarity_rows.append(
            {
                "axis_id": axis,
                "ad_effect_size": "" if math.isnan(ad_effect) else f"{ad_effect:.8g}",
                "pd_effect_size": "" if math.isnan(pd_effect) else f"{pd_effect:.8g}",
                "same_direction": str(sign(ad_effect) == sign(pd_effect) and sign(ad_effect) != 0).lower(),
                "ad_n": str(ad_n),
                "pd_n": str(pd_n),
                "axis_classification": category,
            }
        )
        shared_rows.append(
            {
                "axis_id": axis,
                "classification": category,
                "interpretation": "candidate only; requires independent replication" if "candidate" in category or "enriched" in category else "insufficient or inconclusive evidence",
            }
        )
        claim_rows.append(
            {
                "axis_id": axis,
                "axis_classification": category,
                "claim_strength": claim_strength(category),
                "allowed_claim_text": f"NeuroFate identifies {axis} as a {category} at donor/sample level." if category not in {"insufficient_coverage", "inconclusive_axis"} else f"Evidence is insufficient to classify {axis}.",
                "disallowed_claim_text": "Do not claim causal axis, disease mechanism proven, clinical biomarker, definitive shared mechanism, or validated across diseases.",
                "reviewer_risk": "medium" if "candidate" in category or "enriched" in category else "high",
            }
        )
    return similarity_rows, shared_rows, claim_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare AD and PD NeuroFate axis patterns.")
    parser.add_argument("--axis-stats", type=Path, default=Path("results/tables/phase21_axis_association_statistics.tsv"))
    parser.add_argument("--axis-effects", type=Path, default=Path("results/tables/phase21_axis_effect_sizes.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/108_compare_ad_pd_axis_patterns.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    stats = read_tsv(args.axis_effects) or read_tsv(args.axis_stats)
    similarity_rows, shared_rows, claim_rows = build_rows(stats)
    write_tsv(
        args.outdir / "phase21_ad_pd_axis_similarity.tsv",
        similarity_rows,
        ["axis_id", "ad_effect_size", "pd_effect_size", "same_direction", "ad_n", "pd_n", "axis_classification"],
    )
    write_tsv(
        args.outdir / "phase21_shared_vs_disease_specific_axes.tsv",
        shared_rows,
        ["axis_id", "classification", "interpretation"],
    )
    write_tsv(
        args.outdir / "phase21_axis_claim_strength.tsv",
        claim_rows,
        ["axis_id", "axis_classification", "claim_strength", "allowed_claim_text", "disallowed_claim_text", "reviewer_risk"],
    )
    logging.info("Compared AD/PD axis patterns for axes=%d", len(claim_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
