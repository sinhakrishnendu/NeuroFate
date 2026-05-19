#!/usr/bin/env python3
"""Build conservative cross-cohort NeuroFate-Axis evidence summaries."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


SUMMARY_COLUMNS = [
    "axis_id",
    "sea_ad_effect",
    "sea_ad_p",
    "sea_ad_fdr",
    "gse174367_effect",
    "gse174367_p",
    "gse174367_fdr",
    "gse174367_direction_consistency",
    "gse184950_effect",
    "gse184950_p",
    "gse184950_fdr",
    "gse184950_direction_consistency",
    "gse243639_axis_effect",
    "gse243639_axis_p",
    "gse243639_axis_fdr",
    "gse243639_axis_empirical_p",
    "gse243639_axis_direction_consistency",
    "phase34_best_pd_cohort",
    "phase34_best_pd_effect",
    "phase34_best_pd_p",
    "phase34_best_pd_fdr",
    "phase34_best_pd_label",
    "phase34_best_pd_direction_consistency",
    "phase34_best_pd_divergence_status",
    "pd_gse243639_support",
    "crosscohort_evidence_class",
    "safe_claim",
    "unsafe_claim",
    "next_validation_needed",
]


CLASS_RANK = {
    "strong_ad_axis_with_nominal_external_replication": 0,
    "strong_ad_axis_without_external_replication": 1,
    "preliminary_shared_ad_pd_axis_candidate": 2,
    "pd_preliminary_only": 3,
    "pd_divergent_axis_candidate": 4,
    "inconclusive_axis": 5,
    "insufficient_data": 6,
}


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
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


def by_axis(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row.get("axis_id", ""): row for row in rows if row.get("axis_id")}


def sea_ad_primary(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return by_axis(
        [
            row
            for row in rows
            if row.get("cohort") == "sea_ad" and row.get("endpoint_id") == "sea_ad_cognitive_dementia"
        ]
    )


def pd_gse243639_primary(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return by_axis(
        [
            row
            for row in rows
            if row.get("cohort") == "gse243639_pd_snpc" and row.get("endpoint_id") == "gse243639_pd_diagnosis"
        ]
    )


def pd_gse243639_status(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "unavailable"
    preferred = next(
        (
            row
            for row in rows
            if row.get("model") == "logistic_regression"
            and row.get("validation_mode") == "repeated_stratified_split"
        ),
        rows[0],
    )
    return preferred.get("reliability_flag", "unavailable")


def classify_axis(
    axis: str,
    sea: dict[str, str] | None,
    ad_rep: dict[str, str] | None,
    pd_rep: dict[str, str] | None,
    pd_axis: dict[str, str] | None,
    phase34_pd: dict[str, str] | None,
    pd_support: str,
) -> tuple[str, str, str, str]:
    sea_effect = to_float(sea.get("effect_size") if sea else None)
    sea_p = to_float(sea.get("pvalue") if sea else None, 1.0)
    sea_fdr = to_float(sea.get("fdr") if sea else None, 1.0)
    ad_effect = to_float(ad_rep.get("effect_size") if ad_rep else None)
    ad_p = to_float(ad_rep.get("pvalue") if ad_rep else None, 1.0)
    ad_fdr = to_float(ad_rep.get("fdr") if ad_rep else None, 1.0)
    pd_effect = to_float(pd_rep.get("effect_size") if pd_rep else None)
    pd_p = to_float(pd_rep.get("pvalue") if pd_rep else None, 1.0)
    pd_fdr = to_float(pd_rep.get("fdr") if pd_rep else None, 1.0)
    pd_axis_effect = to_float(pd_axis.get("effect_size") if pd_axis else None)
    pd_axis_p = to_float(pd_axis.get("pvalue") if pd_axis else None, 1.0)
    pd_axis_fdr = to_float(pd_axis.get("fdr") if pd_axis else None, 1.0)
    pd_axis_empirical = to_float(pd_axis.get("empirical_pvalue") if pd_axis else None, 1.0)
    phase34_effect = to_float(phase34_pd.get("effect_size") if phase34_pd else None)
    phase34_p = to_float(phase34_pd.get("pvalue") if phase34_pd else None, 1.0)
    phase34_fdr = to_float(phase34_pd.get("fdr") if phase34_pd else None, 1.0)
    phase34_label = phase34_pd.get("evidence_label", "") if phase34_pd else ""
    ad_consistent = sign(sea_effect) == sign(ad_effect) and sign(sea_effect) != 0
    pd_consistent = bool(pd_rep) and pd_rep.get("directional_consistency") == "consistent"
    pd_axis_consistent = sign(sea_effect) == sign(pd_axis_effect) and sign(sea_effect) != 0
    phase34_consistent = sign(sea_effect) == sign(phase34_effect) and sign(sea_effect) != 0
    sea_supported = sea_p < 0.05 or sea_fdr < 0.10
    ad_nominal = ad_consistent and ad_p < 0.05
    ad_fdr_supported = ad_consistent and ad_fdr < 0.10
    pd_stat_supported = pd_consistent and (pd_p < 0.05 or pd_fdr < 0.10)
    pd_axis_preliminary = pd_axis_consistent and (pd_axis_p < 0.25 or pd_axis_empirical < 0.05)
    pd_axis_supported = pd_axis_consistent and (pd_axis_p < 0.05 or pd_axis_fdr < 0.10)
    phase34_pd_supported = phase34_consistent and (phase34_p < 0.05 or phase34_fdr < 0.10)
    phase34_pd_divergent = phase34_label == "opposite_direction" and (phase34_p < 0.05 or phase34_fdr < 0.10)
    if sea_supported and ad_nominal:
        klass = "strong_ad_axis_with_nominal_external_replication"
        safe = (
            f"{axis} has endpoint-locked SEA-AD support and nominal independent AD replication in GSE174367; "
            "describe it as a directionally consistent AD replication candidate, not a definitive mechanism."
        )
        next_step = "Seek FDR-robust AD replication and stronger independent PD axis replication."
    elif sea_supported:
        klass = "strong_ad_axis_without_external_replication"
        safe = f"{axis} is supported internally in SEA-AD but lacks nominal independent AD replication."
        next_step = "Replicate the axis in an independent AD cohort."
    elif pd_support == "preliminary_pd_internal_signal" or pd_stat_supported:
        klass = "pd_preliminary_only"
        safe = f"{axis} has at most preliminary PD-context support; no shared AD/PD mechanism claim is supported."
        next_step = "Run axis-level PD replication with stronger statistical support."
    elif sea or ad_rep or pd_rep:
        klass = "inconclusive_axis"
        safe = f"{axis} remains inconclusive across the current endpoint-locked evidence."
        next_step = "Improve axis coverage and independent replication."
    else:
        klass = "insufficient_data"
        safe = f"{axis} has insufficient cross-cohort data."
        next_step = "Generate endpoint-locked association statistics."
    if klass == "strong_ad_axis_with_nominal_external_replication" and pd_axis_preliminary and not pd_axis_supported:
        safe = (
            f"{axis} has endpoint-locked SEA-AD support, nominal independent AD replication in GSE174367, "
            "and preliminary same-direction GSE243639 PD axis convergence; describe this as AD-replicated "
            "with exploratory PD convergence, not as a shared mechanism."
        )
        next_step = "Seek statistically supported independent PD axis replication."
    if phase34_pd_divergent:
        klass = "pd_divergent_axis_candidate"
        safe = (
            f"{axis} has statistically supported opposite-direction PD evidence in a Phase 34/37 cohort. "
            "Describe this as a candidate PD-divergent axis requiring direction/probe audit and independent confirmation, "
            "not as shared AD/PD replication."
        )
        next_step = "Validate the opposite-direction PD signal in another PD cohort and audit platform/probe behavior."
    if klass == "strong_ad_axis_with_nominal_external_replication" and phase34_pd_supported:
        klass = "preliminary_shared_ad_pd_axis_candidate"
        safe = (
            f"{axis} has endpoint-locked SEA-AD support, nominal independent AD replication, "
            "and statistically supported Phase 34 PD replication. Describe it as a replicated candidate axis "
            "requiring additional disease cohorts, not as final shared-mechanism evidence."
        )
        next_step = "Add another independent PD cohort and pathway-level validation."
    elif klass == "strong_ad_axis_with_nominal_external_replication" and pd_stat_supported and ad_fdr_supported:
        klass = "preliminary_shared_ad_pd_axis_candidate"
        safe = (
            f"{axis} has nominal-to-statistical support across AD and PD contexts, but still requires additional cohorts "
            "before any definitive shared-mechanism language."
        )
        next_step = "Add independent PD replication and pathway-level validation."
    unsafe = "Do not claim clinical utility, diagnostic use, causality, definitive mechanism, or validated shared AD/PD biology."
    return klass, safe, unsafe, next_step


def build_summary(
    phase22_rows: list[dict[str, str]],
    gse174367_rows: list[dict[str, str]],
    gse184950_rows: list[dict[str, str]],
    gse243639_rows: list[dict[str, str]],
    phase34_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    sea = sea_ad_primary(phase22_rows)
    pd_axis_primary = pd_gse243639_primary(phase22_rows)
    ad = by_axis(gse174367_rows)
    pd_rep = by_axis(gse184950_rows)
    phase34 = best_phase34_by_axis(phase34_rows or [])
    pd_support = pd_gse243639_status(gse243639_rows)
    axes = sorted(set(sea) | set(ad) | set(pd_rep) | set(pd_axis_primary) | set(phase34))
    rows: list[dict[str, str]] = []
    for axis in axes:
        sea_row = sea.get(axis, {})
        ad_row = ad.get(axis, {})
        pd_row = pd_rep.get(axis, {})
        pd_axis_row = pd_axis_primary.get(axis, {})
        phase34_row = phase34.get(axis, {})
        klass, safe, unsafe, next_step = classify_axis(axis, sea_row, ad_row, pd_row, pd_axis_row, phase34_row, pd_support)
        pd_axis_consistency = "not_available"
        if sea_row and pd_axis_row:
            pd_axis_consistency = (
                "consistent"
                if sign(to_float(sea_row.get("effect_size"))) == sign(to_float(pd_axis_row.get("effect_size")))
                and sign(to_float(sea_row.get("effect_size"))) != 0
                else "opposite_or_zero"
            )
        rows.append(
            {
                "axis_id": axis,
                "sea_ad_effect": sea_row.get("effect_size", ""),
                "sea_ad_p": sea_row.get("pvalue", ""),
                "sea_ad_fdr": sea_row.get("fdr", ""),
                "gse174367_effect": ad_row.get("effect_size", ""),
                "gse174367_p": ad_row.get("pvalue", ""),
                "gse174367_fdr": ad_row.get("fdr", ""),
                "gse174367_direction_consistency": ad_row.get("directional_consistency", ""),
                "gse184950_effect": pd_row.get("effect_size", ""),
                "gse184950_p": pd_row.get("pvalue", ""),
                "gse184950_fdr": pd_row.get("fdr", ""),
                "gse184950_direction_consistency": pd_row.get("directional_consistency", ""),
                "gse243639_axis_effect": pd_axis_row.get("effect_size", ""),
                "gse243639_axis_p": pd_axis_row.get("pvalue", ""),
                "gse243639_axis_fdr": pd_axis_row.get("fdr", ""),
                "gse243639_axis_empirical_p": pd_axis_row.get("empirical_pvalue", ""),
                "gse243639_axis_direction_consistency": pd_axis_consistency,
                "phase34_best_pd_cohort": phase34_row.get("cohort_id", ""),
                "phase34_best_pd_effect": phase34_row.get("effect_size", ""),
                "phase34_best_pd_p": phase34_row.get("pvalue", ""),
                "phase34_best_pd_fdr": phase34_row.get("fdr", ""),
                "phase34_best_pd_label": phase34_row.get("evidence_label", ""),
                "phase34_best_pd_direction_consistency": phase34_row.get("phase34_direction_consistency", ""),
                "phase34_best_pd_divergence_status": (
                    "statistically_supported_pd_divergence_candidate"
                    if phase34_row.get("evidence_label") == "opposite_direction"
                    and (to_float(phase34_row.get("pvalue"), 1.0) < 0.05 or to_float(phase34_row.get("fdr"), 1.0) < 0.1)
                    else ""
                ),
                "pd_gse243639_support": pd_support,
                "crosscohort_evidence_class": klass,
                "safe_claim": safe,
                "unsafe_claim": unsafe,
                "next_validation_needed": next_step,
            }
        )
    return rows


def rank_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        rows,
        key=lambda row: (
            CLASS_RANK.get(row.get("crosscohort_evidence_class", ""), 99),
            to_float(row.get("gse174367_p"), 1.0),
            to_float(row.get("sea_ad_fdr"), 1.0),
            -abs(to_float(row.get("gse174367_effect"), 0.0)),
        ),
    )


def best_phase34_by_axis(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    for row in rows:
        axis = row.get("axis_id", "")
        if not axis:
            continue
        current = best.get(axis)
        if current is None or to_float(row.get("pvalue"), 1.0) < to_float(current.get("pvalue"), 1.0):
            best[axis] = dict(row)
    return best


def read_phase34_defaults(paths: list[Path] | None = None) -> list[dict[str, str]]:
    selected = paths if paths is not None else [
        *sorted(Path("results/tables").glob("phase34_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase35_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase36_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase37_*_pd_axis_replication_statistics.tsv")),
    ]
    rows: list[dict[str, str]] = []
    for path in selected:
        for row in read_tsv(path):
            enriched = dict(row)
            enriched.setdefault("cohort_id", path.name.replace("phase34_", "").replace("phase35_", "").replace("phase36_", "").replace("phase37_", "").replace("_pd_axis_replication_statistics.tsv", ""))
            rows.append(enriched)
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Phase 32 cross-cohort axis evidence summary.")
    parser.add_argument("--phase22", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv"))
    parser.add_argument("--gse174367", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv"))
    parser.add_argument("--gse184950", type=Path, default=Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"))
    parser.add_argument("--gse243639", type=Path, default=Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--phase34-pd-stat", action="append", type=Path, default=[], help="Optional Phase 34 PD replication statistics table")
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/159_build_crosscohort_axis_evidence_summary.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = build_summary(read_tsv(args.phase22), read_tsv(args.gse174367), read_tsv(args.gse184950), read_tsv(args.gse243639), read_phase34_defaults(args.phase34_pd_stat or None))
    ranked = rank_rows(rows)
    write_tsv(args.outdir / "phase32_crosscohort_axis_evidence_summary.tsv", rows, SUMMARY_COLUMNS)
    write_tsv(args.outdir / "phase32_axis_evidence_ranked.tsv", ranked, SUMMARY_COLUMNS)
    logging.info("Wrote Phase 32 cross-cohort evidence rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
