#!/usr/bin/env python3
"""Run random-axis negative controls from donor/sample-level feature universes only."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import random
from pathlib import Path


FEATURE_PREFIXES = (
    "gene_mean__",
    "gene_detection__",
    "cell_fraction__",
    "celltype_gene_mean__",
    "celltype_gene_detection__",
    "index__",
)
LABEL_HINTS = ("label__", "diagnosis", "cognitive", "braak", "cerad", "lewy", "pathology", "status")


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


def is_feature_column(column: str) -> bool:
    lowered = column.lower()
    return column.startswith(FEATURE_PREFIXES) and not any(hint in lowered for hint in LABEL_HINTS)


def feature_universe(rows: list[dict[str, str]]) -> list[str]:
    if not rows:
        return []
    return [
        column
        for column in rows[0]
        if is_feature_column(column)
        and any(not math.isnan(to_float(row.get(column))) for row in rows)
    ]


def axis_size(axis: dict[str, str]) -> int:
    return max(1, len([gene for gene in axis.get("gene_members", "").replace(",", ";").split(";") if gene.strip()]))


def mean_abs_feature_value(rows: list[dict[str, str]], features: list[str]) -> float:
    values: list[float] = []
    for row in rows:
        row_values = [to_float(row.get(feature)) for feature in features]
        observed = [value for value in row_values if not math.isnan(value)]
        if observed:
            values.append(abs(sum(observed) / len(observed)))
    return sum(values) / len(values) if values else math.nan


def run_controls_for_cohort(
    rows: list[dict[str, str]],
    cohort: str,
    axes: list[dict[str, str]],
    n_random: int,
    seed: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rng = random.Random(seed)
    universe = feature_universe(rows)
    control_rows: list[dict[str, str]] = []
    empirical_rows: list[dict[str, str]] = []
    for axis in axes:
        size = min(axis_size(axis), len(universe))
        if size == 0:
            empirical_rows.append(
                {
                    "cohort": cohort,
                    "axis_id": axis["axis_id"],
                    "curated_axis_score": "",
                    "random_axis_mean": "",
                    "empirical_pvalue": "",
                    "status": "insufficient_feature_universe",
                }
            )
            continue
        # Curated score uses features whose names contain member gene symbols when available.
        genes = [gene.strip().upper() for gene in axis.get("gene_members", "").replace(",", ";").split(";") if gene.strip()]
        curated_features = [feature for feature in universe if any(gene in feature.upper() for gene in genes)]
        curated_features = curated_features or rng.sample(universe, size)
        curated_score = mean_abs_feature_value(rows, curated_features)
        random_scores = []
        for index in range(n_random):
            sampled = rng.sample(universe, size)
            score = mean_abs_feature_value(rows, sampled)
            random_scores.append(score)
            control_rows.append(
                {
                    "cohort": cohort,
                    "axis_id": axis["axis_id"],
                    "random_id": str(index + 1),
                    "random_feature_count": str(size),
                    "random_axis_score": "" if math.isnan(score) else f"{score:.8g}",
                }
            )
        observed_random = [score for score in random_scores if not math.isnan(score)]
        empirical_p = (
            (sum(score >= curated_score for score in observed_random) + 1) / (len(observed_random) + 1)
            if observed_random and not math.isnan(curated_score)
            else math.nan
        )
        empirical_rows.append(
            {
                "cohort": cohort,
                "axis_id": axis["axis_id"],
                "curated_axis_score": "" if math.isnan(curated_score) else f"{curated_score:.8g}",
                "random_axis_mean": "" if not observed_random else f"{sum(observed_random) / len(observed_random):.8g}",
                "empirical_pvalue": "" if math.isnan(empirical_p) else f"{empirical_p:.8g}",
                "status": "random_control_complete",
            }
        )
    return control_rows, empirical_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run NeuroFate random-axis controls.")
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--pd-features", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--n-random", type=int, default=500)
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/109_run_axis_randomization_controls.log"))
    parser.add_argument("--seed", type=int, default=2109)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    axes = read_tsv(args.axis_registry)
    sea_control, sea_empirical = run_controls_for_cohort(read_tsv(args.sea_ad_features), "sea_ad", axes, args.n_random, args.seed)
    pd_control, pd_empirical = run_controls_for_cohort(read_tsv(args.pd_features), "gse243639_pd_snpc", axes, args.n_random, args.seed + 1)
    write_tsv(
        args.outdir / "phase21_random_axis_controls.tsv",
        [*sea_control, *pd_control],
        ["cohort", "axis_id", "random_id", "random_feature_count", "random_axis_score"],
    )
    write_tsv(
        args.outdir / "phase21_axis_empirical_pvalues.tsv",
        [*sea_empirical, *pd_empirical],
        ["cohort", "axis_id", "curated_axis_score", "random_axis_mean", "empirical_pvalue", "status"],
    )
    logging.info("Random-axis controls generated with n_random=%d", args.n_random)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
