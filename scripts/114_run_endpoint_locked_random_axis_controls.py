#!/usr/bin/env python3
"""Matched random-axis controls using the same locked endpoint statistic as curated axes."""

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
LABEL_HINTS = ("label__", "diagnosis", "cognitive", "braak", "cerad", "lewy", "pathology", "status", "sample_id", "donor_id")


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


def norm(value: str | None) -> str:
    return str(value or "").strip().casefold()


def split_values(value: str) -> list[str]:
    return [part.strip() for part in value.replace(",", ";").split(";") if part.strip()]


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


def field_matches_gene(field: str, gene: str) -> bool:
    upper_field = field.upper()
    token = gene.upper()
    return upper_field.endswith(f"__{token}") or f"__{token}__" in upper_field or upper_field.endswith(f"_{token}")


def axis_features(universe: list[str], axis: dict[str, str]) -> list[str]:
    genes = [gene.strip().upper() for gene in axis.get("gene_members", "").replace(",", ";").split(";") if gene.strip()]
    return [feature for feature in universe if any(field_matches_gene(feature, gene) for gene in genes)]


def mean(values: list[float]) -> float:
    observed = [value for value in values if not math.isnan(value)]
    return sum(observed) / len(observed) if observed else math.nan


def score_rows(rows: list[dict[str, str]], features: list[str]) -> list[float]:
    return [mean([to_float(row.get(feature)) for feature in features]) for row in rows]


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def ranks(values: list[float]) -> list[float]:
    ordered = sorted((value, index) for index, value in enumerate(values))
    out = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][0] == ordered[index][0]:
            end += 1
        rank = (index + 1 + end) / 2.0
        for _value, original_index in ordered[index:end]:
            out[original_index] = rank
        index = end
    return out


def rank_biserial(positive: list[float], negative: list[float]) -> float:
    combined = sorted([(value, 1) for value in positive] + [(value, 0) for value in negative])
    rank_sum = 0.0
    index = 0
    while index < len(combined):
        end = index + 1
        while end < len(combined) and combined[end][0] == combined[index][0]:
            end += 1
        rank = (index + 1 + end) / 2.0
        rank_sum += sum(rank for _value, group in combined[index:end] if group == 1)
        index = end
    n_pos = len(positive)
    n_neg = len(negative)
    u_pos = rank_sum - n_pos * (n_pos + 1) / 2.0
    return 2.0 * u_pos / (n_pos * n_neg) - 1.0 if n_pos and n_neg else math.nan


def spearman(x: list[float], y: list[float]) -> float:
    if len(x) < 4:
        return math.nan
    rx = ranks(x)
    ry = ranks(y)
    mx = mean(rx)
    my = mean(ry)
    numerator = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=False))
    denom = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) or 1.0
    return numerator / denom


def endpoint_effect(rows: list[dict[str, str]], scores: list[float], endpoint: dict[str, str]) -> tuple[float, int, int, int]:
    source = endpoint["source_column"]
    exclude = {norm(value) for value in split_values(endpoint.get("exclude_values", ""))}
    if endpoint["endpoint_type"] == "binary":
        positive_class = norm(endpoint["positive_class"])
        negative_class = norm(endpoint["negative_class"])
        positive: list[float] = []
        negative: list[float] = []
        for row, score in zip(rows, scores, strict=False):
            label = norm(row.get(source))
            if math.isnan(score) or label in exclude or label == "":
                continue
            if label == positive_class:
                positive.append(score)
            elif label == negative_class:
                negative.append(score)
        return rank_biserial(positive, negative), len(positive) + len(negative), len(positive), len(negative)
    order = {norm(value): index for index, value in enumerate(split_values(endpoint.get("ordinal_order", "")), start=1)}
    values: list[float] = []
    labels: list[float] = []
    for row, score in zip(rows, scores, strict=False):
        label = norm(row.get(source))
        if math.isnan(score) or label in exclude or label not in order:
            continue
        values.append(score)
        labels.append(float(order[label]))
    return spearman(values, labels), len(values), len(values), 0


def cohort_feature_rows(sea_rows: list[dict[str, str]], pd_rows: list[dict[str, str]], cohort: str) -> list[dict[str, str]]:
    return sea_rows if cohort == "sea_ad" else pd_rows if cohort == "gse243639_pd_snpc" else []


def run_controls(
    sea_rows: list[dict[str, str]],
    pd_rows: list[dict[str, str]],
    axes: list[dict[str, str]],
    endpoints: list[dict[str, str]],
    n_random: int,
    seed: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rng = random.Random(seed)
    control_rows: list[dict[str, str]] = []
    empirical_rows: list[dict[str, str]] = []
    for endpoint in endpoints:
        rows = cohort_feature_rows(sea_rows, pd_rows, endpoint["cohort"])
        universe = feature_universe(rows)
        if not rows or not universe:
            continue
        for axis in axes:
            curated = axis_features(universe, axis)
            if not curated:
                empirical_rows.append(
                    {
                        "endpoint_id": endpoint["endpoint_id"],
                        "cohort": endpoint["cohort"],
                        "axis_id": axis["axis_id"],
                        "curated_effect_size": "",
                        "random_feature_count": "0",
                        "random_abs_effect_mean": "",
                        "empirical_pvalue": "",
                        "status": "insufficient_curated_axis_features",
                    }
                )
                continue
            curated_effect, n, positive_n, negative_n = endpoint_effect(rows, score_rows(rows, curated), endpoint)
            size = min(len(curated), len(universe))
            random_effects: list[float] = []
            for index in range(n_random):
                sampled = rng.sample(universe, size)
                effect, _n, _pos, _neg = endpoint_effect(rows, score_rows(rows, sampled), endpoint)
                random_effects.append(effect)
                control_rows.append(
                    {
                        "endpoint_id": endpoint["endpoint_id"],
                        "cohort": endpoint["cohort"],
                        "axis_id": axis["axis_id"],
                        "random_id": str(index + 1),
                        "random_feature_count": str(size),
                        "random_effect_size": "" if math.isnan(effect) else f"{effect:.8g}",
                    }
                )
            observed = [effect for effect in random_effects if not math.isnan(effect)]
            empirical = (
                (sum(abs(effect) >= abs(curated_effect) for effect in observed) + 1) / (len(observed) + 1)
                if observed and not math.isnan(curated_effect)
                else math.nan
            )
            empirical_rows.append(
                {
                    "endpoint_id": endpoint["endpoint_id"],
                    "cohort": endpoint["cohort"],
                    "axis_id": axis["axis_id"],
                    "curated_effect_size": "" if math.isnan(curated_effect) else f"{curated_effect:.8g}",
                    "random_feature_count": str(size),
                    "random_abs_effect_mean": "" if not observed else f"{mean([abs(value) for value in observed]):.8g}",
                    "empirical_pvalue": "" if math.isnan(empirical) else f"{empirical:.8g}",
                    "status": "endpoint_locked_random_control_complete",
                    "n": str(n),
                    "positive_n": str(positive_n),
                    "negative_n": str(negative_n),
                }
            )
    return control_rows, empirical_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run endpoint-locked matched random-axis association controls.")
    parser.add_argument("--sea-ad-features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--pd-features", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--endpoint-registry", type=Path, default=Path("metadata/neurofate_axis_endpoint_registry.tsv"))
    parser.add_argument("--n-random", type=int, default=1000)
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/114_run_endpoint_locked_random_axis_controls.log"))
    parser.add_argument("--seed", type=int, default=2214)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    control_rows, empirical_rows = run_controls(
        read_tsv(args.sea_ad_features),
        read_tsv(args.pd_features),
        read_tsv(args.axis_registry),
        read_tsv(args.endpoint_registry),
        args.n_random,
        args.seed,
    )
    write_tsv(
        args.outdir / "phase22_endpoint_locked_random_axis_controls.tsv",
        control_rows,
        ["endpoint_id", "cohort", "axis_id", "random_id", "random_feature_count", "random_effect_size"],
    )
    write_tsv(
        args.outdir / "phase22_endpoint_locked_axis_empirical_pvalues.tsv",
        empirical_rows,
        ["endpoint_id", "cohort", "axis_id", "curated_effect_size", "random_feature_count", "random_abs_effect_mean", "empirical_pvalue", "status", "n", "positive_n", "negative_n"],
    )
    logging.info("Endpoint-locked random controls rows=%d empirical_rows=%d", len(control_rows), len(empirical_rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
