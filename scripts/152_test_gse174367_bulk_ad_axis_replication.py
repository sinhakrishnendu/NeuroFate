#!/usr/bin/env python3
"""Endpoint-locked GSE174367 bulk RNA AD axis replication test."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


AXIS_PREFIX = "axis__"
OUTPUT_COLUMNS = [
    "cohort_id",
    "axis_id",
    "endpoint_column",
    "positive_class",
    "negative_class",
    "effect_size",
    "standardized_mean_difference",
    "pvalue",
    "fdr",
    "n",
    "positive_n",
    "negative_n",
    "phase22_ad_direction",
    "directional_consistency",
    "evidence_label",
]


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


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / (len(values) - 1))


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def sign(value: float) -> int:
    if math.isnan(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def rank_biserial(positive: list[float], negative: list[float]) -> tuple[float, float]:
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
    effect = 2.0 * u_pos / (n_pos * n_neg) - 1.0 if n_pos and n_neg else math.nan
    mean_u = n_pos * n_neg / 2.0
    sd_u = math.sqrt(n_pos * n_neg * (n_pos + n_neg + 1) / 12.0) or 1.0
    pvalue = 2.0 * (1.0 - normal_cdf(abs((u_pos - mean_u) / sd_u)))
    return effect, max(0.0, min(1.0, pvalue))


def smd(positive: list[float], negative: list[float]) -> float:
    pooled = math.sqrt((sd(positive) ** 2 + sd(negative) ** 2) / 2.0) or 1.0
    return (mean(positive) - mean(negative)) / pooled


def bh_fdr(pvalues: list[float]) -> list[float]:
    indexed = sorted((pvalue, index) for index, pvalue in enumerate(pvalues))
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank, (pvalue, index) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, pvalue * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def phase22_ad_direction(axis: str, rows: list[dict[str, str]]) -> str:
    row = next((item for item in rows if item.get("axis_id") == axis and item.get("endpoint_id") == "sea_ad_cognitive_dementia"), None)
    if row is None:
        return "not_available"
    effect = to_float(row.get("effect_size"))
    return "positive" if effect > 0 else "negative" if effect < 0 else "zero"


def evidence_label(effect: float, pvalue: float, fdr: float, phase22: str, n: int, positive_n: int, negative_n: int) -> tuple[str, str]:
    consistency = "not_available"
    if phase22 in {"positive", "negative"} and sign(effect) != 0:
        consistency = "consistent" if (effect > 0 and phase22 == "positive") or (effect < 0 and phase22 == "negative") else "opposite"
    if n < 20 or positive_n < 10 or negative_n < 10:
        return "insufficient_data", consistency
    if consistency == "consistent" and (pvalue < 0.05 or fdr < 0.1):
        return "statistically_supported_ad_replication", consistency
    if consistency == "consistent":
        return "directionally_consistent_but_not_significant", consistency
    if consistency == "opposite":
        return "opposite_direction", consistency
    return "weak_or_no_replication", consistency


def test_axes(rows: list[dict[str, str]], phase22_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    axes = [column for column in (rows[0] if rows else {}) if column.startswith(AXIS_PREFIX)]
    output: list[dict[str, str]] = []
    for axis_col in axes:
        positive: list[float] = []
        negative: list[float] = []
        for row in rows:
            value = to_float(row.get(axis_col))
            if math.isnan(value):
                continue
            if row.get("label__ad_vs_control") == "1":
                positive.append(value)
            elif row.get("label__ad_vs_control") == "0":
                negative.append(value)
        if not positive or not negative:
            continue
        effect, pvalue = rank_biserial(positive, negative)
        axis_id = axis_col.replace(AXIS_PREFIX, "")
        phase22 = phase22_ad_direction(axis_id, phase22_rows)
        output.append(
            {
                "cohort_id": "gse174367_ad_multiomics_bulk",
                "axis_id": axis_id,
                "endpoint_column": "label__ad_vs_control",
                "positive_class": "1",
                "negative_class": "0",
                "effect_size": f"{effect:.8g}",
                "standardized_mean_difference": f"{smd(positive, negative):.8g}",
                "pvalue": f"{pvalue:.8g}",
                "fdr": "",
                "n": str(len(positive) + len(negative)),
                "positive_n": str(len(positive)),
                "negative_n": str(len(negative)),
                "phase22_ad_direction": phase22,
                "directional_consistency": "",
                "evidence_label": "",
            }
        )
    fdrs = bh_fdr([to_float(row["pvalue"], 1.0) for row in output])
    for row, fdr in zip(output, fdrs, strict=False):
        label, consistency = evidence_label(
            to_float(row["effect_size"]),
            to_float(row["pvalue"], 1.0),
            fdr,
            row.get("phase22_ad_direction", ""),
            int(to_float(row.get("n"), 0)),
            int(to_float(row.get("positive_n"), 0)),
            int(to_float(row.get("negative_n"), 0)),
        )
        row["fdr"] = f"{fdr:.8g}"
        row["directional_consistency"] = consistency
        row["evidence_label"] = label
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GSE174367 bulk RNA endpoint-locked AD replication.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase29_gse174367_bulk_axis_scores.tsv"))
    parser.add_argument("--phase22-evidence", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv"))
    parser.add_argument("--fdr-output", type=Path, default=Path("results/tables/phase29_gse174367_bulk_axis_replication_fdr.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/152_test_gse174367_bulk_ad_axis_replication.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = test_axes(read_tsv(args.axis_scores), read_tsv(args.phase22_evidence))
    write_tsv(args.output, rows)
    write_tsv(args.fdr_output, rows)
    logging.info("GSE174367 bulk AD replication rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
