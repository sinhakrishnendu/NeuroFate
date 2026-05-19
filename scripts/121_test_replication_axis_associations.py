#!/usr/bin/env python3
"""Test endpoint-locked binary axis associations in a replication cohort."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


AXIS_PREFIX = "axis__"


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
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


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / (len(values) - 1))


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


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


def test_axes(rows: list[dict[str, str]], cohort_id: str, endpoint_column: str, positive_class: str, negative_class: str) -> list[dict[str, str]]:
    axes = [column for column in (rows[0] if rows else {}) if column.startswith(AXIS_PREFIX)]
    output: list[dict[str, str]] = []
    for axis in axes:
        positive: list[float] = []
        negative: list[float] = []
        for row in rows:
            label = norm(row.get(endpoint_column))
            value = to_float(row.get(axis))
            if math.isnan(value):
                continue
            if label == norm(positive_class):
                positive.append(value)
            elif label == norm(negative_class):
                negative.append(value)
        if not positive or not negative:
            continue
        effect, pvalue = rank_biserial(positive, negative)
        output.append(
            {
                "cohort_id": cohort_id,
                "axis_id": axis.replace(AXIS_PREFIX, ""),
                "endpoint_column": endpoint_column,
                "positive_class": positive_class,
                "negative_class": negative_class,
                "effect_size": f"{effect:.8g}",
                "standardized_mean_difference": f"{smd(positive, negative):.8g}",
                "pvalue": f"{pvalue:.8g}",
                "n": str(len(positive) + len(negative)),
                "positive_n": str(len(positive)),
                "negative_n": str(len(negative)),
            }
        )
    fdr = bh_fdr([to_float(row["pvalue"], 1.0) for row in output])
    for row, value in zip(output, fdr, strict=False):
        row["fdr"] = f"{value:.8g}"
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Test binary endpoint axis associations in a replication cohort.")
    parser.add_argument("--axis-scores", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--endpoint-column", required=True)
    parser.add_argument("--positive-class", required=True)
    parser.add_argument("--negative-class", required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--log-file", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_file = args.log_file or Path(f"results/logs/121_{args.cohort_id}_axis_associations.log")
    configure_logging(log_file)
    rows = test_axes(read_tsv(args.axis_scores), args.cohort_id, args.endpoint_column, args.positive_class, args.negative_class)
    output = args.output or Path(f"results/tables/phase23_{args.cohort_id}_axis_association_statistics.tsv")
    write_tsv(output, rows, ["cohort_id", "axis_id", "endpoint_column", "positive_class", "negative_class", "effect_size", "standardized_mean_difference", "pvalue", "fdr", "n", "positive_n", "negative_n"])
    logging.info("Replication associations cohort=%s rows=%d", args.cohort_id, len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
