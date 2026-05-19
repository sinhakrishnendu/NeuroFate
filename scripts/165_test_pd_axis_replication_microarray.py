#!/usr/bin/env python3
"""Test endpoint-locked PD axis replication in sample-level microarray/LCM cohorts."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path
from statistics import NormalDist


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


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def variance(values: list[float]) -> float:
    if len(values) < 2:
        return math.nan
    mu = mean(values)
    return sum((value - mu) ** 2 for value in values) / (len(values) - 1)


def smd(positive: list[float], negative: list[float]) -> float:
    pooled = math.sqrt(((len(positive) - 1) * variance(positive) + (len(negative) - 1) * variance(negative)) / max(1, len(positive) + len(negative) - 2))
    return (mean(positive) - mean(negative)) / pooled if pooled and not math.isnan(pooled) else math.nan


def rank_biserial(positive: list[float], negative: list[float]) -> float:
    wins = ties = 0.0
    total = len(positive) * len(negative)
    if total == 0:
        return math.nan
    for pos in positive:
        for neg in negative:
            if pos > neg:
                wins += 1
            elif pos == neg:
                ties += 1
    u = wins + 0.5 * ties
    return (2 * u / total) - 1


def mann_whitney_pvalue(positive: list[float], negative: list[float]) -> float:
    n1, n0 = len(positive), len(negative)
    if n1 < 2 or n0 < 2:
        return math.nan
    combined = [(value, 1) for value in positive] + [(value, 0) for value in negative]
    combined.sort(key=lambda item: item[0])
    ranks: list[tuple[float, int]] = []
    index = 0
    while index < len(combined):
        end = index + 1
        while end < len(combined) and combined[end][0] == combined[index][0]:
            end += 1
        rank = (index + 1 + end) / 2
        ranks.extend((rank, group) for _value, group in combined[index:end])
        index = end
    rank_sum_pos = sum(rank for rank, group in ranks if group == 1)
    u1 = rank_sum_pos - n1 * (n1 + 1) / 2
    mean_u = n1 * n0 / 2
    sd_u = math.sqrt(n1 * n0 * (n1 + n0 + 1) / 12)
    if sd_u == 0:
        return math.nan
    z = (u1 - mean_u) / sd_u
    return 2 * (1 - NormalDist().cdf(abs(z)))


def bh_fdr(pvalues: list[float]) -> list[float]:
    indexed = sorted((p if not math.isnan(p) else 1.0, index) for index, p in enumerate(pvalues))
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    m = len(pvalues)
    for rank, (pvalue, index) in enumerate(reversed(indexed), start=1):
        original_rank = m - rank + 1
        running = min(running, pvalue * m / original_rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def sign(value: float) -> int:
    if math.isnan(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def phase_pd_direction(axis_id: str, phase22: list[dict[str, str]], phase32: list[dict[str, str]]) -> int:
    for row in phase32:
        if row.get("axis_id") == axis_id:
            direction = sign(to_float(row.get("gse243639_axis_effect")))
            if direction:
                return direction
    for row in phase22:
        if row.get("axis_id") == axis_id and row.get("cohort") == "gse243639_pd_snpc":
            direction = sign(to_float(row.get("effect_size")))
            if direction:
                return direction
    return 0


def evidence_label(effect: float, reference_direction: int, pvalue: float, fdr: float, n: int, positive_n: int, negative_n: int) -> str:
    if n < 10 or positive_n < 2 or negative_n < 2 or sign(effect) == 0:
        return "insufficient_data"
    consistent = reference_direction != 0 and sign(effect) == reference_direction
    if consistent and (pvalue < 0.05 or fdr < 0.1):
        return "statistically_supported_pd_replication"
    if consistent:
        return "directionally_consistent_but_not_significant"
    if reference_direction != 0 and sign(effect) != reference_direction:
        return "opposite_direction"
    return "weak_or_no_replication"


def test_axes(axis_scores: list[dict[str, str]], cohort_id: str, phase22: list[dict[str, str]], phase32: list[dict[str, str]]) -> list[dict[str, str]]:
    axis_cols = [column for column in (axis_scores[0].keys() if axis_scores else []) if column.startswith("axis__")]
    rows: list[dict[str, str]] = []
    for column in axis_cols:
        axis_id = column.replace("axis__", "", 1)
        positive = [to_float(row.get(column)) for row in axis_scores if row.get("label__pd_vs_control") == "1"]
        negative = [to_float(row.get(column)) for row in axis_scores if row.get("label__pd_vs_control") == "0"]
        positive = [value for value in positive if not math.isnan(value)]
        negative = [value for value in negative if not math.isnan(value)]
        effect = rank_biserial(positive, negative)
        pvalue = mann_whitney_pvalue(positive, negative)
        rows.append(
            {
                "cohort_id": cohort_id,
                "axis_id": axis_id,
                "effect_size": "" if math.isnan(effect) else f"{effect:.8g}",
                "standardized_mean_difference": "" if math.isnan(smd(positive, negative)) else f"{smd(positive, negative):.8g}",
                "pvalue": "" if math.isnan(pvalue) else f"{pvalue:.8g}",
                "fdr": "",
                "n": str(len(positive) + len(negative)),
                "positive_n": str(len(positive)),
                "negative_n": str(len(negative)),
                "phase_pd_reference_direction": str(phase_pd_direction(axis_id, phase22, phase32)),
            }
        )
    fdrs = bh_fdr([to_float(row.get("pvalue"), 1.0) for row in rows])
    for row, fdr in zip(rows, fdrs, strict=False):
        effect = to_float(row.get("effect_size"))
        row["fdr"] = f"{fdr:.8g}"
        row["evidence_label"] = evidence_label(effect, int(row["phase_pd_reference_direction"]), to_float(row.get("pvalue"), 1.0), fdr, int(row["n"]), int(row["positive_n"]), int(row["negative_n"]))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run endpoint-locked PD axis replication in a sample-level cohort.")
    parser.add_argument("--axis-scores", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--phase22-evidence", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv"))
    parser.add_argument("--phase32-evidence", type=Path, default=Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fdr-output", type=Path, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = test_axes(read_tsv(args.axis_scores), args.cohort_id, read_tsv(args.phase22_evidence), read_tsv(args.phase32_evidence))
    columns = ["cohort_id", "axis_id", "effect_size", "standardized_mean_difference", "pvalue", "fdr", "n", "positive_n", "negative_n", "phase_pd_reference_direction", "evidence_label"]
    write_tsv(args.output, rows, columns)
    write_tsv(args.fdr_output, rows, columns)
    logging.info("Wrote Phase 33 PD replication statistics rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
