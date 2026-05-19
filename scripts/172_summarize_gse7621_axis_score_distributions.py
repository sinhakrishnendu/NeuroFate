#!/usr/bin/env python3
"""Summarize GSE7621 axis-score distributions by PD/control label."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")])


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


def sd(values: list[float]) -> float:
    if len(values) < 2:
        return math.nan
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / (len(values) - 1))


def median(values: list[float]) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def fmt(value: float) -> str:
    return "" if math.isnan(value) else f"{value:.8g}"


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    if not rows:
        return []
    axis_cols = [col for col in rows[0] if col.startswith("axis__")]
    out: list[dict[str, str]] = []
    for col in axis_cols:
        axis_id = col.replace("axis__", "", 1)
        control = [to_float(row.get(col)) for row in rows if row.get("label__pd_vs_control") == "0"]
        pd = [to_float(row.get(col)) for row in rows if row.get("label__pd_vs_control") == "1"]
        control = [value for value in control if not math.isnan(value)]
        pd = [value for value in pd if not math.isnan(value)]
        delta = mean(pd) - mean(control)
        med_delta = median(pd) - median(control)
        direction = "pd_lower_than_control" if delta < 0 else "pd_higher_than_control" if delta > 0 else "no_difference"
        out.append(
            {
                "axis_id": axis_id,
                "control_n": str(len(control)),
                "pd_n": str(len(pd)),
                "control_mean": fmt(mean(control)),
                "pd_mean": fmt(mean(pd)),
                "control_sd": fmt(sd(control)),
                "pd_sd": fmt(sd(pd)),
                "pd_minus_control": fmt(delta),
                "median_difference": fmt(med_delta),
                "direction": direction,
            }
        )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize GSE7621 axis score distributions.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_axis_scores.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase38_gse7621_axis_score_distributions.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/172_summarize_gse7621_axis_score_distributions.log"))
    args = parser.parse_args()
    configure_logging(args.log_file)
    rows = summarize(read_tsv(args.axis_scores))
    write_tsv(args.output, rows, ["axis_id", "control_n", "pd_n", "control_mean", "pd_mean", "control_sd", "pd_sd", "pd_minus_control", "median_difference", "direction"])
    logging.info("Wrote GSE7621 axis-score distribution rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
