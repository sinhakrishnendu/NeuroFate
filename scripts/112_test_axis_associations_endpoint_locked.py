#!/usr/bin/env python3
"""Endpoint-locked donor/sample-level NeuroFate axis association tests."""

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


def axis_columns(rows: list[dict[str, str]]) -> list[str]:
    return [column for column in (rows[0] if rows else {}) if column.startswith(AXIS_PREFIX)]


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / (len(values) - 1))


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


def mann_whitney_rank_biserial(positive: list[float], negative: list[float]) -> tuple[float, float]:
    combined = sorted([(value, 1) for value in positive] + [(value, 0) for value in negative])
    ranked: list[tuple[float, int]] = []
    index = 0
    while index < len(combined):
        end = index + 1
        while end < len(combined) and combined[end][0] == combined[index][0]:
            end += 1
        rank = (index + 1 + end) / 2.0
        ranked.extend((rank, group) for _value, group in combined[index:end])
        index = end
    n_pos = len(positive)
    n_neg = len(negative)
    rank_pos = sum(rank for rank, group in ranked if group == 1)
    u_pos = rank_pos - n_pos * (n_pos + 1) / 2.0
    rank_biserial = 2.0 * u_pos / (n_pos * n_neg) - 1.0 if n_pos and n_neg else math.nan
    mean_u = n_pos * n_neg / 2.0
    sd_u = math.sqrt(n_pos * n_neg * (n_pos + n_neg + 1) / 12.0) or 1.0
    z = (u_pos - mean_u) / sd_u
    pvalue = 2.0 * (1.0 - normal_cdf(abs(z)))
    return rank_biserial, max(0.0, min(1.0, pvalue))


def standardized_mean_difference(positive: list[float], negative: list[float]) -> float:
    pooled = math.sqrt((sd(positive) ** 2 + sd(negative) ** 2) / 2.0) or 1.0
    return (mean(positive) - mean(negative)) / pooled


def spearman(x: list[float], y: list[float]) -> tuple[float, float]:
    if len(x) < 4:
        return math.nan, math.nan
    rx = ranks(x)
    ry = ranks(y)
    mx = mean(rx)
    my = mean(ry)
    numerator = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=False))
    denom = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) or 1.0
    rho = numerator / denom
    z = rho * math.sqrt(max(1, len(x) - 3))
    pvalue = 2.0 * (1.0 - normal_cdf(abs(z)))
    return rho, max(0.0, min(1.0, pvalue))


def bh_fdr(pvalues: list[float]) -> list[float]:
    indexed = sorted((pvalue, index) for index, pvalue in enumerate(pvalues))
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    total = len(pvalues)
    for rank, (pvalue, index) in reversed(list(enumerate(indexed, start=1))):
        running = min(running, pvalue * total / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def cohort_rows(sea_rows: list[dict[str, str]], pd_rows: list[dict[str, str]], cohort: str) -> list[dict[str, str]]:
    return sea_rows if cohort == "sea_ad" else pd_rows if cohort == "gse243639_pd_snpc" else []


def test_binary(rows: list[dict[str, str]], endpoint: dict[str, str], axis: str) -> dict[str, str] | None:
    source = endpoint["source_column"]
    positive_class = norm(endpoint["positive_class"])
    negative_class = norm(endpoint["negative_class"])
    exclude = {norm(value) for value in split_values(endpoint.get("exclude_values", ""))}
    positive: list[float] = []
    negative: list[float] = []
    missing = 0
    for row in rows:
        label = norm(row.get(source))
        value = to_float(row.get(axis))
        if math.isnan(value) or label in exclude or label == "":
            missing += 1
            continue
        if label == positive_class:
            positive.append(value)
        elif label == negative_class:
            negative.append(value)
        else:
            missing += 1
    minimum_n = int(to_float(endpoint.get("minimum_n"), 0))
    minimum_per_class = int(to_float(endpoint.get("minimum_per_class"), 0))
    if len(positive) + len(negative) < minimum_n or min(len(positive), len(negative)) < minimum_per_class:
        return None
    effect, pvalue = mann_whitney_rank_biserial(positive, negative)
    smd = standardized_mean_difference(positive, negative)
    axis_id = axis.replace(AXIS_PREFIX, "")
    return {
        "endpoint_id": endpoint["endpoint_id"],
        "cohort": endpoint["cohort"],
        "endpoint_role": endpoint["endpoint_role"],
        "axis_id": axis_id,
        "source_column": source,
        "endpoint_type": "binary",
        "test": "mann_whitney_rank_biserial",
        "effect_size": f"{effect:.8g}",
        "standardized_mean_difference": f"{smd:.8g}",
        "pvalue": f"{pvalue:.8g}",
        "n": str(len(positive) + len(negative)),
        "positive_n": str(len(positive)),
        "negative_n": str(len(negative)),
        "valid_n": str(len(positive) + len(negative)),
        "missingness": str(missing),
        "direction": "positive_minus_negative",
        "status": "tested",
    }


def test_ordinal(rows: list[dict[str, str]], endpoint: dict[str, str], axis: str) -> dict[str, str] | None:
    source = endpoint["source_column"]
    exclude = {norm(value) for value in split_values(endpoint.get("exclude_values", ""))}
    order = {norm(value): index for index, value in enumerate(split_values(endpoint.get("ordinal_order", "")), start=1)}
    values: list[float] = []
    labels: list[float] = []
    missing = 0
    for row in rows:
        label = norm(row.get(source))
        value = to_float(row.get(axis))
        if math.isnan(value) or label in exclude or label not in order:
            missing += 1
            continue
        values.append(value)
        labels.append(float(order[label]))
    minimum_n = int(to_float(endpoint.get("minimum_n"), 0))
    if len(values) < minimum_n or len(set(labels)) < 2:
        return None
    effect, pvalue = spearman(values, labels)
    axis_id = axis.replace(AXIS_PREFIX, "")
    return {
        "endpoint_id": endpoint["endpoint_id"],
        "cohort": endpoint["cohort"],
        "endpoint_role": endpoint["endpoint_role"],
        "axis_id": axis_id,
        "source_column": source,
        "endpoint_type": "ordinal",
        "test": "spearman_ordinal",
        "effect_size": f"{effect:.8g}",
        "standardized_mean_difference": "",
        "pvalue": f"{pvalue:.8g}",
        "n": str(len(values)),
        "positive_n": "",
        "negative_n": "",
        "valid_n": str(len(values)),
        "missingness": str(missing),
        "direction": "increasing_ordinal",
        "status": "tested",
    }


def run_endpoint_tests(
    sea_rows: list[dict[str, str]],
    pd_rows: list[dict[str, str]],
    endpoints: list[dict[str, str]],
) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for endpoint in endpoints:
        rows = cohort_rows(sea_rows, pd_rows, endpoint["cohort"])
        for axis in axis_columns(rows):
            result = (
                test_binary(rows, endpoint, axis)
                if endpoint["endpoint_type"] == "binary"
                else test_ordinal(rows, endpoint, axis)
            )
            if result is not None:
                output.append(result)
    for endpoint_id in sorted({row["endpoint_id"] for row in output}):
        selected = [row for row in output if row["endpoint_id"] == endpoint_id]
        fdrs = bh_fdr([to_float(row.get("pvalue"), 1.0) for row in selected])
        for row, fdr in zip(selected, fdrs, strict=False):
            row["fdr"] = f"{fdr:.8g}"
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run endpoint-locked NeuroFate axis association tests.")
    parser.add_argument("--sea-ad-axis", type=Path, default=Path("results/tables/phase21_sea_ad_axis_scores.tsv"))
    parser.add_argument("--pd-axis", type=Path, default=Path("results/tables/phase21_gse243639_axis_scores.tsv"))
    parser.add_argument("--endpoint-registry", type=Path, default=Path("metadata/neurofate_axis_endpoint_registry.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/112_test_axis_associations_endpoint_locked.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = run_endpoint_tests(read_tsv(args.sea_ad_axis), read_tsv(args.pd_axis), read_tsv(args.endpoint_registry))
    columns = [
        "endpoint_id",
        "cohort",
        "endpoint_role",
        "axis_id",
        "source_column",
        "endpoint_type",
        "test",
        "effect_size",
        "standardized_mean_difference",
        "pvalue",
        "fdr",
        "n",
        "positive_n",
        "negative_n",
        "valid_n",
        "missingness",
        "direction",
        "status",
    ]
    write_tsv(args.outdir / "phase22_endpoint_locked_axis_statistics.tsv", rows, columns)
    write_tsv(args.outdir / "phase22_endpoint_locked_axis_effects.tsv", rows, columns)
    fdr_rows = [
        {
            "endpoint_id": row["endpoint_id"],
            "cohort": row["cohort"],
            "axis_id": row["axis_id"],
            "pvalue": row.get("pvalue", ""),
            "fdr": row.get("fdr", ""),
            "status": "fdr_below_0.10" if to_float(row.get("fdr"), 1.0) <= 0.10 else "not_fdr_supported",
        }
        for row in rows
    ]
    write_tsv(args.outdir / "phase22_endpoint_locked_axis_fdr.tsv", fdr_rows, ["endpoint_id", "cohort", "axis_id", "pvalue", "fdr", "status"])
    logging.info("Endpoint-locked axis statistics rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
