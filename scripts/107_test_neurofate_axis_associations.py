#!/usr/bin/env python3
"""Test donor/sample-level NeuroFate axis associations within AD and PD cohorts."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


AXIS_PREFIX = "axis__"
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


def to_float(value: str | None, default: float = math.nan) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def axis_columns(rows: list[dict[str, str]]) -> list[str]:
    return [column for column in (rows[0] if rows else {}) if column.startswith(AXIS_PREFIX)]


def label_columns(rows: list[dict[str, str]]) -> list[str]:
    return [
        column
        for column in (rows[0] if rows else {})
        if any(hint in column.lower() for hint in LABEL_HINTS)
        and not column.startswith(AXIS_PREFIX)
    ]


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def mann_whitney(values_a: list[float], values_b: list[float]) -> tuple[float, float]:
    combined = sorted([(value, 0) for value in values_a] + [(value, 1) for value in values_b])
    ranks: list[tuple[float, int]] = []
    index = 0
    while index < len(combined):
        end = index + 1
        while end < len(combined) and combined[end][0] == combined[index][0]:
            end += 1
        rank = (index + 1 + end) / 2.0
        ranks.extend((rank, group) for _value, group in combined[index:end])
        index = end
    rank_a = sum(rank for rank, group in ranks if group == 0)
    n_a = len(values_a)
    n_b = len(values_b)
    u_a = rank_a - n_a * (n_a + 1) / 2.0
    mean_u = n_a * n_b / 2.0
    sd_u = math.sqrt(n_a * n_b * (n_a + n_b + 1) / 12.0) or 1.0
    z = (u_a - mean_u) / sd_u
    pvalue = 2.0 * (1.0 - normal_cdf(abs(z)))
    rank_biserial = 2.0 * u_a / (n_a * n_b) - 1.0 if n_a and n_b else math.nan
    return rank_biserial, max(0.0, min(1.0, pvalue))


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / (len(values) - 1))


def standardized_mean_difference(values_a: list[float], values_b: list[float]) -> float:
    pooled = math.sqrt((sd(values_a) ** 2 + sd(values_b) ** 2) / 2.0) or 1.0
    return (mean(values_a) - mean(values_b)) / pooled


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


def test_axis_against_label(rows: list[dict[str, str]], cohort: str, axis: str, label: str) -> dict[str, str] | None:
    pairs = [(to_float(row.get(axis)), row.get(label, "")) for row in rows]
    pairs = [(value, label_value) for value, label_value in pairs if not math.isnan(value) and label_value not in {"", "NA", "nan"}]
    if len(pairs) < 4:
        return None
    numeric_labels = [to_float(label_value) for _value, label_value in pairs]
    axis_values = [value for value, _label_value in pairs]
    axis_id = axis.replace(AXIS_PREFIX, "")
    if all(not math.isnan(value) for value in numeric_labels) and len(set(numeric_labels)) > 2:
        effect, pvalue = spearman(axis_values, numeric_labels)
        test_name = "spearman_ordinal"
        n_a = n_b = ""
        smd = ""
    else:
        groups: dict[str, list[float]] = {}
        for value, label_value in pairs:
            groups.setdefault(label_value, []).append(value)
        if len(groups) != 2:
            return None
        ordered = sorted(groups)
        values_a = groups[ordered[1]]
        values_b = groups[ordered[0]]
        effect, pvalue = mann_whitney(values_a, values_b)
        smd = standardized_mean_difference(values_a, values_b)
        test_name = "mann_whitney_rank_biserial"
        n_a = str(len(values_a))
        n_b = str(len(values_b))
    return {
        "cohort": cohort,
        "axis_id": axis_id,
        "label": label,
        "test": test_name,
        "n": str(len(pairs)),
        "group_a_n": n_a,
        "group_b_n": n_b,
        "effect_size": "" if math.isnan(effect) else f"{effect:.8g}",
        "standardized_mean_difference": "" if smd == "" or math.isnan(smd) else f"{smd:.8g}",
        "pvalue": "" if math.isnan(pvalue) else f"{pvalue:.8g}",
        "missingness": str(len(rows) - len(pairs)),
    }


def run_tests(rows: list[dict[str, str]], cohort: str) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for axis in axis_columns(rows):
        for label in label_columns(rows):
            result = test_axis_against_label(rows, cohort, axis, label)
            if result is not None:
                output.append(result)
    pvalues = [to_float(row.get("pvalue"), 1.0) for row in output]
    fdr = bh_fdr(pvalues) if pvalues else []
    for row, adjusted in zip(output, fdr, strict=False):
        row["fdr"] = f"{adjusted:.8g}"
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Test NeuroFate axis associations at donor/sample level.")
    parser.add_argument("--sea-ad-axis", type=Path, default=Path("results/tables/phase21_sea_ad_axis_scores.tsv"))
    parser.add_argument("--pd-axis", type=Path, default=Path("results/tables/phase21_gse243639_axis_scores.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/107_test_neurofate_axis_associations.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = [
        *run_tests(read_tsv(args.sea_ad_axis), "sea_ad"),
        *run_tests(read_tsv(args.pd_axis), "gse243639_pd_snpc"),
    ]
    columns = [
        "cohort",
        "axis_id",
        "label",
        "test",
        "n",
        "group_a_n",
        "group_b_n",
        "effect_size",
        "standardized_mean_difference",
        "pvalue",
        "fdr",
        "missingness",
    ]
    write_tsv(args.outdir / "phase21_axis_association_statistics.tsv", rows, columns)
    write_tsv(args.outdir / "phase21_axis_effect_sizes.tsv", rows, columns)
    fdr_rows = [
        {
            "cohort": row["cohort"],
            "axis_id": row["axis_id"],
            "label": row["label"],
            "pvalue": row.get("pvalue", ""),
            "fdr": row.get("fdr", ""),
            "status": "fdr_below_0.10" if to_float(row.get("fdr"), 1.0) <= 0.10 else "not_fdr_supported",
        }
        for row in rows
    ]
    write_tsv(args.outdir / "phase21_axis_fdr_summary.tsv", fdr_rows, ["cohort", "axis_id", "label", "pvalue", "fdr", "status"])
    logging.info("Axis association tests generated rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
