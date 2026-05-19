#!/usr/bin/env python3
"""Audit GSE7621 axis direction and GPL570 probe mapping."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


FOCUS_AXES = {
    "synuclein_mitochondrial_axis",
    "neuronal_vulnerability_axis",
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


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def probe_counts(probe_map: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, set[str]] = {}
    for row in probe_map:
        gene = row.get("gene_symbol", "").upper()
        probe = row.get("probe_id", "")
        if gene and probe:
            counts.setdefault(gene, set()).add(probe)
    return {gene: len(probes) for gene, probes in counts.items()}


def classify_axis(axis_id: str, effect: float, pvalue: float, fdr: float, label: str) -> str:
    if axis_id == "synuclein_mitochondrial_axis" and effect < 0 and (pvalue < 0.05 or fdr < 0.1) and label == "opposite_direction":
        return "statistically_significant_opposite_direction"
    if axis_id == "neuronal_vulnerability_axis" and label == "directionally_consistent_but_not_significant":
        return "directionally_consistent_not_significant"
    if label == "opposite_direction" and (pvalue < 0.05 or fdr < 0.1):
        return "statistically_significant_opposite_direction"
    return label or "unclassified"


def summarize_direction(axis_scores: list[dict[str, str]], axis_id: str) -> tuple[str, str, str]:
    col = f"axis__{axis_id}"
    control = [to_float(row.get(col)) for row in axis_scores if row.get("label__pd_vs_control") == "0"]
    pd = [to_float(row.get(col)) for row in axis_scores if row.get("label__pd_vs_control") == "1"]
    control = [value for value in control if not math.isnan(value)]
    pd = [value for value in pd if not math.isnan(value)]
    delta = mean(pd) - mean(control)
    direction = "pd_lower_than_control" if delta < 0 else "pd_higher_than_control" if delta > 0 else "no_difference"
    return f"{mean(control):.8g}", f"{mean(pd):.8g}", direction


def build_audit(
    axis_scores: list[dict[str, str]],
    probe_map_rows: list[dict[str, str]],
    coverage_rows: list[dict[str, str]],
    stats_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    counts = probe_counts(probe_map_rows)
    coverage = {row.get("axis_id", ""): row for row in coverage_rows}
    stats = {row.get("axis_id", ""): row for row in stats_rows}
    rows: list[dict[str, str]] = []
    for axis_id in sorted(set(coverage) | set(stats)):
        cov = coverage.get(axis_id, {})
        stat = stats.get(axis_id, {})
        found = [gene for gene in cov.get("found_gene_members", "").split(";") if gene]
        missing = [gene for gene in cov.get("missing_gene_members", "").split(";") if gene]
        focus = axis_id in FOCUS_AXES
        probe_summary = ";".join(f"{gene}:{counts.get(gene, 0)}" for gene in found) if focus else ""
        multi_probe = any(counts.get(gene, 0) > 1 for gene in found)
        control_mean, pd_mean, direction = summarize_direction(axis_scores, axis_id)
        effect = to_float(stat.get("effect_size"))
        pvalue = to_float(stat.get("pvalue"), 1.0)
        fdr = to_float(stat.get("fdr"), 1.0)
        flag = classify_axis(axis_id, effect, pvalue, fdr, stat.get("evidence_label", ""))
        rows.append(
            {
                "axis_id": axis_id,
                "genes_requested": cov.get("genes_requested", ""),
                "genes_found": cov.get("genes_found", ""),
                "genes_missing": cov.get("genes_missing", ""),
                "found_gene_members": cov.get("found_gene_members", ""),
                "missing_gene_members": cov.get("missing_gene_members", ""),
                "prkn_missing": "true" if "PRKN" in missing else "false",
                "focus_axis_probe_counts": probe_summary,
                "multi_probe_genes_present": "true" if multi_probe else "false",
                "control_axis_mean": control_mean,
                "pd_axis_mean": pd_mean,
                "axis_score_direction": direction,
                "effect_size": stat.get("effect_size", ""),
                "pvalue": stat.get("pvalue", ""),
                "fdr": stat.get("fdr", ""),
                "evidence_label": stat.get("evidence_label", ""),
                "phase38_direction_flag": flag,
                "safe_interpretation": safe_interpretation(axis_id, flag),
            }
        )
    return rows


def safe_interpretation(axis_id: str, flag: str) -> str:
    if flag == "statistically_significant_opposite_direction":
        return f"{axis_id} is a candidate PD-divergent signal in GSE7621, not shared AD/PD replication."
    if flag == "directionally_consistent_not_significant":
        return f"{axis_id} is directionally consistent but not statistically supported in GSE7621."
    return f"{axis_id} remains exploratory in GSE7621."


def write_preview(path: Path, rows: list[dict[str, str]]) -> None:
    syn = next((row for row in rows if row["axis_id"] == "synuclein_mitochondrial_axis"), {})
    neuronal = next((row for row in rows if row["axis_id"] == "neuronal_vulnerability_axis"), {})
    lines = [
        "# Phase 38 GSE7621 Direction/Probe Audit",
        "",
        "GSE7621 is technically valid for sample-level axis testing, but directionality must be interpreted conservatively.",
        "",
        "## Synuclein Mitochondrial Axis",
        f"- Direction flag: {syn.get('phase38_direction_flag', 'missing')}",
        f"- Effect: {syn.get('effect_size', '')}; p={syn.get('pvalue', '')}; FDR={syn.get('fdr', '')}",
        f"- Axis score direction: {syn.get('axis_score_direction', '')}",
        f"- Missing genes: {syn.get('missing_gene_members', '') or 'none'}",
        f"- Probe counts: {syn.get('focus_axis_probe_counts', '')}",
        "",
        "## Neuronal Vulnerability Axis",
        f"- Direction flag: {neuronal.get('phase38_direction_flag', 'missing')}",
        f"- Effect: {neuronal.get('effect_size', '')}; p={neuronal.get('pvalue', '')}; FDR={neuronal.get('fdr', '')}",
        f"- Axis score direction: {neuronal.get('axis_score_direction', '')}",
        f"- Missing genes: {neuronal.get('missing_gene_members', '') or 'none'}",
        "",
        "The significant opposite-direction synuclein-mitochondrial result is not shared AD/PD replication.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit GSE7621 direction and NeuroFate probe mapping.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_axis_scores.tsv"))
    parser.add_argument("--probe-map", type=Path, default=Path("results/tables/phase34_GPL570_axis_probe_mapping.tsv"))
    parser.add_argument("--axis-feature-coverage", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_axis_feature_coverage.tsv"))
    parser.add_argument("--replication-stats", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_pd_axis_replication_statistics.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase38_gse7621_axis_direction_probe_audit.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("results/reports/phase38_gse7621_axis_direction_probe_audit.md"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/171_audit_gse7621_axis_direction_and_probe_mapping.log"))
    args = parser.parse_args()
    configure_logging(args.log_file)
    rows = build_audit(read_tsv(args.axis_scores), read_tsv(args.probe_map), read_tsv(args.axis_feature_coverage), read_tsv(args.replication_stats))
    columns = [
        "axis_id",
        "genes_requested",
        "genes_found",
        "genes_missing",
        "found_gene_members",
        "missing_gene_members",
        "prkn_missing",
        "focus_axis_probe_counts",
        "multi_probe_genes_present",
        "control_axis_mean",
        "pd_axis_mean",
        "axis_score_direction",
        "effect_size",
        "pvalue",
        "fdr",
        "evidence_label",
        "phase38_direction_flag",
        "safe_interpretation",
    ]
    write_tsv(args.output, rows, columns)
    write_preview(args.preview_output, rows)
    logging.info("Wrote GSE7621 direction/probe audit rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
