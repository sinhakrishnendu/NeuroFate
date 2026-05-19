#!/usr/bin/env python3
"""Integrate endpoint-locked discovery evidence with independent replication statistics."""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path


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


def sign(value: float) -> int:
    if math.isnan(value) or value == 0:
        return 0
    return 1 if value > 0 else -1


def parse_replication(items: list[str]) -> dict[str, list[dict[str, str]]]:
    out: dict[str, list[dict[str, str]]] = {}
    if not items:
        defaults = []
        if Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv").exists():
            defaults.append("gse184950_pd_sn=results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv")
        elif Path("results/tables/phase25_gse184950_axis_replication_statistics.tsv").exists():
            defaults.append("gse184950_pd_sn=results/tables/phase25_gse184950_axis_replication_statistics.tsv")
        elif Path("results/tables/phase24_gse184950_axis_replication_statistics.tsv").exists():
            defaults.append("gse184950_pd_sn=results/tables/phase24_gse184950_axis_replication_statistics.tsv")
        for path in sorted(Path("results/tables").glob("phase28_*_axis_replication_statistics.tsv")):
            cohort_id = path.name.replace("phase28_", "").replace("_axis_replication_statistics.tsv", "")
            defaults.append(f"{cohort_id}={path}")
        phase31_ad = Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv")
        phase30_ad = Path("results/tables/phase30_gse174367_bulk_axis_replication_statistics.tsv")
        phase29_ad = Path("results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv")
        if phase31_ad.exists():
            defaults.append(f"gse174367_ad_multiomics_bulk={phase31_ad}")
        elif phase30_ad.exists():
            defaults.append(f"gse174367_ad_multiomics_bulk={phase30_ad}")
        elif phase29_ad.exists():
            defaults.append(f"gse174367_ad_multiomics_bulk={phase29_ad}")
        for path in sorted(Path("results/tables").glob("phase33_*_pd_axis_replication_statistics.tsv")):
            cohort_id = path.name.replace("phase33_", "").replace("_pd_axis_replication_statistics.tsv", "")
            defaults.append(f"{cohort_id}={path}")
        for path in sorted(Path("results/tables").glob("phase34_*_pd_axis_replication_statistics.tsv")):
            cohort_id = path.name.replace("phase34_", "").replace("_pd_axis_replication_statistics.tsv", "")
            defaults.append(f"{cohort_id}={path}")
        for path in sorted(Path("results/tables").glob("phase35_*_pd_axis_replication_statistics.tsv")):
            cohort_id = path.name.replace("phase35_", "").replace("_pd_axis_replication_statistics.tsv", "")
            defaults.append(f"{cohort_id}={path}")
        for path in sorted(Path("results/tables").glob("phase36_*_pd_axis_replication_statistics.tsv")):
            cohort_id = path.name.replace("phase36_", "").replace("_pd_axis_replication_statistics.tsv", "")
            defaults.append(f"{cohort_id}={path}")
        for path in sorted(Path("results/tables").glob("phase37_*_pd_axis_replication_statistics.tsv")):
            cohort_id = path.name.replace("phase37_", "").replace("_pd_axis_replication_statistics.tsv", "")
            defaults.append(f"{cohort_id}={path}")
        items = defaults
    for item in items:
        if "=" not in item:
            continue
        cohort_id, path = item.split("=", 1)
        out[cohort_id] = read_tsv(Path(path))
    return out


def integrate(discovery: list[dict[str, str]], replication: dict[str, list[dict[str, str]]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for evidence in discovery:
        axis = evidence.get("axis_id", "")
        discovery_effect = to_float(evidence.get("effect_size"))
        matches = []
        for cohort_id, rep_rows in replication.items():
            rep = next((row for row in rep_rows if row.get("axis_id") == axis), None)
            if not rep:
                continue
            rep_effect = to_float(rep.get("effect_size"))
            consistent = sign(discovery_effect) == sign(rep_effect) and sign(discovery_effect) != 0
            pvalue = to_float(rep.get("pvalue"), 1.0)
            fdr = to_float(rep.get("fdr"), 1.0)
            n = int(to_float(rep.get("n"), 0))
            positive_n = int(to_float(rep.get("positive_n"), 0))
            negative_n = int(to_float(rep.get("negative_n"), 0))
            label = rep.get("evidence_label", "")
            statistically_supported = consistent and n >= 20 and positive_n >= 10 and negative_n >= 10 and (
                pvalue < 0.05
                or fdr < 0.1
                or label in {"replicated_statistically_supported", "statistically_supported_ad_replication", "statistically_supported_pd_replication"}
            )
            matches.append((cohort_id, rep, consistent, statistically_supported))
        if not matches:
            status = "insufficient_data"
            replicated = "false"
            reason = "No independent replication statistics supplied."
            upgrade = "false"
        elif any(statistically_supported for _cohort, _rep, _consistent, statistically_supported in matches):
            status = "statistically_supported_replication"
            replicated = "true"
            reason = "At least one replication cohort is directionally consistent and passes p/FDR support thresholds."
            upgrade = "true"
        elif any(consistent for _cohort, _rep, consistent, _supported in matches):
            status = "directionally_consistent_preliminary_signal"
            replicated = "false"
            reason = "Direction is consistent, but p/FDR support is insufficient; do not call this replicated."
            upgrade = "false"
        elif any(not consistent and sign(to_float(rep.get("effect_size"))) != 0 for _cohort, rep, consistent, _supported in matches):
            status = "opposite_direction"
            replicated = "false"
            reason = "Replication cohort effect is directionally opposite or incompatible."
            upgrade = "false"
        else:
            status = "weak_or_no_replication"
            replicated = "false"
            reason = "Replication cohort effect is weak, unavailable, or not directionally informative."
            upgrade = "false"
        rows.append(
            {
                "axis_id": axis,
                "discovery_endpoint_id": evidence.get("endpoint_id", ""),
                "discovery_cohort": evidence.get("cohort", ""),
                "discovery_effect_size": evidence.get("effect_size", ""),
                "discovery_axis_claim_class": evidence.get("axis_claim_class", ""),
                "replication_cohorts": ";".join(cohort for cohort, _rep, _consistent, _supported in matches),
                "directionally_consistent_replication": replicated,
                "replication_status": status,
                "claim_upgrade_allowed": upgrade,
                "safe_interpretation": reason,
                "disallowed_interpretation": "Do not claim confirmed cross-disease mechanism, clinical biomarker, diagnostic utility, or causality.",
            }
        )
    return rows


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["replication_status"]] = counts.get(row["replication_status"], 0) + 1
    lines = ["# Phase 23 Replication Readiness Report", "", "Replication integration does not upgrade claims from direction-only evidence. Statistical p/FDR support is required.", ""]
    for key, value in sorted(counts.items()):
        lines.append(f"- {key}: {value}")
    lines.extend(["", "Strong shared AD/PD mechanism language remains disallowed unless statistically supported independent replication exists."])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Integrate endpoint-locked replication evidence.")
    parser.add_argument("--phase22-evidence", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv"))
    parser.add_argument("--replication-stat", action="append", default=[], help="cohort_id=path")
    parser.add_argument("--outdir", type=Path, default=Path("results/tables"))
    parser.add_argument("--report-dir", type=Path, default=Path("results/reports"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/122_integrate_endpoint_locked_replication.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = integrate(read_tsv(args.phase22_evidence), parse_replication(args.replication_stat))
    write_tsv(
        args.outdir / "phase23_replication_integrated_axis_evidence.tsv",
        rows,
        ["axis_id", "discovery_endpoint_id", "discovery_cohort", "discovery_effect_size", "discovery_axis_claim_class", "replication_cohorts", "directionally_consistent_replication", "replication_status", "claim_upgrade_allowed", "safe_interpretation", "disallowed_interpretation"],
    )
    write_report(args.report_dir / "phase23_replication_readiness_report.md", rows)
    logging.info("Integrated replication evidence rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
