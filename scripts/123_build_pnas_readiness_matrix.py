#!/usr/bin/env python3
"""Build a conservative PNAS readiness matrix for NeuroFate-Axis."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


CRITERIA = [
    ("software_reproducibility", "NeuroFate package, docs, demo, CI, release metadata"),
    ("sea_ad_internal_evidence", "SEA-AD donor-level internal AD evidence"),
    ("pd_gse243639_evidence", "GSE243639 Phase 20 preliminary PD internal signal"),
    ("endpoint_locked_axis_evidence", "Phase 22 endpoint-locked axis evidence"),
    ("matched_random_controls", "Endpoint-matched random-axis controls"),
    ("independent_ad_replication", "At least one independent AD replication cohort"),
    ("independent_pd_replication", "At least one independent PD replication cohort"),
    ("pd_divergent_axis_candidate", "Direction-audited statistically supported opposite-direction PD axis candidate"),
    ("gse184950_series_metadata_parsed", "Phase 25 GSE184950 full series-matrix metadata parsed"),
    ("gse184950_archive_reconciled", "Phase 25 GSE184950 RAW archive inventory reconciled with series manifest"),
    ("gse184950_processed_matrices_available", "Phase 25 GSE184950 processed matrix availability established"),
    ("gse184950_axis_replication_pending_or_complete", "Phase 25 GSE184950 axis replication pending or complete"),
    ("shared_ad_pd_axis_claim", "Shared AD/PD axis claim readiness"),
    ("pnas_biological_claim", "Overall PNAS biological claim readiness"),
    ("network_pathway_interpretation", "Network/pathway interpretation beyond target-gene panel"),
    ("no_overclaiming_audit", "No high-severity unallowed overclaiming flags"),
    ("public_reproducibility_package", "Source release and results-review package without raw data"),
]


def exists(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not exists(path):
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def gse184950_statistical_support() -> bool:
    rows = read_tsv(Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"))
    return any(row.get("evidence_label") == "replicated_statistically_supported" for row in rows)


def gse184950_complete_processed_count() -> int:
    nested = read_tsv(Path("results/tables/phase26_gse184950_nested_archive_inventory.tsv"))
    complete_samples = {row.get("sample_id", "") for row in nested if row.get("complete_processed_matrix_set") == "true" and row.get("sample_id")}
    if complete_samples:
        return len(complete_samples)
    audit = read_tsv(Path("results/tables/phase26_gse184950_selected_extraction_audit.tsv"))
    return len({row.get("sample_id", "") for row in audit if row.get("status") == "extracted_processed_matrix_file_only" and row.get("sample_id")})


def phase33_pd_paths() -> list[Path]:
    return sorted(Path("results/tables").glob("phase33_*_pd_axis_replication_statistics.tsv"))


def phase34_pd_paths() -> list[Path]:
    return [
        *sorted(Path("results/tables").glob("phase34_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase35_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase36_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase37_*_pd_axis_replication_statistics.tsv")),
    ]


def phase33_pd_statistical_support() -> bool:
    for path in [*phase33_pd_paths(), *phase34_pd_paths()]:
        rows = read_tsv(path)
        if any(row.get("evidence_label") == "statistically_supported_pd_replication" for row in rows):
            return True
    return False


def phase34_pd_divergence_candidate() -> bool:
    for path in [*phase33_pd_paths(), *phase34_pd_paths()]:
        for row in read_tsv(path):
            if row.get("evidence_label") != "opposite_direction":
                continue
            try:
                pvalue = float(row.get("pvalue", "1") or 1)
                fdr = float(row.get("fdr", "1") or 1)
            except ValueError:
                continue
            if pvalue < 0.05 or fdr < 0.1:
                return True
    if exists(Path("results/tables/phase38_gse7621_axis_direction_probe_audit.tsv")):
        return any(
            row.get("phase38_direction_flag") == "statistically_significant_opposite_direction"
            for row in read_tsv(Path("results/tables/phase38_gse7621_axis_direction_probe_audit.tsv"))
        )
    return False


def phase33_pd_available() -> bool:
    return any(exists(path) for path in [*phase33_pd_paths(), *phase34_pd_paths()])


def phase33_shared_axis_converged() -> bool:
    nominal_ad_axes = {
        row.get("axis_id", "")
        for row in phase32_rows()
        if row.get("crosscohort_evidence_class") == "strong_ad_axis_with_nominal_external_replication"
    }
    if not nominal_ad_axes:
        return False
    for path in [*phase33_pd_paths(), *phase34_pd_paths()]:
        for row in read_tsv(path):
            if row.get("axis_id") in nominal_ad_axes and row.get("evidence_label") == "statistically_supported_pd_replication":
                return True
    return False


def phase28_ad_statistical_support() -> bool:
    paths = list(Path("results/tables").glob("phase28_*_axis_replication_statistics.tsv"))
    paths.append(Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv"))
    paths.append(Path("results/tables/phase30_gse174367_bulk_axis_replication_statistics.tsv"))
    paths.append(Path("results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv"))
    for path in paths:
        rows = read_tsv(path)
        if any(row.get("evidence_label") == "statistically_supported_ad_replication" for row in rows):
            return True
    return False


def phase28_ad_replication_available() -> bool:
    paths = list(Path("results/tables").glob("phase28_*_axis_replication_statistics.tsv"))
    paths.append(Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv"))
    paths.append(Path("results/tables/phase30_gse174367_bulk_axis_replication_statistics.tsv"))
    paths.append(Path("results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv"))
    return any(path.exists() and path.stat().st_size > 0 for path in paths)


def phase32_rows() -> list[dict[str, str]]:
    return read_tsv(Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv"))


def phase32_nominal_ad_replication() -> bool:
    return any(row.get("crosscohort_evidence_class") == "strong_ad_axis_with_nominal_external_replication" for row in phase32_rows())


def phase32_fdr_ad_replication() -> bool:
    for row in phase32_rows():
        if row.get("crosscohort_evidence_class") == "strong_ad_axis_with_nominal_external_replication":
            try:
                if float(row.get("gse174367_fdr", "1") or 1) < 0.1:
                    return True
            except ValueError:
                continue
    return False


def phase32_shared_axis_ready() -> bool:
    return any(row.get("crosscohort_evidence_class") == "preliminary_shared_ad_pd_axis_candidate" for row in phase32_rows())


def status_for(criterion: str) -> tuple[str, str]:
    checks = {
        "software_reproducibility": [Path("README.md"), Path("pyproject.toml"), Path("examples/tiny_demo")],
        "sea_ad_internal_evidence": [Path("results/tables/phase5_donor_feature_table.tsv")],
        "pd_gse243639_evidence": [Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv")],
        "endpoint_locked_axis_evidence": [Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv")],
        "matched_random_controls": [Path("results/tables/phase22_endpoint_locked_axis_empirical_pvalues.tsv")],
        "independent_ad_replication": [
            Path("results/tables/phase23_gse174367_ad_multiomics_axis_association_statistics.tsv"),
            Path("results/tables/phase23_gse147528_ad_progression_axis_association_statistics.tsv"),
            Path("results/tables/phase28_gse174367_ad_multiomics_axis_replication_statistics.tsv"),
            Path("results/tables/phase28_gse147528_ad_progression_axis_replication_statistics.tsv"),
            Path("results/tables/phase28_gse157827_ad_snuc_optional_axis_replication_statistics.tsv"),
            Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv"),
            Path("results/tables/phase30_gse174367_bulk_axis_replication_statistics.tsv"),
            Path("results/tables/phase29_gse174367_bulk_axis_replication_statistics.tsv"),
        ],
        "independent_pd_replication": [
            Path("results/tables/phase23_gse184950_pd_sn_axis_association_statistics.tsv"),
            Path("results/tables/phase24_gse184950_axis_replication_statistics.tsv"),
            Path("results/tables/phase25_gse184950_axis_replication_statistics.tsv"),
            Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"),
            *phase33_pd_paths(),
            *phase34_pd_paths(),
        ],
        "pd_divergent_axis_candidate": [
            Path("results/tables/phase38_gse7621_axis_direction_probe_audit.tsv"),
            *phase34_pd_paths(),
        ],
        "gse184950_series_metadata_parsed": [Path("results/tables/phase25_gse184950_series_sample_metadata.tsv")],
        "gse184950_archive_reconciled": [
            Path("results/tables/phase26_gse184950_nested_archive_inventory.tsv"),
            Path("results/tables/phase26_gse184950_selected_extraction_audit.tsv"),
        ],
        "gse184950_processed_matrices_available": [
            Path("results/tables/phase26_gse184950_nested_archive_inventory.tsv"),
            Path("results/tables/phase26_gse184950_selected_extraction_audit.tsv"),
        ],
        "gse184950_axis_replication_pending_or_complete": [
            Path("results/tables/phase27_gse184950_axis_scores_clean.tsv"),
            Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"),
        ],
        "shared_ad_pd_axis_claim": [Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv")],
        "pnas_biological_claim": [Path("results/tables/phase32_crosscohort_axis_evidence_summary.tsv")],
        "network_pathway_interpretation": [Path("metadata/neurofate_axis_registry.tsv")],
        "no_overclaiming_audit": [Path("results/reports/no_overclaiming_audit.tsv")],
        "public_reproducibility_package": [Path("release_artifacts")],
    }
    paths = checks[criterion]
    if criterion == "independent_ad_replication":
        if phase32_fdr_ad_replication():
            return "statistically_supported", ";".join(str(path) for path in paths)
        if phase32_nominal_ad_replication() or phase28_ad_statistical_support():
            return "nominally_supported", ";".join(str(path) for path in paths)
        if phase28_ad_replication_available() or any(exists(path) for path in paths):
            return "available_but_preliminary", ";".join(str(path) for path in paths)
        return "missing", ";".join(str(path) for path in paths)
    if criterion == "independent_pd_replication":
        if gse184950_statistical_support() or phase33_pd_statistical_support():
            return "statistically_supported", ";".join(str(path) for path in paths)
        if phase34_pd_divergence_candidate():
            return "mixed_pd_evidence", ";".join(str(path) for path in paths)
        if phase33_pd_available():
            return "available_but_preliminary", ";".join(str(path) for path in paths)
        if exists(Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv")):
            return "available_but_preliminary", ";".join(str(path) for path in paths)
    if criterion == "pd_divergent_axis_candidate":
        if phase34_pd_divergence_candidate():
            return "present", ";".join(str(path) for path in paths)
        return "missing_or_pending", ";".join(str(path) for path in paths)
    if criterion == "gse184950_processed_matrices_available":
        count = gse184950_complete_processed_count()
        if count >= 34:
            return "satisfied_34_processed_matrices", ";".join(str(path) for path in paths)
        if count:
            return "partial_processed_matrix_support", ";".join(str(path) for path in paths)
    if criterion == "gse184950_axis_replication_pending_or_complete":
        if gse184950_statistical_support():
            return "available_statistically_supported", ";".join(str(path) for path in paths)
        if exists(Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv")):
            return "available_but_weak", ";".join(str(path) for path in paths)
    if criterion == "shared_ad_pd_axis_claim":
        if phase33_shared_axis_converged():
            return "candidate_ready_for_review_not_definitive", "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
        if phase32_shared_axis_ready():
            return "preliminary_not_definitive", "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
        return "not_ready", "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
    if criterion == "pnas_biological_claim":
        if phase34_pd_divergence_candidate():
            return "promising_but_requires_pd_resolution", "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
        if phase32_nominal_ad_replication() and exists(Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv")):
            return "promising_but_not_ready", "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
        return "not_ready", "results/tables/phase32_crosscohort_axis_evidence_summary.tsv"
    if any(exists(path) for path in paths):
        return "available_or_partially_available", ";".join(str(path) for path in paths)
    return "missing_or_pending", ";".join(str(path) for path in paths)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["criterion", "description", "status", "evidence_paths", "pnas_gap"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, rows: list[dict[str, str]]) -> None:
    missing = [row for row in rows if row["status"] == "missing_or_pending"]
    lines = ["# Phase 23 PNAS Gap Report", "", "PNAS-level biological claims remain pending until independent replication and endpoint-locked controls are complete.", ""]
    for row in missing:
        lines.append(f"- `{row['criterion']}`: {row['pnas_gap']}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build NeuroFate-Axis PNAS readiness matrix.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--gap-report", type=Path, default=Path("results/reports/phase23_pnas_gap_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = []
    for criterion, description in CRITERIA:
        status, paths = status_for(criterion)
        pnas_gap = "satisfied_or_in_progress"
        if status in {"missing_or_pending", "available_but_preliminary", "available_but_weak", "partial_processed_matrix_support", "nominally_supported", "not_ready", "promising_but_not_ready", "mixed_pd_evidence", "present", "promising_but_requires_pd_resolution"}:
            pnas_gap = "required_before_strong_pnas_claims"
        rows.append(
            {
                "criterion": criterion,
                "description": description,
                "status": status,
                "evidence_paths": paths,
                "pnas_gap": pnas_gap,
            }
        )
    write_tsv(args.output, rows)
    write_report(args.gap_report, rows)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.gap_report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
