#!/usr/bin/env python3
"""Generate Phase 38 GSE7621 direction-aware interpretation report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find(rows: list[dict[str, str]], axis_id: str) -> dict[str, str]:
    return next((row for row in rows if row.get("axis_id") == axis_id), {})


def count(rows: list[dict[str, str]], key: str, value: str) -> int:
    return sum(1 for row in rows if row.get(key) == value)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Phase 38 GSE7621 interpretation report.")
    parser.add_argument("--debug", type=Path, default=Path("results/tables/phase37_gse7621_builder_metadata_debug.tsv"))
    parser.add_argument("--coverage", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_axis_feature_coverage.tsv"))
    parser.add_argument("--stats", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_pd_axis_replication_statistics.tsv"))
    parser.add_argument("--audit", type=Path, default=Path("results/tables/phase38_gse7621_axis_direction_probe_audit.tsv"))
    parser.add_argument("--distributions", type=Path, default=Path("results/tables/phase38_gse7621_axis_score_distributions.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase38_gse7621_interpretation_report.md"))
    args = parser.parse_args()

    debug = read_tsv(args.debug)
    coverage = read_tsv(args.coverage)
    stats = read_tsv(args.stats)
    audit = read_tsv(args.audit)
    distributions = read_tsv(args.distributions)
    readiness = read_tsv(args.readiness)
    syn = find(stats, "synuclein_mitochondrial_axis")
    neuronal = find(stats, "neuronal_vulnerability_axis")
    syn_audit = find(audit, "synuclein_mitochondrial_axis")
    neuronal_audit = find(audit, "neuronal_vulnerability_axis")
    syn_dist = find(distributions, "synuclein_mitochondrial_axis")
    neuronal_dist = find(distributions, "neuronal_vulnerability_axis")
    debug_row = debug[0] if debug else {}
    pd_readiness = next((row for row in readiness if row.get("criterion") == "independent_pd_replication"), {})
    divergent = next((row for row in readiness if row.get("criterion") == "pd_divergent_axis_candidate"), {})
    shared = next((row for row in readiness if row.get("criterion") == "shared_ad_pd_axis_claim"), {})
    lines = [
        "# Phase 38 GSE7621 Direction-Aware Interpretation Report",
        "",
        "## 1. GSE7621 Technical Validity",
        f"- Selected join key: {debug_row.get('selected_join_key', 'not_available')}",
        f"- Matched samples: {debug_row.get('matched_sample_count', 'not_available')}/{debug_row.get('expression_sample_count', 'not_available')}",
        f"- Final labelled samples: {debug_row.get('final_labeled_matched_sample_count', 'not_available')}",
        "",
        "## 2. Sample and Label Integrity",
        f"- Label counts: {debug_row.get('label_counts', 'not_available')}",
        "- Endpoint is sample-level PD versus Control.",
        "",
        "## 3. Axis Coverage",
        f"- Axes covered: {count(coverage, 'status', 'ok')}/{len(coverage)}",
        "- NeuroFate gene coverage: 29/30 mapped; PRKN is missing from the GPL570 map used here.",
        "",
        "## 4. Neuronal Vulnerability Result",
        f"- Effect={neuronal.get('effect_size', '')}; SMD={neuronal.get('standardized_mean_difference', '')}; p={neuronal.get('pvalue', '')}; FDR={neuronal.get('fdr', '')}.",
        f"- Evidence label: {neuronal.get('evidence_label', '')}.",
        f"- Distribution direction: {neuronal_dist.get('direction', '')}; PD-minus-Control={neuronal_dist.get('pd_minus_control', '')}.",
        f"- Phase 38 flag: {neuronal_audit.get('phase38_direction_flag', '')}.",
        "",
        "## 5. Synuclein-Mitochondrial Opposite-Direction Result",
        f"- Effect={syn.get('effect_size', '')}; SMD={syn.get('standardized_mean_difference', '')}; p={syn.get('pvalue', '')}; FDR={syn.get('fdr', '')}.",
        f"- Evidence label: {syn.get('evidence_label', '')}.",
        f"- Distribution direction: {syn_dist.get('direction', '')}; PD-minus-Control={syn_dist.get('pd_minus_control', '')}.",
        f"- Probe counts: {syn_audit.get('focus_axis_probe_counts', '')}.",
        f"- Missing genes: {syn_audit.get('missing_gene_members', '') or 'none'}.",
        "",
        "## 6. Why This Does Not Validate a Shared AD/PD Axis",
        "The synuclein-mitochondrial axis is statistically strong in GSE7621 but opposite in direction to the current cross-cohort reference. It must not be counted as shared AD/PD replication.",
        "",
        "## 7. Why This May Indicate PD-Divergent Biology",
        "The result is compatible with a candidate PD-divergent axis, especially because it is strong after FDR correction and involves a PD-relevant synuclein/mitochondrial theme. This remains a candidate interpretation until it is confirmed in another PD cohort and platform/probe behavior is audited.",
        "",
        "## 8. What Requires Further Validation",
        "- Independent PD replication of the synuclein-mitochondrial direction.",
        "- Confirmation that the opposite direction is not driven by platform, probe aggregation, tissue composition, or cohort-specific preprocessing.",
        "- Stronger PD support aligned with an AD-replicated axis before shared-axis claims are made.",
        "",
        "## 9. Journal-Readiness Implication",
        f"- independent_pd_replication: {pd_readiness.get('status', 'not_available')}",
        f"- pd_divergent_axis_candidate: {divergent.get('status', 'not_available')}",
        f"- shared_ad_pd_axis_claim: {shared.get('status', 'not_available')}",
        "PNAS readiness improves in biological nuance but not in definitive shared-axis support.",
        "",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
