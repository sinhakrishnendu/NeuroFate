#!/usr/bin/env python3
"""Generate the Phase 29 GSE174367 bulk AD replication report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_labels(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = row.get(column, "missing") or "missing"
        counts[value] = counts.get(value, 0) + 1
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 29 GSE174367 bulk AD replication report.")
    parser.add_argument("--structure", type=Path, default=Path("results/reports/phase29_gse174367_bulk_rda_structure.tsv"))
    parser.add_argument("--sample-map", type=Path, default=Path("results/tables/phase31_gse174367_bulk_sample_map.tsv"))
    parser.add_argument("--coverage", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_feature_coverage.tsv"))
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_scores.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv"))
    parser.add_argument("--mapping-audit", type=Path, default=Path("results/tables/phase30_gse174367_bulk_sample_mapping_audit.tsv"))
    parser.add_argument("--gene-mapping", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_mapping.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase29_gse174367_ad_replication_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    structure = read_tsv(args.structure)
    sample_map = read_tsv(args.sample_map)
    mapping_audit = read_tsv(args.mapping_audit)
    gene_mapping = read_tsv(args.gene_mapping)
    coverage = read_tsv(args.coverage)
    axis_scores = read_tsv(args.axis_scores)
    replication = read_tsv(args.replication)
    matched = sum(1 for row in sample_map if row.get("match_status") == "matched")
    labels = count_labels(axis_scores, "label__ad_vs_control")
    supported = [row for row in replication if row.get("evidence_label") == "statistically_supported_ad_replication"]
    directional = [row for row in replication if row.get("evidence_label") == "directionally_consistent_but_not_significant"]
    lines = [
        "# Phase 29 GSE174367 Bulk RNA Independent AD Replication",
        "",
        "## 1. Dataset Summary",
        "GSE174367 is used here through the processed bulk RNA resource because NeuroFate-Axis tests donor/sample-level axis associations.",
        "",
        "## 2. Why Bulk RNA Was Used First",
        "The processed bulk RNA file is small and sample-level, making it the appropriate first independent AD replication route. Single-nucleus expression processing remains optional and was not used in Phase 29.",
        "",
        "## 3. RDA Structure",
        f"RDA structure rows available: {len(structure)}.",
        "",
        "## 4. Sample Mapping",
        "Phase 30 uses the RDA-internal `targets` table as the primary sample map when it overlaps expression colnames.",
        f"Expression samples matched to metadata: {matched} of {len(sample_map)}.",
        f"Mapping audit rows: {len(mapping_audit)}.",
        "",
        "## 5. Axis-Gene Coverage",
        f"Axis coverage rows: {len(coverage)}.",
        f"Mapped axis-gene rows: {len(gene_mapping)}.",
        "",
        "## 6. Endpoint Definition",
        "`label__ad_vs_control` is endpoint-locked: AD is positive (1) and Control is negative (0).",
        f"Observed label counts: {labels}.",
        "",
        "## 7. AD vs Control Replication Statistics",
        f"Replication rows: {len(replication)}. Statistically supported AD replication rows: {len(supported)}. Direction-only preliminary rows: {len(directional)}.",
        "",
        "## 8. Directional Agreement With SEA-AD Phase 22",
        "Each axis is compared against the SEA-AD Phase 22 dementia endpoint direction. Direction alone is not treated as replication.",
        "",
        "## 9. PNAS Implication",
        "Independent AD replication can strengthen NeuroFate-Axis only when endpoint-locked effects are directionally consistent and supported by p-value or FDR thresholds.",
        "",
        "## 10. Remaining Limitations",
        "Bulk tissue composition, sample mapping ambiguity, axis-gene coverage, and cohort-specific processing can limit interpretation.",
        "",
        "## 11. Whether snRNA Route Is Still Needed",
        "A single-nucleus route may still be useful if bulk RNA is insufficient, but it should remain processed-matrix based and sample-level.",
        "",
        "No clinical, diagnostic, causal, or definitive cross-disease claim is made from Phase 29 alone.",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
