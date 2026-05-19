#!/usr/bin/env python3
"""Generate the Phase 37 GSE7621 PD replication report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def counts(rows: list[dict[str, str]], key: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for row in rows:
        value = row.get(key, "") or "missing"
        out[value] = out.get(value, 0) + 1
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Phase 37 GSE7621 PD replication report.")
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_sample_metadata.tsv"))
    parser.add_argument("--platform-summary", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_platform_summary.tsv"))
    parser.add_argument("--mapping-audit", type=Path, default=Path("results/tables/phase37_gse7621_sample_mapping_audit.tsv"))
    parser.add_argument("--coverage", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_axis_feature_coverage.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("results/tables/phase37_gse7621_pd_sn_bulk_pd_axis_replication_statistics.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase37_gse7621_pd_replication_report.md"))
    args = parser.parse_args()

    metadata = read_tsv(args.metadata)
    platforms = read_tsv(args.platform_summary)
    mapping = read_tsv(args.mapping_audit)
    coverage = read_tsv(args.coverage)
    stats = read_tsv(args.replication)
    readiness = read_tsv(args.readiness)
    best_join = next((row.get("best_join_key", "") for row in mapping if row.get("is_best_join_key") == "true"), "")
    neuronal = next((row for row in stats if row.get("axis_id") == "neuronal_vulnerability_axis"), {})
    independent_pd = next((row for row in readiness if row.get("criterion") == "independent_pd_replication"), {})
    shared_claim = next((row for row in readiness if row.get("criterion") == "shared_ad_pd_axis_claim"), {})

    lines = [
        "# Phase 37 GSE7621 PD Replication Report",
        "",
        "## 1. Why GSE7621 Was Added",
        "GSE20141 produced technically successful but statistically weak PD axis evidence. GSE7621 is the next donor/sample-level substantia nigra cohort intended to test whether PD axis replication strengthens.",
        "",
        "## 2. Metadata and Endpoint",
    ]
    if metadata:
        lines.append(f"- Samples parsed: {len(metadata)}")
        lines.append(f"- Endpoint counts: {counts(metadata, 'label__pd_vs_control')}")
    else:
        lines.append("- Metadata has not been parsed yet.")
    lines.extend(["", "## 3. Platform/Probe Mapping"])
    if platforms:
        for row in platforms:
            lines.append(f"- Platform `{row.get('platform_id', 'missing')}`: {row.get('sample_count', '')} samples.")
    else:
        lines.append("- Platform summary is not available yet.")
    lines.extend(["", "## 4. Axis Coverage"])
    if coverage:
        lines.append(f"- Axes with mapped genes: {sum(row.get('status') == 'ok' for row in coverage)}/{len(coverage)}")
    else:
        lines.append("- Axis coverage has not been generated.")
    lines.extend(["", "## 5. PD Replication Statistics"])
    if stats:
        lines.append(f"- Evidence categories: {counts(stats, 'evidence_label')}")
    else:
        lines.append("- PD replication statistics are not available yet.")
    lines.extend(["", "## 6. Neuronal Vulnerability Result"])
    if neuronal:
        lines.append(f"- effect={neuronal.get('effect_size', '')}; p={neuronal.get('pvalue', '')}; FDR={neuronal.get('fdr', '')}; label={neuronal.get('evidence_label', '')}.")
    else:
        lines.append("- The neuronal vulnerability axis has not been tested in GSE7621 yet.")
    lines.extend(
        [
            "",
            "## 7. Whether PNAS Bottleneck Improved",
            f"- Best sample join key: {best_join or 'not_available'}",
            f"- independent_pd_replication: {independent_pd.get('status', 'not_available')}",
            f"- shared_ad_pd_axis_claim: {shared_claim.get('status', 'not_available')}",
            "Direction-only support remains preliminary. Shared AD/PD mechanism language remains blocked unless statistically supported independent PD replication aligns with an AD-supported axis.",
            "",
            "## 8. Remaining Limitations",
            "- This route uses sample-level GEO expression values and platform probe mappings only.",
            "- Raw CEL/CHP, FASTQ/SRA, single-cell embedding, clustering, and deep-learning routes are not used.",
            "- Small PD cohorts may remain underpowered even when biologically plausible effects are present.",
            "",
            "## 9. Next Dataset If Needed",
            "If GSE7621 does not produce statistically supported PD replication, proceed to region-aware GSE8397 or clean subseries from GSE20186.",
            "",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
