#!/usr/bin/env python3
"""Generate Phase 28 independent AD replication onboarding report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 28 AD replication report.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase28_ad_replication_registry.tsv"))
    parser.add_argument("--inventory", type=Path, default=Path("results/reports/phase28_ad_replication_file_inventory.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase28_ad_replication_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    registry = read_tsv(args.registry)
    inventory = read_tsv(args.inventory)
    readiness = read_tsv(args.readiness)
    ad_ready = next((row for row in readiness if row.get("criterion") == "independent_ad_replication"), {})
    lines = [
        "# Phase 28 AD Replication Report",
        "",
        "## 1. Why AD Replication Is Required",
        "SEA-AD provides strong internal AD axis evidence, but PNAS-level biological claims need independent AD replication.",
        "",
        "## 2. Candidate Cohorts",
    ]
    for row in registry:
        lines.append(f"- `{row['cohort_id']}` ({row['geo_accession']}): priority {row['priority']}, {row['first_use_strategy']}")
    lines.extend(
        [
            "",
            "## 3. Local File Status",
            f"- Inventory rows: {len(inventory)}",
            "",
            "## 4. Metadata Availability",
            "Use `scripts/143_parse_geo_series_matrix_generic.py` when a GEO series matrix is present locally.",
            "",
            "## 5. Matrix Availability",
            "Use bulk/sample-level matrices first when available. snRNA routes remain manual planning only.",
            "",
            "## 6. Endpoint Definition",
            "The first AD replication endpoint should be endpoint-locked AD versus Control, or high AD pathology versus control/low pathology when diagnosis labels are unavailable.",
            "",
            "## 7. Axis Coverage If Available",
            "Coverage is written by `scripts/145_build_ad_replication_axis_scores_from_matrix.py` after a local sample-level expression matrix is selected.",
            "",
            "## 8. Replication Status",
            f"- Independent AD replication readiness: {ad_ready.get('status', 'missing_or_pending')}",
            "",
            "## 9. PNAS Implication",
            "PNAS-level shared-axis claims remain premature until at least one independent AD cohort shows statistically supported endpoint-locked replication.",
            "",
            "## 10. Next Manual Command",
            "`RUN_MANUAL_DOWNLOAD=YES bash scripts/manual_downloads/download_gse174367_manual.sh` after reviewing GEO GSE174367 files and data-use terms.",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
