#!/usr/bin/env python3
"""Generate the Phase 25 GSE184950 replication-readiness report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_by(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        value = row.get(column, "") or "missing"
        counts[value] = counts.get(value, 0) + 1
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 25 GSE184950 readiness report.")
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--manifest", type=Path, default=Path("results/tables/phase25_gse184950_series_processed_file_manifest.tsv"))
    parser.add_argument("--reconciliation", type=Path, default=Path("results/tables/phase25_gse184950_archive_series_reconciliation.tsv"))
    parser.add_argument("--plan", type=Path, default=Path("results/tables/phase25_gse184950_selective_extraction_plan.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase25_gse184950_replication_readiness_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metadata = read_tsv(args.metadata)
    manifest = read_tsv(args.manifest)
    reconciliation = read_tsv(args.reconciliation)
    plan = read_tsv(args.plan)
    labels = count_by(metadata, "disease_state")
    positives = sum(row.get("label__pd_pdd_vs_control") == "1" for row in metadata)
    negatives = sum(row.get("label__pd_pdd_vs_control") == "0" for row in metadata)
    found_archives = sum(row.get("found_in_archive") == "true" for row in reconciliation)
    processed_ready = sum(row.get("status") == "ready_for_manual_selective_processed_tar_extraction" for row in plan)

    lines = [
        "# Phase 25 GSE184950 Replication Readiness Report",
        "",
        "## 1. Why Series Matrix Replaced Workbook Metadata",
        "The small add2 workbook is incomplete for replication because it exposes only representative sample rows. The GEO series matrix is the primary metadata source for the full 34-sample GSE184950 cohort.",
        "",
        "## 2. Sample Counts and Endpoint",
        f"- Parsed samples: {len(metadata)}",
        f"- PD/PDD positive samples: {positives}",
        f"- Unaffected controls: {negatives}",
        "- Endpoint: `label__pd_pdd_vs_control`, where Parkinson's Disease and Parkinson's Disease Dementia are positive and Unaffected Control is negative.",
        "",
        "Disease-state counts:",
    ]
    lines.extend(f"- {label}: {count}" for label, count in sorted(labels.items()))
    lines.extend(
        [
            "",
            "## 3. Processed File Manifest",
            f"- Expected per-sample supplementary tar entries: {len(manifest)}",
            "",
            "## 4. RAW Archive Reconciliation",
            f"- Expected archives found in current inventory: {found_archives}",
            "- No archive members were extracted by this report.",
            "",
            "## 5. Selective Extraction Readiness",
            f"- Samples ready for manual selective processed-tar extraction: {processed_ready}",
            "- FASTQ/SRA processing remains out of scope unless no processed matrices are available and the user explicitly chooses a separate preprocessing workflow.",
            "",
            "## 6. Whether Processed Matrices Are Available",
            "Use `results/tables/phase25_gse184950_archive_series_reconciliation.tsv` and `results/tables/phase25_gse184950_selective_extraction_plan.tsv` to decide whether processed 10x matrices can be selectively extracted.",
            "",
            "## 7. Replication Endpoint Definition",
            "The replication endpoint is donor/sample-level PD/PDD versus Unaffected Control. This is a biological replication endpoint, not clinical validation.",
            "",
            "## 8. Remaining Blockers",
            "- Manual download of `GSE184950_RAW.tar` if not already present.",
            "- Safe archive listing and reconciliation.",
            "- Manual selective extraction of processed matrix files only.",
            "- Axis-gene sample-level extraction and endpoint-locked replication statistics.",
            "",
            "## 9. Next Manual Command",
            "`python scripts/132_reconcile_gse184950_archive_with_series_metadata.py` after `python scripts/125_list_gse184950_raw_archive.py` has created the archive inventory.",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
