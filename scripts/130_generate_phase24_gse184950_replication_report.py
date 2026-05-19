#!/usr/bin/env python3
"""Generate a conservative Phase 24 GSE184950 replication onboarding report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_rows(path: Path) -> int:
    return len(read_tsv(path))


def build_report(args: argparse.Namespace) -> str:
    metadata_n = count_rows(args.sample_metadata)
    inventory_n = count_rows(args.archive_inventory)
    plan_rows = read_tsv(args.extraction_plan)
    processed_status = "available_in_inventory" if any(row.get("matrix_status") != "processed_10x_matrix_not_detected" for row in plan_rows) else "not_yet_confirmed"
    return f"""# Phase 24 GSE184950 Replication Report

## 1. Dataset Overview

GSE184950 is the priority independent PD replication cohort for endpoint-locked NeuroFate-Axis validation. It is treated as a donor/sample-level replication target, not a clinical validation dataset.

## 2. Metadata Workbook Parsing

Parsed sample metadata rows: {metadata_n}. The workbook is expected to define disease state, donor ID, age, gender, postmortem interval, Braak stage, processed data files, and raw files.

## 3. RAW Archive Inventory

Archive inventory rows: {inventory_n}. The archive is listed only; no members are extracted by the lister.

## 4. Processed Matrix Availability

Processed matrix status: {processed_status}. Prefer processed 10x matrices if present. FASTQ/SRA processing is intentionally avoided in NeuroFate Phase 24.

## 5. Axis Gene Extraction Feasibility

Selective extraction should use manually prepared processed 10x directories and `scripts/127_extract_gse184950_axis_genes_from_10x.py` with the manual extraction guard enabled.

## 6. Replication Endpoint

The intended endpoint is PD/PDD versus Control using workbook disease state mapped to `label__pd_vs_control`.

## 7. Replication Status

Replication is pending until axis scores and endpoint-locked statistics are generated.

## 8. Remaining Blockers

- Manual download of `GSE184950_RAW.tar`.
- Archive inventory.
- Manual selective extraction of processed matrix members.
- Axis-score construction.
- Endpoint-locked replication statistics.

## 9. Next Manual Command

```bash
python scripts/124_parse_gse184950_geo_metadata_workbook.py
```
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 24 GSE184950 replication report.")
    parser.add_argument("--sample-metadata", type=Path, default=Path("results/tables/phase24_gse184950_sample_metadata.tsv"))
    parser.add_argument("--archive-inventory", type=Path, default=Path("results/tables/phase24_gse184950_raw_archive_inventory.tsv"))
    parser.add_argument("--extraction-plan", type=Path, default=Path("results/tables/phase24_gse184950_axis_extraction_plan.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase24_gse184950_replication_report.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(build_report(args), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
