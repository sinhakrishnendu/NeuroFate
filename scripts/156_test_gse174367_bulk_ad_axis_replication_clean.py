#!/usr/bin/env python3
"""Run clean Phase 30 GSE174367 bulk AD axis replication statistics."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


def load_phase29_tester():
    script = Path(__file__).with_name("152_test_gse174367_bulk_ad_axis_replication.py")
    spec = importlib.util.spec_from_file_location("phase29_gse174367_bulk_tester", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 31 clean GSE174367 bulk AD axis replication.")
    parser.add_argument("--axis-scores", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_scores.tsv"))
    parser.add_argument("--phase22-evidence", type=Path, default=Path("results/tables/phase22_endpoint_locked_axis_evidence_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_replication_statistics.tsv"))
    parser.add_argument("--fdr-output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_replication_fdr.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/156_test_gse174367_bulk_ad_axis_replication_clean_phase31.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tester = load_phase29_tester()
    tester.configure_logging(args.log_file)
    rows = tester.test_axes(tester.read_tsv(args.axis_scores), tester.read_tsv(args.phase22_evidence))
    for row in rows:
        row["cohort_id"] = "gse174367_ad_multiomics_bulk_phase31"
    tester.write_tsv(args.output, rows)
    tester.write_tsv(args.fdr_output, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
