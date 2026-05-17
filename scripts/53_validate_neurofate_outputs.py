#!/usr/bin/env python3
"""Validate expected NeuroFate output files without recomputing analyses."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_OUTPUTS = {
    "core_sea_ad": [
        "results/tables/table1_sea_ad_publication_ready.tsv",
        "data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv",
    ],
    "phase3": [
        "results/tables/gene_by_celltype_summary.tsv",
        "results/tables/microglial_activation_signature.tsv",
        "results/figures/figure1_celltype_composition.png",
    ],
    "phase4": [
        "results/tables/phase4_gene_statistics.tsv",
        "results/tables/phase4_composite_indices.tsv",
        "results/figures/figure5_braak_associations.png",
    ],
    "phase5": [
        "results/tables/phase5_donor_feature_table.tsv",
        "results/tables/phase5_model_metrics.tsv",
        "results/tables/phase5_neurofate_scores.tsv",
    ],
    "phase6": [
        "results/tables/phase6_mps_model_metrics.tsv",
        "results/tables/phase6_mps_training_log.tsv",
        "results/tables/phase6_mps_predictions.tsv",
    ],
}

OPTIONAL_OUTPUTS = {
    "phase9_mathys": [
        "results/tables/mathys_2019_phase5_donor_feature_table.tsv",
        "results/tables/phase9_mathys_external_validation_metrics.tsv",
        "results/tables/phase9_mathys_external_predictions.tsv",
    ]
}


def output_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for group, paths in REQUIRED_OUTPUTS.items():
        for path_text in paths:
            path = Path(path_text)
            rows.append(
                {
                    "group": group,
                    "path": path_text,
                    "required": "true",
                    "exists": str(path.exists()).lower(),
                    "size_bytes": str(path.stat().st_size) if path.exists() and path.is_file() else "0",
                    "status": "present" if path.exists() else "missing",
                }
            )
    for group, paths in OPTIONAL_OUTPUTS.items():
        available = any(Path(path).exists() for path in paths)
        for path_text in paths:
            path = Path(path_text)
            rows.append(
                {
                    "group": group,
                    "path": path_text,
                    "required": "false",
                    "exists": str(path.exists()).lower(),
                    "size_bytes": str(path.stat().st_size) if path.exists() and path.is_file() else "0",
                    "status": "present" if path.exists() else "optional_missing" if available else "not_started",
                }
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate NeuroFate expected outputs.")
    parser.add_argument("--output", type=Path, default=Path("results/reports/output_validation_report.tsv"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = output_rows()
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["group", "path", "required", "exists", "size_bytes", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
