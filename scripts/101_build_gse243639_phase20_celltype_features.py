#!/usr/bin/env python3
"""Build Phase 20 GSE243639 cell-type features from the Phase 19 safe map."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


def load_phase85_module():
    script = Path(__file__).resolve().parent / "85_build_gse243639_celltype_feature_table.py"
    spec = importlib.util.spec_from_file_location("phase85_celltype_features", script)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load script 85 feature builder.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Phase 20 safe-map-based GSE243639 cell-type feature table.")
    parser.add_argument("--expression", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_sparse_gene_panel_expression.tsv.gz"))
    parser.add_argument("--annotations", type=Path, default=Path("data/interim/external/gse243639_pd_snpc/gse243639_safe_cell_annotation_map.tsv"))
    parser.add_argument("--clinical", type=Path, default=Path("data/raw/external/gse243639_pd_snpc/GSE243639_Clinical_data.csv.gz"))
    parser.add_argument("--phase5-schema", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--schema-output", type=Path, default=Path("results/tables/phase20_gse243639_celltype_schema_alignment.tsv"))
    parser.add_argument("--label-summary-output", type=Path, default=Path("results/tables/phase20_gse243639_celltype_label_summary.tsv"))
    parser.add_argument("--feature-group-output", type=Path, default=Path("results/tables/phase20_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/101_build_gse243639_phase20_celltype_features.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    phase85 = load_phase85_module()
    phase85.configure_logging(args.log_file)
    clinical = phase85.read_clinical(args.clinical)
    annotations, celltype_counts, sample_counts, join_column, normalization_rule = phase85.read_annotations(args.annotations)
    (
        genes,
        celltypes,
        sample_gene_sum,
        sample_gene_detection,
        sample_celltype_gene_sum,
        sample_celltype_gene_detection,
        unmatched_cells,
        matched_cells,
    ) = phase85.stream_expression(args.expression, annotations, join_column, normalization_rule)
    unique_expression_cells = len(matched_cells | unmatched_cells)
    match_rate = len(matched_cells) / unique_expression_cells if unique_expression_cells else 0.0
    rows, fieldnames = phase85.build_rows(
        clinical,
        genes,
        celltypes,
        sample_gene_sum,
        sample_gene_detection,
        sample_celltype_gene_sum,
        sample_celltype_gene_detection,
        celltype_counts,
        sample_counts,
    )
    phase85.write_tsv(args.output, rows, fieldnames)
    phase85.write_tsv(args.schema_output, phase85.schema_alignment(fieldnames, args.phase5_schema), ["feature", "in_gse243639_celltype", "in_sea_ad_phase5", "status"])
    phase85.write_tsv(args.label_summary_output, phase85.label_summary(rows), ["label_field", "label", "sample_count"])
    phase85.write_tsv(
        args.feature_group_output,
        phase85.feature_group_counts(fieldnames, rows, match_rate, unmatched_cells),
        ["feature_group", "feature_count", "sample_rows", "annotation_match_rate", "unmatched_unique_expression_cells", "warning"],
    )
    print(f"Wrote {args.output}")
    print(f"Wrote {args.feature_group_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
