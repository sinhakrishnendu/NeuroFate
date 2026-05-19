from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FEATURE_BUILDER = ROOT / "scripts/85_build_gse243639_celltype_feature_table.py"
PHASE20_WRAPPER = ROOT / "scripts/101_build_gse243639_phase20_celltype_features.py"


def load_feature_module():
    spec = importlib.util.spec_from_file_location("phase20_feature_builder", FEATURE_BUILDER)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_script85_joins_expression_cell_id_to_safe_map_cell_id_expression(tmp_path: Path) -> None:
    module = load_feature_module()
    annotations = tmp_path / "safe.tsv"
    fields = [
        "cell_id_expression",
        "cell_id_annotation_original",
        "cell_id_annotation_normalized",
        "normalized_cell_id",
        "barcode_core",
        "sample_id",
        "cell_type",
        "cluster_id",
        "annotation_source_sheet",
        "annotation_column_used",
        "match_status",
        "normalization_rule",
    ]
    with annotations.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "cell_id_expression": "s.1_AAA.1",
                "cell_id_annotation_original": "s-1_AAA-1",
                "cell_id_annotation_normalized": "s.1_AAA.1",
                "normalized_cell_id": "s.1_AAA.1",
                "barcode_core": "AAA",
                "sample_id": "s.1",
                "cell_type": "cluster_1",
                "cluster_id": "1",
                "annotation_source_sheet": "Sheet1",
                "annotation_column_used": "CLUSTER",
                "match_status": "matched",
                "normalization_rule": "replace_dash_with_dot",
            }
        )
    expression = tmp_path / "expr.tsv.gz"
    with gzip.open(expression, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_id", "sample_id", "gene_symbol", "expression_value"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"cell_id": "s.1_AAA.1", "sample_id": "s.1", "gene_symbol": "SNCA", "expression_value": "1"})
        writer.writerow({"cell_id": "s.1_AAA.1", "sample_id": "s.1", "gene_symbol": "APOE", "expression_value": "2"})
    annotations_by_id, celltype_counts, sample_counts, join_column, rule = module.read_annotations(annotations)
    result = module.stream_expression(expression, annotations_by_id, join_column, rule)
    unmatched_cells = result[-2]
    matched_cells = result[-1]
    assert join_column == "cell_id_expression"
    assert rule == "replace_dash_with_dot"
    assert matched_cells == {"s.1_AAA.1"}
    assert unmatched_cells == set()
    assert celltype_counts[("s.1", "cluster_1")] == 1
    assert sample_counts["s.1"] == 1


def test_annotation_match_rate_is_unique_cell_based() -> None:
    text = FEATURE_BUILDER.read_text(encoding="utf-8")
    assert "matched_cells" in text
    assert "unmatched_cells" in text
    assert "unique_expression_cells = len(matched_cells | unmatched_cells)" in text
    assert "len(matched_cells) / unique_expression_cells" in text


def test_phase20_wrapper_uses_safe_map_outputs() -> None:
    text = PHASE20_WRAPPER.read_text(encoding="utf-8")
    assert "gse243639_safe_cell_annotation_map.tsv" in text
    assert "phase20_gse243639_celltype_feature_table.tsv" in text
    assert "phase20_gse243639_feature_group_counts.tsv" in text
