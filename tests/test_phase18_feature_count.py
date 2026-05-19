from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FEATURE_SCRIPT = ROOT / "scripts/85_build_gse243639_celltype_feature_table.py"
VALIDATION_SCRIPT = ROOT / "scripts/91_run_gse243639_repaired_celltype_pd_validation.py"


def test_repaired_feature_table_preserves_global_and_celltype_features() -> None:
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    for token in [
        "gene_mean__",
        "gene_detection__",
        "cell_fraction__",
        "celltype_gene_mean__",
        "celltype_gene_detection__",
        "phase18_gse243639_feature_group_counts.tsv",
    ]:
        assert token in text


def test_repaired_feature_table_counts_unmatched_unique_cells_not_sparse_rows() -> None:
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    assert "unmatched_cells: set[str]" in text
    assert "unmatched_cells.add(cell_id)" in text
    assert "unmatched_unique_expression_cells" in text


def test_repaired_feature_table_is_sample_level_only() -> None:
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    assert '"unit_type": "sample"' in text
    assert "dataset_unit_id" in text
    assert "sample_rows" in text


def test_repaired_validation_marks_technical_failure_when_features_remain_too_small() -> None:
    text = VALIDATION_SCRIPT.read_text(encoding="utf-8")
    assert "MIN_REPAIRED_FEATURES" in text
    assert "MIN_MATCH_RATE" in text
    assert "technical_failure_annotation_join" in text
    assert "phase18_gse243639_celltype_validation_metrics.tsv" in text


def test_repaired_scripts_avoid_forbidden_workloads() -> None:
    combined = "\n".join(
        [
            FEATURE_SCRIPT.read_text(encoding="utf-8").lower(),
            VALIDATION_SCRIPT.read_text(encoding="utf-8").lower(),
        ]
    )
    for forbidden in ["scanpy", "anndata", "read_h5ad", "fit_transform", "umap", "leiden", "tensorflow", "keras"]:
        assert forbidden not in combined
