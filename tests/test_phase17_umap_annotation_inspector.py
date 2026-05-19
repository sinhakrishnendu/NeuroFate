from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INSPECTOR = ROOT / "scripts/83_inspect_gse243639_umap_annotations.py"
MAPPER = ROOT / "scripts/84_build_gse243639_cell_annotation_map.py"


def test_umap_inspector_reads_annotations_without_scanpy() -> None:
    text = INSPECTOR.read_text(encoding="utf-8")
    assert "openpyxl" in text
    for forbidden in ["scanpy", "anndata", "read_h5ad", "fit_transform", "neighbors", "leiden"]:
        assert forbidden not in text.lower()


def test_umap_inspector_detects_annotation_and_coordinate_columns() -> None:
    text = INSPECTOR.read_text(encoding="utf-8")
    for token in ["candidate_cell_id", "candidate_sample_id", "candidate_celltype_or_cluster_annotation", "candidate_umap_coordinate"]:
        assert token in text
    assert "preview-output" in text


def test_annotation_mapper_preserves_sample_and_cell_type_fields() -> None:
    text = MAPPER.read_text(encoding="utf-8")
    for token in ["cell_id", "sample_id", "cell_type", "chosen_annotation_column", "annotation_candidate_columns"]:
        assert token in text
    assert "unmatched_cells_total" in text


def test_annotation_mapper_does_not_recompute_embeddings() -> None:
    text = MAPPER.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "fit_transform", "umap(", "pca", "leiden", "neighbors"]:
        assert forbidden not in text
