from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/85_build_gse243639_celltype_feature_table.py"


def test_celltype_feature_table_is_sample_level() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert '"unit_type": "sample"' in text
    assert "dataset_unit_id" in text
    assert "sample_id" in text
    assert "label__diagnosis_binary" in text


def test_celltype_feature_definitions_are_present() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for token in [
        "gene_mean__",
        "gene_detection__",
        "cell_fraction__",
        "celltype_gene_mean__",
        "celltype_gene_detection__",
        "microglial_activation_index",
        "astrocyte_stress_index",
        "neuronal_vulnerability_index",
        "myelin_oligodendrocyte_index",
        "synuclein_axis_index",
    ]:
        assert token in text


def test_celltype_feature_table_uses_existing_annotation_map() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "gse243639_cell_annotation_map.tsv" in text
    assert "read_annotations" in text
    assert "stream_expression" in text


def test_celltype_feature_table_avoids_scanpy_and_h5ad() -> None:
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "h5py", "fit_transform", "leiden"]:
        assert forbidden not in text
