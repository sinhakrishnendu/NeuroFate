from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/119_build_axis_scores_from_sample_matrix.py"


def test_sample_matrix_builder_excludes_labels_and_metadata():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "LABEL_HINTS" in text
    assert "is_label_or_metadata_column" in text
    assert "diagnosis" in text
    assert "sample_id" in text
    assert "donor_id" in text
    assert "not is_label_or_metadata_column" in text


def test_sample_matrix_builder_supports_orientations_and_size_limit():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "genes_rows" in text
    assert "samples_rows" in text
    assert "--max-file-size-mb" in text
    assert "--allow-large-matrix" in text
    assert "read_genes_rows" in text
    assert "read_samples_rows" in text


def test_sample_matrix_builder_avoids_single_cell_tools():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "read_h5ad", ".h5ad", "umap", "leiden", "cluster", "torch", "sra"]:
        assert forbidden not in text
