from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/70_inspect_external_dataset_files.py"


def test_external_file_inspector_avoids_scanpy_and_matrix_readers() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    forbidden = ["scanpy", "read_h5ad", "anndata", "h5py.File", "pandas.read", "polars.read"]
    for token in forbidden:
        assert token not in text


def test_external_file_inspector_detects_required_formats() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for token in [".h5ad", ".h5", ".mtx", ".csv", ".tsv", ".txt", ".rds", ".loom"]:
        assert token in text
    for token in ["count_matrix_candidate", "metadata_or_covariates", "gene_or_feature_file"]:
        assert token in text


def test_external_file_inspector_flags_large_files_without_opening_contents() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "LARGE_FILE_MB" in text
    assert "large_file_flag" in text
    assert ".stat().st_size" in text
