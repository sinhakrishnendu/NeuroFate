from pathlib import Path


SCRIPT = Path("scripts/11_extract_sea_ad_metadata_only.py")


def test_extractor_does_not_import_scanpy_or_anndata_reader():
    text = SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "anndata" not in lowered
    assert "AnnData" not in text


def test_extractor_has_x_access_guard_without_direct_x_indexing():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "FORBIDDEN_ROOT_KEY" in text
    assert "RuntimeError" in text
    assert '["X"]' not in text
    assert "['X']" not in text


def test_extractor_declares_expected_outputs():
    text = SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "sea_ad_obs_metadata_minimal.tsv",
        "sea_ad_var_genes.tsv",
        "sea_ad_metadata_summary.tsv",
        "table1_sea_ad_cohort_cell_summary.tsv",
        "11_extract_sea_ad_metadata_only.log",
    ]:
        assert filename in text


def test_extractor_mentions_no_matrix_workflows():
    text = SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["normalize_total", "pca", "umap", "leiden", "neighbors", "fit("]:
        assert forbidden not in text
