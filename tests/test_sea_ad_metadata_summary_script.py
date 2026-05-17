from pathlib import Path


SCRIPT = Path("scripts/12_summarize_sea_ad_metadata.py")


def test_summary_script_reads_tsv_only():
    text = SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import h5py" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "anndata" not in lowered
    assert "AnnData" not in text


def test_summary_script_declares_expected_outputs():
    text = SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "sea_ad_obs_metadata_minimal.tsv",
        "sea_ad_donor_summary.tsv",
        "sea_ad_celltype_by_ad_pathology.tsv",
        "sea_ad_celltype_by_cognitive_status.tsv",
        "12_summarize_sea_ad_metadata.log",
    ]:
        assert filename in text


def test_summary_script_has_no_matrix_workflows():
    text = SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["normalize_total", "pca", "umap", "leiden", "neighbors", "fit("]:
        assert forbidden not in text
