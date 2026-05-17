from pathlib import Path


FEATURE_SCRIPT = Path("scripts/22_build_donor_feature_table.py")


def test_phase5_feature_table_uses_allowed_inputs_only():
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "import h5py" not in lowered
    assert "import anndata" not in lowered
    assert "import torch" not in lowered
    assert "sparse_gene_panel_expression.tsv.gz" in text
    assert "sea_ad_obs_metadata_decoded.tsv" in text


def test_phase5_feature_table_declares_donor_aggregation_outputs():
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "phase5_donor_feature_table.tsv",
        "donor_gene_sum",
        "donor_gene_nonzero",
        "donor_index_sum",
        "donor_celltype_index_sum",
        "cell_fraction__",
        "celltype_index__",
        "max-row-chunk",
    ]:
        assert expected in text


def test_phase5_feature_table_avoids_dense_and_pipeline_terms():
    text = FEATURE_SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in [
        "toarray(",
        "todense(",
        "np.array",
        "pivot(",
        "pivot_table",
        "umap",
        "leiden",
        "neighbors",
        "scvi",
        "scvelo",
        "cellrank",
    ]:
        assert forbidden not in text
