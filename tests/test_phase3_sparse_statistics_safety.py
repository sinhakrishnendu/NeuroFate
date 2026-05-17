from pathlib import Path


STAT_SCRIPT = Path("scripts/16_compute_sparse_expression_statistics.py")
TEXT_SCRIPT = Path("scripts/18_generate_phase3_results_text.py")


def test_phase3_statistics_uses_sparse_tsv_inputs_only():
    text = STAT_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "import h5py" not in lowered
    assert "import pandas" not in lowered
    assert "sparse_gene_panel_expression.tsv.gz" in text
    assert "sea_ad_obs_metadata_decoded.tsv" in text


def test_phase3_statistics_prohibits_dense_and_pipeline_terms():
    text = STAT_SCRIPT.read_text(encoding="utf-8").lower()
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
        "cellrank",
        "scvelo",
    ]:
        assert forbidden not in text


def test_phase3_statistics_declares_expected_outputs_and_chunking():
    text = STAT_SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "gene_by_celltype_summary.tsv",
        "gene_by_ad_pathology.tsv",
        "gene_by_cognitive_status.tsv",
        "microglial_activation_signature.tsv",
        "astrocyte_stress_signature.tsv",
        "neuronal_signature_summary.tsv",
        "neurodegeneration_signature_summary.tsv",
    ]:
        assert filename in text
    assert "max-row-chunk" in text
    assert "DEFAULT_MAX_ROW_CHUNK" in text


def test_phase3_results_text_reads_summary_tables_only():
    text = TEXT_SCRIPT.read_text(encoding="utf-8").lower()
    assert "import scanpy" not in text
    assert "read_h5ad" not in text
    assert "sparse_gene_panel_expression" not in text
    assert "phase3_results_summary.txt" in text
