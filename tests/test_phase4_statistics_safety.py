from pathlib import Path


STAT_SCRIPT = Path("scripts/19_compute_phase4_statistics.py")
TEXT_SCRIPT = Path("scripts/21_generate_phase4_results_text.py")


def test_phase4_statistics_uses_allowed_inputs_only():
    text = STAT_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "import h5py" not in lowered
    assert "import anndata" not in lowered
    assert "import torch" not in lowered
    assert "sparse_gene_panel_expression.tsv.gz" in text
    assert "sea_ad_obs_metadata_decoded.tsv" in text


def test_phase4_statistics_prohibits_dense_and_pipeline_terms():
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
        "scvelo",
        "cellrank",
    ]:
        assert forbidden not in text


def test_phase4_statistics_declares_methods_and_outputs():
    text = STAT_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "benjamini_hochberg",
        "spearman_rank_trend",
        "kruskal_rank_association",
        "donor_values_for_gene",
        "phase4_gene_statistics.tsv",
        "phase4_celltype_vulnerability.tsv",
        "phase4_apoe_analysis.tsv",
        "phase4_mixed_pathology.tsv",
        "phase4_composite_indices.tsv",
        "max-row-chunk",
    ]:
        assert expected in text


def test_phase4_results_text_reads_summary_tables_only():
    text = TEXT_SCRIPT.read_text(encoding="utf-8").lower()
    assert "read_h5ad" not in text
    assert "sparse_gene_panel_expression" not in text
    assert "phase4_results_summary.txt" in text
    for filename in [
        "phase4_gene_statistics.tsv",
        "phase4_celltype_vulnerability.tsv",
        "phase4_apoe_analysis.tsv",
        "phase4_mixed_pathology.tsv",
    ]:
        assert filename in text
