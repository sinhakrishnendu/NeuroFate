from pathlib import Path


FIGURE_SCRIPT = Path("scripts/20_generate_phase4_figures.py")


def test_phase4_figure_script_uses_matplotlib_summary_tables_only():
    text = FIGURE_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "matplotlib" in lowered
    assert "import pandas" not in lowered
    assert "import seaborn" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "import h5py" not in lowered
    for filename in [
        "phase4_gene_statistics.tsv",
        "phase4_apoe_analysis.tsv",
        "phase4_celltype_vulnerability.tsv",
        "phase4_composite_indices.tsv",
    ]:
        assert filename in text


def test_phase4_figure_script_declares_expected_outputs():
    text = FIGURE_SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "figure5_braak_associations.png",
        "figure6_apoe_microglia.png",
        "figure7_celltype_vulnerability_heatmap.png",
        "figure8_composite_indices.png",
    ]:
        assert filename in text


def test_phase4_figure_script_avoids_forbidden_workflows():
    text = FIGURE_SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in [
        "umap",
        "leiden",
        "neighbors",
        "scvi",
        "scvelo",
        "cellrank",
        "fit(",
        "toarray(",
        "todense(",
    ]:
        assert forbidden not in text
