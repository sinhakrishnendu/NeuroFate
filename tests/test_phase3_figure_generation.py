from pathlib import Path


FIGURE_SCRIPT = Path("scripts/17_generate_phase3_figures.py")


def test_phase3_figure_script_uses_matplotlib_only_for_figures():
    text = FIGURE_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "import matplotlib" in lowered
    assert "import seaborn" not in lowered
    assert "import scanpy" not in lowered
    assert "read_h5ad" not in lowered
    assert "import h5py" not in lowered
    assert "import pandas" not in lowered


def test_phase3_figure_script_declares_expected_outputs():
    text = FIGURE_SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "figure1_celltype_composition.png",
        "figure2_microglial_activation.png",
        "figure3_neurodegeneration_signatures.png",
        "figure4_ad_pathology_gene_trends.png",
    ]:
        assert filename in text


def test_phase3_figure_script_avoids_pipeline_terms():
    text = FIGURE_SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["umap", "leiden", "neighbors", "scvi", "cellrank", "scvelo", "fit("]:
        assert forbidden not in text
