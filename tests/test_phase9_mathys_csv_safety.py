from pathlib import Path


PHASE9_SCRIPTS = [
    Path("scripts/40_inspect_mathys_csv_structure.py"),
    Path("scripts/41_extract_mathys_target_gene_panel.py"),
    Path("scripts/42_build_mathys_feature_table.py"),
    Path("scripts/43_run_mathys_external_validation.py"),
    Path("scripts/44_generate_phase9_figures.py"),
    Path("scripts/45_generate_phase9_results_text.py"),
]


def test_phase9_scripts_avoid_single_cell_workflows():
    for script in PHASE9_SCRIPTS:
        text = script.read_text(encoding="utf-8").lower()
        assert "import scanpy" not in text
        assert "read_h5ad" not in text
        assert "import anndata" not in text
        assert "import h5py" not in text
        assert "scvi" not in text
        assert "scvelo" not in text
        assert "umap" not in text
        assert "toarray(" not in text
        assert "todense(" not in text


def test_mathys_csv_inspector_declares_real_inputs_and_outputs():
    text = Path("scripts/40_inspect_mathys_csv_structure.py").read_text(encoding="utf-8")
    for expected in [
        "GSE138852_counts.csv.gz",
        "GSE138852_covariates.csv.gz",
        "mathys_csv_structure_summary.tsv",
        "mathys_covariates_preview.tsv",
        "oupSample.batchCond",
        "oupSample.cellType",
        "oupSample.cellType_batchCond",
        "oupSample.subclustID",
        "oupSample.subclustCond",
    ]:
        assert expected in text


def test_mathys_target_gene_extractor_supports_both_orientations():
    text = Path("scripts/41_extract_mathys_target_gene_panel.py").read_text(encoding="utf-8")
    for expected in [
        "genes_as_rows",
        "cells_as_rows",
        "infer_orientation",
        "extract_genes_as_rows",
        "extract_cells_as_rows",
        "mathys_sparse_gene_panel_expression.tsv.gz",
        "cell_id",
        "gene_symbol",
        "expression_value",
        "--dry-run",
    ]:
        assert expected in text


def test_phase9_outputs_declared():
    combined = "\n".join(script.read_text(encoding="utf-8") for script in PHASE9_SCRIPTS)
    for expected in [
        "phase9_mathys_external_validation_metrics.tsv",
        "phase9_mathys_external_predictions.tsv",
        "figure20_mathys_gene_overlap.png",
        "figure21_mathys_external_validation.png",
        "figure22_mathys_celltype_composition.png",
        "phase9_results_summary.txt",
    ]:
        assert expected in combined
