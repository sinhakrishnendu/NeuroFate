from pathlib import Path


MODEL_SCRIPT = Path("scripts/23_run_phase5_models.py")
FIGURE_SCRIPT = Path("scripts/24_generate_phase5_figures.py")
TEXT_SCRIPT = Path("scripts/25_generate_phase5_results_text.py")


def test_phase5_model_script_uses_donor_table_only():
    text = MODEL_SCRIPT.read_text(encoding="utf-8")
    lowered = text.lower()
    assert "read_h5ad" not in lowered
    assert "import h5py" not in lowered
    assert "import scanpy" not in lowered
    assert "import anndata" not in lowered
    assert "import torch" not in lowered
    assert "phase5_donor_feature_table.tsv" in text
    assert "sparse_gene_panel_expression" not in lowered


def test_phase5_model_script_declares_tasks_models_and_outputs():
    text = MODEL_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "dementia_vs_reference",
        "high_vs_low_ad_neuropathology",
        "apoe_risk_prediction",
        "mixed_pathology_burden",
        "logistic_regression",
        "elastic_net",
        "random_forest_baseline",
        "gradient_boosting_baseline",
        "train_test_split",
        "StratifiedKFold",
        "phase5_model_metrics.tsv",
        "phase5_feature_importance.tsv",
        "phase5_neurofate_scores.tsv",
        "neurofate_neurodegeneration_risk_score",
    ]:
        assert expected in text


def test_phase5_model_script_excludes_labels_from_features():
    text = MODEL_SCRIPT.read_text(encoding="utf-8")
    assert 'LABEL_PREFIX = "label__"' in text
    assert "FEATURE_PREFIXES" in text
    assert "labels are derived" in text.lower()


def test_phase5_figure_and_text_scripts_are_summary_only():
    for script in [FIGURE_SCRIPT, TEXT_SCRIPT]:
        text = script.read_text(encoding="utf-8").lower()
        assert "read_h5ad" not in text
        assert "sparse_gene_panel_expression" not in text
        assert "import h5py" not in text
        assert "import scanpy" not in text
    figure_text = FIGURE_SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "figure9_model_performance.png",
        "figure10_feature_importance.png",
        "figure11_neurofate_score_distribution.png",
        "figure12_donor_risk_heatmap.png",
    ]:
        assert filename in figure_text
    assert "phase5_results_summary.txt" in TEXT_SCRIPT.read_text(encoding="utf-8")
