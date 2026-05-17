from pathlib import Path


CHECK_SCRIPT = Path("scripts/26_check_mps_device.py")
TRAIN_SCRIPT = Path("scripts/27_train_neurofate_mps_model.py")
FIGURE_SCRIPT = Path("scripts/28_generate_phase6_figures.py")
TEXT_SCRIPT = Path("scripts/29_generate_phase6_results_text.py")


def test_phase6_scripts_do_not_access_single_cell_files_or_workflows():
    for script in [CHECK_SCRIPT, TRAIN_SCRIPT, FIGURE_SCRIPT, TEXT_SCRIPT]:
        text = script.read_text(encoding="utf-8").lower()
        assert "read_h5ad" not in text
        assert "import h5py" not in text
        assert "import scanpy" not in text
        assert "import anndata" not in text
        assert ".h5ad" not in text
        assert "sparse_gene_panel_expression" not in text
        assert "scvi" not in text
        assert "scvelo" not in text


def test_mps_device_logic_exists():
    text = CHECK_SCRIPT.read_text(encoding="utf-8")
    train_text = TRAIN_SCRIPT.read_text(encoding="utf-8")
    assert "torch.backends.mps.is_available()" in text
    assert 'device = "mps" if mps_available else "cpu"' in text
    assert "torch.backends.mps.is_available()" in train_text
    assert 'torch.device("mps")' in train_text
    assert 'torch.device("cpu")' in train_text


def test_phase6_training_uses_donor_feature_table_and_saves_models():
    text = TRAIN_SCRIPT.read_text(encoding="utf-8")
    assert "results/tables/phase5_donor_feature_table.tsv" in text
    assert "results/models" in text
    assert "neurofate_mps_{task_id}.pt" in text
    assert "torch.save" in text
    for expected in [
        "phase6_mps_model_metrics.tsv",
        "phase6_mps_training_log.tsv",
        "phase6_mps_predictions.tsv",
        "BCEWithLogitsLoss",
        "AdamW",
        "Early",
    ]:
        assert expected.lower() in text.lower()


def test_phase6_figures_and_text_outputs_declared():
    figure_text = FIGURE_SCRIPT.read_text(encoding="utf-8")
    for filename in [
        "figure13_mps_model_performance.png",
        "figure14_mps_training_curves.png",
        "figure15_mps_prediction_distribution.png",
    ]:
        assert filename in figure_text
    assert "phase6_results_summary.txt" in TEXT_SCRIPT.read_text(encoding="utf-8")
