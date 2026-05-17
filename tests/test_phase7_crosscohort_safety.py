from pathlib import Path


PHASE7_SCRIPTS = [
    Path("scripts/30_prepare_external_validation_plan.py"),
    Path("scripts/31_sparse_external_gene_extraction.py"),
    Path("scripts/32_build_crosscohort_feature_tables.py"),
    Path("scripts/33_run_crosscohort_validation.py"),
    Path("scripts/34_generate_phase7_figures.py"),
    Path("scripts/35_generate_phase7_results_text.py"),
]


def test_phase7_scripts_avoid_disallowed_single_cell_workflows():
    for script in PHASE7_SCRIPTS:
        text = script.read_text(encoding="utf-8").lower()
        assert "import scanpy" not in text
        assert "read_h5ad" not in text
        assert "import anndata" not in text
        assert "scvi" not in text
        assert "scvelo" not in text
        assert "toarray(" not in text
        assert "todense(" not in text
        assert "np.array(" not in text
        assert "pivot_table" not in text


def test_external_sparse_extraction_is_manual_guarded_and_chunk_limited():
    text = Path("scripts/31_sparse_external_gene_extraction.py").read_text(encoding="utf-8")
    assert "RUN_MANUAL_EXTERNAL_EXTRACTION" in text
    assert "--execute" in text
    assert "--dry-run" in text
    assert "MAX_GENES = 64" in text
    assert "MAX_CHUNK_SIZE = 50000" in text
    assert "memory_limit_mb" in text
    assert "Dry run only. No external expression file was opened." in text


def test_phase7_expected_outputs_declared():
    combined = "\n".join(script.read_text(encoding="utf-8") for script in PHASE7_SCRIPTS)
    for expected in [
        "external_validation_gene_overlap.tsv",
        "external_validation_metadata_overlap.tsv",
        "crosscohort_donor_feature_table.tsv",
        "crosscohort_feature_overlap.tsv",
        "phase7_crosscohort_metrics.tsv",
        "phase7_generalization_summary.tsv",
        "phase7_results_summary.txt",
        "figure16_crosscohort_generalization.png",
        "figure17_cohort_transfer_performance.png",
        "figure18_feature_stability.png",
        "figure19_multicohort_neurofate_scores.png",
    ]:
        assert expected in combined


def test_crosscohort_validation_modes_are_declared():
    text = Path("scripts/33_run_crosscohort_validation.py").read_text(encoding="utf-8")
    for expected in [
        "train_sea_ad_test_external",
        "leave_one_cohort_out",
        "pooled_multicohort_training",
        "roc_auc_score",
        "average_precision_score",
        "balanced_accuracy_score",
        "brier_score_loss",
    ]:
        assert expected in text
