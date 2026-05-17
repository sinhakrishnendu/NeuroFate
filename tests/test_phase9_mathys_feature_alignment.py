from pathlib import Path


FEATURE_SCRIPT = Path("scripts/42_build_mathys_feature_table.py")
VALIDATION_SCRIPT = Path("scripts/43_run_mathys_external_validation.py")


def test_mathys_feature_builder_maps_covariates_to_neurofate_fields():
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "MATHYS_CELLTYPE_FIELD",
        "oupSample.cellType",
        "oupSample.batchCond",
        "oupSample.subclustID",
        "oupSample.subclustCond",
        "diagnosis_from_covariates",
        "infer_sample_id",
        "label__diagnosis",
        "label__Overall_AD_neuropathological_Change",
        "label__subcluster",
    ]:
        assert expected in text


def test_mathys_feature_builder_aligns_to_phase5_schema():
    text = FEATURE_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "phase5_donor_feature_table.tsv",
        "mathys_2019_phase5_donor_feature_table.tsv",
        "mathys_2019_feature_schema_alignment.tsv",
        "mathys_2019_label_summary.tsv",
        "gene_mean__",
        "gene_detection__",
        "index__",
        "cell_fraction__",
        "celltype_index__",
    ]:
        assert expected in text


def test_mathys_validation_uses_shared_feature_columns_and_baseline_model():
    text = VALIDATION_SCRIPT.read_text(encoding="utf-8")
    for expected in [
        "phase5_donor_feature_table.tsv",
        "mathys_2019_phase5_donor_feature_table.tsv",
        "feature_columns",
        "set.intersection",
        "train_sea_ad_test_mathys",
        "mathys_internal_diagnostic",
        "LogisticRegression",
        "roc_auc_score",
        "average_precision_score",
        "balanced_accuracy_score",
        "brier_score_loss",
    ]:
        assert expected in text


def test_phase9_validation_does_not_use_deep_training_or_h5ad():
    combined = FEATURE_SCRIPT.read_text(encoding="utf-8").lower() + VALIDATION_SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in [
        "torch",
        "tensorflow",
        "keras",
        "read_h5ad",
        "import h5py",
        "import scanpy",
        "pca",
        "umap",
        "leiden",
    ]:
        assert forbidden not in combined
