from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/114_run_endpoint_locked_random_axis_controls.py"


def test_random_controls_use_same_endpoint_statistic():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "endpoint_effect" in text
    assert "same locked endpoint statistic" in text
    assert "curated_effect" in text
    assert "random_effect_size" in text
    assert "empirical_pvalue" in text
    assert "abs(effect) >= abs(curated_effect)" in text


def test_random_controls_exclude_label_metadata_columns():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "FEATURE_PREFIXES" in text
    assert "LABEL_HINTS" in text
    assert "is_feature_column" in text
    assert "not any(hint in lowered for hint in LABEL_HINTS)" in text
    assert "sample_id" in text
    assert "donor_id" in text


def test_random_controls_do_not_use_single_cell_frameworks():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "read_h5ad", ".h5ad", "umap", "leiden", "pca", "scvi", "torch"]:
        assert forbidden not in text
