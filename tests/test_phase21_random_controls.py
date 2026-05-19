from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/109_run_axis_randomization_controls.py"


def test_random_controls_use_feature_universe_only():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "feature_universe" in text
    assert "is_feature_column" in text
    assert "LABEL_HINTS" in text
    assert "rng.sample(universe" in text
    assert "phase21_random_axis_controls.tsv" in text
    assert "phase21_axis_empirical_pvalues.tsv" in text


def test_random_controls_are_not_heavy_or_single_cell():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    forbidden = ["scanpy", "read_h5ad", ".h5ad", "umap", "leiden", "pca", "scvi", "torch"]
    for phrase in forbidden:
        assert phrase not in text
