from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/86_run_gse243639_celltype_pd_validation.py"
COMPARE = ROOT / "scripts/87_compare_phase16_phase17_pd_validation.py"


def test_pd_validation_is_sample_level_and_interpretable() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "sample-level" in text
    assert "LogisticRegression" in text
    assert "RandomForestClassifier" in text
    assert "LeaveOneOut" in text
    assert "RepeatedStratifiedKFold" in text


def test_pd_validation_has_robustness_controls() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for token in ["permutation_control", "bootstrap_ci", "empirical_permutation_pvalue", "auroc_ci_low", "auroc_ci_high"]:
        assert token in text


def test_pd_validation_reliability_categories_are_conservative() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for token in ["moderate_pd_internal_validation", "preliminary_pd_internal_signal", "weak_pd_signal"]:
        assert token in text
    assert "observed_auroc >= 0.70" in text
    assert "empirical_p <= 0.05" in text


def test_pd_validation_avoids_forbidden_workloads() -> None:
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "read_h5ad", "anndata", "tensorflow", "keras", "umap", "leiden"]:
        assert forbidden not in text


def test_phase16_phase17_comparison_exists() -> None:
    text = COMPARE.read_text(encoding="utf-8")
    assert "phase16_global_sample_features" in text
    assert "phase17_celltype_aware_features" in text
    assert "phase17_pd_validation_comparison.tsv" in text
