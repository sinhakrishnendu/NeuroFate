from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATION_SCRIPT = ROOT / "scripts/80_run_gse243639_pd_external_validation.py"
REPORT_SCRIPT = ROOT / "scripts/82_generate_phase16_gse243639_report.py"
MULTI_EXTERNAL_SCRIPT = ROOT / "scripts/75_run_multi_external_validation.py"
CLAIM_SCRIPT = ROOT / "scripts/65_build_claim_strength_table.py"


def test_validation_script_uses_conservative_pd_reliability_categories() -> None:
    text = VALIDATION_SCRIPT.read_text(encoding="utf-8")
    assert "moderate_pd_internal_validation" in text
    assert "preliminary_cross_disease_feature_transfer" in text
    assert "insufficient_cross_disease_validation" in text
    assert "not_biologically_equivalent" in text


def test_validation_script_does_not_use_deep_learning_or_scanpy() -> None:
    text = VALIDATION_SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["torch", "tensorflow", "keras", "scanpy", "read_h5ad", "UMAP", "Leiden"]:
        assert forbidden not in text


def test_multi_external_and_claim_strength_support_phase16() -> None:
    multi_text = MULTI_EXTERNAL_SCRIPT.read_text(encoding="utf-8")
    claim_text = CLAIM_SCRIPT.read_text(encoding="utf-8")
    assert "phase16_gse243639_feature_table.tsv" in multi_text
    assert "phase16_gse243639_external_validation_metrics.tsv" in claim_text
    assert "moderate_pd_internal_validation" in multi_text
    assert "moderate_pd_internal_validation" in claim_text


def test_phase16_report_avoids_forbidden_overclaiming_phrases() -> None:
    text = REPORT_SCRIPT.read_text(encoding="utf-8").lower()
    forbidden = [
        "clinical-grade",
        "diagnostic tool",
        "validated across diseases",
        "causal",
        "foundation model",
    ]
    for phrase in forbidden:
        assert phrase not in text
