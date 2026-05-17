from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/63_classify_evidence_strength.py"


def test_evidence_categories_are_declared() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for category in [
        "strong_internal",
        "moderate_internal",
        "preliminary_external",
        "insufficient",
        "failed_or_unstable",
    ]:
        assert category in text


def test_evidence_classifier_uses_required_signals() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for token in [
        "sample_size",
        "auroc_mean",
        "auroc_sd",
        "empirical_pvalue",
        "ablation_consistency",
        "external_validation_available",
        "no_overclaiming_high_flags",
    ]:
        assert token in text


def test_evidence_classifier_strong_internal_path() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "sample_size >= 40" in text
    assert "auroc_mean >= 0.75" in text
    assert "auroc_sd <= 0.05" in text
    assert "empirical_pvalue <= 0.05" in text
    assert 'category = "strong_internal"' in text
