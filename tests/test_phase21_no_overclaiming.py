from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_no_overclaiming_audit_knows_axis_phrases():
    text = (ROOT / "scripts/54_no_overclaiming_audit.py").read_text(encoding="utf-8")
    for phrase in [
        "causal axis",
        "disease mechanism proven",
        "clinical biomarker",
        "definitive shared mechanism",
        "validated across diseases",
    ]:
        assert phrase in text
    for allowed in [
        "candidate shared axis",
        "preliminary disease-specific axis",
        "donor-level association",
        "exploratory cross-disease convergence",
    ]:
        assert allowed in text


def test_phase21_documents_are_conservative():
    strategy = (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8").lower()
    pnas_doc = (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8").lower()
    combined = strategy + "\n" + pnas_doc
    assert "not a clinical" in combined or "not clinical" in combined
    assert "cannot yet claim" in combined
    assert "candidate shared axis" in combined
    assert "preliminary disease-specific axis" in combined
    assert "validated across diseases" in combined
    assert "do not describe" in combined or "do not claim" in combined
