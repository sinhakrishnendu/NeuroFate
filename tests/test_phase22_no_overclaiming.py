from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_no_overclaiming_audit_blocks_phase22_forbidden_terms():
    text = (ROOT / "scripts/54_no_overclaiming_audit.py").read_text(encoding="utf-8")
    for phrase in [
        "causal axis",
        "proven mechanism",
        "clinical biomarker",
        "diagnostic",
        "validated across diseases",
        "definitive shared mechanism",
    ]:
        assert phrase in text
    assert "endpoint-locked" in text


def test_docs_say_phase22_supersedes_phase21_for_pnas_claims():
    docs = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 22 supersedes phase 21" in docs or "phase 22 replaces phase 21" in docs
    assert "largest-effect-across-label" in docs
    assert "clinical" in docs
    assert "causal" in docs
    assert "validated across diseases" in docs
    assert "do not claim" in docs
