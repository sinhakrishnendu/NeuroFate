from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase27_docs_are_conservative():
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 27" in combined
    assert "direction-only" in combined or "directionally consistent" in combined
    assert "not enough" in combined or "not treated as replication" in combined
    for phrase in [
        "phase 27 proves",
        "phase 27 validates",
        "phase 27 establishes a clinical biomarker",
        "phase 27 validates a diagnostic axis",
        "phase 27 validates a shared neurodegeneration axis",
        "phase 27 proves a replicated ad/pd mechanism",
    ]:
        assert phrase not in combined


def test_no_overclaiming_audit_contains_phase27_forbidden_claims():
    text = (ROOT / "scripts/54_no_overclaiming_audit.py").read_text(encoding="utf-8").lower()
    assert "replicated ad/pd mechanism" in text
    assert "validated shared neurodegeneration axis" in text
    assert "diagnostic axis" in text
    assert "causal mechanism" in text
