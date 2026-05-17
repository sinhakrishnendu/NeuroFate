from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/54_no_overclaiming_audit.py"


def test_no_overclaiming_audit_flags_required_phrases() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for phrase in [
        "validated across cohorts",
        "externally validated",
        "clinical-grade",
        "diagnostic",
        "causal",
        "foundation model",
        "state-of-the-art",
        "generalizable across diseases",
        "biomarker",
        "patient-level diagnosis",
    ]:
        assert phrase in text


def test_no_overclaiming_audit_allows_cautious_contexts() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for phrase in [
        "preliminary external feasibility",
        "research software",
        "candidate biomarker-like signal",
    ]:
        assert phrase in text


def test_validated_across_cohorts_not_used_as_unqualified_claim() -> None:
    interpretation = (ROOT / "RESULTS_INTERPRETATION.md").read_text(encoding="utf-8").lower()
    readme = (ROOT / "README.md").read_text(encoding="utf-8").lower()
    assert "validated across cohorts" not in interpretation
    assert "validated across cohorts" not in readme


def test_claim_language_matrix_exists() -> None:
    matrix = ROOT / "results/reports/claim_language_matrix.tsv"
    assert matrix.exists()
    text = matrix.read_text(encoding="utf-8")
    assert "internal prediction" in text
    assert "Apple Silicon optimization" in text
