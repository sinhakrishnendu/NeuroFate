from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase26_docs_keep_replication_claims_conservative():
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 26" in combined
    assert "pd/pdd" in combined
    assert "claims remain" in combined or "candidate or preliminary" in combined
    for phrase in [
        "phase 26 proves",
        "phase 26 is clinical-grade",
        "phase 26 is a diagnostic tool",
        "phase 26 proves causality",
        "phase 26 is validated across diseases",
    ]:
        assert phrase not in combined


def test_phase26_report_script_uses_limitations_language():
    text = (ROOT / "scripts/137_generate_phase26_gse184950_replication_report.py").read_text(encoding="utf-8").lower()
    assert "not a clinical validation" in text
    assert "no causal mechanism" in text
    assert "candidate or preliminary" in text
