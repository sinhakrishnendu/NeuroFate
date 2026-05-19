from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "scripts/89_generate_phase17_gse243639_report.py"
AUDIT = ROOT / "scripts/54_no_overclaiming_audit.py"
CLAIMS = ROOT / "scripts/65_build_claim_strength_table.py"


def test_phase17_report_avoids_forbidden_claim_phrases() -> None:
    text = REPORT.read_text(encoding="utf-8").lower()
    forbidden = [
        "clinical pd prediction",
        "diagnostic pd classifier",
        "validated pd biomarker",
        "cross-disease diagnostic transfer",
        "clinical-grade",
        "foundation model",
    ]
    for phrase in forbidden:
        assert phrase not in text


def test_no_overclaiming_audit_knows_phase17_context() -> None:
    text = AUDIT.read_text(encoding="utf-8")
    assert "cell-type-aware pd signal" in text
    assert "preliminary pd internal signal" in text
    assert "clinical pd prediction" in text
    assert "diagnostic pd classifier" in text
    assert "validated pd biomarker" in text
    assert "cross-disease diagnostic transfer" in text


def test_claim_strength_reads_phase17_metrics_without_ad_upgrade() -> None:
    text = CLAIMS.read_text(encoding="utf-8")
    assert "phase17_gse243639_celltype_validation_metrics.tsv" in text
    assert "phase17_claim_strength_delta.tsv" in text
    assert "no existing AD claim was modified" in text
    assert "moderate_pd_internal_validation" in text
    assert "preliminary_pd_internal_signal" in text


def test_docs_phrase_phase17_conservatively() -> None:
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "does not recompute umap" in combined
    assert "sample-level" in combined
    assert "exploratory pd signal" in combined
    assert "direct ad-to-pd disease-label transfer" in combined
