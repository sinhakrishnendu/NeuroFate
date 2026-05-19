from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "scripts/94_generate_phase18_gse243639_repair_report.py"
CLAIMS = ROOT / "scripts/65_build_claim_strength_table.py"
README = ROOT / "README.md"
DOC = ROOT / "docs/external_validation_expansion.md"


def test_phase18_report_avoids_forbidden_claim_language() -> None:
    text = REPORT.read_text(encoding="utf-8").lower()
    forbidden = [
        "clinical-grade",
        "diagnostic tool",
        "validated across diseases",
        "causal",
        "foundation model",
    ]
    for phrase in forbidden:
        assert phrase not in text


def test_phase18_claim_strength_supersedes_phase17_when_available() -> None:
    text = CLAIMS.read_text(encoding="utf-8")
    assert "phase18_gse243639_celltype_validation_metrics.tsv" in text
    assert "phase18_claim_strength_delta.tsv" in text
    assert "superseded_by_phase18" in text
    assert "technical_failure_annotation_join" in text
    assert "no existing AD claim was modified" in text


def test_phase18_docs_are_conservative() -> None:
    combined = "\n".join([README.read_text(encoding="utf-8"), DOC.read_text(encoding="utf-8")]).lower()
    assert "phase 17 should be treated as a technical audit result" in combined
    assert "final pd interpretation should use phase 18" in combined
    assert "technical_failure_annotation_join" in combined
    assert "not a biological conclusion" in combined


def test_phase18_report_and_docs_do_not_claim_recomputed_umap_or_h5ad() -> None:
    combined = "\n".join([REPORT.read_text(encoding="utf-8"), README.read_text(encoding="utf-8"), DOC.read_text(encoding="utf-8")]).lower()
    for forbidden in ["umap.fit", "fit_transform", "scanpy workflow", "annData object creation"]:
        assert forbidden.lower() not in combined
