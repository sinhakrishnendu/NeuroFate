from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
DOC = ROOT / "docs/external_validation_expansion.md"
REPORT = ROOT / "scripts/105_generate_phase20_gse243639_report.py"
FIGURES = ROOT / "scripts/104_generate_phase20_gse243639_figures.py"


def test_phase20_docs_are_conservative() -> None:
    combined = "\n".join([README.read_text(encoding="utf-8"), DOC.read_text(encoding="utf-8")]).lower()
    assert "phase 19 found safe normalized id linkage" in combined
    assert "phase 20 repairs" in combined
    assert "only phase 20 should be used" in combined
    assert "phase 16 remains the valid global sample-level pd extension" in combined


def test_phase20_report_avoids_forbidden_claims() -> None:
    text = REPORT.read_text(encoding="utf-8").lower()
    for forbidden in [
        "clinical-grade",
        "diagnostic tool",
        "validated across diseases",
        "causal",
        "foundation model",
        "clinical pd prediction",
        "diagnostic pd classifier",
    ]:
        assert forbidden not in text


def test_phase20_reporting_scripts_are_lightweight_and_no_scanpy() -> None:
    combined = "\n".join([REPORT.read_text(encoding="utf-8"), FIGURES.read_text(encoding="utf-8")]).lower()
    assert "matplotlib" in combined
    for forbidden in ["scanpy", "anndata", "read_h5ad", "umap.fit", "fit_transform", "leiden", "neighbors"]:
        assert forbidden not in combined
