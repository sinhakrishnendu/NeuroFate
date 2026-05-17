from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_no_overclaiming_audit_tracks_unsafe_phrases() -> None:
    source = (ROOT / "scripts/54_no_overclaiming_audit.py").read_text(encoding="utf-8")
    expected = [
        "validated across cohorts",
        "foundation model",
        "causal",
        "clinical-grade",
        "diagnostic tool",
        "no_overclaiming_audit.tsv",
    ]
    for token in expected:
        assert token in source


def test_no_overclaiming_audit_is_text_only() -> None:
    source = (ROOT / "scripts/54_no_overclaiming_audit.py").read_text(encoding="utf-8").lower()
    forbidden = ["scanpy", "read_h5ad", "anndata", "h5py", "torch", "pca", "umap", "clustering"]
    for token in forbidden:
        assert token not in source


def test_phase9_text_keeps_external_feasibility_warning() -> None:
    source = (ROOT / "scripts/45_generate_phase9_results_text.py").read_text(encoding="utf-8").lower()
    assert "preliminary external feasibility" in source
    assert "n=6" in source or "small" in source
