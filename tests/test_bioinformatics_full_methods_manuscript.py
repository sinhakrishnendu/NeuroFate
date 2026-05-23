from __future__ import annotations

from pathlib import Path


MANUSCRIPT = Path("manuscript/bioinformatics/neurofate_bioinformatics_full_methods_paper.tex")


def test_full_methods_manuscript_mentions_release_candidate_workflow() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8").lower()
    required = [
        "format-aware",
        "geo series matrix",
        "neurofate run",
        "research-use risk",
        "gse20141",
        "0.3.0",
    ]
    for phrase in required:
        assert phrase in text


def test_full_methods_manuscript_avoids_clinical_claims() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8").lower()
    forbidden = [
        "clinically validated diagnosis",
        "medical device",
        "patient diagnosis",
        "diagnostic biomarker",
        "fda/ce-ready",
        "validated shared ad/pd mechanism",
        "clinical decision-making tool",
    ]
    for phrase in forbidden:
        assert phrase not in text


def test_full_methods_manuscript_uses_single_mechanistic_figure() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    assert text.count("\\includegraphics") == 1
    assert "figures/figure1_full_workflow.pdf" in text
    assert "figure2_ingestion_engine" not in text
    assert "figure3_cli_outputs" not in text
    assert "figure4_multicohort_demonstration" not in text
