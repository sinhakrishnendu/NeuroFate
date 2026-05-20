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
