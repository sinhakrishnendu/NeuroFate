from __future__ import annotations

from pathlib import Path


def _manuscript_text() -> str:
    paths = [
        "manuscript/bioinformatics/neurofate_bioinformatics_application_note.tex",
        "manuscript/bioinformatics/neurofate_bioinformatics_full_methods_paper.tex",
    ]
    return "\n".join(Path(path).read_text().lower() for path in paths)


def test_bioinformatics_manuscript_is_software_framed() -> None:
    text = _manuscript_text()
    assert "command-line" in text
    assert "pypi" in text
    assert "research software" in text
    assert "biological discovery manuscript" not in text


def test_bioinformatics_manuscript_avoids_forbidden_claims() -> None:
    text = _manuscript_text()
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
    assert "not clinical diagnostic software" not in text
