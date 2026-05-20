from __future__ import annotations

from pathlib import Path


def test_research_use_only_docs_exist_and_forbid_deployment_language() -> None:
    text = Path("docs/research_use_only.md").read_text().lower()
    assert "research use only" in text
    assert "not validated for" in text
    assert "treatment selection" in text


def test_cli_reports_do_not_print_deployment_claims() -> None:
    combined = "\n".join(
        [
            Path("neurofate/axis.py").read_text(),
            Path("neurofate/demo.py").read_text(),
            Path("docs/cli_reference.md").read_text(),
        ]
    ).lower()
    forbidden_positive = [
        "diagnosis confirmed",
        "patient diagnosis",
        "clinical recommendation",
        "fda/ce-ready",
    ]
    for phrase in forbidden_positive:
        assert phrase not in combined

