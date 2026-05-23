from __future__ import annotations

from neurofate.cli import build_parser


def test_public_commands_are_available() -> None:
    parser = build_parser()
    actions = [action for action in parser._actions if getattr(action, "choices", None)]
    choices = set(actions[0].choices)
    expected = {
        "check-system",
        "doctor",
        "run-demo",
        "ingest",
        "run",
        "build-axis-scores",
        "score-risk",
        "adapt-endpoint",
        "audit-leakage",
        "train-baseline",
        "train-mps",
        "validate-external",
        "make-report",
        "benchmark",
        "benchmark-report",
    }
    assert expected <= choices


def test_help_text_is_research_software_oriented() -> None:
    parser = build_parser()
    help_text = parser.format_help().lower()
    assert "command-line research software" in help_text
    assert "clinical diagnostic" not in help_text
