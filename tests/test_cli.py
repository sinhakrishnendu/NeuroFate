from __future__ import annotations

from neurofate import cli


def test_cli_exposes_required_subcommands() -> None:
    parser = cli.build_parser()
    subcommands = parser._subparsers._group_actions[0].choices
    expected = {
        "check-system",
        "inspect-sea-ad",
        "extract-metadata",
        "extract-panel",
        "summarize",
        "train-baseline",
        "train-mps",
        "validate-external",
        "make-report",
        "doctor",
    }
    assert expected.issubset(set(subcommands))


def test_cli_safe_commands_return_success() -> None:
    assert cli.main(["check-system"]) == 0
    assert cli.main(["make-report"]) == 0
    assert cli.main(["doctor"]) == 0


def test_cli_entrypoint_declared() -> None:
    pyproject = (cli.PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'neurofate = "neurofate.cli:main"' in pyproject
