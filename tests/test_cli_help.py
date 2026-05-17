from __future__ import annotations

from neurofate import cli


def test_cli_help_mentions_demo_and_doctor() -> None:
    help_text = cli.build_parser().format_help()
    assert "NeuroFate" in help_text
    subcommands = cli.build_parser()._subparsers._group_actions[0].choices
    assert "run-demo" in subcommands
    assert "doctor" in subcommands


def test_workflow_help_text_is_clear() -> None:
    for name, workflow in cli.WORKFLOWS.items():
        assert workflow["help"], name
        assert workflow["command"].startswith("python scripts/"), name
        assert "hint" in workflow, name


def test_cli_run_demo_returns_success() -> None:
    assert cli.main(["run-demo"]) == 0


def test_python_module_entrypoint_exists() -> None:
    assert (cli.PROJECT_ROOT / "neurofate/__main__.py").exists()
