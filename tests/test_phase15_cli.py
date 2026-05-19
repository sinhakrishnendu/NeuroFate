from __future__ import annotations

from neurofate import cli


def test_phase15_cli_commands_exist() -> None:
    parser = cli.build_parser()
    subcommands = parser._subparsers._group_actions[0].choices
    for command in [
        "external-triage",
        "inspect-external",
        "plan-external-extraction",
        "validate-multi-external",
        "external-report",
    ]:
        assert command in subcommands


def test_phase15_cli_workflows_are_safe_by_default() -> None:
    for command in [
        "external-triage",
        "inspect-external",
        "plan-external-extraction",
        "validate-multi-external",
        "external-report",
    ]:
        assert cli.main([command]) == 0


def test_phase15_cli_commands_point_to_phase15_scripts() -> None:
    workflows = cli.WORKFLOWS
    assert "scripts/69_triage_external_validation_candidates.py" in workflows["external-triage"]["command"]
    assert "scripts/70_inspect_external_dataset_files.py" in workflows["inspect-external"]["command"]
    assert "scripts/73_prepare_external_sparse_extraction_plan.py" in workflows["plan-external-extraction"]["command"]
    assert "scripts/75_run_multi_external_validation.py" in workflows["validate-multi-external"]["command"]
    assert "scripts/77_generate_phase15_external_validation_report.py" in workflows["external-report"]["command"]
