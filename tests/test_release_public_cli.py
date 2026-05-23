from __future__ import annotations

import subprocess
import sys


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "neurofate", *args],
        check=False,
        text=True,
        capture_output=True,
    )


def test_release_public_cli_help_commands() -> None:
    commands = [
        ("--help",),
        ("ingest", "--help"),
        ("run", "--help"),
        ("adapt-endpoint", "--help"),
        ("build-axis-scores", "--help"),
        ("score-risk", "--help"),
    ]
    for command in commands:
        completed = run_cli(*command)
        assert completed.returncode == 0, completed.stderr
        assert "clinical diagnostic" not in (completed.stdout + completed.stderr).lower()


def test_release_public_cli_smoke_commands() -> None:
    for command in [("check-system",), ("doctor",), ("run-demo",)]:
        completed = run_cli(*command)
        assert completed.returncode == 0, completed.stderr

