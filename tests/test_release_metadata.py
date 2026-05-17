from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_release_metadata_files_exist() -> None:
    expected = [
        "CITATION.cff",
        "codemeta.json",
        "CHANGELOG.md",
        "RELEASE_CHECKLIST.md",
        "CONTRIBUTING.md",
        "CODE_OF_CONDUCT.md",
        "docs/license_and_data_use.md",
        "docs/runtime_benchmarks.md",
        ".github/workflows/ci.yml",
    ]
    for path in expected:
        assert (ROOT / path).exists(), path


def test_codemeta_is_valid_json() -> None:
    metadata = json.loads((ROOT / "codemeta.json").read_text(encoding="utf-8"))
    assert metadata["name"] == "NeuroFate"
    assert metadata["version"] == "0.1.0"


def test_ci_runs_release_smoke_checks() -> None:
    ci = (ROOT / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    assert "python -m pip install -e ." in ci
    assert "python -m py_compile scripts/*.py neurofate/*.py" in ci
    assert "python -m pytest tests/" in ci
    assert "neurofate check-system" in ci
    assert "neurofate doctor" in ci
    assert "neurofate run-demo" in ci
