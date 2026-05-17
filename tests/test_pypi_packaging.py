from __future__ import annotations

import tomllib
from importlib import resources
from pathlib import Path

from neurofate import cli


ROOT = Path(__file__).resolve().parents[1]


def test_pyproject_has_pypi_metadata() -> None:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = pyproject["project"]
    assert project["name"] == "neurofate"
    assert project["version"] == "0.1.0"
    assert project["requires-python"] == ">=3.11,<3.13"
    assert "neurofate" in pyproject["project"]["scripts"]
    assert pyproject["project"]["scripts"]["neurofate"] == "neurofate.cli:main"
    assert "classifiers" in project
    assert "urls" in project


def test_package_data_includes_tiny_demo_resources() -> None:
    demo_files = resources.files("neurofate.resources.tiny_demo")
    assert demo_files.joinpath("tiny_metadata.tsv").is_file()
    assert demo_files.joinpath("tiny_gene_panel.tsv").is_file()
    assert demo_files.joinpath("tiny_sparse_expression.tsv").is_file()


def test_manifest_and_license_files_exist() -> None:
    expected = ["LICENSE", "MANIFEST.in", "PYPI_RELEASE_CHECKLIST.md", "docs/pypi_release.md"]
    for path in expected:
        assert (ROOT / path).exists(), path


def test_doctor_succeeds_in_repository_layout() -> None:
    assert cli.main(["doctor"]) == 0
