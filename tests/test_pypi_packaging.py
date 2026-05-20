from __future__ import annotations

import tomllib
from pathlib import Path


def test_pyproject_has_console_entry_point_and_lightweight_metadata() -> None:
    pyproject = tomllib.loads(Path("pyproject.toml").read_text())
    project = pyproject["project"]
    assert project["name"] == "neurofate"
    assert project["version"] == "0.3.0"
    assert project["scripts"]["neurofate"] == "neurofate.cli:main"
    assert "scanpy" not in "\n".join(project["dependencies"]).lower()
    assert "anndata" not in "\n".join(project["dependencies"]).lower()
    assert "mps" in project["optional-dependencies"]
    assert "plotting" in project["optional-dependencies"]


def test_release_metadata_points_to_public_repository() -> None:
    citation = Path("CITATION.cff").read_text()
    codemeta = Path("codemeta.json").read_text()
    assert "https://github.com/sinhakrishnendu/NeuroFate" in citation
    assert "https://github.com/sinhakrishnendu/NeuroFate" in codemeta
    assert "research software" in citation.lower()


def test_package_data_manifest_includes_demo_and_metadata() -> None:
    manifest = Path("MANIFEST.in").read_text()
    pyproject = Path("pyproject.toml").read_text()
    assert "recursive-include examples/tiny_demo" in manifest
    assert "recursive-include examples/format_examples" in manifest
    assert "recursive-include metadata *.tsv" in manifest
    assert "resources/tiny_demo/*.tsv" in pyproject
    assert "resources/*.tsv" in pyproject
