from __future__ import annotations

import json
import tomllib
from pathlib import Path

import neurofate


def test_release_version_and_authors_are_consistent() -> None:
    pyproject = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))
    citation = Path("CITATION.cff").read_text(encoding="utf-8")
    codemeta = json.loads(Path("codemeta.json").read_text(encoding="utf-8"))
    readme = Path("README.md").read_text(encoding="utf-8")

    assert pyproject["project"]["version"] == "0.3.0"
    assert neurofate.__version__ == "0.3.0"
    assert "version: 0.3.0" in citation
    assert codemeta["version"] == "0.3.0"
    assert "Current release-candidate version: **0.3.0**" in readme

    authors = pyproject["project"]["authors"]
    assert authors[0]["name"] == "Nabanita Ghosh"
    assert authors[1]["name"] == "Krishnendu Sinha"


def test_release_pyproject_has_pypi_ready_metadata() -> None:
    pyproject = tomllib.loads(Path("pyproject.toml").read_text(encoding="utf-8"))
    project = pyproject["project"]
    classifiers = set(project["classifiers"])
    dependencies = "\n".join(project["dependencies"]).lower()

    assert project["name"] == "neurofate"
    assert project["scripts"]["neurofate"] == "neurofate.cli:main"
    assert "License :: OSI Approved :: MIT License" in classifiers
    assert "Topic :: Scientific/Engineering :: Bio-Informatics" in classifiers
    assert "scanpy" not in dependencies
    assert "anndata" not in dependencies
    for extra in ["mps", "plotting", "dev", "docs"]:
        assert extra in project["optional-dependencies"]


def test_release_manifest_keeps_package_small_and_reviewable() -> None:
    manifest = Path("MANIFEST.in").read_text(encoding="utf-8")
    pyproject = Path("pyproject.toml").read_text(encoding="utf-8")
    assert "recursive-include examples/tiny_demo" in manifest
    assert "recursive-include examples/format_examples" in manifest
    assert "recursive-include metadata *.tsv" in manifest
    assert "resources/tiny_demo/*.tsv" in pyproject
    assert "resources/*.tsv" in pyproject

