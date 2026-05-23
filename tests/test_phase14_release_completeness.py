from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_license_is_mit_and_metadata_agrees() -> None:
    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    codemeta = (ROOT / "codemeta.json").read_text(encoding="utf-8")
    assert license_text.startswith("MIT License")
    assert 'license = { text = "MIT" }' in pyproject
    assert "License :: OSI Approved :: MIT License" in pyproject
    assert "license: MIT" in citation
    assert '"license": "MIT"' in codemeta


def test_source_release_includes_release_completeness_files() -> None:
    text = (ROOT / "scripts/66_create_source_release_package.py").read_text(encoding="utf-8")
    for path in ["LICENSE", "PYPI_RELEASE_CHECKLIST.md", "environment.yml"]:
        assert f'"{path}"' in text


def test_results_review_package_includes_license_and_pypi_checklist() -> None:
    text = (ROOT / "scripts/67_create_results_review_package.py").read_text(encoding="utf-8")
    assert '"LICENSE"' in text
    assert '"PYPI_RELEASE_CHECKLIST.md"' in text


def test_manifest_includes_conda_and_pypi_release_files() -> None:
    text = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "include LICENSE" in text
    assert "include environment.yml" in text
    assert "include PYPI_RELEASE_CHECKLIST.md" in text


def test_review_archives_do_not_use_pypi_dist_directory() -> None:
    source_builder = (ROOT / "scripts/66_create_source_release_package.py").read_text(encoding="utf-8")
    results_builder = (ROOT / "scripts/67_create_results_review_package.py").read_text(encoding="utf-8")
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    checklist = (ROOT / "RELEASE_CHECKLIST.md").read_text(encoding="utf-8")
    assert 'default=Path("release_artifacts")' in source_builder
    assert 'default=Path("release_artifacts")' in results_builder
    assert "release_artifacts/neurofate_source_release_<timestamp>.zip" in readme
    assert "release_artifacts/neurofate_results_review_<timestamp>.zip" in readme
    assert "dist/` should contain only `.whl` and `.tar.gz`" in checklist
