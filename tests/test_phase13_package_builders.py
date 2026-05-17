from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_source_package_builder_excludes_large_paths() -> None:
    text = (ROOT / "scripts/66_create_source_release_package.py").read_text(encoding="utf-8")
    assert '"data/"' in text
    assert '"results/models/"' in text
    assert '"results/figures/"' in text
    assert '"results/tables/"' in text
    assert '".h5ad"' in text
    assert '".pt"' in text
    assert "neurofate_source_release_" in text
    assert "source_release_manifest.tsv" in text
    assert 'default=Path("release_artifacts")' in text
    assert 'release_artifacts/"' in text


def test_results_package_builder_excludes_raw_data_and_models() -> None:
    text = (ROOT / "scripts/67_create_results_review_package.py").read_text(encoding="utf-8")
    assert '"data/"' in text
    assert '"results/models/"' in text
    assert '".h5ad"' in text
    assert '".pt"' in text
    assert "neurofate_results_review_" in text
    assert "results_review_manifest.tsv" in text
    assert 'default=Path("release_artifacts")' in text
    assert 'release_artifacts/"' in text


def test_results_interpretation_exists() -> None:
    text = (ROOT / "RESULTS_INTERPRETATION.md").read_text(encoding="utf-8")
    assert "Executive Summary" in text
    assert "Mathys" in text
    assert "not clinical-grade" in text
    assert "not a foundation model" in text
