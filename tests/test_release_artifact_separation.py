from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_review_zip_builders_write_to_release_artifacts() -> None:
    for script_name in ["66_create_source_release_package.py", "67_create_results_review_package.py"]:
        text = (ROOT / "scripts" / script_name).read_text(encoding="utf-8")
        assert 'default=Path("release_artifacts")' in text
        assert "zip_path = args.output_dir" in text
        assert "manifest_path = args.output_dir" in text
        assert 'default=Path("dist")' not in text


def test_dist_is_documented_as_pypi_only() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    checklist = (ROOT / "RELEASE_CHECKLIST.md").read_text(encoding="utf-8")
    assert "`dist/` is reserved for PyPI artifacts" in readme
    assert "dist/` should contain only `.whl` and `.tar.gz`" in checklist
