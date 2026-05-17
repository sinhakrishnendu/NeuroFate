from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_end_user_report_generator_is_lightweight() -> None:
    source = (ROOT / "scripts/51_generate_end_user_report.py").read_text(encoding="utf-8").lower()
    forbidden = ["scanpy", "read_h5ad", "anndata", "h5py", "torch", "fit("]
    for token in forbidden:
        assert token not in source
    assert "does not run analysis" in source
    assert "neurofate_analysis_report.html" in source
    assert "neurofate_analysis_report.md" in source


def test_readme_mentions_platform_safety() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "What NeuroFate Does" in readme
    assert "Safety And Memory Design" in readme
    assert "Current Validation Status" in readme
    assert "python scripts/51_generate_end_user_report.py" in readme
