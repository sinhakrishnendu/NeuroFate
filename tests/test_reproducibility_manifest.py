from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_reproducibility_manifest_captures_expected_metadata() -> None:
    source = (ROOT / "scripts/52_generate_reproducibility_manifest.py").read_text(encoding="utf-8")
    expected = [
        "python_version",
        "platform",
        "package_versions",
        "git_commit",
        "input_file_inventory",
        "output_inventory",
        "sha256",
        "reproducibility_manifest.json",
    ]
    for token in expected:
        assert token in source


def test_reproducibility_manifest_avoids_analysis_engines() -> None:
    source = (ROOT / "scripts/52_generate_reproducibility_manifest.py").read_text(encoding="utf-8").lower()
    forbidden = ["scanpy", "read_h5ad", "umap", "leiden", "fit("]
    for token in forbidden:
        assert token not in source
