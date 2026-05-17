from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_documentation_pages_exist() -> None:
    expected = [
        "docs/installation.md",
        "docs/quickstart.md",
        "docs/input_formats.md",
        "docs/tutorial_sea_ad.md",
        "docs/tutorial_mathys_external_validation.md",
        "docs/apple_silicon_mps.md",
        "docs/interpretation_guide.md",
        "docs/troubleshooting.md",
        "docs/reproducibility.md",
        "docs/limitations.md",
    ]
    for path in expected:
        full_path = ROOT / path
        assert full_path.exists(), path
        assert full_path.read_text(encoding="utf-8").strip(), path


def test_configuration_templates_exist() -> None:
    expected = [
        "configs/templates/sea_ad_minimal.yaml",
        "configs/templates/mathys_csv.yaml",
        "configs/templates/target_gene_panel.yaml",
        "configs/templates/mps_model.yaml",
        "configs/templates/external_validation.yaml",
    ]
    for path in expected:
        assert (ROOT / path).exists(), path
