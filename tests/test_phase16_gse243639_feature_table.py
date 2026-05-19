from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/79_build_gse243639_feature_table.py"


def load_feature_module():
    spec = importlib.util.spec_from_file_location("gse243639_features", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_feature_builder_parses_clinical_header_line_6_and_semicolon(tmp_path: Path) -> None:
    module = load_feature_module()
    clinical = tmp_path / "clinical.csv.gz"
    lines = [
        "preamble 1",
        "preamble 2",
        "preamble 3",
        "preamble 4",
        "",
        "Sample ID;Clinical diagnosis;Age;Sex;PMI hours;RIN measure",
        "s.0096;Parkinson's;77;M;12;7.8",
        "s.0100;Control;69;F;10;8.1",
    ]
    with gzip.open(clinical, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(lines) + "\n")
    rows = module.read_clinical(clinical, header_line=6, delimiter=";")
    assert rows["s.0096"]["diagnosis"] == "Parkinson's"
    assert rows["s.0096"]["label__diagnosis_binary"] == "1"
    assert rows["s.0100"]["label__diagnosis_binary"] == "0"


def test_feature_builder_is_sample_level() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert '"unit_type": "sample"' in text
    assert "sample_id" in text
    assert "dataset_unit_id" in text
    assert "label__diagnosis_binary" in text


def test_feature_builder_creates_gene_mean_and_detection_features() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "gene_mean__" in text
    assert "gene_detection__" in text
    assert "phase5_schema" in text
