from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/71_inspect_external_metadata_safe.py"


def load_metadata_module():
    spec = importlib.util.spec_from_file_location("metadata_safe", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_metadata_parser_supports_semicolon_and_header_line_6(tmp_path: Path) -> None:
    module = load_metadata_module()
    clinical = tmp_path / "GSE243639_Clinical_data.csv.gz"
    lines = [
        "prose line 1",
        "prose line 2",
        "prose line 3",
        "prose line 4",
        "",
        "Sample ID;Clinical diagnosis;Age;Sex;PMI hours;RIN measure",
        "s.0096;Parkinson's;77;M;12;7.8",
        "s.0100;Control;69;F;10;8.1",
    ]
    with gzip.open(clinical, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\n".join(lines) + "\n")
    header, preview, delimiter, header_index = module.inspect_text_header(clinical, ";", header_line=6)
    assert delimiter == ";"
    assert header_index == 5
    assert header[:2] == ["Sample ID", "Clinical diagnosis"]
    assert preview[0][0] == "s.0096"


def test_metadata_parser_canonical_gse243639_fields() -> None:
    module = load_metadata_module()
    expected = {
        "Sample ID": "sample_id",
        "Clinical diagnosis": "diagnosis",
        "Age": "age",
        "Sex": "sex",
        "PMI hours": "pmi",
        "RIN measure": "rin",
        "Lewy bodies presence in midbrain": "lewy_body_midbrain",
        "Lewy bodies presence in limbic regions (amygdala)": "lewy_body_limbic",
        "Lewy bodies presence in neocortical regions (frontal cortex)": "lewy_body_neocortical",
        "CERAD score for neuritic plaques": "cerad",
        "Braak stage for neurofibrillary tangles": "braak",
    }
    for source, canonical in expected.items():
        assert module.suggest_mapping(source)[0] == canonical


def test_metadata_inspector_has_cli_header_and_delimiter_arguments() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "--header-line" in text
    assert "--delimiter" in text
    assert "metadata_preview" in text
