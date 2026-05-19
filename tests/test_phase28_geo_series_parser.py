import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/143_parse_geo_series_matrix_generic.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase28_series", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_generic_geo_series_parser_extracts_ad_metadata(tmp_path):
    module = load_module()
    series = tmp_path / "series.txt.gz"
    rows = [
        ["!Sample_title", "AD_1", "CTRL_1"],
        ["!Sample_geo_accession", "GSM1", "GSM2"],
        ["!Sample_source_name_ch1", "brain", "brain"],
        ["!Sample_characteristics_ch1", "diagnosis: Alzheimer's disease", "diagnosis: Control"],
        ["!Sample_characteristics_ch1", "brain region: frontal cortex", "brain region: frontal cortex"],
        ["!Sample_characteristics_ch1", "age: 80", "age: 77"],
        ["!Sample_characteristics_ch1", "sex: female", "sex: male"],
        ["!Sample_characteristics_ch1", "postmortem interval: 5", "postmortem interval: 7"],
        ["!Sample_characteristics_ch1", "Braak stage: VI", "Braak stage: I"],
        ["!Sample_supplementary_file_1", "ftp://example/AD_1.csv.gz", "ftp://example/CTRL_1.csv.gz"],
    ]
    with gzip.open(series, "wt", encoding="utf-8") as handle:
        for row in rows:
            handle.write("\t".join(f'"{value}"' for value in row) + "\n")
    metadata, manifest, labels = module.build_metadata(series, "mock_ad")
    assert len(metadata) == 2
    assert metadata[0]["diagnosis"] == "Alzheimer's disease"
    assert metadata[0]["inferred_ad_endpoint"] == "ad"
    assert metadata[1]["inferred_ad_endpoint"] == "control"
    assert len(manifest) == 2
    assert any(row["label"] == "ad" and row["count"] == "1" for row in labels)


def test_geo_series_parser_is_metadata_only():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "tarfile.open", "extractall", "umap", "leiden"]:
        assert forbidden not in text
