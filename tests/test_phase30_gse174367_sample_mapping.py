import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/155_audit_gse174367_bulk_sample_mapping.py"
CONVERTER = ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_audit_script_prefers_rda_targets_candidates():
    text = AUDIT.read_text(encoding="utf-8")
    assert "normExpr.reg" in text
    assert "targets" in text
    assert "SampleID" in text
    assert "geo_accession" in text


def test_converter_prefers_targets_over_series_metadata():
    module = load(CONVERTER, "phase30_converter")
    r_code = module.r_axis_extraction_code(["SNCA", "MAPT"])
    assert "normExpr.reg" in r_code
    assert "targets" in r_code
    assert "SampleID" in r_code
    assert "series_metadata_role <- 'secondary_annotation_only'" in r_code


def test_converter_refuses_zero_overlap_mapping():
    module = load(CONVERTER, "phase30_converter_zero")
    r_code = module.r_axis_extraction_code(["SNCA"])
    assert "best_overlap == 0" in r_code
    assert "stop rather than mapping to GEO series metadata" in r_code
