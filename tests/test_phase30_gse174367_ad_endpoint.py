import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py"
CLEAN_TESTER = ROOT / "scripts/156_test_gse174367_bulk_ad_axis_replication_clean.py"
PHASE29_TESTER = ROOT / "scripts/152_test_gse174367_bulk_ad_axis_replication.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_ad_control_label_mapping():
    module = load(CONVERTER, "phase30_endpoint_converter")
    assert module.label_for_endpoint("AD") == "1"
    assert module.label_for_endpoint("Control") == "0"
    assert module.label_for_endpoint("Parkinson's Disease") == ""


def test_clean_tester_defaults_to_phase30_outputs():
    text = CLEAN_TESTER.read_text(encoding="utf-8")
    assert "phase30_gse174367_bulk_axis_scores.tsv" in text
    assert "phase30_gse174367_bulk_axis_replication_statistics.tsv" in text
    assert "phase30_gse174367_bulk_axis_replication_fdr.tsv" in text


def test_expected_sample_support_is_sufficient_for_endpoint():
    module = load(PHASE29_TESTER, "phase30_endpoint_tester")
    label, consistency = module.evidence_label(
        effect=0.2,
        pvalue=0.2,
        fdr=0.3,
        phase22="positive",
        n=90,
        positive_n=44,
        negative_n=46,
    )
    assert consistency == "consistent"
    assert label == "directionally_consistent_but_not_significant"
