import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PARSER = ROOT / "scripts/131_parse_gse184950_series_matrix.py"
SCORER = ROOT / "scripts/128_build_gse184950_axis_scores.py"
TESTER = ROOT / "scripts/129_test_gse184950_axis_replication.py"


def load_parser():
    spec = importlib.util.spec_from_file_location("phase25_series_endpoint", PARSER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_pd_pdd_endpoint_mapping():
    module = load_parser()
    assert module.disease_label("Unaffected Control") == "0"
    assert module.disease_label("Parkinson's Disease") == "1"
    assert module.disease_label("Parkinson's Disease Dementia") == "1"


def test_axis_scorer_and_tester_use_phase25_endpoint():
    scorer = SCORER.read_text(encoding="utf-8")
    tester = TESTER.read_text(encoding="utf-8")
    assert "phase25_gse184950_series_sample_metadata.tsv" in scorer
    assert "label__pd_pdd_vs_control" in scorer
    assert "phase25_gse184950_axis_scores.tsv" in scorer
    assert "label__pd_pdd_vs_control" in tester
    assert "phase25_gse184950_axis_replication_statistics.tsv" in tester
    assert "phase25_gse184950_axis_replication_fdr.tsv" in tester
