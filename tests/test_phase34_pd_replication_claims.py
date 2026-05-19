import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TESTER = ROOT / "scripts/166_test_phase34_pd_axis_replication.py"
READINESS = ROOT / "scripts/123_build_pnas_readiness_matrix.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_phase34_direction_only_does_not_upgrade_pd_claim():
    module = load(TESTER, "phase34_tester")
    assert module.evidence_label(effect=-0.2, ref_direction=-1, pvalue=0.4, fdr=0.8, n=18, positive_n=10, negative_n=8) == "directionally_consistent_but_not_significant"
    assert module.evidence_label(effect=-0.2, ref_direction=-1, pvalue=0.03, fdr=0.2, n=18, positive_n=10, negative_n=8) == "statistically_supported_pd_replication"


def test_readiness_matrix_knows_phase34_outputs():
    text = READINESS.read_text(encoding="utf-8")
    assert "phase34_*_pd_axis_replication_statistics.tsv" in text
    assert "statistically_supported_pd_replication" in text
