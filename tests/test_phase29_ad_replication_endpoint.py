import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/152_test_gse174367_bulk_ad_axis_replication.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase29_ad_replication", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_endpoint_is_locked_to_ad_vs_control():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "label__ad_vs_control" in text
    assert "positive_class" in text
    assert "negative_class" in text


def test_evidence_label_requires_statistical_support():
    module = load_module()
    label, consistency = module.evidence_label(
        effect=0.3,
        pvalue=0.4,
        fdr=0.8,
        phase22="positive",
        n=230,
        positive_n=118,
        negative_n=112,
    )
    assert consistency == "consistent"
    assert label == "directionally_consistent_but_not_significant"
    label, consistency = module.evidence_label(
        effect=0.3,
        pvalue=0.01,
        fdr=0.05,
        phase22="positive",
        n=230,
        positive_n=118,
        negative_n=112,
    )
    assert consistency == "consistent"
    assert label == "statistically_supported_ad_replication"
