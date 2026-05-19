import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/139_audit_gse184950_sample_integrity.py"


def load_module():
    spec = importlib.util.spec_from_file_location("phase27_integrity", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_sample_integrity_flags_processed_matrices_pseudo_sample():
    module = load_module()
    metadata = [{"sample_name": "A10"}, {"sample_name": "A15"}]
    scores = [{"sample_id": "A10"}, {"sample_id": "A15"}, {"sample_id": "processed_matrices"}]
    audit = [{"sample_id": "A10"}, {"sample_id": "processed_matrices"}]
    row = module.audit(scores, metadata, audit)
    assert row["expected_samples"] == "2"
    assert row["observed_axis_score_samples"] == "3"
    assert row["valid_axis_score_samples"] == "2"
    assert row["invalid_axis_score_samples"] == "1"
    assert row["invalid_sample_ids"] == "processed_matrices"
    assert row["invalid_axis_gene_audit_ids"] == "processed_matrices"
