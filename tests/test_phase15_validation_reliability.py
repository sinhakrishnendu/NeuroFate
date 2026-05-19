from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/75_run_multi_external_validation.py"


def load_phase15_validation_module():
    spec = importlib.util.spec_from_file_location("phase15_validation", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reliability_categories_are_declared() -> None:
    module = load_phase15_validation_module()
    for category in [
        "reliable_external_validation",
        "preliminary_external_feasibility",
        "insufficient_sample_size",
        "insufficient_feature_overlap",
        "failed_label_mapping",
    ]:
        assert category in module.RELIABILITY_CATEGORIES


def test_small_n_is_not_reliable_external_validation() -> None:
    module = load_phase15_validation_module()
    assert module.reliability(n_test=6, positives=3, negatives=3, overlap=30) == "preliminary_external_feasibility"
    assert module.reliability(n_test=5, positives=2, negatives=3, overlap=30) == "insufficient_sample_size"


def test_missing_class_or_low_overlap_prevents_overclaiming() -> None:
    module = load_phase15_validation_module()
    assert module.reliability(n_test=30, positives=30, negatives=0, overlap=30) == "failed_label_mapping"
    assert module.reliability(n_test=30, positives=15, negatives=15, overlap=5) == "insufficient_feature_overlap"
    assert module.reliability(n_test=30, positives=15, negatives=15, overlap=20) == "reliable_external_validation"


def test_validation_runner_does_not_train_deep_models() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for forbidden in ["torch", "tensorflow", "keras", "scanpy", "read_h5ad", "PCA", "UMAP", "Leiden"]:
        assert forbidden not in text
