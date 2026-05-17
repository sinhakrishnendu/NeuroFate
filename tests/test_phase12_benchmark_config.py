from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_benchmark_config_has_required_values() -> None:
    text = (ROOT / "configs/benchmark_config.yaml").read_text(encoding="utf-8")
    assert "manual_execution_only: true" in text
    assert "seed_list: [1, 7, 11, 19, 42, 101, 202]" in text
    assert "n_permutations: 100" in text
    assert "n_bootstrap: 200" in text
    assert "test_size: 0.25" in text
    assert "min_class_count: 5" in text
    assert "max_runtime_minutes: 30" in text


def test_benchmark_config_models_and_tasks() -> None:
    text = (ROOT / "configs/benchmark_config.yaml").read_text(encoding="utf-8")
    for model in {
        "logistic_regression",
        "elastic_net",
        "random_forest",
        "gradient_boosting",
    }:
        assert f"  - {model}" in text
    for task in {
        "dementia_vs_reference",
        "high_vs_low_ad_neuropathology",
        "apoe_risk_prediction",
        "mixed_pathology_burden",
    }:
        assert f"  - {task}" in text
