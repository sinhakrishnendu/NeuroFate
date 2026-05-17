from pathlib import Path

import yaml


def test_project_config_safety_defaults():
    config = yaml.safe_load(Path("configs/project_config.yaml").read_text(encoding="utf-8"))
    assert config["safety"]["allow_downloads"] is False
    assert config["safety"]["allow_training"] is False
    assert config["safety"]["require_manual_heavy_commands"] is True


def test_dataset_config_contains_placeholders():
    config = yaml.safe_load(Path("configs/datasets.yaml").read_text(encoding="utf-8"))
    datasets = config["datasets"]
    assert datasets
    assert all(item["status"] == "placeholder" for item in datasets)
