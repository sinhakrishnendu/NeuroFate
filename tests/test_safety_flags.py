from pathlib import Path

import yaml


def test_phase_1b_safety_flags_are_strict():
    config = yaml.safe_load(Path("configs/project_config.yaml").read_text(encoding="utf-8"))
    safety = config["safety"]
    assert safety["allow_downloads"] is False
    assert safety["allow_heavy_processing"] is False
    assert safety["allow_model_training"] is False
    assert safety["allow_h5ad_processing"] is False
    assert safety["manual_execution_only"] is True


def test_legacy_safety_flags_remain_strict():
    config = yaml.safe_load(Path("configs/project_config.yaml").read_text(encoding="utf-8"))
    safety = config["safety"]
    assert safety["allow_large_hdf5_processing"] is False
    assert safety["allow_training"] is False
    assert safety["require_manual_heavy_commands"] is True
