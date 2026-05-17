from pathlib import Path

import yaml


def test_intake_config_safety_flags_are_strict():
    config = yaml.safe_load(Path("configs/intake_config.yaml").read_text(encoding="utf-8"))
    intake = config["intake"]
    assert intake["manual_execution_only"] is True
    assert intake["allow_downloads"] is False
    assert intake["allow_remote_fetch"] is False
    assert intake["allow_checksum_only"] is True
    assert intake["allow_file_opening"] is False
    assert intake["allow_h5ad_processing"] is False
    assert intake["allow_large_file_validation"] is False
    assert intake["max_lightweight_file_size_mb"] == 5
