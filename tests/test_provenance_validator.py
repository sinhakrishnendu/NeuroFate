import importlib.util
from pathlib import Path


def load_validator_module():
    path = Path("scripts/06_validate_provenance.py")
    spec = importlib.util.spec_from_file_location("validate_provenance", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_provenance_validator_constants_include_required_columns():
    module = load_validator_module()
    assert "dataset_id" in module.PROVENANCE_COLUMNS
    assert "checksum_value" in module.PROVENANCE_COLUMNS
    assert "download_command_manual_only" in module.PROVENANCE_COLUMNS


def test_provenance_validator_complete_status_policy():
    module = load_validator_module()
    assert "complete" in module.COMPLETE_STATUSES
    assert "pending" in module.MISSING_CHECKSUM_VALUES
