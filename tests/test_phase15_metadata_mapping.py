from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/71_inspect_external_metadata_safe.py"


def test_metadata_inspector_does_not_import_scanpy_or_anndata() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "scanpy" not in text
    assert "read_h5ad" not in text
    assert "anndata" not in text


def test_h5_metadata_inspector_marks_x_forbidden_without_accessing_it() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert 'key == "X"' in text
    assert "FORBIDDEN_X_PRESENT_NOT_ACCESSED" in text
    assert 'handle["X"]' not in text
    assert "continue" in text


def test_metadata_mapping_covers_required_canonical_fields() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for field in [
        "donor_id",
        "sample_id",
        "cell_id",
        "diagnosis",
        "disease_status",
        "pathology",
        "brain_region",
        "age",
        "sex",
        "cell_type",
        "batch",
        "apoe_genotype",
        "sequencing_platform",
    ]:
        assert field in text
