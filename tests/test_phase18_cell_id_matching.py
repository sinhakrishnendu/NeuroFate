from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/90_audit_gse243639_cell_id_matching.py"


def load_audit_module():
    spec = importlib.util.spec_from_file_location("phase18_cell_id_audit", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_cell_id_normalization_handles_prefixed_and_unprefixed_ids() -> None:
    module = load_audit_module()
    prefixed = module.cell_id_variants("s.0096_AAACCCAAGTACGAGC.1")
    unprefixed = module.cell_id_variants("AAACCCAAGTACGAGC.1")
    core = module.cell_id_variants("AAACCCAAGTACGAGC")
    assert prefixed["sample_id"] == "s.0096"
    assert prefixed["barcode_core"] == "AAACCCAAGTACGAGC"
    assert unprefixed["barcode_core"] == "AAACCCAAGTACGAGC"
    assert core["barcode_core"] == "AAACCCAAGTACGAGC"


def test_cell_id_audit_reports_required_metrics() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for metric in [
        "unique_expression_cell_ids",
        "unique_cell_sample_map_cell_ids",
        "unique_workbook_cell_ids",
        "direct_overlap_count",
        "overlap_after_removing_sample_prefix",
        "overlap_after_adding_sample_prefix",
        "overlap_after_removing_trailing_dot_suffix",
        "overlap_after_normalizing_punctuation",
        "recommended_normalization_rule",
    ]:
        assert metric in text


def test_cell_id_audit_avoids_forbidden_workloads() -> None:
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "h5py", "fit_transform", "leiden", "neighbors"]:
        assert forbidden not in text
