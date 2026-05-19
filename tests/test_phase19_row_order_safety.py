from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/97_audit_gse243639_row_order_annotation_link.py"


def load_row_order_module():
    spec = importlib.util.spec_from_file_location("phase19_row_order", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_row_order_is_not_accepted_without_independent_checks() -> None:
    module = load_row_order_module()
    decision, reason = module.decide_row_order(
        count_equal=True,
        sample_column_present=False,
        sample_values_match_clinical=False,
        first_last_consistent=False,
        sheet_consistent=True,
    )
    assert decision == "inconclusive_row_order_linkage"
    assert "sample" in reason.lower()


def test_row_order_requires_all_audit_checks() -> None:
    module = load_row_order_module()
    decision, _ = module.decide_row_order(
        count_equal=True,
        sample_column_present=True,
        sample_values_match_clinical=True,
        first_last_consistent=True,
        sheet_consistent=True,
    )
    assert decision == "safe_row_order_linkage"
    decision, _ = module.decide_row_order(
        count_equal=False,
        sample_column_present=True,
        sample_values_match_clinical=True,
        first_last_consistent=True,
        sheet_consistent=True,
    )
    assert decision == "unsafe_row_order_linkage"


def test_row_order_script_states_hypothesis_not_assumption() -> None:
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "hypothesis" in text
    assert "safe_row_order_linkage" in text
    assert "inconclusive_row_order_linkage" in text
    assert "unsafe_row_order_linkage" in text
    for forbidden in ["scanpy", "anndata", "read_h5ad", "fit_transform", "leiden", "neighbors", "model.fit("]:
        assert forbidden not in text
