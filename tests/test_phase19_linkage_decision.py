from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/98_decide_gse243639_annotation_linkage.py"
REPORT = ROOT / "scripts/94_generate_phase18_gse243639_repair_report.py"


def load_decision_module():
    spec = importlib.util.spec_from_file_location("phase19_decision", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_decision_marks_zero_overlap_as_unsafe() -> None:
    module = load_decision_module()
    decision = module.decide(
        [{"rule_name": "raw_id", "overlap_rate": "0", "overlap_count": "0"}],
        [{"decision": "unsafe_row_order_linkage"}],
        [{"match_rate": "0"}],
    )
    assert decision["decision_category"] == "annotation_linkage_unsafe"
    assert decision["safe_to_build_annotation_map"] == "false"
    assert "Phase 16" in decision["recommended_action"] or "Retire" in decision["recommended_action"]


def test_decision_requires_high_overlap_for_id_linkage() -> None:
    module = load_decision_module()
    decision = module.decide(
        [{"rule_name": "collapse_punctuation", "overlap_rate": "0.949", "overlap_count": "79000"}],
        [{"decision": "inconclusive_row_order_linkage"}],
        [{"match_rate": "0"}],
    )
    assert decision["decision_category"] == "annotation_linkage_inconclusive"
    decision = module.decide(
        [{"rule_name": "collapse_punctuation", "overlap_rate": "0.951", "overlap_count": "79300"}],
        [{"decision": "unsafe_row_order_linkage"}],
        [{"match_rate": "0"}],
    )
    assert decision["decision_category"] == "normalized_id_linkage_safe"


def test_phase18_report_reads_phase19_status() -> None:
    text = REPORT.read_text(encoding="utf-8")
    assert "phase19_gse243639_annotation_linkage_decision.tsv" in text
    assert "Cell-type-aware GSE243639 validation is not currently supported" in text
