from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCHMARKING = ROOT / "neurofate/benchmarking.py"


def test_label_and_identifier_columns_are_excluded_in_source() -> None:
    text = BENCHMARKING.read_text(encoding="utf-8")
    assert 'LABEL_PREFIX = "label__"' in text
    assert 'EXCLUDED_COLUMNS = {"donor_id", "sample_id", "cell_id", "cohort_id", "n_cells"}' in text
    assert "column.startswith(LABEL_PREFIX)" in text
    assert "column in EXCLUDED_COLUMNS" in text
    assert "column.startswith(FEATURE_PREFIXES)" in text


def test_task_labels_are_derived_from_label_columns_in_source() -> None:
    text = BENCHMARKING.read_text(encoding="utf-8")
    for token in [
        "label__Cognitive_Status",
        "label__Overall_AD_neuropathological_Change",
        "label__APOE_Genotype",
        "label__Highest_Lewy_Body_Disease",
        "label__LATE",
        "label__Overall_CAA_Score",
        "dementia_vs_reference",
        "high_vs_low_ad_neuropathology",
        "apoe_risk_prediction",
        "mixed_pathology_burden",
    ]:
        assert token in text
