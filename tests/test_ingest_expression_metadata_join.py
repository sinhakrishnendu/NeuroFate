from __future__ import annotations

from neurofate.ingest import validate_expression_metadata_join


def test_expression_metadata_join_uses_safe_normalization() -> None:
    audit = validate_expression_metadata_join(["GSM-001", "GSM-002"], ["GSM001", "GSM002", "GSM003"])
    assert audit["expression_sample_count"] == 2
    assert audit["metadata_sample_count"] == 3
    assert audit["matched_sample_count"] == 2
    assert "GSM003" in audit["unmatched_metadata_samples"]

