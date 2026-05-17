"""Quality-control placeholders."""

from __future__ import annotations


def describe_qc_plan() -> dict[str, str]:
    """Return the intended QC scope without executing data processing."""
    return {
        "status": "placeholder_only",
        "single_cell": "future Scanpy QC workflow, manual execution required",
        "microbiome": "future table-level QC workflow, manual execution required",
    }
