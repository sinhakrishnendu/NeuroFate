"""Reporting placeholders."""

from __future__ import annotations


def report_sections() -> list[str]:
    """List planned reporting sections without generating files."""
    return [
        "dataset_registry",
        "quality_control_summary",
        "multimodal_features",
        "model_performance",
        "interpretability",
    ]
