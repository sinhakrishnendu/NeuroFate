"""Modeling placeholders with training disabled by default."""

from __future__ import annotations


def training_enabled() -> bool:
    """Training is disabled in the skeleton phase."""
    return False


def model_workflow_status() -> str:
    """Return a status string without initializing ML frameworks."""
    return "placeholder_only_training_disabled"
