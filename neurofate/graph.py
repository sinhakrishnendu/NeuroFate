"""Network feature placeholders for GRN, PPI, and pathway graphs."""

from __future__ import annotations


def graph_workflow_status() -> str:
    """Return a status string without constructing large graphs."""
    return "placeholder_only_manual_network_inputs_required"
