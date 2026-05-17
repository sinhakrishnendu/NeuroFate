"""Single-cell workflow placeholders.

HEAVY future work in this module may involve Scanpy, AnnData, scVI, RNA velocity,
or large HDF5 processing. Nothing here launches those workflows automatically.
"""

from __future__ import annotations


def singlecell_workflow_status() -> str:
    """Return a status string without importing Scanpy or reading data."""
    return "placeholder_only_manual_execution_required"
