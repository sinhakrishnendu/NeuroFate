"""Feature registry and feature engineering placeholders."""

from __future__ import annotations


def planned_feature_groups() -> list[str]:
    """List planned multimodal feature groups."""
    return [
        "gene_expression",
        "cell_type_abundance",
        "microbiome_signatures",
        "metabolite_signatures",
        "regulatory_network_metrics",
        "protein_interaction_metrics",
        "evolutionary_conservation",
    ]
