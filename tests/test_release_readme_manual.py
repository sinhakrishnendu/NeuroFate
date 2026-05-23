from __future__ import annotations

from pathlib import Path


def test_release_readme_is_exhaustive_user_manual() -> None:
    readme = Path("README.md").read_text(encoding="utf-8")
    required_sections = [
        "## Research-Use-Only Notice",
        "## Key Features",
        "## Installation",
        "## Quick Start",
        "## Public CLI Overview",
        "## Input Formats",
        "## Metadata Requirements",
        "## Output File Dictionary",
        "## Real-World Example: GSE20141",
        "## NeuroFate Axes",
        "## Reproducibility",
        "## Testing",
        "## Packaging and Release",
        "## Troubleshooting",
        "## Citation",
        "## License",
        "## Contributing",
    ]
    for section in required_sections:
        assert section in readme


def test_release_readme_documents_public_workflow_and_outputs() -> None:
    readme = Path("README.md").read_text(encoding="utf-8")
    required_tokens = [
        "neurofate run-demo",
        "neurofate ingest",
        "neurofate run",
        "neurofate adapt-endpoint",
        "standardized_expression.tsv.gz",
        "standardized_metadata.tsv",
        "gene_mapping_report.tsv",
        "neurofate_risk_scores.tsv",
        "run_config.yaml",
        "GSE20141",
        "18/18",
        "29/30",
        "10/10",
    ]
    for token in required_tokens:
        assert token in readme


def test_release_readme_preserves_review_archive_separation() -> None:
    readme = Path("README.md").read_text(encoding="utf-8")
    checklist = Path("RELEASE_CHECKLIST.md").read_text(encoding="utf-8")
    assert "release_artifacts/neurofate_source_release_<timestamp>.zip" in readme
    assert "release_artifacts/neurofate_results_review_<timestamp>.zip" in readme
    assert "dist/` should contain only `.whl` and `.tar.gz`" in checklist

