from __future__ import annotations

from pathlib import Path

from neurofate.cli import main


def test_public_run_writes_expected_reports_and_configs(tmp_path: Path) -> None:
    rc = main(
        [
            "run",
            "--expression",
            "examples/format_examples/genes_by_samples/expression.tsv",
            "--metadata",
            "examples/format_examples/genes_by_samples/metadata.tsv",
            "--outdir",
            str(tmp_path / "run"),
        ]
    )
    assert rc == 0
    expected = [
        "ingest/standardized_expression.tsv.gz",
        "ingest/standardized_metadata.tsv",
        "ingest/input_schema_detected.tsv",
        "ingest/expression_metadata_join.tsv",
        "ingest/gene_mapping_report.tsv",
        "ingest/ingest_warnings.tsv",
        "ingest/run_config.yaml",
        "axis/axis_scores.tsv",
        "axis/axis_feature_coverage.tsv",
        "axis/label_summary.tsv",
        "axis/run_config.yaml",
        "risk/neurofate_risk_scores.tsv",
        "risk/risk_score_report.md",
        "neurofate_run_report.md",
        "run_config.yaml",
    ]
    for relative in expected:
        assert (tmp_path / "run" / relative).exists()


def test_public_reports_avoid_unsafe_claims(tmp_path: Path) -> None:
    rc = main(
        [
            "run",
            "--expression",
            "examples/format_examples/genes_by_samples/expression.tsv",
            "--metadata",
            "examples/format_examples/genes_by_samples/metadata.tsv",
            "--outdir",
            str(tmp_path / "run"),
        ]
    )
    assert rc == 0
    text = "\n".join(
        path.read_text(encoding="utf-8")
        for path in [
            tmp_path / "run/ingest/ingest_report.md",
            tmp_path / "run/risk/risk_score_report.md",
            tmp_path / "run/neurofate_run_report.md",
        ]
    ).lower()
    assert "research use only" in text
    forbidden = [
        "clinical diagnostic software",
        "patient diagnosis",
        "medical device",
        "treatment recommendation",
        "clinically validated diagnostic biomarker",
    ]
    assert not any(phrase in text for phrase in forbidden)
