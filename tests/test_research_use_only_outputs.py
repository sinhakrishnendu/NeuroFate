from __future__ import annotations

from pathlib import Path

from neurofate.cli import main


def test_reports_use_research_language_without_clinical_claims(tmp_path: Path) -> None:
    rc = main(
        [
            "run",
            "--expression",
            "examples/format_examples/microarray_probe_map/expression.tsv",
            "--metadata",
            "examples/format_examples/microarray_probe_map/metadata.tsv",
            "--gene-map",
            "examples/format_examples/microarray_probe_map/probe_map.tsv",
            "--outdir",
            str(tmp_path / "run"),
        ]
    )
    assert rc == 0
    text = "\n".join(
        [
            (tmp_path / "run/ingest/ingest_report.md").read_text(encoding="utf-8"),
            (tmp_path / "run/risk/risk_score_report.md").read_text(encoding="utf-8"),
            (tmp_path / "run/neurofate_run_report.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "research use only" in text
    forbidden_claims = [
        "clinical diagnostic software",
        "patient diagnosis",
        "medical device",
        "treatment recommendation",
        "clinically validated diagnostic biomarker",
    ]
    assert not any(claim in text for claim in forbidden_claims)
