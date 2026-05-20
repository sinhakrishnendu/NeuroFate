from __future__ import annotations

from pathlib import Path

from neurofate.cli import main


def test_neurofate_run_completes_on_gene_row_example(tmp_path: Path) -> None:
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
    assert (tmp_path / "run/ingest/standardized_expression.tsv.gz").exists()
    assert (tmp_path / "run/axis/axis_scores.tsv").exists()
    assert (tmp_path / "run/risk/neurofate_risk_scores.tsv").exists()
    report = (tmp_path / "run/neurofate_run_report.md").read_text(encoding="utf-8")
    assert "research use only" in report.lower()

