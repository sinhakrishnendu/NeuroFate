from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from neurofate.cli import main
from neurofate.ingest import detect_delimiter, infer_expression_orientation, is_geo_series_matrix


def test_release_run_workflow_creates_expected_outputs(tmp_path: Path) -> None:
    outdir = tmp_path / "run"
    rc = main(
        [
            "run",
            "--expression",
            "examples/format_examples/genes_by_samples/expression.tsv",
            "--metadata",
            "examples/format_examples/genes_by_samples/metadata.tsv",
            "--outdir",
            str(outdir),
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
        "axis/axis_scores.tsv",
        "axis/axis_feature_coverage.tsv",
        "axis/label_summary.tsv",
        "risk/neurofate_risk_scores.tsv",
        "risk/risk_score_report.md",
        "neurofate_run_report.md",
        "run_config.yaml",
    ]
    for relpath in expected:
        assert (outdir / relpath).exists(), relpath
    assert "research use only" in (outdir / "neurofate_run_report.md").read_text(encoding="utf-8").lower()


def test_release_ingest_handles_geo_series_matrix_preamble(tmp_path: Path) -> None:
    series = tmp_path / "GSE_release_series_matrix.txt.gz"
    with gzip.open(series, "wt", encoding="utf-8") as handle:
        handle.write('!Series_title\t"Release test"\n')
        handle.write("!series_matrix_table_begin\n")
        handle.write('"ID_REF"\t"S1"\t"S2"\n')
        handle.write('"SNCA"\t1.0\t2.0\n')
        handle.write('"MAPT"\t3.0\t4.0\n')
        handle.write('"APP"\t5.0\t6.0\n')
        handle.write("!series_matrix_table_end\n")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample_id\tdiagnosis\nS1\tControl\nS2\tPD\n", encoding="utf-8")

    assert is_geo_series_matrix(series)
    assert detect_delimiter(series) == "\t"
    assert infer_expression_orientation(series) == "genes_rows"

    outdir = tmp_path / "ingest"
    rc = main(
        [
            "ingest",
            "--expression",
            str(series),
            "--metadata",
            str(metadata),
            "--outdir",
            str(outdir),
            "--endpoint-column",
            "diagnosis",
            "--positive-class",
            "PD",
            "--negative-class",
            "Control",
            "--min-axis-genes",
            "2",
        ]
    )
    assert rc == 0
    expression = pd.read_csv(outdir / "standardized_expression.tsv.gz", sep="\t")
    assert set(expression["gene_symbol"]) <= {"SNCA", "MAPT", "APP"}
    join = pd.read_csv(outdir / "expression_metadata_join.tsv", sep="\t")
    assert int(join.loc[0, "matched_sample_count"]) == 2
