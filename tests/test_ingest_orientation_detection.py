from __future__ import annotations

from pathlib import Path

from neurofate.ingest import infer_expression_orientation


def test_orientation_detection_for_common_formats(tmp_path: Path) -> None:
    genes_rows = tmp_path / "genes.tsv"
    genes_rows.write_text("gene_symbol\tS1\tS2\nSNCA\t1\t2\n", encoding="utf-8")
    samples_rows = tmp_path / "samples.tsv"
    samples_rows.write_text("sample_id\tSNCA\tGFAP\nS1\t1\t2\n", encoding="utf-8")
    long_format = tmp_path / "long.tsv"
    long_format.write_text(
        "sample_id\tgene_symbol\texpression_value\nS1\tSNCA\t1\n",
        encoding="utf-8",
    )

    assert infer_expression_orientation(genes_rows) == "genes_rows"
    assert infer_expression_orientation(samples_rows) == "samples_rows"
    assert infer_expression_orientation(long_format) == "long"

