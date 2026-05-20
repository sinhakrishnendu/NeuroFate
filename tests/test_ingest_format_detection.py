from __future__ import annotations

import gzip
from pathlib import Path

from neurofate.ingest import detect_compression, detect_delimiter, inspect_table_shape


def test_detects_csv_tsv_and_gzip(tmp_path: Path) -> None:
    csv_path = tmp_path / "expression.csv"
    csv_path.write_text("gene_symbol,S1,S2\nSNCA,1,2\n", encoding="utf-8")
    tsv_gz = tmp_path / "expression.tsv.gz"
    with gzip.open(tsv_gz, "wt", encoding="utf-8") as handle:
        handle.write("gene_symbol\tS1\tS2\nGFAP\t1\t2\n")

    assert detect_delimiter(csv_path) == ","
    assert detect_delimiter(tsv_gz) == "\t"
    assert detect_compression(tsv_gz) == "gzip"
    shape = inspect_table_shape(tsv_gz)
    assert shape["columns"] == 3
    assert shape["compression"] == "gzip"

