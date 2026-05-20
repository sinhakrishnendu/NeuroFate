from __future__ import annotations

import gzip
from pathlib import Path

from neurofate.ingest import detect_delimiter, infer_expression_orientation, inspect_table_shape, is_geo_series_matrix


def test_geo_series_matrix_table_is_detected_after_preamble(tmp_path: Path) -> None:
    series = tmp_path / "GSE00000_series_matrix.txt.gz"
    with gzip.open(series, "wt", encoding="utf-8") as handle:
        handle.write('!Series_title\t"Example GEO dataset"\n')
        handle.write('!Sample_title\t"S1"\t"S2"\n')
        handle.write("!series_matrix_table_begin\n")
        handle.write('"ID_REF"\t"GSM1"\t"GSM2"\n')
        handle.write('"1007_s_at"\t1.0\t2.0\n')
        handle.write('"1053_at"\t3.0\t4.0\n')
        handle.write("!series_matrix_table_end\n")

    assert is_geo_series_matrix(series)
    assert detect_delimiter(series) == "\t"
    assert infer_expression_orientation(series) == "genes_rows"
    shape = inspect_table_shape(series)
    assert shape["columns"] == 3
    assert shape["preview_rows"] == 2
    assert "ID_REF" in shape["column_names"]
