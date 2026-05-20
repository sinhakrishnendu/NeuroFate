from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from neurofate.ingest import infer_endpoint_column, infer_positive_negative_classes, standardize_metadata_table


def test_endpoint_and_classes_are_inferred() -> None:
    metadata = pd.DataFrame(
        {"sample_id": ["S1", "S2", "S3", "S4"], "diagnosis": ["Control", "Control", "AD", "AD"]}
    )
    endpoint = infer_endpoint_column(metadata)
    assert endpoint == "diagnosis"
    assert infer_positive_negative_classes(metadata, endpoint) == ("AD", "Control")


def test_ambiguous_endpoint_fails_safely(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "sample_id\tgroup\nS1\tA\nS2\tB\nS3\tC\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="Could not infer"):
        standardize_metadata_table(metadata)

