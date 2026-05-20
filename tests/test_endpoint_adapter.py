from __future__ import annotations

from pathlib import Path

import pandas as pd

from neurofate.adapters import adapt_endpoint_metadata, ensure_label_endpoint_aliases
from neurofate.cli import main


def test_endpoint_adapter_maps_public_label_to_pd_alias(tmp_path: Path) -> None:
    metadata = tmp_path / "standardized_metadata.tsv"
    metadata.write_text(
        "sample_id\tlabel__endpoint\tresearch_use_only\nS1\t0\ttrue\nS2\t1\ttrue\n",
        encoding="utf-8",
    )

    result = adapt_endpoint_metadata(metadata, tmp_path / "adapted", task="pd_vs_control")
    adapted = pd.read_csv(result.adapted_metadata, sep="\t", dtype=str)
    assert adapted["label__pd_vs_control"].tolist() == ["0", "1"]
    assert adapted["endpoint_label"].tolist() == ["0", "1"]
    assert adapted["label"].tolist() == ["0", "1"]
    report = result.endpoint_adapter_report.read_text(encoding="utf-8").lower()
    assert "research use only" in report
    assert "no biological label direction was changed" in report


def test_endpoint_adapter_refuses_ambiguous_labels() -> None:
    frame = pd.DataFrame({"sample_id": ["S1", "S2"], "label__endpoint": ["Control", "maybe"]})
    try:
        ensure_label_endpoint_aliases(frame, task="ad_vs_control")
    except ValueError as exc:
        assert "not unambiguous binary labels" in str(exc)
    else:  # pragma: no cover - explicit failure branch for clarity
        raise AssertionError("adapter accepted an ambiguous label")


def test_adapt_endpoint_cli_creates_outputs(tmp_path: Path) -> None:
    metadata = tmp_path / "standardized_metadata.tsv"
    metadata.write_text("sample_id\tlabel__endpoint\nS1\t0\nS2\t1\n", encoding="utf-8")
    rc = main(
        [
            "adapt-endpoint",
            "--metadata",
            str(metadata),
            "--outdir",
            str(tmp_path / "adapted"),
            "--task",
            "pd_vs_control",
        ]
    )
    assert rc == 0
    assert (tmp_path / "adapted/adapted_metadata.tsv").exists()
    assert (tmp_path / "adapted/endpoint_adapter_report.md").exists()
    assert (tmp_path / "adapted/endpoint_aliases.tsv").exists()
