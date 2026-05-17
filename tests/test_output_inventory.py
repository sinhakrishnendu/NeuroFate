from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_output_inventory_script_has_expected_columns() -> None:
    source = (ROOT / "scripts/56_inventory_outputs.py").read_text(encoding="utf-8")
    expected = ["path", "size_bytes", "modified_time", "artifact_type", "phase", "user_facing"]
    for token in expected:
        assert token in source
    forbidden = ["scanpy", "read_h5ad", "anndata", "h5py", "torch", "fit("]
    for token in forbidden:
        assert token not in source.lower()


def test_output_inventory_writes_report(tmp_path: Path) -> None:
    output = tmp_path / "output_inventory.tsv"
    result = subprocess.run(
        [sys.executable, "scripts/56_inventory_outputs.py", "--output", str(output)],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with output.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == [
            "path",
            "size_bytes",
            "modified_time",
            "artifact_type",
            "phase",
            "user_facing",
        ]
