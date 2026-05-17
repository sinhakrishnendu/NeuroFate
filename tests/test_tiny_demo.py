from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_tiny_demo_inputs_are_small_and_complete() -> None:
    demo = ROOT / "examples/tiny_demo"
    assert (demo / "tiny_metadata.tsv").exists()
    assert (demo / "tiny_gene_panel.tsv").exists()
    assert (demo / "tiny_sparse_expression.tsv.gz").exists()
    assert (demo / "tiny_sparse_expression.tsv.gz").stat().st_size < 10_000

    with (demo / "tiny_metadata.tsv").open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len({row["donor_id"] for row in rows}) == 4
    assert len({row["cell_type"] for row in rows}) == 3

    with gzip.open(demo / "tiny_sparse_expression.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
        expression_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len({row["gene_symbol"] for row in expression_rows}) == 5


def test_tiny_demo_script_writes_expected_outputs(tmp_path: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            "scripts/55_run_tiny_demo.py",
            "--outdir",
            str(tmp_path),
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "demo_donor_feature_table.tsv").exists()
    assert (tmp_path / "demo_model_metrics.tsv").exists()
    assert (tmp_path / "demo_report.md").exists()
