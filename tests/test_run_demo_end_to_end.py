from __future__ import annotations

import subprocess
import sys
from pathlib import Path
import os


REPO_ROOT = Path(__file__).resolve().parents[1]


def python_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    return env


def test_run_demo_end_to_end(tmp_path: Path) -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "neurofate", "run-demo"],
        cwd=tmp_path,
        env=python_env(),
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    outdir = tmp_path / "results" / "demo"
    expected = [
        "demo_donor_feature_table.tsv",
        "demo_model_metrics.tsv",
        "axis_scores.tsv",
        "neurofate_risk_scores.tsv",
        "risk_score_report.md",
        "demo_report.md",
    ]
    for name in expected:
        assert (outdir / name).is_file()
    assert "research use only" in (outdir / "demo_report.md").read_text().lower()


def test_build_axis_scores_and_score_risk_cli(tmp_path: Path) -> None:
    expression = tmp_path / "expression.tsv"
    metadata = tmp_path / "metadata.tsv"
    registry = tmp_path / "axis.tsv"
    expression.write_text(
        "sample_id\tAPOE\tGFAP\tSLC17A7\tSNCA\n"
        "S1\t1.0\t0.5\t2.0\t0.2\n"
        "S2\t2.0\t1.2\t1.0\t0.8\n",
        encoding="utf-8",
    )
    metadata.write_text("sample_id\tdiagnosis\nS1\tControl\nS2\tAD\n", encoding="utf-8")
    registry.write_text(
        "axis_id\taxis_name\tbiological_theme\tgene_members\n"
        "demo_axis\tDemo axis\tdemo\tAPOE;GFAP;SNCA\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "axis_out"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "neurofate",
            "build-axis-scores",
            "--expression",
            str(expression),
            "--metadata",
            str(metadata),
            "--axis-registry",
            str(registry),
            "--sample-id-column",
            "sample_id",
            "--endpoint-column",
            "diagnosis",
            "--positive-class",
            "AD",
            "--negative-class",
            "Control",
            "--outdir",
            str(outdir),
        ],
        text=True,
        env=python_env(),
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert (outdir / "axis_scores.tsv").is_file()
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "neurofate",
            "score-risk",
            "--axis-scores",
            str(outdir / "axis_scores.tsv"),
            "--outdir",
            str(outdir),
        ],
        text=True,
        env=python_env(),
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert (outdir / "neurofate_risk_scores.tsv").is_file()
