from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
BLOCKER = ROOT / "scripts/99_build_gse243639_safe_annotation_map_if_valid.py"
README = ROOT / "README.md"
DOC = ROOT / "docs/external_validation_expansion.md"
DECISION = ROOT / "scripts/98_decide_gse243639_annotation_linkage.py"


def test_unsafe_linkage_blocks_annotation_map_creation(tmp_path: Path) -> None:
    decision_path = tmp_path / "decision.tsv"
    blocked_path = tmp_path / "blocked.md"
    output_path = tmp_path / "safe_map.tsv"
    with decision_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "dataset_id",
                "decision_category",
                "safe_to_build_annotation_map",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "dataset_id": "gse243639_pd_snpc",
                "decision_category": "annotation_linkage_unsafe",
                "safe_to_build_annotation_map": "false",
            }
        )
    subprocess.run(
        [
            sys.executable,
            str(BLOCKER),
            "--decision",
            str(decision_path),
            "--output",
            str(output_path),
            "--blocked-output",
            str(blocked_path),
        ],
        check=True,
        cwd=ROOT,
    )
    assert blocked_path.exists()
    assert not output_path.exists()
    assert "not currently supported" in blocked_path.read_text(encoding="utf-8")


def test_docs_do_not_overclaim_celltype_pd_validation() -> None:
    combined = "\n".join([README.read_text(encoding="utf-8"), DOC.read_text(encoding="utf-8")]).lower()
    assert "phase 16 remains the valid global sample-level pd extension" in combined
    assert "should not be interpreted biologically" in combined
    assert "not to force it" in combined
    for forbidden in ["clinical-grade", "diagnostic tool", "validated across diseases", "foundation model"]:
        assert forbidden not in combined


def test_decision_and_blocker_avoid_forbidden_workloads() -> None:
    combined = "\n".join([DECISION.read_text(encoding="utf-8"), BLOCKER.read_text(encoding="utf-8")]).lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "fit_transform", "leiden", "neighbors", "model.fit("]:
        assert forbidden not in combined
