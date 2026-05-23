from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


def test_release_overclaiming_audit_tracks_required_unsafe_phrases() -> None:
    source = Path("scripts/54_no_overclaiming_audit.py").read_text(encoding="utf-8")
    required_phrases = [
        "diagnostic tool",
        "clinical diagnosis",
        "patient diagnosis",
        "medical device",
        "causal",
        "biomarker",
        "treatment recommendation",
        "validated shared mechanism",
    ]
    for phrase in required_phrases:
        assert f'"phrase": "{phrase}"' in source or f'"{phrase}"' in source


def test_release_no_overclaiming_audit_has_zero_high_flags(tmp_path: Path) -> None:
    output = tmp_path / "no_overclaiming.tsv"
    completed = subprocess.run(
        [sys.executable, "scripts/54_no_overclaiming_audit.py", "--output", str(output)],
        check=False,
        text=True,
        capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr
    rows = list(csv.DictReader(output.open(encoding="utf-8"), delimiter="\t"))
    high = [row for row in rows if row["severity"] == "high"]
    assert high == []


def test_release_docs_use_research_use_only_boundary() -> None:
    docs = "\n".join(
        [
            Path("README.md").read_text(encoding="utf-8"),
            Path("docs/research_use_only.md").read_text(encoding="utf-8"),
            Path("manuscript/bioinformatics/neurofate_bioinformatics_full_methods_paper.tex").read_text(
                encoding="utf-8"
            ),
        ]
    ).lower()
    assert "research use only" in docs
    forbidden_affirmative = [
        "is a clinical diagnostic software",
        "is a diagnostic tool",
        "is a medical device",
        "provides patient diagnosis",
        "provides treatment recommendation",
        "is a validated shared mechanism",
    ]
    for phrase in forbidden_affirmative:
        assert phrase not in docs

