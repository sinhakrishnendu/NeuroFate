from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/reports/phase39_elife_readiness_report.md"
DOSSIER = ROOT / "MANUSCRIPT_DECISION_DOSSIER.md"
MANUSCRIPT = ROOT / "manuscript/neurofate_elife_manuscript.tex"


def test_phase39_readiness_report_targets_elife_not_pnas_readiness():
    text = REPORT.read_text(encoding="utf-8").lower()
    assert "elife" in text
    assert "pnas is premature" in text
    assert "pnas-ready" not in text
    assert "definitive shared" not in text
    assert "candidate pd-divergent" in text


def test_phase39_dossier_and_manuscript_preserve_conservative_main_claim():
    dossier = DOSSIER.read_text(encoding="utf-8")
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    assert "Primary target: eLife" in dossier
    assert "NeuroFate-Axis is an endpoint-locked" in dossier
    assert "\\documentclass[9pt,lineno]{eLife_LaTeX_template/elife}" in manuscript
    assert "AD-replicated neuronal vulnerability" in manuscript
    assert "candidate PD-divergent" in manuscript
    assert "not shared AD/PD replication" in manuscript
