from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FILES = [
    ROOT / "MANUSCRIPT_DECISION_DOSSIER.md",
    ROOT / "results/reports/phase39_elife_readiness_report.md",
    ROOT / "manuscript/elife_neurofate_axis_outline.md",
    ROOT / "manuscript/neurofate_elife_manuscript.tex",
]


def test_phase39_manuscript_materials_do_not_make_affirmative_forbidden_claims():
    combined = "\n".join(path.read_text(encoding="utf-8").lower() for path in FILES)
    forbidden_affirmative = [
        "is a clinical diagnostic",
        "is a diagnostic tool",
        "is a clinical biomarker",
        "establishes causality",
        "proves causality",
        "supports a validated shared ad/pd mechanism",
        "establishes a definitive shared ad/pd mechanism",
        "pnas-ready",
    ]
    for phrase in forbidden_affirmative:
        assert phrase not in combined


def test_phase39_pd_divergent_axis_is_not_shared_replication():
    manuscript = (ROOT / "manuscript/neurofate_elife_manuscript.tex").read_text(encoding="utf-8").lower()
    assert "candidate pd-divergent axis" in manuscript
    assert "not shared ad/pd replication" in manuscript
    assert "not a clinical diagnostic tool" in manuscript
    assert "not a causal inference framework" in manuscript
