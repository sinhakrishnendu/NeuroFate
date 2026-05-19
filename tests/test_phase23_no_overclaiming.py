from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase23_docs_do_not_overclaim():
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 23" in combined
    assert "gse184950" in combined
    assert "gse174367" in combined
    assert "gse147528" in combined
    assert "does not upgrade claims" in combined or "do not upgrade claims" in combined
    assert "is clinical validation" not in combined
    clinical_context = combined.replace("not as clinical validation", "").replace("must not be described as clinical validation", "")
    assert "as clinical validation" not in clinical_context
    assert "is a diagnostic tool" not in combined
    assert "is causal" not in combined
    assert "proves causal" not in combined
    assert "validated across diseases" in combined


def test_phase23_scripts_do_not_download_or_run_sra():
    scripts = [
        ROOT / "scripts/118_triage_replication_cohort_files.py",
        ROOT / "scripts/120_prepare_replication_snrna_axis_plan.py",
        ROOT / "scripts/122_integrate_endpoint_locked_replication.py",
    ]
    for script in scripts:
        text = script.read_text(encoding="utf-8").lower()
        for forbidden in ["wget ", "curl ", "prefetch(", "fasterq", "scanpy", "read_h5ad"]:
            assert forbidden not in text
