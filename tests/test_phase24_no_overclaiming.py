from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase24_docs_are_conservative():
    combined = "\n".join(
        [
            (ROOT / "README.md").read_text(encoding="utf-8"),
            (ROOT / "PNAS_DISCOVERY_STRATEGY.md").read_text(encoding="utf-8"),
            (ROOT / "docs/pnas_validation_strategy.md").read_text(encoding="utf-8"),
            (ROOT / "docs/external_validation_expansion.md").read_text(encoding="utf-8"),
        ]
    ).lower()
    assert "phase 24" in combined
    assert "excel file alone is insufficient" in combined or "workbook alone is insufficient" in combined
    assert "prefer processed 10x matrices" in combined
    assert "fastq/sra processing is avoided" in combined or "avoid fastq/sra processing" in combined
    affirmative_overclaims = [
        "neurofate is clinical-grade",
        "neurofate is a diagnostic tool",
        "phase 24 proves causality",
        "phase 24 is validated across diseases",
        "gse184950 validates neurofate across diseases",
    ]
    for overclaim in affirmative_overclaims:
        assert overclaim not in combined


def test_phase24_report_uses_pending_replication_language():
    text = (ROOT / "scripts/130_generate_phase24_gse184950_replication_report.py").read_text(encoding="utf-8").lower()
    assert "replication is pending" in text
    assert "not a clinical validation dataset" in text
    assert "fastq/sra processing is intentionally avoided" in text
