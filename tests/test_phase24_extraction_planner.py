from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/126_plan_gse184950_processed_matrix_extraction.py"


def test_extraction_planner_prefers_processed_10x_and_manual_template():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "per_sample_processed_archive" in text
    assert "tenx_matrix" in text
    assert "manual_extract_per_sample_tar_then_run_127_on_processed_matrices" in text
    assert "manual_phase24_gse184950_axis_extraction_template.sh" in text
    assert "RUN_MANUAL_EXTRACTION" in text


def test_extraction_planner_does_not_process_fastq():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "fastq_processing_status" in text
    assert "not_used" in text
    assert "blocked" in text
    for forbidden in ["fasterq", "prefetch ", "cellranger", "STARsolo".lower()]:
        assert forbidden not in text
