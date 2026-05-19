import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "scripts/150_convert_gse174367_bulk_rda_to_axis_matrix.py"
BUILDER = ROOT / "scripts/151_build_gse174367_bulk_axis_scores.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_axis_gene_coverage_reports_only_axis_genes():
    module = load(CONVERTER, "phase30_axis_converter")
    coverage = module.coverage_rows(
        [{"axis_id": "synuclein_mitochondrial_axis", "gene_members": "SNCA;PINK1;PRKN"}],
        {"SNCA", "PRKN"},
    )
    assert coverage[0]["genes_found"] == "2"
    assert coverage[0]["genes_missing"] == "1"
    assert coverage[0]["missing_gene_members"] == "PINK1"


def test_axis_score_builder_uses_ad_control_labels_only():
    module = load(BUILDER, "phase30_axis_builder")
    matrix_values = {
        "S1": {"SNCA": 1.0, "MAPT": 2.0},
        "S2": {"SNCA": 2.0, "MAPT": 3.0},
        "S3": {"SNCA": 4.0, "MAPT": 5.0},
    }
    sample_map = [
        {"expression_sample_id": "S1", "diagnosis": "Control", "inferred_ad_endpoint": "Control", "label__ad_vs_control": "0", "match_status": "matched"},
        {"expression_sample_id": "S2", "diagnosis": "AD", "inferred_ad_endpoint": "AD", "label__ad_vs_control": "1", "match_status": "matched"},
        {"expression_sample_id": "S3", "diagnosis": "Other", "inferred_ad_endpoint": "Other", "label__ad_vs_control": "", "match_status": "matched"},
    ]
    axes = [{"axis_id": "amyloid_tau_axis", "gene_members": "SNCA;MAPT"}]
    rows, _coverage, labels = module.build_scores(matrix_values, sample_map, axes)
    assert [row["sample_id"] for row in rows] == ["S1", "S2"]
    assert labels == [{"label__ad_vs_control": "0", "count": "1"}, {"label__ad_vs_control": "1", "count": "1"}]


def test_converter_source_avoids_genome_wide_matrix_claim():
    text = CONVERTER.read_text(encoding="utf-8").lower()
    assert "axis_genes" in text
    assert "axis_df" in text
    assert "write.table(axis_df" in text
    assert "write.table(expr" not in text
