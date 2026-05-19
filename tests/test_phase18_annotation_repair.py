from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/84_build_gse243639_cell_annotation_map.py"


def load_annotation_module():
    spec = importlib.util.spec_from_file_location("phase18_annotation_repair", SCRIPT)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_annotation_repair_normalizes_expression_and_workbook_ids() -> None:
    module = load_annotation_module()
    expression = module.cell_id_parts("s.0096_AAACCCAAGTACGAGC.1")
    workbook = module.cell_id_parts("AAACCCAAGTACGAGC.1", fallback_sample="s.0096")
    core = module.cell_id_parts("AAACCCAAGTACGAGC", fallback_sample="s.0096")
    assert expression["normalized_cell_id"] == "s.0096_AAACCCAAGTACGAGC"
    assert workbook["normalized_cell_id"] == "s.0096_AAACCCAAGTACGAGC"
    assert core["barcode_core"] == "AAACCCAAGTACGAGC"


def test_annotation_repair_outputs_required_columns_and_summary() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for column in [
        "cell_id_expression",
        "cell_id_annotation_original",
        "normalized_cell_id",
        "barcode_core",
        "sample_id",
        "cell_type",
        "annotation_source_sheet",
        "annotation_column_used",
        "match_status",
        "phase18_gse243639_annotation_match_summary.tsv",
        "phase18_gse243639_annotation_column_candidates.tsv",
    ]:
        assert column in text


def test_annotation_repair_celltype_selection_is_conservative() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    for label in ["astrocyte", "microglia", "oligodendrocyte", "dopaminergic", "excitatory", "inhibitory"]:
        assert label in text
    assert "cluster_" in text
    assert "biological_celltype_confidence" in text


def test_annotation_repair_avoids_forbidden_workloads() -> None:
    text = SCRIPT.read_text(encoding="utf-8").lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "h5py", "fit_transform", "leiden", "neighbors"]:
        assert forbidden not in text
