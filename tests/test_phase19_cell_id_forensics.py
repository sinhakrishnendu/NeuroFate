from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORENSIC = ROOT / "scripts/95_forensic_gse243639_workbook_cell_ids.py"
NORMALIZATION = ROOT / "scripts/96_deep_gse243639_cell_id_normalization_audit.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_forensic_id_helpers_classify_expected_patterns() -> None:
    module = load_module(FORENSIC, "phase19_forensic")
    assert module.id_type_guess("s.0096_AAACCCAAGTACGAGC.1") == "sample_prefixed_barcode"
    assert module.id_type_guess("AAACCCAAGTACGAGC.1") == "barcode_like"
    assert module.id_type_guess("12345") == "numeric_id"
    assert module.barcode_core("s.0096_AAACCCAAGTACGAGC.1") == "AAACCCAAGTACGAGC"


def test_deep_normalization_rules_are_declared() -> None:
    text = NORMALIZATION.read_text(encoding="utf-8")
    for rule in [
        "raw_id",
        "lowercase",
        "remove_quotes",
        "replace_dash_with_dot",
        "replace_dot_with_dash",
        "remove_trailing_dot_or_dash_one",
        "keep_barcode_only",
        "keep_sample_prefix_only",
        "remove_sample_prefix",
        "collapse_punctuation",
        "seurat_sample_barcode_dot_one",
        "seurat_sample_barcode_dash_one",
        "seurat_barcode_dash_one_sample",
        "seurat_barcode_dot_one_sample",
        "workbook_row_number_mapping_candidate",
        "count_header_column_order_mapping_candidate",
    ]:
        assert rule in text


def test_phase19_forensic_scripts_avoid_forbidden_workloads() -> None:
    combined = "\n".join([FORENSIC.read_text(encoding="utf-8"), NORMALIZATION.read_text(encoding="utf-8")]).lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "fit_transform", "leiden", "neighbors", "model.fit("]:
        assert forbidden not in combined
