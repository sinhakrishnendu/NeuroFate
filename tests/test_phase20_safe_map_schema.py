from __future__ import annotations

import csv
import gzip
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts/100_audit_phase19_safe_annotation_map_schema.py"
BUILDER99 = ROOT / "scripts/99_build_gse243639_safe_annotation_map_if_valid.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_safe_map_required_columns_are_declared() -> None:
    text = BUILDER99.read_text(encoding="utf-8")
    for column in [
        "cell_id_expression",
        "cell_id_annotation_original",
        "cell_id_annotation_normalized",
        "normalized_cell_id",
        "barcode_core",
        "sample_id",
        "cell_type",
        "cluster_id",
        "annotation_source_sheet",
        "annotation_column_used",
        "match_status",
        "normalization_rule",
    ]:
        assert column in text


def test_schema_audit_recommends_direct_expression_column(tmp_path: Path) -> None:
    module = load_module(AUDIT, "phase20_schema_audit")
    expression = tmp_path / "expr.tsv.gz"
    with gzip.open(expression, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_id", "sample_id", "gene_symbol", "expression_value"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"cell_id": "s.1_AAA.1", "sample_id": "s.1", "gene_symbol": "SNCA", "expression_value": "1"})
    safe_map = tmp_path / "safe.tsv"
    with safe_map.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=module.REQUIRED_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "cell_id_expression": "s.1_AAA.1",
                "cell_id_annotation_original": "s-1_AAA-1",
                "cell_id_annotation_normalized": "s.1_AAA.1",
                "normalized_cell_id": "s.1_AAA.1",
                "barcode_core": "AAA",
                "sample_id": "s.1",
                "cell_type": "cluster_1",
                "cluster_id": "1",
                "annotation_source_sheet": "Sheet1",
                "annotation_column_used": "CLUSTER",
                "match_status": "matched",
                "normalization_rule": "replace_dash_with_dot",
            }
        )
    rows, preview, _fields = module.build_audit(safe_map, expression)
    assert preview[0]["cell_id_expression"] == "s.1_AAA.1"
    assert any(row["recommended_join_column"] == "cell_id_expression" for row in rows)
    assert any(row["audit_item"] == "missing_required_columns" and row["value"] == "" for row in rows)


def test_phase20_schema_scripts_avoid_forbidden_workloads() -> None:
    combined = "\n".join([AUDIT.read_text(encoding="utf-8"), BUILDER99.read_text(encoding="utf-8")]).lower()
    for forbidden in ["scanpy", "anndata", "read_h5ad", "fit_transform", "leiden", "neighbors"]:
        assert forbidden not in combined
