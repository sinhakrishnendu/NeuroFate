from __future__ import annotations

import gzip
from pathlib import Path

from neurofate.ingest import IngestConfig, detect_gene_identifier_type, run_ingest


def test_ensembl_ids_map_to_neurofate_symbols(tmp_path: Path) -> None:
    expression = tmp_path / "expression.tsv"
    expression.write_text(
        "ensembl_gene_id\tS1\tS2\n"
        "ENSG00000145335\t1\t2\n"
        "ENSG00000131095\t3\t4\n"
        "ENSG00000104725\t5\t6\n",
        encoding="utf-8",
    )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample_id\tdiagnosis\nS1\tControl\nS2\tAD\n", encoding="utf-8")
    result = run_ingest(
        IngestConfig(
            expression=expression,
            metadata=metadata,
            outdir=tmp_path / "out",
            axis_registry=Path("metadata/neurofate_axis_registry.tsv"),
            alias_table=Path("metadata/neurofate_axis_gene_aliases.tsv"),
            min_axis_genes=2,
        )
    )
    with gzip.open(result.standardized_expression, "rt", encoding="utf-8") as handle:
        standardized = handle.read()
    mapping = result.gene_mapping_report.read_text(encoding="utf-8")
    assert detect_gene_identifier_type(["ENSG00000145335", "ENSG00000131095"]) == "ensembl_gene_id"
    assert "SNCA" in mapping
    assert "GFAP" in mapping
    assert "ENSG" not in standardized
