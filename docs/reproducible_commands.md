# Reproducible Commands

Run commands from the repository root unless noted otherwise.

## Install From Source

```bash
python -m pip install -e .
neurofate check-system
neurofate doctor
```

## Tiny Demo

```bash
neurofate run-demo
```

## GSE20141 Public CLI Smoke Test

The real-world smoke test uses the public GSE20141 GEO series matrix and GPL570
platform annotation. Prepare metadata and the NeuroFate-restricted GPL570 probe
map, then run:

```bash
neurofate run \
  --expression data/raw/end_user_smoke/gse20141/GSE20141_series_matrix.txt.gz \
  --metadata results/end_user_smoke/gse20141/sample_metadata.tsv \
  --gene-map results/end_user_smoke/gse20141/gpl570_axis_probe_mapping.tsv \
  --outdir results/end_user_smoke/gse20141/neurofate_public_run \
  --sample-id-column geo_accession \
  --endpoint-column label__pd_vs_control \
  --positive-class 1 \
  --negative-class 0 \
  --orientation auto \
  --min-axis-genes 10
```

## Endpoint Adapter

```bash
neurofate adapt-endpoint \
  --metadata results/end_user_smoke/gse20141/neurofate_public_run/ingest/standardized_metadata.tsv \
  --endpoint-column label__endpoint \
  --task pd_vs_control \
  --outdir results/end_user_smoke/gse20141/neurofate_public_run/adapted
```

## Tests

```bash
python -m py_compile scripts/*.py neurofate/*.py
python -m pytest \
  tests/test_ingest_geo_series_matrix.py \
  tests/test_ingest_format_detection.py \
  tests/test_ingest_orientation_detection.py \
  tests/test_ingest_gene_identifier_mapping.py \
  tests/test_ingest_expression_metadata_join.py \
  tests/test_neurofate_run_end_to_end.py \
  tests/test_endpoint_adapter.py \
  tests/test_public_cli_reports.py \
  tests/test_research_use_only_outputs.py \
  tests/test_bioinformatics_full_methods_manuscript.py
```

## Package Build

```bash
python -m build --outdir dist_phase42
python -m twine check dist_phase42/*
```

Do not upload to PyPI until release metadata, repository visibility, and final
tagging are confirmed.

## Manuscript Compile

```bash
cd manuscript/bioinformatics
latexmk -pdf neurofate_bioinformatics_full_methods_paper.tex
```

NeuroFate is intended for research use only. It is not validated for clinical
diagnosis, patient-level decision-making, or treatment selection.
