# Quickstart

Install in editable mode from a checkout:

```bash
python -m pip install -e .
```

Start with the lightweight checks:

```bash
neurofate check-system
neurofate doctor
neurofate run-demo
```

The demo writes `results/demo/axis_scores.tsv`, `results/demo/neurofate_risk_scores.tsv`, and Markdown reports without downloads.

For user-supplied compact expression and metadata tables, the easiest route is
the full public workflow:

```bash
neurofate run \
  --expression examples/format_examples/genes_by_samples/expression.tsv \
  --metadata examples/format_examples/genes_by_samples/metadata.tsv \
  --outdir results/neurofate_run
```

This runs format-aware ingestion, endpoint locking, axis scoring, research-use
risk scoring, and report generation.

GEO series matrix files are supported directly when they contain a
`!series_matrix_table_begin` expression section:

```bash
neurofate run \
  --expression GSE20141_series_matrix.txt.gz \
  --metadata sample_metadata.tsv \
  --gene-map gpl570_axis_probe_mapping.tsv \
  --outdir results/neurofate_run \
  --sample-id-column geo_accession \
  --endpoint-column label__pd_vs_control \
  --positive-class 1 \
  --negative-class 0
```

To inspect and standardize inputs without scoring yet:

```bash
neurofate ingest \
  --expression examples/format_examples/samples_by_genes/expression.csv \
  --metadata examples/format_examples/samples_by_genes/metadata.csv \
  --outdir results/neurofate_ingest
```

Build axis scores for a compact sample-level matrix:

```bash
neurofate build-axis-scores \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --axis-registry metadata/neurofate_axis_registry.tsv \
  --sample-id-column sample_id \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control \
  --outdir results/neurofate_axis
```

Compute a research-use risk score:

```bash
neurofate score-risk \
  --axis-scores results/neurofate_axis/axis_scores.tsv \
  --outdir results/neurofate_axis
```

Create endpoint aliases for validation-script compatibility:

```bash
neurofate adapt-endpoint \
  --metadata results/neurofate_run/ingest/standardized_metadata.tsv \
  --endpoint-column label__endpoint \
  --task pd_vs_control \
  --outdir results/neurofate_run/adapted
```

Generate an end-user report from existing outputs:

```bash
python scripts/51_generate_end_user_report.py
python scripts/52_generate_reproducibility_manifest.py
python scripts/54_no_overclaiming_audit.py
```

NeuroFate does not download datasets or run heavy analysis automatically.

NeuroFate is intended for research use only. It is not validated for clinical diagnosis, patient-level decision-making, or treatment selection.
