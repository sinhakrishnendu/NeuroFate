# Real-World GEO Smoke Test: GSE20141

This smoke test validates the public `neurofate run` workflow on downloaded
public GEO data. It is a software usability test, not a clinical or mechanistic
validation claim.

## Dataset

- GEO accession: GSE20141
- Study: laser-dissected substantia nigra pars compacta neurons in Parkinson's
  disease and control postmortem brain
- Expression input: `GSE20141_series_matrix.txt.gz`
- Platform annotation: `GPL570.annot.gz`
- Endpoint: PD vs control

Downloaded files used in the local smoke test:

- `data/raw/end_user_smoke/gse20141/GSE20141_series_matrix.txt.gz`
- `data/raw/end_user_smoke/gpl570/GPL570.annot.gz`

Checksums:

- GSE20141 series matrix:
  `8975344b5a4715032bd07e08a7a94a68b811fddc59b1fbc53dcf204d1005cf4b`
- GPL570 annotation:
  `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`

## Preparation

Sample metadata were parsed from the GEO series matrix, and GPL570 probe
annotations were restricted to probes mapping to NeuroFate axis genes.

- Samples: 18
- Controls: 8
- PD: 10
- Platform: GPL570
- Retained GPL570 probes: 79

## Public Command

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

## Result

- Run status: passed
- Standardized metadata rows: 18
- Expression samples: 18
- Matched samples: 18/18
- Ambiguous labels: 0
- Retained NeuroFate genes: 29/30
- Axes scored: 10/10
- Research-use risk scores: 18
- Informative warnings: 29/30 axis genes retained, 54,596 input probes not
  mapped to NeuroFate genes, and 20 retained genes represented by multiple
  probes

Main outputs:

- `ingest/standardized_expression.tsv.gz`
- `ingest/standardized_metadata.tsv`
- `ingest/input_schema_detected.tsv`
- `ingest/expression_metadata_join.tsv`
- `ingest/gene_mapping_report.tsv`
- `ingest/ingest_warnings.tsv`
- `axis/axis_scores.tsv`
- `axis/axis_feature_coverage.tsv`
- `axis/label_summary.tsv`
- `risk/neurofate_risk_scores.tsv`
- `risk/risk_score_report.md`
- `neurofate_run_report.md`

## Safety Boundary

No CEL/CHP files, FASTQ/SRA files, H5AD/H5 files, Scanpy, AnnData, UMAP,
clustering, dense genome-wide converted expression output, or model training
were used.

NeuroFate is intended for research use only. It is not validated for clinical
diagnosis, patient-level decision-making, or treatment selection.
