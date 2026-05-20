# Tutorial: NeuroFate Axis Scoring

This tutorial demonstrates the public software workflow for user-supplied donor/sample-level expression data.

## 1. Prepare Inputs

Expression can be provided as sample rows:

```text
sample_id  APOE  GFAP  SNCA  SLC17A7
S1         1.2   0.8   0.4   2.1
S2         2.0   1.5   0.9   1.2
```

or gene rows:

```text
gene_symbol  S1   S2
APOE         1.2  2.0
GFAP         0.8  1.5
```

Metadata must contain a sample identifier and endpoint column:

```text
sample_id  diagnosis
S1         Control
S2         AD
```

## 2. Build Axis Scores

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

## 3. Score Research Risk

```bash
neurofate score-risk \
  --axis-scores results/neurofate_axis/axis_scores.tsv \
  --outdir results/neurofate_axis
```

## 4. Interpret Conservatively

Axis scores and risk scores are exploratory research outputs. They can support cohort-level stratification, validation, and methods development. They must not be interpreted as clinical diagnosis, patient-level decision-making, or treatment selection.

