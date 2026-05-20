# NeuroFate CLI Reference

NeuroFate exposes a small stable CLI for research-use transcriptomic neurodegeneration-axis analysis. Historical phase scripts remain in `scripts/` for reproducibility, but the commands below are the public surface used by the Bioinformatics manuscript.

## Stable Commands

### `neurofate check-system`

Reports Python, platform, and optional dependency availability.

### `neurofate doctor`

Checks packaged demo resources and, in a repository checkout, core workflow files.

### `neurofate run-demo`

Runs the bundled synthetic tiny demo without downloads and writes outputs to `results/demo/`.

### `neurofate ingest`

Inspects user-supplied expression and metadata tables, infers common schemas,
validates sample overlap and endpoint labels, maps genes/probes where possible,
and writes standardized NeuroFate-ready inputs.

```bash
neurofate ingest \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_ingest \
  --endpoint-column auto \
  --positive-class auto \
  --negative-class auto
```

Outputs:

- `standardized_expression.tsv.gz`
- `standardized_metadata.tsv`
- `input_schema_detected.tsv`
- `expression_metadata_join.tsv`
- `gene_mapping_report.tsv`
- `ingest_warnings.tsv`
- `ingest_report.md`
- `run_config.yaml`

### `neurofate run`

Runs the complete public workflow:

```text
ingest -> build-axis-scores -> score-risk -> report
```

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_run \
  --endpoint-column auto \
  --positive-class auto \
  --negative-class auto
```

This is the recommended first command for new user-supplied compact expression
tables.

GEO series matrices with a `!series_matrix_table_begin` expression section are
supported directly; users do not need to manually strip the GEO preamble.

The complete workflow writes `run_config.yaml` at the run root and preserves
the ingest-specific configuration under `ingest/run_config.yaml`.

### `neurofate build-axis-scores`

Builds sample-level axis scores from a compact expression table and metadata.

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

Outputs:

- `axis_scores.tsv`
- `axis_feature_coverage.tsv`
- `label_summary.tsv`
- `run_config.yaml`
- `warnings.tsv`

### `neurofate score-risk`

Computes a donor/sample-level research-use risk score from an axis-score table.

```bash
neurofate score-risk \
  --axis-scores results/neurofate_axis/axis_scores.tsv \
  --outdir results/neurofate_axis
```

Outputs:

- `neurofate_risk_scores.tsv`
- `risk_score_report.md`

### `neurofate adapt-endpoint`

Creates explicit endpoint-label aliases for downstream validation scripts that
expect task-specific column names.

```bash
neurofate adapt-endpoint \
  --metadata results/neurofate_run/ingest/standardized_metadata.tsv \
  --endpoint-column label__endpoint \
  --task pd_vs_control \
  --outdir results/neurofate_run/adapted
```

Outputs:

- `adapted_metadata.tsv`
- `endpoint_adapter_report.md`
- `endpoint_aliases.tsv`

The adapter copies binary 0/1 labels only; it does not reinterpret biological
class direction.

## Guarded Research Wrappers

The CLI also exposes guarded wrappers such as `audit-leakage`, `train-baseline`, `train-mps`, `validate-external`, `make-report`, `benchmark`, and `benchmark-report`. These wrappers print dry-run commands by default and require `--run` to execute.

## Research-Use-Only Notice

NeuroFate is intended for research use only. It is not validated for clinical diagnosis, patient-level decision-making, or treatment selection.
