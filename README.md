# NeuroFate

**NeuroFate: format-aware command-line software for endpoint-locked transcriptomic neurodegeneration risk scoring.**

NeuroFate is a Python command-line research software package for reproducible donor/sample-level transcriptomic neurodegeneration-axis analysis. It inspects user-supplied expression and metadata tables, detects common table layouts, harmonizes gene/probe identifiers, locks endpoints before scoring, builds curated neurodegeneration-axis scores, writes research-use risk scores, and creates auditable reports.

Repository: <https://github.com/sinhakrishnendu/NeuroFate.git>

Current release-candidate version: **0.3.0**

## Research-Use-Only Notice

NeuroFate is intended for research use only. It is not a clinical biomarker and is not validated for clinical diagnosis, patient-level decision-making, treatment selection, or care-delivery use. NeuroFate outputs are intended for cohort-level transcriptomic research, diagnosis-oriented research, endpoint-locked disease-state modelling, and reproducible software demonstrations.

## Key Features

- CLI/PyPI-ready package with the console command `neurofate`.
- Format-aware ingestion through `neurofate ingest`.
- Complete public workflow through `neurofate run`.
- GEO series matrix support through direct parsing of `!series_matrix_table_begin` expression sections.
- CSV, TSV, TXT, and `.gz` input support.
- Genes-by-samples, samples-by-genes, and long-format expression support.
- Ensembl ID, gene-symbol, and microarray probe mapping support.
- Endpoint locking with explicit positive and negative classes.
- Curated NeuroFate axis scoring.
- Research-use risk scoring and Markdown reports.
- Leakage-audit and no-overclaiming audit scripts for repository-level validation.
- Endpoint adapter for compatibility between public CLI outputs and validation scripts.
- Real-world public GEO smoke test using GSE20141 and GPL570.
- Buildable wheel/sdist artifacts and reviewer-facing manuscript assets.

## Installation

### Install From PyPI

After public release:

```bash
python -m pip install neurofate
```

### Install From GitHub

```bash
python -m pip install git+https://github.com/sinhakrishnendu/NeuroFate.git
```

### Developer Install

```bash
git clone https://github.com/sinhakrishnendu/NeuroFate.git
cd NeuroFate/NeuroFate
python -m pip install -e ".[dev]"
```

### Optional Extras

```bash
python -m pip install -e ".[plotting]"
python -m pip install -e ".[docs]"
python -m pip install -e ".[mps]"
python -m pip install -e ".[dev]"
```

The default package does not require Scanpy, AnnData, PyTorch, or matplotlib. PyTorch/MPS and plotting dependencies are optional.

## Quick Start

Check the installation:

```bash
neurofate check-system
neurofate doctor
```

Run the bundled no-download demo:

```bash
neurofate run-demo
```

Run the full public workflow on a compact user dataset:

```bash
neurofate run \
  --expression examples/format_examples/genes_by_samples/expression.tsv \
  --metadata examples/format_examples/genes_by_samples/metadata.tsv \
  --outdir results/neurofate_run
```

Expected top-level outputs include:

- `ingest/standardized_expression.tsv.gz`
- `ingest/standardized_metadata.tsv`
- `axis/axis_scores.tsv`
- `axis/axis_feature_coverage.tsv`
- `axis/label_summary.tsv`
- `risk/neurofate_risk_scores.tsv`
- `risk/risk_score_report.md`
- `neurofate_run_report.md`
- `run_config.yaml`

## Public CLI Overview

Stable user-facing commands:

```bash
neurofate check-system
neurofate doctor
neurofate run-demo
neurofate ingest
neurofate build-axis-scores
neurofate score-risk
neurofate run
neurofate adapt-endpoint
```

### `neurofate check-system`

Reports Python version, platform, and optional dependency availability.

### `neurofate doctor`

Checks packaged resources and, in a repository checkout, core project files.

### `neurofate run-demo`

Runs a small synthetic dataset without downloads and writes demo outputs under `results/demo/`.

### `neurofate ingest`

Inspects expression and metadata tables, infers format, validates sample overlap and endpoint labels, maps genes/probes, writes standardized inputs, and reports warnings.

```bash
neurofate ingest \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_ingest
```

### `neurofate build-axis-scores`

Builds sample-level NeuroFate axis scores from compact or standardized inputs.

```bash
neurofate build-axis-scores \
  --expression results/neurofate_ingest/standardized_expression.tsv.gz \
  --metadata results/neurofate_ingest/standardized_metadata.tsv \
  --axis-registry metadata/neurofate_axis_registry.tsv \
  --sample-id-column sample_id \
  --endpoint-column label__endpoint \
  --positive-class 1 \
  --negative-class 0 \
  --outdir results/neurofate_axis
```

### `neurofate score-risk`

Computes an exploratory research-use score from axis scores.

```bash
neurofate score-risk \
  --axis-scores results/neurofate_axis/axis_scores.tsv \
  --outdir results/neurofate_axis
```

### `neurofate run`

Runs the complete public workflow:

```text
ingest -> build-axis-scores -> score-risk -> report
```

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_run
```

### `neurofate adapt-endpoint`

Creates explicit endpoint aliases for validation scripts that expect task-specific label columns.

```bash
neurofate adapt-endpoint \
  --metadata results/neurofate_run/ingest/standardized_metadata.tsv \
  --endpoint-column label__endpoint \
  --task pd_vs_control \
  --outdir results/neurofate_run/adapted
```

Outputs:

- `adapted_metadata.tsv`
- `endpoint_aliases.tsv`
- `endpoint_adapter_report.md`

The adapter copies binary 0/1 labels only. It does not reinterpret biological class direction.

### `neurofate make-report`

`make-report` is a guarded repository workflow for generating reports from existing project outputs. It is useful in the full repository checkout but is not required for the public ingest/run workflow.

Advanced or experimental commands such as `train-baseline`, `train-mps`, `validate-external`, `benchmark`, and historical phase wrappers are retained for reproducibility. They are not the recommended first commands for new users.

## Input Formats

NeuroFate public ingestion accepts compact text tables. It does not process raw FASTQ/FQ, SRA, CEL/CHP, H5AD/AnnData, or HDF5 single-cell containers.

### Genes-by-Samples Matrix

```text
gene_symbol    S01    S02    S03
SNCA           0.2    0.4    0.8
GFAP           0.1    0.3    1.1
NEFL           1.2    1.0    0.7
```

### Samples-by-Genes Matrix

```text
sample_id    SNCA    GFAP    NEFL
S01          0.2     0.1     1.2
S02          0.4     0.3     1.0
S03          0.8     1.1     0.7
```

### Long Format

```text
sample_id    gene_symbol    expression_value
S01          SNCA           0.2
S01          GFAP           0.1
S02          SNCA           0.4
```

### GEO Series Matrix

```text
!Series_title    "Example GEO dataset"
!Sample_geo_accession    "GSM1"    "GSM2"
!series_matrix_table_begin
"ID_REF"    "GSM1"    "GSM2"
"1007_s_at"    1.2    1.5
!series_matrix_table_end
```

NeuroFate reads the expression table between `!series_matrix_table_begin` and `!series_matrix_table_end`. Supply a separate metadata table with sample identifiers matching the expression columns.

### Ensembl-ID Matrix

```text
ensembl_gene_id    S01    S02
ENSG00000145335    0.2    0.4
ENSG00000131095    0.1    0.3
```

NeuroFate maps curated axis genes using `metadata/neurofate_axis_gene_aliases.tsv`.

### Microarray Probe Matrix With Gene Map

Expression:

```text
ID_REF       GSM1    GSM2
probe_SNCA   0.2     0.4
probe_GFAP   0.1     0.3
```

Probe map:

```text
probe_id     gene_symbol
probe_SNCA   SNCA
probe_GFAP   GFAP
```

Command:

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --gene-map probe_map.tsv \
  --outdir results/neurofate_run
```

Compressed `.gz` files are supported for CSV, TSV, TXT, and GEO series matrix inputs.

## Metadata Requirements

Metadata must contain:

- A sample identifier column such as `sample_id`, `geo_accession`, `donor_id`, `subject_id`, or `participant_id`.
- An endpoint column such as `diagnosis`, `disease_state`, `condition`, `group`, `status`, `phenotype`, or `label`.
- Positive and negative classes, either inferred or passed explicitly.

Example:

```text
sample_id    diagnosis
S01          Control
S02          AD
S03          AD
```

Explicit endpoint locking:

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control \
  --outdir results/neurofate_run
```

Endpoint locking ensures the disease-state contrast is defined before score interpretation. NeuroFate does not scan all metadata labels to choose the strongest result.

Optional covariates such as age, sex, postmortem interval, brain region, and batch can be retained in source metadata, but the public axis-scoring workflow uses only the locked endpoint label and expression values.

## Output File Dictionary

`neurofate ingest` writes:

- `standardized_expression.tsv.gz`: NeuroFate axis-gene expression matrix with genes as rows and samples as columns.
- `standardized_metadata.tsv`: standardized sample metadata with `sample_id`, `endpoint`, `label__endpoint`, and `research_use_only`.
- `input_schema_detected.tsv`: detected delimiter, orientation, endpoint settings, feature counts, and retained genes.
- `expression_metadata_join.tsv`: expression/metadata sample-overlap audit.
- `gene_mapping_report.tsv`: input feature mapping and retention status.
- `ingest_warnings.tsv`: non-fatal warnings.
- `ingest_report.md`: human-readable ingest report.
- `run_config.yaml`: reproducibility settings for ingestion.

`neurofate run` additionally writes:

- `axis/axis_scores.tsv`: sample-level axis scores.
- `axis/axis_feature_coverage.tsv`: mapped and missing genes per axis.
- `axis/label_summary.tsv`: locked endpoint label counts.
- `axis/warnings.tsv`: scoring warnings.
- `risk/neurofate_risk_scores.tsv`: exploratory research-use sample scores.
- `risk/risk_score_report.md`: risk-score report with research-use-only notice.
- `neurofate_run_report.md`: complete workflow report.
- `run_config.yaml`: top-level workflow configuration.

`neurofate adapt-endpoint` writes:

- `adapted_metadata.tsv`: standardized metadata plus endpoint aliases.
- `endpoint_aliases.tsv`: alias mapping audit.
- `endpoint_adapter_report.md`: human-readable adapter report.

## Real-World Example: GSE20141

GSE20141 is a public GEO laser-dissected substantia nigra pars compacta microarray cohort for Parkinson's disease versus control research. The final public CLI smoke test used:

- `GSE20141_series_matrix.txt.gz`
- `GPL570.annot.gz`
- parsed sample metadata
- GPL570 NeuroFate axis probe map

Command:

```bash
neurofate run \
  --expression data/raw/end_user_smoke/gse20141/GSE20141_series_matrix.txt.gz \
  --metadata results/end_user_smoke/gse20141/sample_metadata.tsv \
  --gene-map results/end_user_smoke/gse20141/gpl570_axis_probe_mapping.tsv \
  --outdir results/end_user_smoke/gse20141/neurofate_public_run_final \
  --sample-id-column geo_accession \
  --endpoint-column label__pd_vs_control \
  --positive-class 1 \
  --negative-class 0 \
  --orientation auto \
  --min-axis-genes 10
```

Result:

- Run status: passed.
- Samples matched: 18/18.
- Label counts: 10 PD and 8 controls.
- Retained NeuroFate genes: 29/30.
- Retained GPL570 probes: 79.
- Axes scored: 10/10.
- Research-use risk scores generated for 18 samples.
- No fatal ingest errors.
- Informative warnings: incomplete axis-gene coverage (29/30), unmapped non-axis probes, and multiple probes per retained gene.

Outputs are written under:

```text
results/end_user_smoke/gse20141/neurofate_public_run_final/
```

Detailed smoke-test documentation:

```text
docs/real_world_geo_smoke_test_gse20141.md
results/reports/final_gse20141_public_cli_smoke_test.md
```

## NeuroFate Axes

The default axis registry is stored in `metadata/neurofate_axis_registry.tsv` and bundled as package data.

- `neuronal_vulnerability_axis`: inhibitory/excitatory neuronal vulnerability markers and neurofilament genes.
- `synuclein_mitochondrial_axis`: synuclein, mitochondrial stress, and PD-relevant genes.
- `astrocyte_stress_axis`: astrocyte activation and stress-associated genes.
- `inflammatory_microglial_axis`: microglial and inflammatory response genes.
- `myelin_oligodendrocyte_axis`: myelin and oligodendrocyte genes.
- `proteostasis_autophagy_axis`: proteostasis, autophagy, and lysosomal/mitochondrial stress genes.
- `amyloid_tau_axis`: amyloid, presenilin, tau, and APOE-related genes.
- `immune_antigen_presentation_axis`: immune and antigen-presentation genes.
- `vascular_barrier_axis`: vascular, barrier, and inflammatory interaction genes.
- `global_neurodegeneration_axis`: broad neurodegeneration-associated axis.

Axes are research summaries of available expression features. They are not by themselves causal mechanisms or care-delivery tools.

## Reproducibility

Install from source:

```bash
python -m pip install -e .
```

Run the demo:

```bash
neurofate run-demo
```

Run the real GEO smoke test after acquiring the public files:

```bash
neurofate run \
  --expression data/raw/end_user_smoke/gse20141/GSE20141_series_matrix.txt.gz \
  --metadata results/end_user_smoke/gse20141/sample_metadata.tsv \
  --gene-map results/end_user_smoke/gse20141/gpl570_axis_probe_mapping.tsv \
  --outdir results/end_user_smoke/gse20141/neurofate_public_run_final \
  --sample-id-column geo_accession \
  --endpoint-column label__pd_vs_control \
  --positive-class 1 \
  --negative-class 0 \
  --orientation auto \
  --min-axis-genes 10
```

GSE20141 checksums used in the local smoke test:

- `GSE20141_series_matrix.txt.gz`: `8975344b5a4715032bd07e08a7a94a68b811fddc59b1fbc53dcf204d1005cf4b`
- `GPL570.annot.gz`: `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`

Build the package:

```bash
python -m build --outdir dist_final
python -m twine check dist_final/*
```

Compile the manuscript:

```bash
latexmk -pdf manuscript/bioinformatics/neurofate_bioinformatics_full_methods_paper.tex
```

## Testing

Core checks:

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
  tests/test_pypi_packaging.py \
  tests/test_cli_public_commands.py \
  tests/test_bioinformatics_full_methods_manuscript.py
```

Test coverage includes:

- Public CLI availability.
- GEO series matrix parsing.
- CSV/TSV/GZ detection.
- Expression orientation detection.
- Ensembl and probe mapping.
- Expression/metadata sample joins.
- End-to-end `neurofate run`.
- Endpoint adapter safety.
- Research-use-only report language.
- Bioinformatics manuscript claim-safety checks.

## Packaging and Release

Version: `0.3.0`

Build artifacts:

```bash
python -m build --outdir dist_final
python -m twine check dist_final/*
```

Historical reviewer archive builders remain separate from PyPI artifacts. When
used, they write review ZIPs such as:

- `release_artifacts/neurofate_source_release_<timestamp>.zip`
- `release_artifacts/neurofate_results_review_<timestamp>.zip`

Before release:

1. Confirm version consistency in `pyproject.toml`, `neurofate/__init__.py`, `CITATION.cff`, `codemeta.json`, `CHANGELOG.md`, README, docs, and manuscript.
2. Confirm tests pass.
3. Confirm wheel and source distribution pass `twine check`.
4. Confirm GitHub repository visibility.
5. Create release tag `v0.3.0`.
6. Optionally dry-run TestPyPI.
7. Publish to PyPI.
8. Archive a GitHub release on Zenodo and update citation metadata with DOI.

Do not bundle large public datasets, controlled data, raw matrices, trained real-data models, or generated heavy outputs in the PyPI package.

## Troubleshooting

### Sample IDs Do Not Match

Inspect:

```text
ingest/expression_metadata_join.tsv
```

Common causes include whitespace, punctuation differences, using sample titles instead of accessions, or choosing the wrong sample ID column. Rerun with `--sample-id-column`.

### Ambiguous Endpoint Column

Rerun with explicit endpoint settings:

```bash
--endpoint-column diagnosis --positive-class AD --negative-class Control
```

### Too Few Axis Genes

Check:

```text
ingest/gene_mapping_report.tsv
axis/axis_feature_coverage.tsv
```

Use `--gene-map` for microarray probes or an alias table for Ensembl IDs.

### Unsupported Raw Formats

The public CLI rejects FASTQ/FQ, SRA, CEL/CHP, H5AD/AnnData, and HDF5 containers. Convert outside NeuroFate to compact sample-level or target-gene tables first.

### Missing Gene Map for Microarray

Prepare a table with at least:

```text
probe_id    gene_symbol
```

Then pass:

```bash
--gene-map probe_map.tsv
```

### GEO File Not Parsed

Confirm the file contains:

```text
!series_matrix_table_begin
```

If the file is a SOFT/MINiML/platform annotation rather than a series matrix expression table, prepare the expression table separately.

### Low Coverage Warnings

Low axis-gene coverage does not necessarily mean the run failed. It means interpretation should be cautious and platform coverage should be reported.

## Citation

Use `CITATION.cff` for the software citation. Cite the Bioinformatics manuscript after publication and cite each external dataset according to its source-specific instructions.

Manuscript citation placeholder:

```text
Ghosh N, Sinha K. NeuroFate: format-aware command-line software for endpoint-locked transcriptomic neurodegeneration risk scoring. Bioinformatics. In preparation.
```

Zenodo DOI placeholder: add after archiving the release.

## License

NeuroFate is released under the MIT License. See `LICENSE`.

## Contributing

See:

- `CONTRIBUTING.md`
- `CODE_OF_CONDUCT.md`

Contributions should preserve the research-use-only safety boundary, avoid care-delivery claims, and keep public commands reproducible on compact donor/sample-level inputs.
