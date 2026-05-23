# NeuroFate

[![DOI](https://zenodo.org/badge/1241603569.svg)](https://doi.org/10.5281/zenodo.20319379)

**NeuroFate: format-aware command-line software for endpoint-locked transcriptomic neurodegeneration risk scoring.**

NeuroFate is a Python CLI package for reproducible donor/sample-level transcriptomic analysis in neurodegeneration research. It takes compact expression and metadata tables, harmonizes common gene/probe identifiers, locks the endpoint before scoring, computes curated NeuroFate axis scores, creates research-use risk summaries, and writes audit-friendly reports.

Repository: <https://github.com/sinhakrishnendu/NeuroFate.git>

## What NeuroFate Does

NeuroFate turns a user-supplied transcriptomic cohort into a compact, auditable set of research outputs:

1. It inspects expression and metadata files.
2. It standardizes sample IDs, endpoint labels, and gene/probe identifiers.
3. It retains NeuroFate axis genes rather than exporting whole-genome matrices.
4. It computes donor/sample-level axis scores.
5. It writes a research-use risk score and plain-language report.
6. It records every important choice in join audits, mapping reports, warnings, and run configuration files.

The usual public workflow is:

```text
expression table + metadata
  -> neurofate ingest
  -> standardized expression + standardized metadata
  -> neurofate build-axis-scores
  -> axis scores + feature coverage
  -> neurofate score-risk
  -> research-use score report
```

For most users, the single command `neurofate run` performs the complete workflow.

## Research-Use-Only Notice

NeuroFate is intended for research use only. It is not for care delivery, individual-level disease calls, therapeutic guidance, or regulated health use. NeuroFate outputs are designed for cohort-level transcriptomic research, endpoint-locked disease-state modelling, reproducible software demonstrations, and transparent evidence grading.

Every public report includes the same boundary statement:

```text
Research use only. NeuroFate is not validated for care-delivery use,
individual-level disease calls, or therapeutic guidance.
```

## Key Features

- CLI/PyPI-ready package with the console command `neurofate`.
- One-command public workflow with `neurofate run`.
- Separate inspection/standardization workflow with `neurofate ingest`.
- GEO series matrix support, including embedded `!series_matrix_table_begin` expression tables.
- CSV, TSV, TXT, and `.gz` input support.
- Genes-by-samples, samples-by-genes, and long-format expression support.
- Ensembl ID, gene-symbol, and microarray probe mapping support.
- Endpoint locking with explicit positive and negative classes.
- Curated NeuroFate axis scoring at donor/sample level.
- Research-use risk scoring and Markdown reports.
- Endpoint adapter for compatibility between public CLI outputs and validation scripts.
- Leakage-audit and no-overclaiming audit scripts for repository-level checks.
- Bundled no-download tiny demo.
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

### 1. Check The Installation

```bash
neurofate --help
neurofate check-system
neurofate doctor
```

### 2. Run The Bundled Demo

The demo needs no downloads and finishes in seconds:

```bash
neurofate run-demo
```

Demo outputs are written under:

```text
results/demo/
```

Inspect:

```text
results/demo/axis_scores.tsv
results/demo/neurofate_risk_scores.tsv
results/demo/risk_score_report.md
```

### 3. Run Your Own Dataset

Minimal command:

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_run
```

Recommended explicit command:

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_run \
  --sample-id-column sample_id \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control \
  --orientation auto
```

Expected top-level outputs:

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
- `run_config.yaml`

## Choose The Right Workflow

Use this quick decision guide.

| Goal | Command |
| --- | --- |
| Check install and packaged resources | `neurofate check-system` and `neurofate doctor` |
| Run a no-download example | `neurofate run-demo` |
| Inspect and standardize data before scoring | `neurofate ingest` |
| Run the complete public workflow | `neurofate run` |
| Score an already standardized expression table | `neurofate build-axis-scores` |
| Create a research-use score from axis scores | `neurofate score-risk` |
| Create endpoint aliases for validation scripts | `neurofate adapt-endpoint` |

Most users should start with `neurofate run-demo`, then run `neurofate run` on their own expression and metadata files.

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

Reports Python version, platform, executable path, and optional dependency availability.

```bash
neurofate check-system
```

### `neurofate doctor`

Checks packaged resources and, in a repository checkout, core project files.

```bash
neurofate doctor
```

### `neurofate run-demo`

Runs a small synthetic dataset without downloads and writes demo outputs under `results/demo/`.

```bash
neurofate run-demo
```

### `neurofate ingest`

Inspects expression and metadata tables, infers format, validates sample overlap and endpoint labels, maps genes/probes, writes standardized inputs, and reports warnings.

```bash
neurofate ingest \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/neurofate_ingest \
  --sample-id-column sample_id \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control
```

Useful when you want to examine the detected input schema before score calculation.

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
  --outdir results/neurofate_run \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control
```

Use `--gene-map` for microarray probe inputs:

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --gene-map probe_map.tsv \
  --outdir results/neurofate_run \
  --sample-id-column sample_id \
  --endpoint-column group \
  --positive-class PD \
  --negative-class Control
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

The adapter copies binary 0/1 labels only. It does not reinterpret class direction.

### Advanced Repository Commands

`make-report` and historical phase scripts are retained for reproducibility in the full repository checkout. They are not required for normal public use.

Experimental or internal research commands such as `train-baseline`, `train-mps`, `validate-external`, `benchmark`, and phase-specific wrappers should be treated as advanced workflows unless documented otherwise.

## Prepare Your Files

Before running NeuroFate, confirm four things.

1. Your expression table contains sample-level values.
2. Your metadata table contains one row per sample.
3. Sample identifiers can be matched between expression and metadata.
4. You know the endpoint column and the two classes to compare.

Recommended metadata columns:

```text
sample_id    diagnosis
S01          Control
S02          AD
S03          AD
```

For GEO studies, accessions such as `GSM503950` are usually the safest sample IDs.

For microarray data, prepare a probe map if the expression rows are probes rather than gene symbols:

```text
probe_id      gene_symbol
1007_s_at     DDR1
207827_x_at   SNCA
```

For Ensembl matrices, NeuroFate can use the bundled curated alias table:

```text
metadata/neurofate_axis_gene_aliases.tsv
```

## Input Formats

NeuroFate public ingestion accepts compact text tables. It does not process raw FASTQ/FQ, SRA, CEL/CHP, H5AD/AnnData, or HDF5 single-cell containers.

Compressed `.gz` files are supported for CSV, TSV, TXT, and GEO series matrix inputs.

### Genes-by-Samples Matrix

```text
gene_symbol    S01    S02    S03
SNCA           0.2    0.4    0.8
GFAP           0.1    0.3    1.1
NEFL           1.2    1.0    0.7
```

Command:

```bash
neurofate run \
  --expression genes_by_samples.tsv \
  --metadata metadata.tsv \
  --outdir results/neurofate_run \
  --orientation genes_rows
```

### Samples-by-Genes Matrix

```text
sample_id    SNCA    GFAP    NEFL
S01          0.2     0.1     1.2
S02          0.4     0.3     1.0
S03          0.8     1.1     0.7
```

Command:

```bash
neurofate run \
  --expression samples_by_genes.tsv \
  --metadata metadata.tsv \
  --outdir results/neurofate_run \
  --orientation samples_rows
```

### Long Format

```text
sample_id    gene_symbol    expression_value
S01          SNCA           0.2
S01          GFAP           0.1
S02          SNCA           0.4
```

Command:

```bash
neurofate run \
  --expression long_expression.tsv \
  --metadata metadata.tsv \
  --outdir results/neurofate_run \
  --orientation long
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

Command:

```bash
neurofate run \
  --expression GSE_series_matrix.txt.gz \
  --metadata sample_metadata.tsv \
  --gene-map platform_axis_probe_mapping.tsv \
  --outdir results/neurofate_run \
  --sample-id-column geo_accession \
  --endpoint-column label__pd_vs_control \
  --positive-class 1 \
  --negative-class 0
```

### Ensembl-ID Matrix

```text
ensembl_gene_id    S01    S02
ENSG00000145335    0.2    0.4
ENSG00000131095    0.1    0.3
```

NeuroFate maps curated axis genes using:

```text
metadata/neurofate_axis_gene_aliases.tsv
```

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
- `ingest_warnings.tsv`: non-fatal warnings. This file is written even when empty.
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

## How To Read The Outputs

Start with these files.

### `neurofate_run_report.md`

This is the human-readable run summary. It tells you whether ingestion, axis scoring, and risk scoring completed.

### `ingest/expression_metadata_join.tsv`

Use this first if samples are missing. Important columns:

- `expression_sample_count`
- `metadata_sample_count`
- `matched_sample_count`
- `unmatched_expression_samples`
- `unmatched_metadata_samples`

### `ingest/gene_mapping_report.tsv`

Use this to see which input genes/probes were retained for NeuroFate axes.

### `axis/axis_feature_coverage.tsv`

Use this to decide how cautious interpretation should be for each axis. Low coverage means the axis score is less complete.

### `axis/axis_scores.tsv`

This is the core sample-level output. Each row is a sample and each axis column is a NeuroFate score.

### `risk/neurofate_risk_scores.tsv`

This is a research-use score summary derived from axis scores. It should be interpreted only as an exploratory cohort-level ranking or stratification output.

### `run_config.yaml`

This records the command settings needed to repeat the run.

## Practical Recipes

### Recipe 1: Let NeuroFate Infer Common Settings

```bash
neurofate run \
  --expression expression.tsv \
  --metadata metadata.tsv \
  --outdir results/my_run
```

Use this for simple files with obvious `sample_id` and `diagnosis` or `group` columns.

### Recipe 2: Lock The Endpoint Explicitly

```bash
neurofate run \
  --expression expression.tsv \
  --metadata metadata.tsv \
  --outdir results/my_run \
  --sample-id-column sample_id \
  --endpoint-column disease_state \
  --positive-class PD \
  --negative-class Control
```

Use this when metadata contains several candidate label columns.

### Recipe 3: Score A GEO Series Matrix With A Platform Map

```bash
neurofate run \
  --expression GSE20141_series_matrix.txt.gz \
  --metadata sample_metadata.tsv \
  --gene-map gpl570_axis_probe_mapping.tsv \
  --outdir results/gse20141_neurofate_run \
  --sample-id-column geo_accession \
  --endpoint-column label__pd_vs_control \
  --positive-class 1 \
  --negative-class 0 \
  --orientation auto
```

Use this when expression columns are GEO sample accessions and rows are platform probes.

### Recipe 4: Inspect First, Score Later

```bash
neurofate ingest \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --outdir results/ingest_check \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control
```

After checking `results/ingest_check/ingest_report.md`, run:

```bash
neurofate build-axis-scores \
  --expression results/ingest_check/standardized_expression.tsv.gz \
  --metadata results/ingest_check/standardized_metadata.tsv \
  --outdir results/axis_scores \
  --sample-id-column sample_id \
  --endpoint-column label__endpoint \
  --positive-class 1 \
  --negative-class 0
```

Then:

```bash
neurofate score-risk \
  --axis-scores results/axis_scores/axis_scores.tsv \
  --outdir results/axis_scores
```

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

Axes are research summaries of available expression features. They should not be read as proof of disease biology by themselves.

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
latexmk -pdf -cd manuscript/bioinformatics/neurofate_bioinformatics_full_methods_paper.tex
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

Release-manual checks:

```bash
python -m pytest \
  tests/test_release_readme_manual.py \
  tests/test_release_research_use_only.py \
  tests/test_release_packaging_metadata.py
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

`dist/` is reserved for PyPI artifacts. Review ZIPs and manuscript/reviewer packages should use `release_artifacts/` or another explicit review directory.

Build artifacts:

```bash
python -m build --outdir dist_final
python -m twine check dist_final/*
```

Historical reviewer archive builders remain separate from PyPI artifacts. When used, they write review ZIPs such as:

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

## Safety And Memory Design

NeuroFate public commands operate on compact donor/sample-level or axis-gene/probe tables. The public ingestion workflow does not process raw FASTQ/SRA, CEL/CHP, H5AD/AnnData, UMAP, clustering, or dense genome-wide converted matrices.

Large study-specific scripts are kept outside the recommended public workflow. They remain in the repository to document analyses, but the PyPI-style interface is centered on compact user-supplied tables.

## Current Validation Status

The current release is validated as research software through public CLI tests, ingestion tests, a bundled tiny demo, a real-world GSE20141 GEO smoke test, package build checks, and no-overclaiming audits. Biological cohort results are demonstration evidence and should not be interpreted as care-delivery validation.

Reviewer report generators remain lightweight and can be run from existing outputs, for example:

```bash
python scripts/51_generate_end_user_report.py --tables-dir results/tables --reports-dir results/reports
```

## Troubleshooting

### Sample IDs Do Not Match

Inspect:

```text
ingest/expression_metadata_join.tsv
```

Common causes:

- Metadata uses titles while expression uses accessions.
- One table has whitespace around IDs.
- Expression sample IDs include punctuation or prefixes not present in metadata.
- The wrong sample ID column was selected.

Fix:

```bash
neurofate run \
  --expression expression.tsv \
  --metadata metadata.tsv \
  --sample-id-column geo_accession \
  --outdir results/neurofate_run
```

### Endpoint Column Is Ambiguous

Rerun with explicit endpoint settings:

```bash
neurofate run \
  --expression expression.tsv.gz \
  --metadata metadata.tsv \
  --endpoint-column diagnosis \
  --positive-class AD \
  --negative-class Control \
  --outdir results/neurofate_run
```

### Too Few Axis Genes

Check:

```text
ingest/gene_mapping_report.tsv
axis/axis_feature_coverage.tsv
```

Common fixes:

- Use `--gene-map` for microarray probes.
- Confirm whether row IDs are gene symbols or Ensembl IDs.
- Lower `--min-axis-genes` only for exploration and report low coverage clearly.

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

### Multiple Probes Map To One Gene

This is common for microarray platforms. NeuroFate reports multi-probe mapping and aggregates retained probes conservatively for axis scoring.

### Output Directory Already Exists

Use a fresh `--outdir` for a clean run:

```bash
neurofate run \
  --expression expression.tsv \
  --metadata metadata.tsv \
  --outdir results/neurofate_run_001
```

## Citation

Use `CITATION.cff` for the software citation. Cite the Bioinformatics manuscript after publication and cite each external dataset according to its source-specific instructions.

Manuscript citation placeholder:

```text
Ghosh N, Sinha K. NeuroFate: command-line research software for endpoint-locked transcriptomic neurodegeneration risk scoring. Bioinformatics. In preparation.
```

Zenodo DOI placeholder: add after archiving the release.

## License

NeuroFate is released under the MIT License. See `LICENSE`.

## Contributing

See:

- `CONTRIBUTING.md`
- `CODE_OF_CONDUCT.md`

Contributions should preserve the research-use-only safety boundary, avoid care-delivery claims, and keep public commands reproducible on compact donor/sample-level inputs.
