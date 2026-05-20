# NeuroFate Input and Output Schema

NeuroFate accepts compact donor/sample-level transcriptomic tables and converts
them into standardized inputs for endpoint-locked axis scoring. The public
`neurofate ingest` command does not process raw FASTQ, SRA, CEL, CHP, H5AD, or
large single-cell containers.

## Supported Input Formats

- CSV, TSV, and TXT tables
- `.gz` compressed CSV/TSV/TXT files
- GEO series matrix files containing `!series_matrix_table_begin` and
  `!series_matrix_table_end`
- Gene-by-sample expression matrices
- Sample-by-gene expression matrices
- Long tables with `sample_id`, `gene_symbol`, and `expression_value`
- Probe-by-sample microarray matrices with `--gene-map`
- Ensembl-ID matrices using `metadata/neurofate_axis_gene_aliases.tsv`

## Required Metadata Fields

At minimum, metadata must contain:

- a sample identifier column, such as `sample_id`, `geo_accession`, `donor_id`,
  `subject_id`, or `participant_id`
- an endpoint/group column, such as `diagnosis`, `disease_state`, `condition`,
  `group`, `status`, or `phenotype`

The endpoint must contain an unambiguous positive and negative class.

## Auto-Detected Fields

`neurofate ingest` can infer:

- table delimiter and gzip compression
- expression orientation
- sample ID column
- endpoint column
- positive and negative endpoint classes when values resemble AD/PD/control
  labels
- gene/probe identifier column
- gene identifier type: symbol, Ensembl ID, versioned Ensembl ID, Entrez-like
  numeric ID, or unknown

If inference is ambiguous, the command stops and prints a suggested explicit
argument rather than guessing.

## Expression Orientation Examples

GEO series matrix:

```text
!Series_title  "Example"
!Sample_geo_accession  "GSM1"  "GSM2"
!series_matrix_table_begin
"ID_REF"  "GSM1"  "GSM2"
"1007_s_at"  1.2  1.4
!series_matrix_table_end
```

NeuroFate reads only the tabular expression section and ignores the GEO
metadata preamble for expression scoring. Supply metadata separately with
sample identifiers matching the expression columns, such as GEO accessions.

Gene-by-sample:

```text
gene_symbol  S01  S02  S03
SNCA         0.2  0.3  0.9
GFAP         0.1  0.2  1.0
```

Sample-by-gene:

```text
sample_id  SNCA  GFAP  NEFL
S01        0.2   0.1   1.0
S02        0.3   0.2   0.9
```

## Long-Format Example

```text
sample_id  gene_symbol  expression_value
S01        SNCA         0.2
S01        GFAP         0.1
S02        SNCA         0.3
```

## Ensembl-ID Example

```text
ensembl_gene_id  S01  S02
ENSG00000145335  0.2  0.3
ENSG00000131095  0.1  0.2
```

The alias table maps curated NeuroFate axis genes to human Ensembl gene IDs.
Only mapped NeuroFate axis genes are retained.

## Microarray Probe-Map Example

Expression:

```text
ID_REF      S01  S02
probe_SNCA  0.2  0.3
probe_GFAP  0.1  0.2
```

Probe map:

```text
probe_id    gene_symbol
probe_SNCA  SNCA
probe_GFAP  GFAP
```

Run with `--gene-map path/to/probe_map.tsv`.

## Unsupported Raw Formats

The public ingestion command intentionally rejects raw or container formats:

- FASTQ/FQ
- SRA
- CEL/CHP
- H5AD/AnnData
- HDF5 single-cell matrices

Prepare compact sample-level or target-gene expression tables before using the
public CLI.

## Output File Dictionary

`neurofate ingest` writes:

- `standardized_expression.tsv.gz`: NeuroFate axis-gene expression matrix with
  genes as rows and samples as columns
- `standardized_metadata.tsv`: sample metadata with `sample_id`,
  `label__endpoint`, and research-use marker
- `input_schema_detected.tsv`: inferred table structure and endpoint settings
- `expression_metadata_join.tsv`: sample overlap audit
- `gene_mapping_report.tsv`: retained and discarded feature mapping audit
- `ingest_warnings.tsv`: non-fatal warnings
- `ingest_report.md`: human-readable ingest report
- `run_config.yaml`: reproducibility settings

`neurofate run` additionally writes:

- `axis/axis_scores.tsv`
- `axis/axis_feature_coverage.tsv`
- `axis/label_summary.tsv`
- `risk/neurofate_risk_scores.tsv`
- `risk/risk_score_report.md`
- `neurofate_run_report.md`
- `run_config.yaml`

`neurofate adapt-endpoint` writes:

- `adapted_metadata.tsv`: standardized metadata plus explicit endpoint aliases
- `endpoint_aliases.tsv`: source-to-alias mapping audit
- `endpoint_adapter_report.md`: human-readable compatibility report

## Troubleshooting

- If no endpoint is inferred, pass `--endpoint-column`, `--positive-class`, and
  `--negative-class`.
- If no genes are retained, check whether the expression table uses Ensembl IDs
  or probes and provide `--gene-map` when needed.
- If samples do not join, inspect `expression_metadata_join.tsv` for whitespace,
  punctuation, or GEO accession mismatches.
- If a raw container is rejected, convert it outside the public NeuroFate CLI to
  a compact sample-level expression table first.

NeuroFate is intended for research use only. It is not validated for clinical
diagnosis, patient-level decision-making, or treatment selection.
