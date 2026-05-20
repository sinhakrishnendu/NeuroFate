# Changelog

## 0.3.0 - 2026-05-20

- Added `neurofate ingest` for format-aware inspection, validation, harmonization, endpoint inference, sample-join auditing, and NeuroFate axis-gene standardization from CSV/TSV/TXT/GZ expression and metadata tables.
- Added `neurofate run` as the main public workflow: ingest -> axis scoring -> research-use risk scoring -> report generation.
- Added GEO series matrix support for expression sections marked by `!series_matrix_table_begin`.
- Added `neurofate adapt-endpoint` to create explicit endpoint aliases for validation-script compatibility.
- Added examples for gene-by-sample, sample-by-gene, long-format, Ensembl-ID, and microarray probe-map inputs.
- Added a real-world GSE20141 public CLI smoke test documenting sample matching, GPL570 probe mapping, axis coverage, and research-use outputs.
- Added input/output schema documentation and expanded Bioinformatics full-methods manuscript assets.
- Kept raw FASTQ/SRA/CEL/H5AD processing out of the public ingestion command and preserved research-use-only safeguards.

## 0.2.0 - 2026-05-19

- Reframed NeuroFate as command-line, PyPI-ready research software for endpoint-locked transcriptomic neurodegeneration risk scoring.
- Added stable public CLI commands for axis-score building and research-use risk scoring.
- Strengthened tiny demo outputs with axis scores, risk-score table, and research-use-only report language.
- Moved heavy plotting/MPS dependencies to optional extras and removed mandatory Scanpy/AnnData-style dependencies from package metadata.
- Added Bioinformatics-oriented documentation, validation summaries, tests, and manuscript assets.

## 0.1.0 - 2026-05-17

- Added guarded NeuroFate CLI.
- Added tiny synthetic demo workflow.
- Added release metadata, documentation, CI, report generation, reproducibility manifest, output inventory, and no-overclaiming audit.
- Preserved real-data workflows as manual, guarded commands.
