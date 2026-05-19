#!/usr/bin/env python3
"""Audit GSE174367 bulk RNA sample-ID mapping candidates without conversion."""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path


AUDIT_COLUMNS = [
    "comparison",
    "expression_object",
    "target_object",
    "candidate_column",
    "expression_sample_count",
    "candidate_sample_count",
    "overlap_count",
    "overlap_rate",
    "diagnosis_column",
    "ad_count",
    "control_count",
    "best_mapping",
    "note",
]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def write_missing_runtime_outputs(output: Path, preview_output: Path, reason: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "comparison": "runtime_unavailable",
                "expression_object": "",
                "target_object": "",
                "candidate_column": "",
                "expression_sample_count": "",
                "candidate_sample_count": "",
                "overlap_count": "0",
                "overlap_rate": "0",
                "diagnosis_column": "",
                "ad_count": "",
                "control_count": "",
                "best_mapping": "false",
                "note": reason,
            }
        )
    preview_output.parent.mkdir(parents=True, exist_ok=True)
    preview_output.write_text(
        "GSE174367 sample mapping audit could not inspect the RDA because Rscript is unavailable.\n"
        f"Reason: {reason}\n"
        "No expression conversion or analysis was attempted.\n",
        encoding="utf-8",
    )


def r_code() -> str:
    return r"""
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
series_metadata <- args[[2]]
output <- args[[3]]
preview_output <- args[[4]]
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(preview_output), recursive=TRUE, showWarnings=FALSE)
tmp <- tempfile(fileext='.rda')
con_in <- gzfile(input, 'rb')
con_out <- file(tmp, 'wb')
repeat {
  chunk <- readBin(con_in, what='raw', n=1048576)
  if (length(chunk) == 0) break
  writeBin(chunk, con_out)
}
close(con_in); close(con_out)
env <- new.env(parent=emptyenv())
load(tmp, envir=env)
objects <- ls(env)
expr_name <- if ('normExpr.reg' %in% objects) 'normExpr.reg' else NA_character_
if (is.na(expr_name)) {
  for (name in objects) {
    obj <- get(name, envir=env)
    if ((is.matrix(obj) || is.data.frame(obj)) && !is.null(colnames(obj)) && length(colnames(obj)) > 10) {
      expr_name <- name
      break
    }
  }
}
target_name <- if ('targets' %in% objects) 'targets' else NA_character_
if (is.na(expr_name) || is.na(target_name)) stop('Required expression object or targets table was not found.')
expr <- get(expr_name, envir=env)
targets <- as.data.frame(get(target_name, envir=env), stringsAsFactors=FALSE)
expr_samples <- as.character(colnames(expr))
candidate_sets <- list()
for (candidate in c('SampleID', 'Sample.ID')) {
  if (candidate %in% colnames(targets)) candidate_sets[[paste0('targets$', candidate)]] <- as.character(targets[[candidate]])
}
candidate_sets[['targets_rownames']] <- as.character(rownames(targets))
if (file.exists(series_metadata)) {
  series <- tryCatch(read.delim(series_metadata, stringsAsFactors=FALSE, check.names=FALSE), error=function(e) NULL)
  if (!is.null(series)) {
    for (candidate in c('sample_id', 'sample_title', 'title', 'geo_accession', 'sample_name', 'source_name')) {
      if (candidate %in% colnames(series)) candidate_sets[[paste0('series$', candidate)]] <- as.character(series[[candidate]])
    }
  }
}
detect_diagnosis <- function(tab) {
  preferred <- c('Diagnosis', 'diagnosis', 'Group', 'condition', 'disease', 'Clinical.Dx', 'Neuropath.Dx')
  for (candidate in preferred) {
    if (candidate %in% colnames(tab)) {
      vals <- trimws(as.character(tab[[candidate]]))
      if (any(vals %in% c('AD', 'Control'))) return(candidate)
    }
  }
  for (candidate in colnames(tab)) {
    vals <- trimws(as.character(tab[[candidate]]))
    if (any(vals == 'AD') && any(vals == 'Control')) return(candidate)
  }
  return('')
}
diagnosis_column <- detect_diagnosis(targets)
diagnosis_values <- if (diagnosis_column != '') trimws(as.character(targets[[diagnosis_column]])) else character()
ad_count <- sum(diagnosis_values == 'AD')
control_count <- sum(diagnosis_values == 'Control')
rows <- list()
best_idx <- 0
best_overlap <- -1
i <- 0
for (candidate_name in names(candidate_sets)) {
  values <- unique(candidate_sets[[candidate_name]][!is.na(candidate_sets[[candidate_name]]) & candidate_sets[[candidate_name]] != ''])
  overlap <- length(intersect(expr_samples, values))
  i <- i + 1
  if (overlap > best_overlap) {
    best_overlap <- overlap
    best_idx <- i
  }
  rows[[i]] <- data.frame(
    comparison='expression_colnames_vs_candidate',
    expression_object=expr_name,
    target_object=target_name,
    candidate_column=candidate_name,
    expression_sample_count=length(expr_samples),
    candidate_sample_count=length(values),
    overlap_count=overlap,
    overlap_rate=ifelse(length(expr_samples) > 0, overlap / length(expr_samples), 0),
    diagnosis_column=diagnosis_column,
    ad_count=ad_count,
    control_count=control_count,
    best_mapping='false',
    note='sample-id overlap audit',
    stringsAsFactors=FALSE
  )
}
audit <- do.call(rbind, rows)
if (nrow(audit) > 0) audit$best_mapping[best_idx] <- 'true'
write.table(audit, output, sep='\t', quote=FALSE, row.names=FALSE)
preview <- c(
  'GSE174367 bulk sample mapping audit',
  paste('expression_object:', expr_name),
  paste('target_object:', target_name),
  paste('expression sample count:', length(expr_samples)),
  paste('expression sample preview:', paste(head(expr_samples, 20), collapse=';')),
  paste('targets columns:', paste(colnames(targets), collapse=';')),
  paste('diagnosis column:', diagnosis_column),
  paste('AD count:', ad_count),
  paste('Control count:', control_count)
)
writeLines(preview, preview_output)
unlink(tmp)
"""


def run_audit(rda_gz: Path, series_metadata: Path, output: Path, preview_output: Path) -> bool:
    rscript = shutil.which("Rscript")
    if not rscript:
        return False
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(r_code())
        script = Path(handle.name)
    try:
        subprocess.run([rscript, str(script), str(rda_gz), str(series_metadata), str(output), str(preview_output)], check=True)
    finally:
        script.unlink(missing_ok=True)
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit GSE174367 RDA-internal sample ID mapping.")
    parser.add_argument("--rda-gz", type=Path, default=Path("data/raw/external/gse174367_ad_multiomics/GSE174367_bulkRNA_processed.rda.gz"))
    parser.add_argument("--series-metadata", type=Path, default=Path("results/tables/phase28_gse174367_ad_multiomics_sample_metadata.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase30_gse174367_bulk_sample_mapping_audit.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("results/reports/phase30_gse174367_bulk_sample_mapping_preview.txt"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/155_audit_gse174367_bulk_sample_mapping.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    if not args.rda_gz.exists():
        raise SystemExit(f"Missing RDA input: {args.rda_gz}")
    if run_audit(args.rda_gz, args.series_metadata, args.output, args.preview_output):
        return 0
    write_missing_runtime_outputs(args.output, args.preview_output, "Rscript is required for RDA-internal sample mapping audit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
