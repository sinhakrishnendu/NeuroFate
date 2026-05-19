#!/usr/bin/env python3
"""Audit gene identifiers in GSE174367 normExpr.reg without conversion."""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path


AUDIT_COLUMNS = [
    "rule_name",
    "identifier_type",
    "expression_gene_count",
    "axis_gene_count",
    "matched_axis_genes",
    "matched_gene_symbols",
    "example_expression_ids",
    "safe_to_use",
    "notes",
]


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def axis_genes(axis_rows: list[dict[str, str]]) -> list[str]:
    genes: set[str] = set()
    for row in axis_rows:
        genes.update(parse_genes(row.get("gene_members", "")))
    return sorted(genes)


def r_string_vector(values: list[str]) -> str:
    escaped = [value.replace("\\", "\\\\").replace('"', '\\"') for value in values]
    return "c(" + ",".join(f'"{value}"' for value in escaped) + ")"


def r_audit_code(genes: list[str]) -> str:
    return f"""
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
alias_table <- args[[2]]
output <- args[[3]]
preview_output <- args[[4]]
axis_genes <- {r_string_vector(genes)}
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(preview_output), recursive=TRUE, showWarnings=FALSE)
tmp <- tempfile(fileext='.rda')
con_in <- gzfile(input, 'rb')
con_out <- file(tmp, 'wb')
repeat {{
  chunk <- readBin(con_in, what='raw', n=1048576)
  if (length(chunk) == 0) break
  writeBin(chunk, con_out)
}}
close(con_in); close(con_out)
env <- new.env(parent=emptyenv())
load(tmp, envir=env)
if (!('normExpr.reg' %in% ls(env))) stop('normExpr.reg not found.')
expr <- get('normExpr.reg', envir=env)
ids <- as.character(rownames(expr))
if (is.null(ids)) ids <- character()
classify <- function(values) {{
  examples <- values[!is.na(values) & values != '']
  if (length(examples) == 0) return('unknown')
  if (mean(grepl('^ENSG[0-9]+\\\\.[0-9]+$', examples)) > 0.5) return('ensembl_gene_id_versioned')
  if (mean(grepl('^ENSG[0-9]+$', examples)) > 0.5) return('ensembl_gene_id')
  if (mean(grepl('^[0-9]+$', examples)) > 0.5) return('entrez_numeric')
  if (mean(grepl('^[A-Za-z0-9.-]+$', examples)) > 0.5) return('gene_symbol')
  return('unknown')
}}
identifier_type <- classify(ids)
alias <- if (file.exists(alias_table)) read.delim(alias_table, stringsAsFactors=FALSE, check.names=FALSE) else data.frame(gene_symbol=character(), ensembl_gene_id=character(), alias_type=character())
alias <- alias[!is.na(alias$ensembl_gene_id) & alias$ensembl_gene_id != '' & !is.na(alias$gene_symbol) & alias$gene_symbol != '', , drop=FALSE]
axis_set <- unique(axis_genes)
collapse_punct <- function(x) gsub('[.-]', '', toupper(x))
strip_version <- function(x) sub('\\\\.[0-9]+$', '', x)
rules <- list(
  direct_symbol_match=setNames(toupper(ids), ids),
  uppercase_match=setNames(toupper(ids), ids),
  strip_version_suffix=setNames(toupper(strip_version(ids)), ids),
  replace_hyphen_dot_variants=setNames(collapse_punct(ids), ids)
)
alias_map <- character()
if (nrow(alias) > 0) {{
  alias_map <- toupper(alias$gene_symbol)
  names(alias_map) <- toupper(alias$ensembl_gene_id)
}}
rows <- list()
i <- 0
for (rule_name in names(rules)) {{
  keys <- rules[[rule_name]]
  matched_symbols <- unique(keys[keys %in% axis_set])
  if (rule_name == 'replace_hyphen_dot_variants') {{
    collapsed_axis <- setNames(axis_set, collapse_punct(axis_set))
    matched_keys <- unique(keys[keys %in% names(collapsed_axis)])
    matched_symbols <- unique(collapsed_axis[matched_keys])
  }}
  i <- i + 1
  rows[[i]] <- data.frame(
    rule_name=rule_name,
    identifier_type=identifier_type,
    expression_gene_count=length(ids),
    axis_gene_count=length(axis_set),
    matched_axis_genes=length(matched_symbols),
    matched_gene_symbols=paste(sort(matched_symbols), collapse=';'),
    example_expression_ids=paste(head(ids, 20), collapse=';'),
    safe_to_use=ifelse(length(matched_symbols) >= 10, 'true', 'false'),
    notes='built-in identifier rule',
    stringsAsFactors=FALSE
  )
}}
if (length(alias_map) > 0) {{
  key_sets <- list(
    ensembl_alias_match=toupper(ids),
    ensembl_alias_version_stripped=toupper(strip_version(ids))
  )
  for (rule_name in names(key_sets)) {{
    keys <- key_sets[[rule_name]]
    hit <- keys[keys %in% names(alias_map)]
    matched_symbols <- unique(alias_map[hit])
    matched_symbols <- matched_symbols[matched_symbols %in% axis_set]
    i <- i + 1
    rows[[i]] <- data.frame(
      rule_name=rule_name,
      identifier_type=identifier_type,
      expression_gene_count=length(ids),
      axis_gene_count=length(axis_set),
      matched_axis_genes=length(matched_symbols),
      matched_gene_symbols=paste(sort(matched_symbols), collapse=';'),
      example_expression_ids=paste(head(ids, 20), collapse=';'),
      safe_to_use=ifelse(length(matched_symbols) >= 10, 'true', 'false'),
      notes='explicit alias table rule',
      stringsAsFactors=FALSE
    )
  }}
}}
audit <- do.call(rbind, rows)
write.table(audit, output, sep='\\t', quote=FALSE, row.names=FALSE)
preview <- c(
  'GSE174367 normExpr.reg gene identifier audit',
  paste('identifier_type:', identifier_type),
  paste('row identifier examples:', paste(head(ids, 30), collapse=';')),
  paste('axis genes requested:', length(axis_set)),
  paste('best matched axis genes:', max(audit$matched_axis_genes))
)
writeLines(preview, preview_output)
unlink(tmp)
"""


def run_audit(rda_gz: Path, alias_table: Path, output: Path, preview_output: Path, genes: list[str]) -> bool:
    rscript = shutil.which("Rscript")
    if not rscript:
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerow(
                {
                    "rule_name": "runtime_unavailable",
                    "identifier_type": "unknown",
                    "expression_gene_count": "",
                    "axis_gene_count": str(len(genes)),
                    "matched_axis_genes": "0",
                    "matched_gene_symbols": "",
                    "example_expression_ids": "",
                    "safe_to_use": "false",
                    "notes": "Rscript is required for RDA gene identifier audit.",
                }
            )
        preview_output.parent.mkdir(parents=True, exist_ok=True)
        preview_output.write_text("Rscript unavailable; no expression conversion was attempted.\n", encoding="utf-8")
        return False
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(r_audit_code(genes))
        script = Path(handle.name)
    try:
        subprocess.run([rscript, str(script), str(rda_gz), str(alias_table), str(output), str(preview_output)], check=True)
    finally:
        script.unlink(missing_ok=True)
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit GSE174367 bulk RNA row identifiers against NeuroFate axis genes.")
    parser.add_argument("--rda-gz", type=Path, default=Path("data/raw/external/gse174367_ad_multiomics/GSE174367_bulkRNA_processed.rda.gz"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--alias-table", type=Path, default=Path("metadata/neurofate_axis_gene_aliases.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_gene_identifier_audit.tsv"))
    parser.add_argument("--preview-output", type=Path, default=Path("results/reports/phase31_gse174367_bulk_gene_identifier_preview.txt"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/157_audit_gse174367_bulk_gene_identifiers.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    genes = axis_genes(read_tsv(args.axis_registry))
    run_audit(args.rda_gz, args.alias_table, args.output, args.preview_output, genes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
