#!/usr/bin/env python3
"""Convert GSE174367 bulk RDA into an RDA-targets-mapped axis-gene matrix."""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def axis_gene_map(axis_rows: list[dict[str, str]]) -> dict[str, list[str]]:
    return {row["axis_id"]: parse_genes(row.get("gene_members", "")) for row in axis_rows}


def all_axis_genes(axis_rows: list[dict[str, str]]) -> list[str]:
    genes: set[str] = set()
    for row in axis_rows:
        genes.update(parse_genes(row.get("gene_members", "")))
    return sorted(genes)


def label_for_endpoint(value: str | None) -> str:
    lowered = str(value or "").strip().casefold()
    if lowered in {"ad", "alzheimer disease", "alzheimers disease", "alzheimer's disease", "case"}:
        return "1"
    if lowered in {"control", "normal", "unaffected control", "non-ad", "non ad"}:
        return "0"
    return ""


def build_alias_lookup(alias_rows: list[dict[str, str]]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for row in alias_rows:
        symbol = row.get("gene_symbol", "").strip().upper()
        ensembl = row.get("ensembl_gene_id", "").strip().upper()
        if symbol and ensembl:
            lookup.setdefault(ensembl, symbol)
    return lookup


def read_axis_matrix(path: Path) -> tuple[list[str], set[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        samples = [field for field in (reader.fieldnames or []) if field != "gene_symbol"]
        genes = {row.get("gene_symbol", "").upper() for row in reader if row.get("gene_symbol")}
    return samples, genes


def coverage_rows(axis_rows: list[dict[str, str]], found_genes: set[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for axis_id, genes in axis_gene_map(axis_rows).items():
        found = sorted(set(genes) & found_genes)
        missing = sorted(set(genes) - set(found))
        rows.append(
            {
                "axis_id": axis_id,
                "genes_requested": str(len(genes)),
                "genes_found": str(len(found)),
                "genes_missing": str(len(missing)),
                "found_gene_members": ";".join(found),
                "missing_gene_members": ";".join(missing),
                "status": "ok" if found else "insufficient_coverage",
            }
        )
    return rows


def r_string_vector(values: list[str]) -> str:
    escaped = [value.replace("\\", "\\\\").replace('"', '\\"') for value in values]
    return "c(" + ",".join(f'"{value}"' for value in escaped) + ")"


def r_axis_extraction_code(axis_genes: list[str], min_mapped_genes: int) -> str:
    return f"""
args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
series_metadata <- args[[2]]
alias_table <- args[[3]]
out <- args[[4]]
sample_map_out <- args[[5]]
candidates_out <- args[[6]]
gene_mapping_out <- args[[7]]
axis_genes <- {r_string_vector(axis_genes)}
min_mapped_genes <- {int(min_mapped_genes)}
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(sample_map_out), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(candidates_out), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(gene_mapping_out), recursive=TRUE, showWarnings=FALSE)
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
objects <- ls(env)
if (!('normExpr.reg' %in% objects)) stop('Required expression object normExpr.reg was not found; stop rather than guessing.')
if (!('targets' %in% objects)) stop('Required metadata object targets was not found; stop rather than guessing.')
expr <- get('normExpr.reg', envir=env)
targets <- as.data.frame(get('targets', envir=env), stringsAsFactors=FALSE)
if (is.null(colnames(expr))) stop('Expression object has no sample column names.')
expr_samples <- as.character(colnames(expr))
sample_candidates <- list()
if ('SampleID' %in% colnames(targets)) sample_candidates[['SampleID']] <- as.character(targets[['SampleID']])
if ('Sample.ID' %in% colnames(targets)) sample_candidates[['Sample.ID']] <- as.character(targets[['Sample.ID']])
sample_candidates[['rownames']] <- as.character(rownames(targets))
overlap_rows <- list()
best_name <- ''
best_overlap <- -1
i <- 0
for (candidate_name in names(sample_candidates)) {{
  values <- sample_candidates[[candidate_name]]
  overlap <- length(intersect(expr_samples, values))
  i <- i + 1
  overlap_rows[[i]] <- data.frame(candidate_column=candidate_name, overlap_count=overlap, candidate_sample_count=length(unique(values)), stringsAsFactors=FALSE)
  if (candidate_name == 'SampleID' && overlap > 0) {{
    best_name <- candidate_name
    best_overlap <- overlap
  }} else if (best_name == '' && overlap > best_overlap) {{
    best_name <- candidate_name
    best_overlap <- overlap
  }}
}}
if (best_name == '' || best_overlap == 0) stop('No overlap between expression colnames and RDA targets sample IDs; stop rather than mapping to GEO series metadata.')
sample_ids <- sample_candidates[[best_name]]
detect_diagnosis <- function(tab) {{
  preferred <- c('Diagnosis', 'diagnosis', 'Group', 'condition', 'disease', 'Clinical.Dx', 'Neuropath.Dx')
  for (candidate in preferred) {{
    if (candidate %in% colnames(tab)) {{
      vals <- trimws(as.character(tab[[candidate]]))
      if (any(vals %in% c('AD', 'Control'))) return(candidate)
    }}
  }}
  for (candidate in colnames(tab)) {{
    vals <- trimws(as.character(tab[[candidate]]))
    if (any(vals == 'AD') && any(vals == 'Control')) return(candidate)
  }}
  return('')
}}
diagnosis_column <- detect_diagnosis(targets)
if (diagnosis_column == '') stop('Could not identify an AD/Control diagnosis column in targets.')
target_lookup <- targets[match(expr_samples, sample_ids), , drop=FALSE]
matched <- !is.na(match(expr_samples, sample_ids))
if (sum(matched) == 0) stop('RDA targets mapping produced zero matched samples.')
diagnosis_values <- trimws(as.character(target_lookup[[diagnosis_column]]))
labels <- ifelse(diagnosis_values == 'AD', '1', ifelse(diagnosis_values == 'Control', '0', ''))
sample_map <- data.frame(
  expression_sample_id=expr_samples,
  sample_id=expr_samples,
  target_sample_id=if (best_name == 'rownames') as.character(rownames(target_lookup)) else as.character(target_lookup[[best_name]]),
  diagnosis=diagnosis_values,
  inferred_ad_endpoint=diagnosis_values,
  label__ad_vs_control=labels,
  match_status=ifelse(matched, 'matched', 'unmatched'),
  mapping_source='rda_targets',
  mapping_column=best_name,
  diagnosis_column=diagnosis_column,
  stringsAsFactors=FALSE
)
sample_map <- sample_map[sample_map$match_status == 'matched' & sample_map$label__ad_vs_control %in% c('0', '1'), , drop=FALSE]
if (nrow(sample_map) == 0) stop('No matched AD/Control samples remained after target mapping.')
strip_version <- function(x) sub('\\\\.[0-9]+$', '', x)
collapse_punct <- function(x) gsub('[.-]', '', toupper(x))
expr_ids <- as.character(rownames(expr))
if (is.null(expr_ids)) stop('normExpr.reg has no row identifiers for gene mapping.')
alias <- if (file.exists(alias_table)) read.delim(alias_table, stringsAsFactors=FALSE, check.names=FALSE) else data.frame(gene_symbol=character(), ensembl_gene_id=character(), alias_type=character())
alias <- alias[!is.na(alias$ensembl_gene_id) & alias$ensembl_gene_id != '' & !is.na(alias$gene_symbol) & alias$gene_symbol != '', , drop=FALSE]
alias_map <- character()
if (nrow(alias) > 0) {{
  alias_map <- toupper(alias$gene_symbol)
  names(alias_map) <- toupper(alias$ensembl_gene_id)
}}
alias_map_stripped <- alias_map
names(alias_map_stripped) <- strip_version(names(alias_map))
map_one <- function(identifier) {{
  id <- as.character(identifier)
  upper <- toupper(id)
  stripped <- toupper(strip_version(id))
  collapsed <- collapse_punct(id)
  axis_set <- unique(axis_genes)
  if (id %in% axis_set) return(c(symbol=id, rule='direct_gene_symbol'))
  if (upper %in% axis_set) return(c(symbol=upper, rule='uppercase_gene_symbol'))
  if (upper %in% names(alias_map)) return(c(symbol=alias_map[[upper]], rule='ensembl_alias_match'))
  if (stripped %in% names(alias_map_stripped)) return(c(symbol=alias_map_stripped[[stripped]], rule='ensembl_alias_version_stripped'))
  if (collapsed %in% collapse_punct(axis_set)) {{
    collapsed_axis <- setNames(axis_set, collapse_punct(axis_set))
    return(c(symbol=collapsed_axis[[collapsed]], rule='replace_hyphen_dot_variant'))
  }}
  return(c(symbol='', rule='unmapped'))
}}
mapped <- t(vapply(expr_ids, map_one, FUN.VALUE=c(symbol='', rule='')))
gene_map <- data.frame(
  expression_gene_id=expr_ids,
  gene_symbol=toupper(mapped[, 'symbol']),
  mapping_rule=mapped[, 'rule'],
  is_neurofate_axis_gene=ifelse(mapped[, 'symbol'] != '', 'true', 'false'),
  stringsAsFactors=FALSE
)
gene_map <- gene_map[gene_map$is_neurofate_axis_gene == 'true', , drop=FALSE]
mapped_gene_count <- length(unique(gene_map$gene_symbol))
write.table(gene_map, gene_mapping_out, sep='\\t', quote=FALSE, row.names=FALSE)
if (mapped_gene_count < min_mapped_genes) stop(paste('Fewer than', min_mapped_genes, 'NeuroFate axis genes mapped in normExpr.reg; stop rather than producing undercovered replication matrix.'))
keep_gene <- match(gene_map$expression_gene_id, rownames(expr))
keep_sample <- match(sample_map$expression_sample_id, colnames(expr))
axis_obj <- expr[keep_gene, keep_sample, drop=FALSE]
axis_df <- as.data.frame(axis_obj, check.names=FALSE)
axis_df$gene_symbol <- gene_map$gene_symbol
axis_df <- stats::aggregate(. ~ gene_symbol, data=axis_df, FUN=mean)
write.table(axis_df, out, sep='\\t', quote=FALSE, row.names=FALSE)
write.table(sample_map, sample_map_out, sep='\\t', quote=FALSE, row.names=FALSE)
candidates <- do.call(rbind, overlap_rows)
candidates$best_mapping <- candidates$candidate_column == best_name
candidates$expression_object <- 'normExpr.reg'
candidates$target_object <- 'targets'
candidates$expression_sample_count <- length(expr_samples)
candidates$diagnosis_column <- diagnosis_column
candidates$ad_count <- sum(sample_map$label__ad_vs_control == '1')
candidates$control_count <- sum(sample_map$label__ad_vs_control == '0')
candidates$mapped_axis_gene_count <- mapped_gene_count
if (file.exists(series_metadata)) candidates$series_metadata_role <- 'secondary_annotation_only' else candidates$series_metadata_role <- 'not_available'
write.table(candidates, candidates_out, sep='\\t', quote=FALSE, row.names=FALSE)
unlink(tmp)
"""


def run_r_axis_extraction(
    rda_gz: Path,
    series_metadata: Path,
    alias_table: Path,
    axis_genes: list[str],
    output_matrix: Path,
    sample_map_output: Path,
    candidate_output: Path,
    gene_mapping_output: Path,
    min_mapped_genes: int,
) -> None:
    rscript = shutil.which("Rscript")
    if not rscript:
        raise SystemExit("Rscript is required for GSE174367 RDA conversion. Run the inspector first; no conversion was attempted.")
    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False, encoding="utf-8") as handle:
        handle.write(r_axis_extraction_code(axis_genes, min_mapped_genes))
        script = Path(handle.name)
    try:
        subprocess.run(
            [
                rscript,
                str(script),
                str(rda_gz),
                str(series_metadata),
                str(alias_table),
                str(output_matrix),
                str(sample_map_output),
                str(candidate_output),
                str(gene_mapping_output),
            ],
            check=True,
        )
    finally:
        script.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract only mapped NeuroFate axis genes from GSE174367 bulk RDA using RDA targets metadata.")
    parser.add_argument("--rda-gz", type=Path, default=Path("data/raw/external/gse174367_ad_multiomics/GSE174367_bulkRNA_processed.rda.gz"))
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/phase28_gse174367_ad_multiomics_sample_metadata.tsv"), help="Optional GEO series metadata for secondary annotation only.")
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--alias-table", type=Path, default=Path("metadata/neurofate_axis_gene_aliases.tsv"))
    parser.add_argument("--min-mapped-genes", type=int, default=10)
    parser.add_argument("--output-matrix", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_matrix.tsv"))
    parser.add_argument("--sample-map-output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_sample_map.tsv"))
    parser.add_argument("--coverage-output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_coverage.tsv"))
    parser.add_argument("--gene-mapping-output", type=Path, default=Path("results/tables/phase31_gse174367_bulk_axis_gene_mapping.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/150_convert_gse174367_bulk_rda_to_axis_matrix_phase31.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    axis_rows = read_tsv(args.axis_registry)
    alias_rows = read_tsv(args.alias_table)
    _alias_lookup = build_alias_lookup(alias_rows)
    genes = all_axis_genes(axis_rows)
    candidate_output = args.coverage_output.with_name(args.coverage_output.stem + "_rda_target_mapping_candidates.tsv")
    run_r_axis_extraction(
        args.rda_gz,
        args.metadata,
        args.alias_table,
        genes,
        args.output_matrix,
        args.sample_map_output,
        candidate_output,
        args.gene_mapping_output,
        args.min_mapped_genes,
    )
    samples, found_genes = read_axis_matrix(args.output_matrix)
    sample_map = read_tsv(args.sample_map_output)
    matched = [row for row in sample_map if row.get("match_status") == "matched" and row.get("label__ad_vs_control") in {"0", "1"}]
    if len(matched) == 0:
        raise SystemExit("No expression sample IDs matched RDA targets; stop rather than guessing.")
    if len(found_genes) < args.min_mapped_genes:
        raise SystemExit(f"Only {len(found_genes)} mapped NeuroFate axis genes; required at least {args.min_mapped_genes}.")
    write_tsv(args.coverage_output, coverage_rows(axis_rows, found_genes), ["axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    logging.info("Wrote Phase 31 axis-gene matrix genes=%d samples=%d mapped_samples=%d", len(found_genes), len(samples), len(matched))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
