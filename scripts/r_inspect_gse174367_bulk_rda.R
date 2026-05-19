#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

input <- arg_value("--input")
output <- arg_value("--output")
preview_output <- arg_value("--preview-output")

if (input == "" || output == "" || preview_output == "") {
  stop("Usage: Rscript r_inspect_gse174367_bulk_rda.R --input FILE --output TSV --preview-output TXT")
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(preview_output), recursive = TRUE, showWarnings = FALSE)

input_to_load <- input
temp_file <- NULL
if (grepl("\\.gz$", input, ignore.case = TRUE)) {
  temp_file <- tempfile(fileext = ".rda")
  in_con <- gzfile(input, "rb")
  out_con <- file(temp_file, "wb")
  repeat {
    chunk <- readBin(in_con, what = "raw", n = 1024 * 1024)
    if (length(chunk) == 0) {
      break
    }
    writeBin(chunk, out_con)
  }
  close(in_con)
  close(out_con)
  input_to_load <- temp_file
}

env <- new.env(parent = emptyenv())
loaded <- load(input_to_load, envir = env)

preview <- c("GSE174367 bulk RDA structure preview", paste("Objects:", paste(loaded, collapse = ", ")), "")
summary_rows <- list()

collapse_preview <- function(values) {
  if (is.null(values)) {
    return("")
  }
  paste(utils::head(as.character(values), 8), collapse = ";")
}

for (object_name in loaded) {
  object <- get(object_name, envir = env)
  cls <- paste(class(object), collapse = ";")
  dims <- dim(object)
  nrow_value <- ifelse(is.null(dims), length(object), dims[[1]])
  ncol_value <- ifelse(is.null(dims) || length(dims) < 2, "", dims[[2]])
  row_preview <- collapse_preview(rownames(object))
  col_preview <- collapse_preview(colnames(object))
  numeric_like <- is.matrix(object) || is.data.frame(object)
  expression_like <- FALSE
  metadata_like <- FALSE
  if (numeric_like && !is.null(dims) && length(dims) >= 2) {
    expression_like <- nrow_value > 100 && ncol_value > 5
    metadata_like <- nrow_value > 5 && ncol_value > 2 && ncol_value < 200
  }
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    object_name = object_name,
    object_class = cls,
    nrow = as.character(nrow_value),
    ncol = as.character(ncol_value),
    row_names_preview = row_preview,
    col_names_preview = col_preview,
    expression_matrix_like = as.character(expression_like),
    metadata_table_like = as.character(metadata_like),
    stringsAsFactors = FALSE
  )
  preview <- c(preview, paste("##", object_name), paste(capture.output(str(object, max.level = 1)), collapse = "\n"), "")
}

if (length(summary_rows) == 0) {
  summary_df <- data.frame(
    object_name = character(),
    object_class = character(),
    nrow = character(),
    ncol = character(),
    row_names_preview = character(),
    col_names_preview = character(),
    expression_matrix_like = character(),
    metadata_table_like = character(),
    stringsAsFactors = FALSE
  )
} else {
  summary_df <- do.call(rbind, summary_rows)
}
utils::write.table(summary_df, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(preview, con = preview_output)

if (!is.null(temp_file) && file.exists(temp_file)) {
  unlink(temp_file)
}
