"""Format-aware data ingestion for the public NeuroFate CLI.

The ingestion layer standardizes compact transcriptomic tables into the
donor/sample-level TSV files consumed by ``neurofate build-axis-scores``.  It is
deliberately conservative: raw FASTQ/SRA/CEL/H5AD containers are reported as
unsupported, endpoint ambiguity stops the run, and only NeuroFate axis genes are
written to the standardized expression output.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections import Counter
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

from neurofate.axis import RESEARCH_USE_NOTICE, build_axis_score_tables, load_axis_registry, score_research_risk


SUPPORTED_TEXT_SUFFIXES = {".csv", ".tsv", ".txt", ".gz"}
UNSUPPORTED_RAW_SUFFIXES = {".fastq", ".fq", ".sra", ".cel", ".chp", ".h5ad", ".h5"}
ENDPOINT_CANDIDATES = (
    "diagnosis",
    "disease_state",
    "disease state",
    "disease",
    "condition",
    "group",
    "status",
    "phenotype",
    "case_control",
    "class",
    "label",
)
SAMPLE_ID_CANDIDATES = (
    "sample_id",
    "sampleid",
    "sample id",
    "geo_accession",
    "gsm",
    "donor_id",
    "donorid",
    "subject_id",
    "participant_id",
    "id",
)
GENE_ID_CANDIDATES = (
    "gene_symbol",
    "gene symbol",
    "symbol",
    "gene",
    "gene_id",
    "ensembl_gene_id",
    "ensembl",
    "probe_id",
    "id_ref",
    "id",
    "feature",
    "feature_id",
)
VALUE_CANDIDATES = ("expression_value", "expression", "value", "count", "counts", "normalized")
POSITIVE_TOKENS = {
    "1",
    "true",
    "case",
    "disease",
    "ad",
    "alzheimer",
    "alzheimers",
    "alzheimer's disease",
    "pd",
    "parkinson",
    "parkinson disease",
    "parkinson's disease",
    "parkinsons disease",
    "dementia",
}
NEGATIVE_TOKENS = {
    "0",
    "false",
    "control",
    "unaffected control",
    "reference",
    "normal",
    "healthy",
    "no dementia",
    "non-disease",
}


@dataclass
class IngestConfig:
    expression: Path
    metadata: Path
    outdir: Path
    axis_registry: Path
    sample_id_column: str = "auto"
    endpoint_column: str = "auto"
    positive_class: str = "auto"
    negative_class: str = "auto"
    gene_id_column: str = "auto"
    orientation: str = "auto"
    gene_map: Path | None = None
    alias_table: Path | None = None
    assist: bool = False
    min_axis_genes: int = 3


@dataclass
class IngestResult:
    standardized_expression: Path
    standardized_metadata: Path
    input_schema_detected: Path
    expression_metadata_join: Path
    gene_mapping_report: Path
    ingest_warnings: Path
    ingest_report: Path
    run_config: Path


def open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def normalize_name(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().casefold())


def normalize_id(value: object) -> str:
    return re.sub(r"[\s._-]+", "", str(value).strip().casefold())


def normalize_class(value: object) -> str:
    return re.sub(r"\s+", " ", str(value).strip().casefold())


def detect_compression(path: Path) -> str:
    return "gzip" if str(path).endswith(".gz") else "none"


def _reject_unsupported_raw(path: Path) -> None:
    suffixes = [suffix.casefold() for suffix in path.suffixes]
    if any(suffix in UNSUPPORTED_RAW_SUFFIXES for suffix in suffixes):
        raise ValueError(
            f"Unsupported raw/container input {path}. NeuroFate ingest expects compact CSV/TSV/TXT tables; "
            "raw FASTQ/SRA/CEL/H5AD/H5 files are intentionally not processed by this command."
        )


def detect_delimiter(path: Path) -> str:
    _reject_unsupported_raw(path)
    if is_geo_series_matrix(path):
        return "\t"
    sample = ""
    with open_text(path, "rt") as handle:
        for _ in range(10):
            line = handle.readline()
            if not line:
                break
            if line.strip():
                sample += line
    if not sample:
        raise ValueError(f"Could not read a table header from {path}")
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        return dialect.delimiter
    except csv.Error:
        first_line = sample.splitlines()[0]
        if first_line.count("\t") >= first_line.count(","):
            return "\t"
        return ","


def is_geo_series_matrix(path: Path, max_lines: int = 5000) -> bool:
    _reject_unsupported_raw(path)
    with open_text(path, "rt") as handle:
        for index, line in enumerate(handle):
            if index >= max_lines:
                break
            if line.strip() == "!series_matrix_table_begin":
                return True
    return False


def _read_geo_series_matrix_table(path: Path, nrows: int | None = None) -> pd.DataFrame:
    _reject_unsupported_raw(path)
    table_lines: list[str] = []
    in_table = False
    with open_text(path, "rt") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            marker = stripped.strip()
            if marker == "!series_matrix_table_begin":
                in_table = True
                continue
            if marker == "!series_matrix_table_end":
                break
            if in_table:
                table_lines.append(stripped)
                if nrows is not None and len(table_lines) > nrows:
                    break
    if not table_lines:
        raise ValueError(
            f"Could not locate a GEO series matrix expression table in {path}; "
            "expected !series_matrix_table_begin and !series_matrix_table_end markers."
        )
    return pd.read_csv(StringIO("\n".join(table_lines)), sep="\t", nrows=nrows, dtype=str)


def _read_table(path: Path, nrows: int | None = None) -> pd.DataFrame:
    if is_geo_series_matrix(path):
        return _read_geo_series_matrix_table(path, nrows=nrows)
    delimiter = detect_delimiter(path)
    return pd.read_csv(path, sep=delimiter, nrows=nrows, dtype=str, compression="infer")


def inspect_table_shape(path: Path, max_rows: int = 50) -> dict[str, object]:
    preview = _read_table(path, nrows=max_rows)
    return {
        "path": str(path),
        "compression": detect_compression(path),
        "delimiter": "\\t" if detect_delimiter(path) == "\t" else detect_delimiter(path),
        "preview_rows": len(preview),
        "columns": len(preview.columns),
        "column_names": ";".join(map(str, preview.columns[:20])),
    }


def _find_column(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    normalized = {normalize_name(col): col for col in columns}
    for candidate in candidates:
        hit = normalized.get(normalize_name(candidate))
        if hit is not None:
            return hit
    return None


def infer_gene_id_column(columns: Iterable[str]) -> str | None:
    return _find_column(columns, GENE_ID_CANDIDATES)


def infer_sample_id_column(columns: Iterable[str]) -> str | None:
    return _find_column(columns, SAMPLE_ID_CANDIDATES)


def infer_expression_orientation(path: Path) -> str:
    preview = _read_table(path, nrows=20)
    columns = list(preview.columns)
    sample_col = infer_sample_id_column(columns)
    gene_col = infer_gene_id_column(columns)
    value_col = _find_column(columns, VALUE_CANDIDATES)
    if sample_col and gene_col and value_col:
        return "long"
    if gene_col and len(columns) > 2:
        return "genes_rows"
    if sample_col and len(columns) > 2:
        return "samples_rows"
    first_col = str(columns[0]) if columns else ""
    if normalize_name(first_col) in {normalize_name(name) for name in GENE_ID_CANDIDATES}:
        return "genes_rows"
    if normalize_name(first_col) in {normalize_name(name) for name in SAMPLE_ID_CANDIDATES}:
        return "samples_rows"
    return "genes_rows"


def infer_sample_columns(expression: pd.DataFrame, metadata: pd.DataFrame) -> tuple[list[str], str | None]:
    metadata_sample_col = infer_sample_id_column(metadata.columns)
    metadata_ids = set(metadata[metadata_sample_col].map(normalize_id)) if metadata_sample_col else set()
    candidates = []
    for col in expression.columns:
        if normalize_name(col) in {normalize_name(name) for name in GENE_ID_CANDIDATES + SAMPLE_ID_CANDIDATES}:
            continue
        overlap = normalize_id(col) in metadata_ids
        candidates.append((col, overlap))
    matched = [col for col, overlap in candidates if overlap]
    return matched or [col for col, _ in candidates], metadata_sample_col


def infer_metadata_columns(metadata: pd.DataFrame) -> dict[str, str | None]:
    return {
        "sample_id_column": infer_sample_id_column(metadata.columns),
        "endpoint_column": infer_endpoint_column(metadata),
    }


def infer_endpoint_column(metadata: pd.DataFrame) -> str | None:
    for col in metadata.columns:
        if normalize_name(col) in {normalize_name(name) for name in ENDPOINT_CANDIDATES}:
            values = [value for value in metadata[col].dropna().astype(str).map(str.strip) if value]
            if len(set(values)) >= 2:
                return col
    # Fallback: choose a low-cardinality text column containing case/control-like values.
    for col in metadata.columns:
        values = [normalize_class(value) for value in metadata[col].dropna().astype(str) if str(value).strip()]
        unique = set(values)
        if 2 <= len(unique) <= 6 and (
            any(value in POSITIVE_TOKENS for value in unique)
            or any(value in NEGATIVE_TOKENS for value in unique)
        ):
            return col
    return None


def infer_positive_negative_classes(metadata: pd.DataFrame, endpoint_column: str) -> tuple[str, str]:
    values = [str(value).strip() for value in metadata[endpoint_column].dropna() if str(value).strip()]
    counts = Counter(values)
    unique = list(counts)
    positive = [value for value in unique if normalize_class(value) in POSITIVE_TOKENS]
    negative = [value for value in unique if normalize_class(value) in NEGATIVE_TOKENS]
    if positive and negative:
        return positive[0], negative[0]
    if len(unique) == 2:
        lowered = {normalize_class(value): value for value in unique}
        for token in POSITIVE_TOKENS:
            if token in lowered:
                pos = lowered[token]
                neg = next(value for value in unique if value != pos)
                return pos, neg
        for token in NEGATIVE_TOKENS:
            if token in lowered:
                neg = lowered[token]
                pos = next(value for value in unique if value != neg)
                return pos, neg
    raise ValueError(
        f"Could not infer positive/negative classes from endpoint {endpoint_column!r}. "
        f"Observed values: {dict(counts)}. Pass --positive-class and --negative-class explicitly."
    )


def detect_gene_identifier_type(gene_ids: Iterable[object]) -> str:
    values = [str(value).strip() for value in gene_ids if str(value).strip()]
    if not values:
        return "unknown"
    upper = [value.upper() for value in values[:200]]
    if sum(value.startswith("ENSG") for value in upper) / len(upper) > 0.8:
        if any("." in value for value in upper):
            return "ensembl_gene_id_versioned"
        return "ensembl_gene_id"
    if sum(value.isdigit() for value in values[:200]) / len(values[:200]) > 0.8:
        return "entrez_numeric"
    if sum(bool(re.match(r"^[A-Za-z0-9_.-]+$", value)) for value in values[:200]) / len(values[:200]) > 0.8:
        return "gene_symbol"
    return "unknown"


def _axis_gene_set(axis_registry: Path) -> set[str]:
    axes = load_axis_registry(axis_registry)
    return {gene.upper() for genes in axes.values() for gene in genes}


def _load_alias_map(alias_table: Path | None) -> dict[str, str]:
    alias_map: dict[str, str] = {}
    if alias_table is None or not alias_table.exists():
        return alias_map
    rows = _read_table(alias_table)
    for _, row in rows.iterrows():
        symbol = str(row.get("gene_symbol", "")).strip().upper()
        if not symbol:
            continue
        for col in rows.columns:
            if col == "gene_symbol":
                continue
            value = str(row.get(col, "")).strip()
            if value and value.lower() != "nan":
                alias_map[value.upper()] = symbol
                alias_map[value.split(".")[0].upper()] = symbol
        alias_map[symbol] = symbol
    return alias_map


def _load_gene_map(gene_map: Path | None) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if gene_map is None:
        return mapping
    table = _read_table(gene_map)
    feature_col = _find_column(table.columns, ("probe_id", "id_ref", "id", "feature", "feature_id"))
    gene_col = _find_column(table.columns, ("gene_symbol", "gene symbol", "symbol", "gene", "gene_assignment"))
    if not feature_col or not gene_col:
        raise ValueError(
            f"Gene/probe map {gene_map} must contain a probe/feature column and a gene symbol column."
        )
    for _, row in table.iterrows():
        feature = str(row.get(feature_col, "")).strip().upper()
        raw_symbols = str(row.get(gene_col, "")).replace("///", ";").replace("//", ";").replace(",", ";")
        symbols = [symbol.strip().upper() for symbol in raw_symbols.split(";") if symbol.strip()]
        if feature and symbols:
            mapping[feature] = symbols[0]
    return mapping


def _map_feature_to_gene(
    feature: object,
    axis_genes: set[str],
    alias_map: dict[str, str],
    gene_map: dict[str, str],
) -> str | None:
    value = str(feature).strip().upper()
    if not value:
        return None
    candidates = [
        value,
        value.split(".")[0],
        value.replace("-", "."),
        value.replace(".", "-"),
        gene_map.get(value, ""),
        alias_map.get(value, ""),
        alias_map.get(value.split(".")[0], ""),
    ]
    for candidate in candidates:
        candidate = str(candidate).strip().upper()
        if candidate in axis_genes:
            return candidate
    return None


def _coerce_numeric(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(pd.to_numeric, errors="coerce")


def validate_expression_metadata_join(
    expression_samples: Iterable[str],
    metadata_samples: Iterable[str],
) -> dict[str, object]:
    expr = list(expression_samples)
    meta = list(metadata_samples)
    expr_norm = {normalize_id(value): value for value in expr}
    meta_norm = {normalize_id(value): value for value in meta}
    matched = sorted(set(expr_norm) & set(meta_norm))
    return {
        "expression_sample_count": len(expr),
        "metadata_sample_count": len(meta),
        "matched_sample_count": len(matched),
        "unmatched_expression_samples": ";".join(expr_norm[key] for key in sorted(set(expr_norm) - set(meta_norm))[:10]),
        "unmatched_metadata_samples": ";".join(meta_norm[key] for key in sorted(set(meta_norm) - set(expr_norm))[:10]),
    }


def _standardize_genes_rows(
    expression: pd.DataFrame,
    gene_col: str,
    sample_columns: list[str],
    axis_genes: set[str],
    alias_map: dict[str, str],
    gene_map: dict[str, str],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    mapping_rows: list[dict[str, object]] = []
    mapped_rows = []
    for _, row in expression.iterrows():
        raw_feature = row.get(gene_col, "")
        gene = _map_feature_to_gene(raw_feature, axis_genes, alias_map, gene_map)
        mapping_rows.append(
            {
                "input_feature_id": raw_feature,
                "mapped_gene_symbol": gene or "",
                "retained": "true" if gene else "false",
            }
        )
        if gene:
            values = pd.to_numeric(row[sample_columns], errors="coerce")
            mapped_rows.append(pd.Series([gene, *values.tolist()], index=["gene_symbol", *sample_columns]))
    if not mapped_rows:
        return pd.DataFrame(columns=["gene_symbol", *sample_columns]), mapping_rows
    mapped = pd.DataFrame(mapped_rows)
    numeric = _coerce_numeric(mapped.drop(columns=["gene_symbol"]))
    numeric.insert(0, "gene_symbol", mapped["gene_symbol"].values)
    grouped = numeric.groupby("gene_symbol", as_index=False).mean(numeric_only=True)
    return grouped, mapping_rows


def _standardize_samples_rows(
    expression: pd.DataFrame,
    sample_col: str,
    axis_genes: set[str],
    alias_map: dict[str, str],
    gene_map: dict[str, str],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    mapping_rows: list[dict[str, object]] = []
    retained: dict[str, str] = {}
    for col in expression.columns:
        if col == sample_col:
            continue
        gene = _map_feature_to_gene(col, axis_genes, alias_map, gene_map)
        mapping_rows.append({"input_feature_id": col, "mapped_gene_symbol": gene or "", "retained": "true" if gene else "false"})
        if gene:
            retained[col] = gene
    if not retained:
        return pd.DataFrame(columns=["gene_symbol"]), mapping_rows
    work = expression[[sample_col, *retained]].copy()
    work = work.rename(columns={sample_col: "sample_id", **retained})
    numeric = _coerce_numeric(work.drop(columns=["sample_id"]))
    numeric.insert(0, "sample_id", work["sample_id"].astype(str).str.strip().values)
    long = numeric.melt(id_vars="sample_id", var_name="gene_symbol", value_name="expression")
    pivot = long.pivot_table(index="gene_symbol", columns="sample_id", values="expression", aggfunc="mean").reset_index()
    return pivot, mapping_rows


def _standardize_long(
    expression: pd.DataFrame,
    sample_col: str,
    gene_col: str,
    value_col: str,
    axis_genes: set[str],
    alias_map: dict[str, str],
    gene_map: dict[str, str],
) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    mapping_rows: list[dict[str, object]] = []
    work = expression[[sample_col, gene_col, value_col]].copy()
    work["sample_id"] = work[sample_col].astype(str).str.strip()
    work["mapped_gene_symbol"] = [
        _map_feature_to_gene(value, axis_genes, alias_map, gene_map) for value in work[gene_col]
    ]
    for raw, mapped in zip(work[gene_col], work["mapped_gene_symbol"], strict=False):
        mapping_rows.append({"input_feature_id": raw, "mapped_gene_symbol": mapped or "", "retained": "true" if mapped else "false"})
    work = work[work["mapped_gene_symbol"].notna()].copy()
    if work.empty:
        return pd.DataFrame(columns=["gene_symbol"]), mapping_rows
    work["expression"] = pd.to_numeric(work[value_col], errors="coerce")
    pivot = work.pivot_table(
        index="mapped_gene_symbol",
        columns="sample_id",
        values="expression",
        aggfunc="mean",
    ).reset_index()
    pivot = pivot.rename(columns={"mapped_gene_symbol": "gene_symbol"})
    return pivot, mapping_rows


def standardize_expression_table(
    expression_path: Path,
    metadata: pd.DataFrame,
    axis_registry: Path,
    sample_id_column: str,
    orientation: str = "auto",
    gene_id_column: str = "auto",
    gene_map: Path | None = None,
    alias_table: Path | None = None,
) -> tuple[pd.DataFrame, list[dict[str, object]], dict[str, object]]:
    expression = _read_table(expression_path)
    orientation = infer_expression_orientation(expression_path) if orientation == "auto" else orientation
    axis_genes = _axis_gene_set(axis_registry)
    alias_map = _load_alias_map(alias_table)
    explicit_gene_map = _load_gene_map(gene_map)
    if orientation == "long":
        sample_col = infer_sample_id_column(expression.columns)
        gene_col = infer_gene_id_column(expression.columns)
        value_col = _find_column(expression.columns, VALUE_CANDIDATES)
        if not sample_col or not gene_col or not value_col:
            raise ValueError("Long expression format requires sample, gene/probe, and expression value columns.")
        standardized, mapping_rows = _standardize_long(
            expression, sample_col, gene_col, value_col, axis_genes, alias_map, explicit_gene_map
        )
    elif orientation == "samples_rows":
        sample_col = infer_sample_id_column(expression.columns)
        if not sample_col:
            raise ValueError("Samples-row expression format requires a sample identifier column.")
        standardized, mapping_rows = _standardize_samples_rows(
            expression, sample_col, axis_genes, alias_map, explicit_gene_map
        )
    elif orientation == "genes_rows":
        gene_col = infer_gene_id_column(expression.columns) if gene_id_column == "auto" else gene_id_column
        gene_col = gene_col or expression.columns[0]
        sample_columns, _ = infer_sample_columns(expression, metadata)
        sample_columns = [col for col in sample_columns if col != gene_col]
        standardized, mapping_rows = _standardize_genes_rows(
            expression, gene_col, sample_columns, axis_genes, alias_map, explicit_gene_map
        )
    else:
        raise ValueError("orientation must be auto, genes_rows, samples_rows, or long")
    retained_genes = sorted(set(standardized.get("gene_symbol", [])))
    mapped_input_features = [
        str(row.get("mapped_gene_symbol", "")).strip().upper()
        for row in mapping_rows
        if str(row.get("retained", "")).strip().casefold() == "true"
    ]
    mapped_feature_counts = Counter(gene for gene in mapped_input_features if gene)
    schema = {
        "expression_format": orientation,
        "gene_identifier_type": detect_gene_identifier_type(expression.iloc[:, 0].dropna().head(200)),
        "input_feature_count": len(mapping_rows),
        "retained_input_feature_count": len(mapped_input_features),
        "unmapped_input_feature_count": len(mapping_rows) - len(mapped_input_features),
        "multi_feature_gene_count": sum(count > 1 for count in mapped_feature_counts.values()),
        "retained_axis_gene_count": len(retained_genes),
        "retained_axis_genes": ";".join(retained_genes),
    }
    return standardized, mapping_rows, schema


def standardize_metadata_table(
    metadata_path: Path,
    sample_id_column: str = "auto",
    endpoint_column: str = "auto",
    positive_class: str = "auto",
    negative_class: str = "auto",
) -> tuple[pd.DataFrame, dict[str, object]]:
    metadata = _read_table(metadata_path)
    metadata.columns = [str(col).strip() for col in metadata.columns]
    metadata = metadata.dropna(axis=1, how="all")
    sample_col = infer_sample_id_column(metadata.columns) if sample_id_column == "auto" else sample_id_column
    endpoint_col = infer_endpoint_column(metadata) if endpoint_column == "auto" else endpoint_column
    if not sample_col or sample_col not in metadata.columns:
        raise ValueError("Could not infer sample identifier column. Pass --sample-id-column explicitly.")
    if not endpoint_col or endpoint_col not in metadata.columns:
        raise ValueError("Could not infer endpoint column. Pass --endpoint-column explicitly.")
    if positive_class == "auto" or negative_class == "auto":
        inferred_positive, inferred_negative = infer_positive_negative_classes(metadata, endpoint_col)
        positive_class = inferred_positive if positive_class == "auto" else positive_class
        negative_class = inferred_negative if negative_class == "auto" else negative_class
    positive_norm = normalize_class(positive_class)
    negative_norm = normalize_class(negative_class)
    rows = []
    ambiguous = []
    for _, row in metadata.iterrows():
        sample_id = str(row.get(sample_col, "")).strip()
        endpoint_value = str(row.get(endpoint_col, "")).strip()
        endpoint_norm = normalize_class(endpoint_value)
        label = None
        if endpoint_norm == positive_norm or endpoint_norm in POSITIVE_TOKENS and positive_class == "auto":
            label = 1
        elif endpoint_norm == negative_norm or endpoint_norm in NEGATIVE_TOKENS and negative_class == "auto":
            label = 0
        elif endpoint_norm == positive_norm:
            label = 1
        elif endpoint_norm == negative_norm:
            label = 0
        if sample_id and label is not None:
            rows.append(
                {
                    "sample_id": sample_id,
                    "endpoint": endpoint_value,
                    "label__endpoint": label,
                    "research_use_only": "true",
                }
            )
        elif sample_id:
            ambiguous.append(sample_id)
    if not rows:
        raise ValueError(
            f"No unambiguous endpoint labels found in {endpoint_col!r}; pass explicit classes or inspect metadata."
        )
    if ambiguous and len(rows) < 2:
        raise ValueError(f"Endpoint labels are ambiguous for samples: {ambiguous[:10]}")
    standardized = pd.DataFrame(rows).drop_duplicates(subset=["sample_id"], keep="first")
    schema = {
        "sample_id_column": sample_col,
        "endpoint_column": endpoint_col,
        "positive_class": positive_class,
        "negative_class": negative_class,
        "metadata_rows": len(metadata),
        "standardized_metadata_rows": len(standardized),
        "ambiguous_label_count": len(ambiguous),
    }
    return standardized, schema


def _write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_frame(path: Path, frame: pd.DataFrame, gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if gzip_output or path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression)


def write_ingest_report(
    path: Path,
    schema: dict[str, object],
    join: dict[str, object],
    warnings: list[str],
) -> None:
    lines = [
        "# NeuroFate Ingest Report",
        "",
        RESEARCH_USE_NOTICE,
        "",
        "## Detected Schema",
    ]
    lines.extend(f"- {key}: {value}" for key, value in schema.items())
    lines.extend(["", "## Expression-Metadata Join"])
    lines.extend(f"- {key}: {value}" for key, value in join.items())
    lines.extend(["", "## Warnings"])
    lines.extend(f"- {warning}" for warning in warnings) if warnings else lines.append("- none")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_ingest(config: IngestConfig) -> IngestResult:
    config.outdir.mkdir(parents=True, exist_ok=True)
    if config.assist:
        # The CLI currently reports ambiguity rather than silently guessing. This
        # placeholder records the user's intent while keeping non-interactive
        # behavior deterministic for tests and batch pipelines.
        pass
    metadata, metadata_schema = standardize_metadata_table(
        config.metadata,
        sample_id_column=config.sample_id_column,
        endpoint_column=config.endpoint_column,
        positive_class=config.positive_class,
        negative_class=config.negative_class,
    )
    expression, mapping_rows, expression_schema = standardize_expression_table(
        config.expression,
        metadata,
        config.axis_registry,
        sample_id_column="sample_id",
        orientation=config.orientation,
        gene_id_column=config.gene_id_column,
        gene_map=config.gene_map,
        alias_table=config.alias_table,
    )
    if expression.empty or expression.shape[0] < config.min_axis_genes:
        raise ValueError(
            f"Only {expression.shape[0]} NeuroFate axis genes were retained; "
            f"minimum required is {config.min_axis_genes}. Check gene symbols, Ensembl IDs, or --gene-map."
        )
    expression_samples = [col for col in expression.columns if col != "gene_symbol"]
    join = validate_expression_metadata_join(expression_samples, metadata["sample_id"])
    if int(join["matched_sample_count"]) == 0:
        raise ValueError("No expression samples matched metadata sample IDs after safe normalization.")
    matched_norm = {normalize_id(col) for col in expression_samples}
    metadata = metadata[metadata["sample_id"].map(normalize_id).isin(matched_norm)].copy()
    if metadata.empty:
        raise ValueError("No standardized metadata rows remained after expression join filtering.")
    matched_samples = set(metadata["sample_id"])
    retained_columns = ["gene_symbol"] + [col for col in expression.columns if col in matched_samples]
    # Preserve normalized matching for cases where expression columns differ only by safe punctuation.
    if len(retained_columns) == 1:
        expression_lookup = {normalize_id(col): col for col in expression_samples}
        retained_columns = ["gene_symbol"] + [expression_lookup[normalize_id(sample)] for sample in metadata["sample_id"]]
        expression = expression.rename(columns={expression_lookup[normalize_id(sample)]: sample for sample in metadata["sample_id"]})
        retained_columns = ["gene_symbol", *metadata["sample_id"].tolist()]
    expression = expression[retained_columns]

    warnings = []
    total_axis_genes = len(_axis_gene_set(config.axis_registry))
    retained_axis_genes = int(expression_schema.get("retained_axis_gene_count", 0))
    if retained_axis_genes < total_axis_genes:
        warnings.append(f"axis_gene_coverage_incomplete={retained_axis_genes}/{total_axis_genes}")
    unmapped_features = int(expression_schema.get("unmapped_input_feature_count", 0))
    if unmapped_features:
        warnings.append(f"unmapped_input_features={unmapped_features}")
    multi_feature_genes = int(expression_schema.get("multi_feature_gene_count", 0))
    if multi_feature_genes:
        warnings.append(f"multiple_input_features_per_gene={multi_feature_genes}")
    if int(join["matched_sample_count"]) < int(join["expression_sample_count"]):
        warnings.append("some_expression_samples_did_not_match_metadata")
    if int(join["matched_sample_count"]) < int(join["metadata_sample_count"]):
        warnings.append("some_metadata_samples_did_not_match_expression")

    paths = {
        "standardized_expression": config.outdir / "standardized_expression.tsv.gz",
        "standardized_metadata": config.outdir / "standardized_metadata.tsv",
        "input_schema_detected": config.outdir / "input_schema_detected.tsv",
        "expression_metadata_join": config.outdir / "expression_metadata_join.tsv",
        "gene_mapping_report": config.outdir / "gene_mapping_report.tsv",
        "ingest_warnings": config.outdir / "ingest_warnings.tsv",
        "ingest_report": config.outdir / "ingest_report.md",
        "run_config": config.outdir / "run_config.yaml",
    }
    _write_frame(paths["standardized_expression"], expression, gzip_output=True)
    _write_frame(paths["standardized_metadata"], metadata)
    schema = {**metadata_schema, **expression_schema}
    _write_rows(
        paths["input_schema_detected"],
        [{"key": key, "value": value} for key, value in schema.items()],
        ["key", "value"],
    )
    _write_rows(paths["expression_metadata_join"], [join], list(join))
    _write_rows(
        paths["gene_mapping_report"],
        mapping_rows,
        ["input_feature_id", "mapped_gene_symbol", "retained"],
    )
    _write_rows(
        paths["ingest_warnings"],
        [{"warning": warning} for warning in warnings] or [{"warning": "none"}],
        ["warning"],
    )
    write_ingest_report(paths["ingest_report"], schema, join, warnings)
    config_payload = {
        "expression": str(config.expression),
        "metadata": str(config.metadata),
        "axis_registry": str(config.axis_registry),
        "gene_map": str(config.gene_map) if config.gene_map else "",
        "alias_table": str(config.alias_table) if config.alias_table else "",
        "research_use_only": RESEARCH_USE_NOTICE,
        **schema,
    }
    paths["run_config"].write_text(yaml.safe_dump(config_payload, sort_keys=False), encoding="utf-8")
    return IngestResult(**paths)


def run_complete_workflow(config: IngestConfig) -> dict[str, Path]:
    ingest_outdir = config.outdir / "ingest"
    scoring_outdir = config.outdir / "axis"
    risk_outdir = config.outdir / "risk"
    ingest_result = run_ingest(IngestConfig(**{**config.__dict__, "outdir": ingest_outdir}))
    score_outputs = build_axis_score_tables(
        expression=ingest_result.standardized_expression,
        metadata=ingest_result.standardized_metadata,
        axis_registry=config.axis_registry,
        outdir=scoring_outdir,
        sample_id_column="sample_id",
        endpoint_column="label__endpoint",
        positive_class="1",
        negative_class="0",
        orientation="genes_rows",
        gene_column="gene_symbol",
    )
    risk_outputs = score_research_risk(score_outputs["axis_scores"], risk_outdir)
    report_path = config.outdir / "neurofate_run_report.md"
    workflow_config_path = config.outdir / "run_config.yaml"
    report_path.write_text(
        "\n".join(
            [
                "# NeuroFate Run Report",
                "",
                RESEARCH_USE_NOTICE,
                "",
                "## Outputs",
                f"- Ingest report: `{ingest_result.ingest_report}`",
                f"- Axis scores: `{score_outputs['axis_scores']}`",
                f"- Axis feature coverage: `{score_outputs['axis_feature_coverage']}`",
                f"- Research-use risk scores: `{risk_outputs['risk_scores']}`",
                f"- Research-use risk report: `{risk_outputs['risk_score_report']}`",
                "",
                "This report summarizes a research workflow for cohort-level transcriptomic analysis.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    workflow_config_path.write_text(
        yaml.safe_dump(
            {
                "expression": str(config.expression),
                "metadata": str(config.metadata),
                "axis_registry": str(config.axis_registry),
                "gene_map": str(config.gene_map) if config.gene_map else "",
                "alias_table": str(config.alias_table) if config.alias_table else "",
                "research_use_only": RESEARCH_USE_NOTICE,
                "ingest_outdir": str(ingest_outdir),
                "axis_outdir": str(scoring_outdir),
                "risk_outdir": str(risk_outdir),
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return {
        "standardized_expression": ingest_result.standardized_expression,
        "standardized_metadata": ingest_result.standardized_metadata,
        "axis_scores": score_outputs["axis_scores"],
        "axis_feature_coverage": score_outputs["axis_feature_coverage"],
        "label_summary": score_outputs["label_summary"],
        "neurofate_risk_scores": risk_outputs["risk_scores"],
        "risk_score_report": risk_outputs["risk_score_report"],
        "neurofate_run_report": report_path,
        "warnings": ingest_result.ingest_warnings,
        "run_config": workflow_config_path,
        "ingest_run_config": ingest_result.run_config,
    }
