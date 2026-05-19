#!/usr/bin/env python3
"""Prepare Phase 34 platform probe mappings restricted to NeuroFate axis genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from collections import Counter
from pathlib import Path


PROBE_COLUMN_CANDIDATES = ("ID", "ID_REF", "probe_id", "ProbeID", "Probe Set ID", "Reporter Identifier")
GENE_COLUMN_CANDIDATES = (
    "Gene Symbol",
    "Gene symbol",
    "gene_symbol",
    "Symbol",
    "Gene Assignment",
    "gene_assignment",
    "GENE_SYMBOL",
    "Ensembl",
    "ENSEMBL",
    "ENTREZ_GENE_ID",
    "Gene ID",
    "GB_ACC",
    "GenBank Accession",
)


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(), logging.FileHandler(log_file, mode="w", encoding="utf-8")],
    )


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") if path.name.endswith(".gz") else path.open("r", encoding="utf-8", errors="replace", newline="")


def delimiter_for(path: Path) -> str:
    return "\t" if path.name.endswith((".txt", ".txt.gz", ".tsv", ".tsv.gz", ".soft", ".soft.gz", ".annot", ".annot.gz")) else ","


def clean(value: str | None) -> str:
    return str(value or "").strip().strip('"')


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def axis_genes(path: Path) -> set[str]:
    genes: set[str] = set()
    for row in read_tsv(path):
        for gene in row.get("gene_members", "").replace(",", ";").split(";"):
            if gene.strip():
                genes.add(gene.strip().upper())
    return genes


def aliases(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for row in read_tsv(path):
        symbol = clean(row.get("gene_symbol", "")).upper()
        if not symbol:
            continue
        for key in ("ensembl_gene_id", "alias"):
            value = clean(row.get(key, ""))
            if value:
                out[value.upper().split(".", 1)[0]] = symbol
    return out


def choose_column(fieldnames: list[str], candidates: tuple[str, ...]) -> str:
    lowered = {field.casefold(): field for field in fieldnames}
    for candidate in candidates:
        if candidate.casefold() in lowered:
            return lowered[candidate.casefold()]
    for field in fieldnames:
        low = field.casefold()
        if any(candidate.casefold() in low for candidate in candidates):
            return field
    return ""


def choose_gene_columns(fieldnames: list[str]) -> list[str]:
    lowered = {field.casefold(): field for field in fieldnames}
    chosen: list[str] = []
    for candidate in GENE_COLUMN_CANDIDATES:
        field = lowered.get(candidate.casefold())
        if field and field not in chosen:
            chosen.append(field)
    if chosen:
        return chosen
    fallback_tokens = ("symbol", "ensembl", "entrez", "gene assignment", "gene id", "gb_acc", "genbank accession")
    for field in fieldnames:
        low = field.casefold()
        if any(token in low for token in fallback_tokens):
            chosen.append(field)
    return chosen


def parse_symbols(value: str, genes: set[str], alias_map: dict[str, str]) -> list[str]:
    found: set[str] = set()
    for token in re.split(r"\s*///\s*|\s*//\s*|[;,|]", clean(value)):
        for part in re.split(r"\s+", token.strip()):
            upper = part.strip().upper().split(".", 1)[0]
            if not upper:
                continue
            if upper in genes:
                found.add(upper)
            elif upper in alias_map and alias_map[upper] in genes:
                found.add(alias_map[upper])
    return sorted(found)


def candidate_header(fields: list[str]) -> bool:
    return bool(choose_column(fields, PROBE_COLUMN_CANDIDATES) and choose_gene_columns(fields))


def iter_platform_records(path: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    header: list[str] = []
    header_line = 0
    table_begin_line = 0
    rows: list[dict[str, str]] = []
    in_platform_table = False
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            stripped = line.strip()
            if not stripped:
                continue
            if stripped == "!platform_table_begin":
                table_begin_line = line_number
                in_platform_table = True
                continue
            if stripped == "!platform_table_end":
                break
            if not header:
                if stripped.startswith(("^", "!", "#")) and not in_platform_table:
                    continue
                fields = [clean(value) for value in next(csv.reader([line], delimiter=delimiter_for(path)))]
                if in_platform_table or candidate_header(fields):
                    if candidate_header(fields):
                        header = fields
                        header_line = line_number
                    continue
                continue
            if stripped.startswith(("^", "!")):
                break
            if stripped.startswith("#"):
                continue
            values = [clean(value) for value in next(csv.reader([line], delimiter=delimiter_for(path)))]
            row = {field: values[index] if index < len(values) else "" for index, field in enumerate(header)}
            rows.append(row)
    audit = {
        "parser_mode": "geo_platform_table" if table_begin_line else "plain_table",
        "table_begin_line": str(table_begin_line),
        "header_line": str(header_line),
        "header_column_count": str(len(header)),
        "header_columns": ";".join(header),
    }
    if not header:
        audit["warning"] = "no_tabular_header_detected"
    return rows, audit


def parse_platform(path: Path, platform_id: str, genes: set[str], alias_map: dict[str, str]) -> tuple[list[dict[str, str]], dict[str, str]]:
    records, audit = iter_platform_records(path)
    fieldnames = list(records[0].keys()) if records else audit.get("header_columns", "").split(";")
    probe_col = choose_column(fieldnames, PROBE_COLUMN_CANDIDATES)
    gene_cols = choose_gene_columns(fieldnames)
    if not probe_col or not gene_cols:
        audit.update({"platform_id": platform_id, "probe_column": probe_col, "gene_columns": ";".join(gene_cols), "annotation_rows_seen": str(len(records)), "retained_probe_rows": "0", "mapped_gene_count": "0", "mapped_genes": "", "warning": "probe_or_gene_columns_not_detected"})
        raise SystemExit(f"Could not identify probe/gene columns in platform annotation. Columns={fieldnames}")
    raw_rows: list[dict[str, str]] = []
    for record in records:
        probe = clean(record.get(probe_col, ""))
        if not probe:
            continue
        symbols: set[str] = set()
        source_fields: list[str] = []
        for gene_col in gene_cols:
            parsed = parse_symbols(record.get(gene_col, ""), genes, alias_map)
            if parsed:
                symbols.update(parsed)
                source_fields.append(gene_col)
        for symbol in sorted(symbols):
            raw_rows.append({"platform_id": platform_id, "probe_id": probe, "gene_symbol": symbol, "source_gene_field": ";".join(source_fields)})
    counts = Counter(row["gene_symbol"] for row in raw_rows)
    for row in raw_rows:
        row["multi_probe_gene"] = "true" if counts[row["gene_symbol"]] > 1 else "false"
        row["mapping_status"] = "neurofate_axis_probe"
    audit.update(
        {
            "platform_id": platform_id,
            "platform_file": str(path),
            "probe_column": probe_col,
            "gene_columns": ";".join(gene_cols),
            "annotation_rows_seen": str(len(records)),
            "retained_probe_rows": str(len(raw_rows)),
            "mapped_gene_count": str(len(counts)),
            "mapped_genes": ";".join(sorted(counts)),
            "warning": "" if raw_rows else "no_neurofate_axis_probes_retained",
        }
    )
    return raw_rows, audit


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare Phase 34 probe mapping restricted to NeuroFate axis genes.")
    parser.add_argument("--platform-file", type=Path, required=True)
    parser.add_argument("--platform-id", required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--alias-table", type=Path, default=Path("metadata/neurofate_axis_gene_aliases.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path)
    parser.add_argument("--log-file", type=Path, required=True)
    args = parser.parse_args()
    configure_logging(args.log_file)
    rows, audit = parse_platform(args.platform_file, args.platform_id, axis_genes(args.axis_registry), aliases(args.alias_table))
    write_tsv(args.output, rows, ["platform_id", "probe_id", "gene_symbol", "source_gene_field", "multi_probe_gene", "mapping_status"])
    audit_output = args.audit_output or Path(f"results/tables/phase34_{args.platform_id}_platform_parse_audit.tsv")
    write_tsv(
        audit_output,
        [audit],
        [
            "platform_id",
            "platform_file",
            "parser_mode",
            "table_begin_line",
            "header_line",
            "header_column_count",
            "probe_column",
            "gene_columns",
            "annotation_rows_seen",
            "retained_probe_rows",
            "mapped_gene_count",
            "mapped_genes",
            "warning",
        ],
    )
    logging.info("Wrote Phase 34 platform probe mappings rows=%d", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
