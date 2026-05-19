#!/usr/bin/env python3
"""Build Phase 34 sample-level PD axis scores from GEO expression files."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
from pathlib import Path


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
    return "\t" if path.name.endswith((".txt", ".txt.gz", ".tsv", ".tsv.gz")) else ","


def clean(value: str | None) -> str:
    text = str(value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1]
    return text.strip()


def normalize(value: str | None) -> str:
    return "".join(ch for ch in clean(value).casefold() if ch.isalnum())


POSITIVE_PD_LABELS = {
    "1",
    "pd",
    "parkinson",
    "parkinsons",
    "parkinsondisease",
    "parkinsonsdisease",
    "parkinsonsdiseasecase",
}

NEGATIVE_PD_LABELS = {
    "0",
    "control",
    "controls",
    "unaffectedcontrol",
    "unaffectedcontrols",
    "normal",
    "healthycontrol",
    "healthycontrols",
}


def canonical_pd_label(value: str | None) -> str:
    raw = clean(value)
    try:
        number = float(raw)
        if number == 1:
            return "1"
        if number == 0:
            return "0"
    except ValueError:
        pass
    token = normalize(value)
    if token in POSITIVE_PD_LABELS:
        return "1"
    if token in NEGATIVE_PD_LABELS:
        return "0"
    return ""


def column_lookup(rows: list[dict[str, str]]) -> dict[str, str]:
    columns: dict[str, str] = {}
    for row in rows:
        for key in row:
            columns.setdefault(key.casefold(), key)
    return columns


def select_label_column(rows: list[dict[str, str]]) -> str:
    columns = column_lookup(rows)
    preferred = columns.get("label__pd_vs_control")
    if preferred and any(canonical_pd_label(row.get(preferred)) for row in rows):
        return preferred
    for candidate in ["disease_state", "diagnosis", "condition", "group", "phenotype"]:
        column = columns.get(candidate)
        if column and any(canonical_pd_label(row.get(column)) for row in rows):
            return column
    return preferred or ""


def select_endpoint_status_column(rows: list[dict[str, str]]) -> str:
    return column_lookup(rows).get("endpoint_status", "")


def endpoint_status_allows(value: str | None) -> bool:
    token = normalize(value)
    return token in {"", "unambiguous"}


def metadata_identifier(row: dict[str, str]) -> str:
    return clean(row.get("geo_accession") or row.get("sample_id") or row.get("sample_title") or "")


def prepare_metadata_for_join(metadata_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, str]]:
    label_column = select_label_column(metadata_rows)
    status_column = select_endpoint_status_column(metadata_rows)
    metadata_columns = sorted({key for row in metadata_rows for key in row})
    status_filtered: list[dict[str, str]] = []
    filtered_examples: list[str] = []
    for row in metadata_rows:
        if status_column and not endpoint_status_allows(row.get(status_column)):
            if len(filtered_examples) < 5:
                filtered_examples.append(f"{metadata_identifier(row)}:endpoint_status={clean(row.get(status_column))}")
            continue
        status_filtered.append(row)
    label_counts = {"0": 0, "1": 0}
    valid_rows: list[dict[str, str]] = []
    observed_values = sorted({clean(row.get(label_column)) for row in status_filtered}) if label_column else []
    for row in status_filtered:
        label = canonical_pd_label(row.get(label_column)) if label_column else ""
        if label not in {"0", "1"}:
            if len(filtered_examples) < 5:
                filtered_examples.append(f"{metadata_identifier(row)}:invalid_label={clean(row.get(label_column)) if label_column else ''}")
            continue
        fixed = dict(row)
        fixed["label__pd_vs_control"] = label
        fixed["_selected_label_column"] = label_column
        fixed["_selected_endpoint_status_column"] = status_column
        valid_rows.append(fixed)
        label_counts[label] += 1
    debug = {
        "metadata_columns_present": ";".join(metadata_columns),
        "selected_label_column": label_column,
        "selected_endpoint_status_column": status_column,
        "metadata_rows_before_filtering": str(len(metadata_rows)),
        "metadata_rows_after_endpoint_status_filtering": str(len(status_filtered)),
        "label_values_observed": ";".join(observed_values),
        "label_counts": f"0={label_counts['0']};1={label_counts['1']}",
        "filtered_out_metadata_examples": ";".join(example for example in filtered_examples if example),
    }
    return valid_rows, debug


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


DEBUG_COLUMNS = [
    "metadata_columns_present",
    "selected_join_key",
    "selected_label_column",
    "selected_endpoint_status_column",
    "metadata_rows_before_filtering",
    "metadata_rows_after_endpoint_status_filtering",
    "label_values_observed",
    "label_counts",
    "expression_sample_count",
    "matched_sample_count",
    "final_labeled_matched_sample_count",
    "unmatched_expression_sample_examples",
    "filtered_out_metadata_examples",
]


class BuilderError(Exception):
    def __init__(self, message: str, debug: dict[str, str] | None = None) -> None:
        super().__init__(message)
        self.debug = debug or {}


def to_float(value: str | None) -> float:
    try:
        return float(clean(value))
    except ValueError:
        return math.nan


def parse_genes(value: str) -> list[str]:
    return [gene.strip().upper() for gene in value.replace(",", ";").split(";") if gene.strip()]


def axis_map(axes: list[dict[str, str]]) -> dict[str, list[str]]:
    return {row["axis_id"]: parse_genes(row.get("gene_members", "")) for row in axes}


def mean(values: list[float]) -> float:
    obs = [value for value in values if not math.isnan(value)]
    return sum(obs) / len(obs) if obs else math.nan


def standardize(values: list[float]) -> list[float]:
    obs = [value for value in values if not math.isnan(value)]
    if not obs:
        return [math.nan for _ in values]
    mu = sum(obs) / len(obs)
    sd = math.sqrt(sum((value - mu) ** 2 for value in obs) / max(1, len(obs) - 1)) or 1.0
    return [(value - mu) / sd if not math.isnan(value) else math.nan for value in values]


def probe_map(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    for row in rows:
        probe = clean(row.get("probe_id", ""))
        gene = clean(row.get("gene_symbol", "")).upper()
        if probe and gene:
            mapping.setdefault(probe, []).append(gene)
            mapping.setdefault(probe.upper(), []).append(gene)
    return {key: sorted(set(value)) for key, value in mapping.items()}


def expression_sample_ids(expression: dict[str, dict[str, list[float]]]) -> list[str]:
    return sorted(expression)


def metadata_lookup(rows: list[dict[str, str]], key: str, normalized: bool = False) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        if canonical_pd_label(row.get("label__pd_vs_control")) not in {"0", "1"}:
            continue
        value = normalize(row.get(key, "")) if normalized else clean(row.get(key, ""))
        if value:
            out[value] = row
    return out


def select_metadata_join(expression_samples: list[str], metadata_rows: list[dict[str, str]]) -> tuple[dict[str, dict[str, str]], list[dict[str, str]], str]:
    candidates = [
        ("geo_accession", False),
        ("sample_id", False),
        ("sample_title", False),
        ("sample_title", True),
    ]
    best_key = ""
    best_matches: dict[str, dict[str, str]] = {}
    best_join_rows: list[dict[str, str]] = []
    for key, normalized in candidates:
        label = f"normalized_{key}" if normalized else key
        lookup = metadata_lookup(metadata_rows, key, normalized)
        matches: dict[str, dict[str, str]] = {}
        join_rows: list[dict[str, str]] = []
        for sample in expression_samples:
            join_value = normalize(sample) if normalized else sample
            metadata = lookup.get(join_value)
            if metadata:
                matches[sample] = metadata
                status = "matched"
            else:
                status = "unmatched_expression_sample"
            join_rows.append(
                {
                    "expression_sample_id": sample,
                    "selected_join_key": label,
                    "join_value": join_value,
                    "metadata_sample_id": metadata.get("sample_id", "") if metadata else "",
                    "metadata_geo_accession": metadata.get("geo_accession", "") if metadata else "",
                    "metadata_sample_title": metadata.get("sample_title", "") if metadata else "",
                    "label__pd_vs_control": metadata.get("label__pd_vs_control", "") if metadata else "",
                    "join_status": status,
                }
            )
        if len(matches) > len(best_matches):
            best_key = label
            best_matches = matches
            best_join_rows = join_rows
    metadata_ids = {metadata_identifier(row) for row in metadata_rows if canonical_pd_label(row.get("label__pd_vs_control")) in {"0", "1"}}
    matched_metadata_ids = {
        row.get("metadata_geo_accession") or row.get("metadata_sample_id") or row.get("metadata_sample_title")
        for row in best_join_rows
        if row.get("join_status") == "matched"
    }
    for sample_id in sorted(metadata_ids - matched_metadata_ids):
        best_join_rows.append(
            {
                "expression_sample_id": "",
                "selected_join_key": best_key,
                "join_value": "",
                "metadata_sample_id": sample_id,
                "metadata_geo_accession": sample_id,
                "metadata_sample_title": "",
                "label__pd_vs_control": "",
                "join_status": "unmatched_metadata_sample",
            }
        )
    return best_matches, best_join_rows, best_key


def detect_orientation(path: Path, requested: str) -> str:
    if requested != "auto":
        return requested
    saw_series_preamble = False
    with open_text(path) as handle:
        for row in csv.reader(handle, delimiter="\t" if "series_matrix" in path.name else delimiter_for(path)):
            if not row:
                continue
            first = clean(row[0]).lower()
            if first == "!series_matrix_table_begin":
                return "probes_rows"
            if first.startswith("!series_") or first.startswith("!sample_") or first.startswith("!platform_"):
                saw_series_preamble = True
                continue
            if first in {"id_ref", "id", "probe_id", "probe", "gene"}:
                return "probes_rows"
            if saw_series_preamble:
                continue
            return "samples_rows"
    return "probes_rows"


def iter_expression_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    rows: list[list[str]] = []
    header: list[str] = []
    in_series_table = False
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t" if "series_matrix" in path.name else delimiter_for(path))
        for row in reader:
            if not row:
                continue
            key = clean(row[0])
            if key == "!series_matrix_table_begin":
                in_series_table = True
                continue
            if key == "!series_matrix_table_end":
                break
            if in_series_table or key.lower() in {"id_ref", "id", "probe_id", "probe", "gene"}:
                if not header:
                    header = [clean(value) for value in row]
                    continue
                rows.append([clean(value) for value in row])
            elif header:
                rows.append([clean(value) for value in row])
    if not header:
        with open_text(path) as handle:
            reader = csv.reader(handle, delimiter=delimiter_for(path))
            header = [clean(value) for value in next(reader)]
            rows = [[clean(value) for value in row] for row in reader]
    return header, rows


def read_probes_rows(path: Path, mapping: dict[str, list[str]]) -> dict[str, dict[str, list[float]]]:
    header, rows = iter_expression_rows(path)
    samples = header[1:]
    values: dict[str, dict[str, list[float]]] = {}
    for row in rows:
        if not row:
            continue
        genes = mapping.get(row[0], []) or mapping.get(row[0].upper(), [])
        if not genes:
            continue
        for sample, raw in zip(samples, row[1:], strict=False):
            value = to_float(raw)
            if math.isnan(value):
                continue
            for gene in genes:
                values.setdefault(sample, {}).setdefault(gene, []).append(value)
    return values


def read_samples_rows(path: Path, mapping: dict[str, list[str]]) -> dict[str, dict[str, list[float]]]:
    values: dict[str, dict[str, list[float]]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(path))
        fields = reader.fieldnames or []
        sample_col = next((field for field in fields if field.lower() in {"sample_id", "sample", "geo_accession", "gsm"}), fields[0])
        probe_cols = [field for field in fields if field in mapping or field.upper() in mapping]
        for row in reader:
            sample = clean(row.get(sample_col, ""))
            if not sample:
                continue
            for probe in probe_cols:
                for gene in mapping.get(probe, []) + mapping.get(probe.upper(), []):
                    value = to_float(row.get(probe))
                    if not math.isnan(value):
                        values.setdefault(sample, {}).setdefault(gene, []).append(value)
    return values


def build_scores(
    expression: dict[str, dict[str, list[float]]],
    metadata_rows: list[dict[str, str]],
    axes: list[dict[str, str]],
    cohort_id: str,
    expected_sample_count: int = 0,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], str, dict[str, str]]:
    sample_gene = {sample: {gene: mean(values) for gene, values in genes.items()} for sample, genes in expression.items()}
    valid_metadata_rows, debug = prepare_metadata_for_join(metadata_rows)
    if not debug.get("selected_label_column") or not valid_metadata_rows:
        raise BuilderError("No valid binary PD/control label column was found in metadata.", debug)
    metadata, join_rows, join_key = select_metadata_join(sorted(sample_gene), valid_metadata_rows)
    samples = sorted(sample for sample in sample_gene if sample in metadata)
    debug.update(
        {
            "selected_join_key": join_key,
            "expression_sample_count": str(len(sample_gene)),
            "matched_sample_count": str(sum(row.get("join_status") == "matched" for row in join_rows)),
            "final_labeled_matched_sample_count": str(len(samples)),
            "unmatched_expression_sample_examples": ";".join(row.get("expression_sample_id", "") for row in join_rows if row.get("join_status") == "unmatched_expression_sample" and row.get("expression_sample_id"))[:500],
        }
    )
    if not samples:
        raise BuilderError("No expression samples matched unambiguous PD/control metadata.", debug)
    if expected_sample_count and len(samples) != expected_sample_count:
        raise BuilderError(f"Expected {expected_sample_count} matched samples but found {len(samples)} using {join_key}.", debug)
    labels = {canonical_pd_label(metadata[sample].get("label__pd_vs_control")) for sample in samples}
    if labels != {"0", "1"}:
        raise BuilderError("PD/control endpoint is ambiguous or one class is missing.", debug)
    amap = axis_map(axes)
    raw_scores = {axis_id: [] for axis_id in amap}
    for sample in samples:
        for axis_id, genes in amap.items():
            raw_scores[axis_id].append(mean([sample_gene[sample].get(gene, math.nan) for gene in genes]))
    standardized = {axis_id: standardize(values) for axis_id, values in raw_scores.items()}
    rows: list[dict[str, str]] = []
    counts: dict[str, int] = {}
    for index, sample in enumerate(samples):
        meta = metadata[sample]
        label = canonical_pd_label(meta.get("label__pd_vs_control"))
        counts[label] = counts.get(label, 0) + 1
        row = {
            "cohort_id": cohort_id,
            "sample_id": meta.get("sample_id", "") or meta.get("geo_accession", "") or sample,
            "expression_sample_id": sample,
            "metadata_join_key": join_key,
            "label__pd_vs_control": label,
            "disease_state": meta.get("disease_state", ""),
            "tissue_or_region": meta.get("tissue_or_region", ""),
        }
        for axis_id in amap:
            value = standardized[axis_id][index]
            row[f"axis__{axis_id}"] = "" if math.isnan(value) else f"{value:.8g}"
        rows.append(row)
    found_universe = {gene for genes in sample_gene.values() for gene in genes}
    coverage = []
    for axis_id, genes in amap.items():
        found = sorted(set(genes) & found_universe)
        missing = sorted(set(genes) - set(found))
        coverage.append({"cohort_id": cohort_id, "axis_id": axis_id, "genes_requested": str(len(genes)), "genes_found": str(len(found)), "genes_missing": str(len(missing)), "found_gene_members": ";".join(found), "missing_gene_members": ";".join(missing), "status": "ok" if found else "insufficient_coverage"})
    label_rows = [{"cohort_id": cohort_id, "label__pd_vs_control": key, "count": str(value)} for key, value in sorted(counts.items())]
    return rows, coverage, label_rows, join_rows, join_key, debug


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Phase 34 PD axis scores from GEO expression files.")
    parser.add_argument("--expression-file", type=Path, required=True)
    parser.add_argument("--sample-metadata", type=Path, required=True)
    parser.add_argument("--probe-map", type=Path, required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--orientation", choices=["auto", "probes_rows", "samples_rows"], default="auto")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coverage-output", type=Path, required=True)
    parser.add_argument("--label-summary-output", type=Path, required=True)
    parser.add_argument("--join-output", type=Path)
    parser.add_argument("--metadata-debug-output", type=Path)
    parser.add_argument("--expected-sample-count", type=int, default=0)
    parser.add_argument("--log-file", type=Path, required=True)
    args = parser.parse_args()
    configure_logging(args.log_file)
    mapping = probe_map(read_tsv(args.probe_map))
    orientation = detect_orientation(args.expression_file, args.orientation)
    expression = read_probes_rows(args.expression_file, mapping) if orientation == "probes_rows" else read_samples_rows(args.expression_file, mapping)
    expected = args.expected_sample_count or (18 if args.cohort_id == "gse20141_pd_snpc_lcm" else 0)
    try:
        rows, coverage, labels, join_rows, join_key, debug = build_scores(expression, read_tsv(args.sample_metadata), read_tsv(args.axis_registry), args.cohort_id, expected)
    except BuilderError as error:
        if args.metadata_debug_output and error.debug:
            write_tsv(args.metadata_debug_output, [error.debug], DEBUG_COLUMNS)
        raise SystemExit(str(error)) from error
    columns = ["cohort_id", "sample_id", "expression_sample_id", "metadata_join_key", "label__pd_vs_control", "disease_state", "tissue_or_region"] + [f"axis__{row['axis_id']}" for row in read_tsv(args.axis_registry)]
    write_tsv(args.output, rows, columns)
    write_tsv(args.coverage_output, coverage, ["cohort_id", "axis_id", "genes_requested", "genes_found", "genes_missing", "found_gene_members", "missing_gene_members", "status"])
    write_tsv(args.label_summary_output, labels, ["cohort_id", "label__pd_vs_control", "count"])
    join_output = args.join_output or args.output.with_name(args.output.stem + "_expression_metadata_join.tsv")
    write_tsv(join_output, join_rows, ["expression_sample_id", "selected_join_key", "join_value", "metadata_sample_id", "metadata_geo_accession", "metadata_sample_title", "label__pd_vs_control", "join_status"])
    if args.metadata_debug_output:
        write_tsv(args.metadata_debug_output, [debug], DEBUG_COLUMNS)
    logging.info("Built Phase 34 axis scores cohort=%s samples=%d orientation=%s join_key=%s", args.cohort_id, len(rows), orientation, join_key)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
