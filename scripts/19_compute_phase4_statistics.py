#!/usr/bin/env python3
"""Compute Phase 4 donor-aware statistics from sparse target-gene expression.

Inputs are the pre-extracted sparse target-gene TSV and decoded metadata TSV.
The implementation streams expression rows, stores only compact per-donor and
per-group aggregates, and writes statistical tables for manuscript review.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import math
from collections import defaultdict
from pathlib import Path
from typing import TextIO

try:
    from scipy import stats as scipy_stats
except ImportError:  # pragma: no cover - the environment file includes scipy.
    scipy_stats = None


DEFAULT_MAX_ROW_CHUNK = 1_000_000
MAX_ALLOWED_ROW_CHUNK = 5_000_000

DONOR_FIELD = "Donor ID"
CELLTYPE_FIELD = "Subclass"
STAT_VARIABLES = [
    "Braak",
    "CERAD score",
    "Cognitive Status",
    "Overall AD neuropathological Change",
]
APOE_VARIABLE = "APOE Genotype"
MIXED_PATHOLOGY_VARIABLES = [
    "Highest Lewy Body Disease",
    "LATE",
    "Overall CAA Score",
]
GROUP_VARIABLES = [CELLTYPE_FIELD, *STAT_VARIABLES, APOE_VARIABLE, *MIXED_PATHOLOGY_VARIABLES]
METADATA_FIELDS = [DONOR_FIELD, *GROUP_VARIABLES]

MICROGLIAL_GENES = ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1"]
ASTROCYTE_GENES = ["GFAP"]
NEURONAL_GENES = ["SLC17A7", "SST", "PVALB", "LAMP5"]
NEURODEGENERATION_GENES = ["PINK1", "PRKN", "SNCA", "MAPT", "APOE"]
INFLAMMATORY_GENES = ["TREM2", "TYROBP", "GPNMB", "HLA-DRA", "AIF1", "IL1B", "TNF", "NFKB1", "B2M"]
MITOCHONDRIAL_GENES = ["PINK1", "PRKN", "LRRK2"]

COMPOSITE_INDICES = {
    "MAI": MICROGLIAL_GENES,
    "ASI": ASTROCYTE_GENES,
    "NVI": NEURODEGENERATION_GENES,
}

CELLTYPE_SIGNATURES = {
    "neurodegeneration_signature": NEURODEGENERATION_GENES,
    "inflammatory_signature": INFLAMMATORY_GENES,
    "mitochondrial_dysfunction_signature": MITOCHONDRIAL_GENES,
}


def configure_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
        ],
    )


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def clean_label(value: str | None) -> str:
    if value is None:
        return "missing"
    value = value.strip()
    return value if value else "missing"


def load_metadata_vectors(
    metadata_path: Path,
) -> tuple[
    dict[str, list[str]],
    dict[tuple[str, str], int],
    dict[tuple[str, str, str], int],
    dict[tuple[str, str], set[str]],
]:
    vectors: dict[str, list[str]] = {field: [] for field in METADATA_FIELDS}
    group_cell_counts: dict[tuple[str, str], int] = defaultdict(int)
    donor_group_cell_counts: dict[tuple[str, str, str], int] = defaultdict(int)
    donors_by_group: dict[tuple[str, str], set[str]] = defaultdict(set)

    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            donor = clean_label(row.get(DONOR_FIELD) or f"missing_donor_{row_number}")
            vectors[DONOR_FIELD].append(donor)
            for field in GROUP_VARIABLES:
                label = clean_label(row.get(field))
                vectors[field].append(label)
                group_cell_counts[(field, label)] += 1
                donor_group_cell_counts[(field, label, donor)] += 1
                donors_by_group[(field, label)].add(donor)

    logging.info("Loaded decoded metadata rows: %d", len(vectors[DONOR_FIELD]))
    return vectors, group_cell_counts, donor_group_cell_counts, donors_by_group


def empty_stat() -> dict[str, float]:
    return {"sum": 0.0, "nonzero": 0.0}


def add_stat(stats: dict[tuple[str, ...], dict[str, float]], key: tuple[str, ...], value: float) -> None:
    stats[key]["sum"] += value
    stats[key]["nonzero"] += 1.0


def build_gene_to_signature_lookup() -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = defaultdict(list)
    for index_id, genes in COMPOSITE_INDICES.items():
        for gene in genes:
            lookup[gene].append(index_id)
    for signature_id, genes in CELLTYPE_SIGNATURES.items():
        for gene in genes:
            lookup[gene].append(signature_id)
    return lookup


def stream_sparse_expression(
    expression_path: Path,
    metadata_vectors: dict[str, list[str]],
    max_row_chunk: int,
) -> tuple[
    set[str],
    dict[tuple[str, str, str, str], dict[str, float]],
    dict[tuple[str, str, str], dict[str, float]],
    dict[tuple[str, str, str, str], dict[str, float]],
    dict[tuple[str, str, str], dict[str, float]],
]:
    if max_row_chunk > MAX_ALLOWED_ROW_CHUNK:
        raise ValueError(f"max_row_chunk may not exceed {MAX_ALLOWED_ROW_CHUNK}")

    gene_to_signatures = build_gene_to_signature_lookup()
    genes_seen: set[str] = set()
    gene_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]] = defaultdict(empty_stat)
    gene_group_stats: dict[tuple[str, str, str], dict[str, float]] = defaultdict(empty_stat)
    signature_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]] = defaultdict(empty_stat)
    signature_group_stats: dict[tuple[str, str, str], dict[str, float]] = defaultdict(empty_stat)
    processed = 0

    with open_text(expression_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            processed += 1
            if processed % max_row_chunk == 0:
                logging.info("Processed sparse expression rows: %d", processed)

            row_index = int(row["row_index"])
            if row_index >= len(metadata_vectors[DONOR_FIELD]):
                logging.warning("Skipping out-of-range row_index=%d", row_index)
                continue

            gene = row.get("gene_symbol", "")
            value = float(row["expression_value"])
            donor = metadata_vectors[DONOR_FIELD][row_index]
            genes_seen.add(gene)

            for variable in GROUP_VARIABLES:
                label = metadata_vectors[variable][row_index]
                add_stat(gene_donor_group_stats, (gene, variable, label, donor), value)
                add_stat(gene_group_stats, (gene, variable, label), value)
                for signature_id in gene_to_signatures.get(gene, []):
                    add_stat(signature_donor_group_stats, (signature_id, variable, label, donor), value)
                    add_stat(signature_group_stats, (signature_id, variable, label), value)

    logging.info("Total sparse expression rows processed: %d", processed)
    logging.info("Genes observed in sparse table: %d", len(genes_seen))
    return (
        genes_seen,
        gene_donor_group_stats,
        gene_group_stats,
        signature_donor_group_stats,
        signature_group_stats,
    )


def to_float_text(value: float) -> str:
    if math.isnan(value):
        return "nan"
    return f"{value:.8g}"


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    avg = mean(values)
    variance = sum((value - avg) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(max(variance, 0.0))


def mean_ci95(values: list[float]) -> tuple[float, float, float]:
    if not values:
        return float("nan"), float("nan"), float("nan")
    avg = mean(values)
    if len(values) < 2:
        return avg, avg, avg
    margin = 1.96 * standard_deviation(values) / math.sqrt(len(values))
    return avg, avg - margin, avg + margin


def label_rank(variable: str, label: str) -> float | None:
    lowered = label.lower()
    if label == "missing" or "missing" in lowered or "unknown" in lowered:
        return None
    if variable == "Braak":
        roman = {"0": 0, "i": 1, "ii": 2, "iii": 3, "iv": 4, "v": 5, "vi": 6}
        tokens = lowered.replace("/", " ").replace("-", " ").split()
        for token in tokens:
            token = token.strip("()[],:;")
            if token in roman:
                return float(roman[token])
            if token.isdigit():
                return float(token)
    if variable == "CERAD score":
        mapping = {"absent": 0, "none": 0, "sparse": 1, "moderate": 2, "frequent": 3}
        for key, rank in mapping.items():
            if key in lowered:
                return float(rank)
    if variable == "Overall AD neuropathological Change":
        mapping = {
            "not": 0,
            "none": 0,
            "low": 1,
            "intermediate": 2,
            "moderate": 2,
            "high": 3,
        }
        for key, rank in mapping.items():
            if key in lowered:
                return float(rank)
    return None


def sort_labels(variable: str, labels: list[str]) -> list[str]:
    return sorted(
        labels,
        key=lambda label: (
            label_rank(variable, label) is None,
            label_rank(variable, label) if label_rank(variable, label) is not None else label.lower(),
        ),
    )


def present_gene_count(signature_genes: list[str], genes_seen: set[str]) -> int:
    return max(1, len([gene for gene in signature_genes if gene in genes_seen]))


def group_mean_strings(
    variable: str,
    labels: list[str],
    donor_values_by_label: dict[str, list[float]],
    group_cell_counts: dict[tuple[str, str], int],
    group_nonzero_counts: dict[str, float],
) -> tuple[str, str, str, str]:
    mean_parts: list[str] = []
    detection_parts: list[str] = []
    ci_parts: list[str] = []
    donor_parts: list[str] = []
    for label in labels:
        values = donor_values_by_label.get(label, [])
        avg, low, high = mean_ci95(values)
        cell_count = group_cell_counts.get((variable, label), 0)
        detection = group_nonzero_counts.get(label, 0.0) / cell_count if cell_count else 0.0
        mean_parts.append(f"{label}={to_float_text(avg)}")
        detection_parts.append(f"{label}={to_float_text(detection)}")
        ci_parts.append(f"{label}={to_float_text(low)}..{to_float_text(high)}")
        donor_parts.append(f"{label}={len(values)}")
    return "; ".join(mean_parts), "; ".join(detection_parts), "; ".join(ci_parts), "; ".join(donor_parts)


def donor_values_for_gene(
    gene: str,
    variable: str,
    label: str,
    donors_by_group: dict[tuple[str, str], set[str]],
    donor_group_cell_counts: dict[tuple[str, str, str], int],
    gene_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]],
) -> list[float]:
    values: list[float] = []
    for donor in sorted(donors_by_group.get((variable, label), set())):
        cells = donor_group_cell_counts[(variable, label, donor)]
        stat = gene_donor_group_stats.get((gene, variable, label, donor), empty_stat())
        values.append(stat["sum"] / cells if cells else 0.0)
    return values


def donor_values_for_signature(
    signature_id: str,
    signature_genes: list[str],
    genes_seen: set[str],
    variable: str,
    label: str,
    donors_by_group: dict[tuple[str, str], set[str]],
    donor_group_cell_counts: dict[tuple[str, str, str], int],
    signature_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]],
) -> list[float]:
    values: list[float] = []
    gene_count = present_gene_count(signature_genes, genes_seen)
    for donor in sorted(donors_by_group.get((variable, label), set())):
        cells = donor_group_cell_counts[(variable, label, donor)]
        stat = signature_donor_group_stats.get((signature_id, variable, label, donor), empty_stat())
        denominator = cells * gene_count
        values.append(stat["sum"] / denominator if denominator else 0.0)
    return values


def rank_based_test(
    variable: str,
    labels: list[str],
    donor_values_by_label: dict[str, list[float]],
) -> tuple[str, float, float]:
    samples = [donor_values_by_label[label] for label in labels if donor_values_by_label.get(label)]
    if len(samples) < 2 or scipy_stats is None:
        return "not_tested", float("nan"), float("nan")

    ranked_labels = [label for label in labels if label_rank(variable, label) is not None]
    if len(ranked_labels) >= 2:
        x_values: list[float] = []
        y_values: list[float] = []
        for label in ranked_labels:
            rank = label_rank(variable, label)
            if rank is None:
                continue
            for value in donor_values_by_label.get(label, []):
                x_values.append(rank)
                y_values.append(value)
        if len(set(x_values)) >= 2 and len(set(y_values)) >= 2:
            statistic, p_value = scipy_stats.spearmanr(x_values, y_values)
            return "spearman_rank_trend", float(statistic), float(p_value)

    non_empty_samples = [sample for sample in samples if sample]
    if len(non_empty_samples) >= 2:
        try:
            statistic, p_value = scipy_stats.kruskal(*non_empty_samples)
        except ValueError as exc:
            logging.warning("Skipping rank association for %s: %s", variable, exc)
            return "not_tested_constant_values", float("nan"), float("nan")
        return "kruskal_rank_association", float(statistic), float(p_value)
    return "not_tested", float("nan"), float("nan")


def effect_size_and_direction(
    variable: str,
    labels: list[str],
    donor_values_by_label: dict[str, list[float]],
) -> tuple[float, str]:
    means = [(label, mean(donor_values_by_label.get(label, []))) for label in labels]
    means = [(label, value) for label, value in means if not math.isnan(value)]
    if not means:
        return float("nan"), "not_available"
    low_label, low_value = min(means, key=lambda item: item[1])
    high_label, high_value = max(means, key=lambda item: item[1])
    effect = high_value - low_value
    if all(label_rank(variable, label) is not None for label, _ in means):
        ordered = sort_labels(variable, [label for label, _ in means])
        first = mean(donor_values_by_label[ordered[0]])
        last = mean(donor_values_by_label[ordered[-1]])
        direction = "increases_with_rank" if last >= first else "decreases_with_rank"
    else:
        direction = f"{high_label}_higher_than_{low_label}"
    return effect, direction


def benjamini_hochberg(rows: list[dict[str, str]], p_key: str = "p_value", out_key: str = "fdr_p_value") -> None:
    valid: list[tuple[int, float]] = []
    for index, row in enumerate(rows):
        try:
            p_value = float(row[p_key])
        except (KeyError, TypeError, ValueError):
            row[out_key] = "nan"
            continue
        if math.isnan(p_value):
            row[out_key] = "nan"
            continue
        valid.append((index, p_value))

    valid.sort(key=lambda item: item[1])
    total = len(valid)
    running = 1.0
    adjusted: dict[int, float] = {}
    for rank_from_end, (index, p_value) in enumerate(reversed(valid), start=1):
        rank = total - rank_from_end + 1
        running = min(running, p_value * total / rank)
        adjusted[index] = min(running, 1.0)
    for index, value in adjusted.items():
        rows[index][out_key] = to_float_text(value)


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote %s rows: %d", path, len(rows))


def build_gene_statistics_rows(
    genes_seen: set[str],
    group_cell_counts: dict[tuple[str, str], int],
    donor_group_cell_counts: dict[tuple[str, str, str], int],
    donors_by_group: dict[tuple[str, str], set[str]],
    gene_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]],
    gene_group_stats: dict[tuple[str, str, str], dict[str, float]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for gene in sorted(genes_seen):
        for variable in STAT_VARIABLES:
            labels = sort_labels(variable, [label for field, label in group_cell_counts if field == variable])
            donor_values_by_label = {
                label: donor_values_for_gene(
                    gene,
                    variable,
                    label,
                    donors_by_group,
                    donor_group_cell_counts,
                    gene_donor_group_stats,
                )
                for label in labels
            }
            group_nonzero = {
                label: gene_group_stats.get((gene, variable, label), empty_stat())["nonzero"]
                for label in labels
            }
            means, detection, intervals, donor_counts = group_mean_strings(
                variable,
                labels,
                donor_values_by_label,
                group_cell_counts,
                group_nonzero,
            )
            method, statistic, p_value = rank_based_test(variable, labels, donor_values_by_label)
            effect, direction = effect_size_and_direction(variable, labels, donor_values_by_label)
            rows.append(
                {
                    "gene_symbol": gene,
                    "test_variable": variable,
                    "test_method": method,
                    "rank_statistic": to_float_text(statistic),
                    "p_value": to_float_text(p_value),
                    "fdr_p_value": "nan",
                    "effect_size_max_minus_min": to_float_text(effect),
                    "direction": direction,
                    "donor_count": str(len({donor for label in labels for donor in donors_by_group.get((variable, label), set())})),
                    "cell_count": str(sum(group_cell_counts.get((variable, label), 0) for label in labels)),
                    "group_mean_expression": means,
                    "group_detection_rate": detection,
                    "group_mean_ci95": intervals,
                    "group_donor_count": donor_counts,
                    "notes": "donor-aware rank statistics from sparse target-gene table",
                }
            )
    benjamini_hochberg(rows)
    return rows


def build_celltype_vulnerability_rows(
    genes_seen: set[str],
    group_cell_counts: dict[tuple[str, str], int],
    donor_group_cell_counts: dict[tuple[str, str, str], int],
    donors_by_group: dict[tuple[str, str], set[str]],
    signature_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]],
    signature_group_stats: dict[tuple[str, str, str], dict[str, float]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    celltypes = sort_labels(CELLTYPE_FIELD, [label for field, label in group_cell_counts if field == CELLTYPE_FIELD])
    for signature_id, signature_genes in CELLTYPE_SIGNATURES.items():
        gene_count = present_gene_count(signature_genes, genes_seen)
        ranked: list[dict[str, str]] = []
        for celltype in celltypes:
            cells = group_cell_counts.get((CELLTYPE_FIELD, celltype), 0)
            stat = signature_group_stats.get((signature_id, CELLTYPE_FIELD, celltype), empty_stat())
            denominator = cells * gene_count
            values = donor_values_for_signature(
                signature_id,
                signature_genes,
                genes_seen,
                CELLTYPE_FIELD,
                celltype,
                donors_by_group,
                donor_group_cell_counts,
                signature_donor_group_stats,
            )
            avg, low, high = mean_ci95(values)
            row = {
                "signature_id": signature_id,
                "cell_subclass": celltype,
                "rank": "0",
                "cell_count": str(cells),
                "donor_count": str(len(values)),
                "genes_present": str(gene_count),
                "mean_index": to_float_text(stat["sum"] / denominator if denominator else 0.0),
                "donor_mean_index": to_float_text(avg),
                "ci95_low": to_float_text(low),
                "ci95_high": to_float_text(high),
                "detection_rate": to_float_text(stat["nonzero"] / denominator if denominator else 0.0),
                "notes": "ranked within signature by mean index",
            }
            ranked.append(row)
        ranked.sort(key=lambda row: float(row["mean_index"]), reverse=True)
        for rank, row in enumerate(ranked, start=1):
            row["rank"] = str(rank)
            rows.append(row)
    return rows


def build_index_group_rows(
    table_kind: str,
    variables: list[str],
    genes_seen: set[str],
    group_cell_counts: dict[tuple[str, str], int],
    donor_group_cell_counts: dict[tuple[str, str, str], int],
    donors_by_group: dict[tuple[str, str], set[str]],
    signature_donor_group_stats: dict[tuple[str, str, str, str], dict[str, float]],
    signature_group_stats: dict[tuple[str, str, str], dict[str, float]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    test_rows: list[dict[str, str]] = []
    for index_id, index_genes in COMPOSITE_INDICES.items():
        gene_count = present_gene_count(index_genes, genes_seen)
        for variable in variables:
            labels = sort_labels(variable, [label for field, label in group_cell_counts if field == variable])
            donor_values_by_label = {
                label: donor_values_for_signature(
                    index_id,
                    index_genes,
                    genes_seen,
                    variable,
                    label,
                    donors_by_group,
                    donor_group_cell_counts,
                    signature_donor_group_stats,
                )
                for label in labels
            }
            method, statistic, p_value = rank_based_test(variable, labels, donor_values_by_label)
            effect, direction = effect_size_and_direction(variable, labels, donor_values_by_label)
            test_row = {
                "index_id": index_id,
                "grouping_variable": variable,
                "test_method": method,
                "rank_statistic": to_float_text(statistic),
                "p_value": to_float_text(p_value),
                "fdr_p_value": "nan",
            }
            test_rows.append(test_row)

            for label in labels:
                cells = group_cell_counts.get((variable, label), 0)
                stat = signature_group_stats.get((index_id, variable, label), empty_stat())
                denominator = cells * gene_count
                values = donor_values_by_label.get(label, [])
                avg, low, high = mean_ci95(values)
                row: dict[str, str] = {
                    "index_id": index_id,
                    "grouping_variable": variable,
                    "group_label": label,
                    "donor_count": str(len(values)),
                    "cell_count": str(cells),
                    "genes_present": str(gene_count),
                    "mean_index": to_float_text(stat["sum"] / denominator if denominator else 0.0),
                    "donor_mean_index": to_float_text(avg),
                    "ci95_low": to_float_text(low),
                    "ci95_high": to_float_text(high),
                    "detection_rate": to_float_text(stat["nonzero"] / denominator if denominator else 0.0),
                    "test_method": method,
                    "rank_statistic": to_float_text(statistic),
                    "p_value": to_float_text(p_value),
                    "fdr_p_value": "nan",
                    "effect_size_max_minus_min": to_float_text(effect),
                    "direction": direction,
                    "notes": table_kind,
                }
                if variable == APOE_VARIABLE:
                    row["apoe_genotype"] = label
                rows.append(row)

    benjamini_hochberg(test_rows)
    fdr_lookup = {
        (row["index_id"], row["grouping_variable"]): row["fdr_p_value"]
        for row in test_rows
    }
    for row in rows:
        row["fdr_p_value"] = fdr_lookup.get((row["index_id"], row["grouping_variable"]), "nan")
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute Phase 4 statistical biology tables.")
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv"),
    )
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/19_compute_phase4_statistics.log"),
    )
    parser.add_argument("--max-row-chunk", type=int, default=DEFAULT_MAX_ROW_CHUNK)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting Phase 4 sparse statistical biology analysis.")
    logging.info("Expression input: %s", args.expression)
    logging.info("Metadata input: %s", args.metadata)

    (
        metadata_vectors,
        group_cell_counts,
        donor_group_cell_counts,
        donors_by_group,
    ) = load_metadata_vectors(args.metadata)
    (
        genes_seen,
        gene_donor_group_stats,
        gene_group_stats,
        signature_donor_group_stats,
        signature_group_stats,
    ) = stream_sparse_expression(args.expression, metadata_vectors, args.max_row_chunk)

    gene_rows = build_gene_statistics_rows(
        genes_seen,
        group_cell_counts,
        donor_group_cell_counts,
        donors_by_group,
        gene_donor_group_stats,
        gene_group_stats,
    )
    write_tsv(
        args.tables_dir / "phase4_gene_statistics.tsv",
        gene_rows,
        [
            "gene_symbol",
            "test_variable",
            "test_method",
            "rank_statistic",
            "p_value",
            "fdr_p_value",
            "effect_size_max_minus_min",
            "direction",
            "donor_count",
            "cell_count",
            "group_mean_expression",
            "group_detection_rate",
            "group_mean_ci95",
            "group_donor_count",
            "notes",
        ],
    )

    vulnerability_rows = build_celltype_vulnerability_rows(
        genes_seen,
        group_cell_counts,
        donor_group_cell_counts,
        donors_by_group,
        signature_donor_group_stats,
        signature_group_stats,
    )
    write_tsv(
        args.tables_dir / "phase4_celltype_vulnerability.tsv",
        vulnerability_rows,
        [
            "signature_id",
            "cell_subclass",
            "rank",
            "cell_count",
            "donor_count",
            "genes_present",
            "mean_index",
            "donor_mean_index",
            "ci95_low",
            "ci95_high",
            "detection_rate",
            "notes",
        ],
    )

    apoe_rows = build_index_group_rows(
        "APOE genotype comparison of composite indices",
        [APOE_VARIABLE],
        genes_seen,
        group_cell_counts,
        donor_group_cell_counts,
        donors_by_group,
        signature_donor_group_stats,
        signature_group_stats,
    )
    write_tsv(
        args.tables_dir / "phase4_apoe_analysis.tsv",
        apoe_rows,
        [
            "index_id",
            "grouping_variable",
            "apoe_genotype",
            "group_label",
            "donor_count",
            "cell_count",
            "genes_present",
            "mean_index",
            "donor_mean_index",
            "ci95_low",
            "ci95_high",
            "detection_rate",
            "test_method",
            "rank_statistic",
            "p_value",
            "fdr_p_value",
            "effect_size_max_minus_min",
            "direction",
            "notes",
        ],
    )

    mixed_rows = build_index_group_rows(
        "mixed pathology comparison of composite indices",
        MIXED_PATHOLOGY_VARIABLES,
        genes_seen,
        group_cell_counts,
        donor_group_cell_counts,
        donors_by_group,
        signature_donor_group_stats,
        signature_group_stats,
    )
    write_tsv(
        args.tables_dir / "phase4_mixed_pathology.tsv",
        mixed_rows,
        [
            "index_id",
            "grouping_variable",
            "group_label",
            "donor_count",
            "cell_count",
            "genes_present",
            "mean_index",
            "donor_mean_index",
            "ci95_low",
            "ci95_high",
            "detection_rate",
            "test_method",
            "rank_statistic",
            "p_value",
            "fdr_p_value",
            "effect_size_max_minus_min",
            "direction",
            "notes",
        ],
    )

    composite_rows = build_index_group_rows(
        "disease-associated composite index comparison",
        [*STAT_VARIABLES, APOE_VARIABLE, *MIXED_PATHOLOGY_VARIABLES],
        genes_seen,
        group_cell_counts,
        donor_group_cell_counts,
        donors_by_group,
        signature_donor_group_stats,
        signature_group_stats,
    )
    write_tsv(
        args.tables_dir / "phase4_composite_indices.tsv",
        composite_rows,
        [
            "index_id",
            "grouping_variable",
            "group_label",
            "donor_count",
            "cell_count",
            "genes_present",
            "mean_index",
            "donor_mean_index",
            "ci95_low",
            "ci95_high",
            "detection_rate",
            "test_method",
            "rank_statistic",
            "p_value",
            "fdr_p_value",
            "effect_size_max_minus_min",
            "direction",
            "notes",
        ],
    )

    logging.info("Phase 4 statistics complete.")
    logging.info("No full matrix, dense matrix, or workflow-engine analysis was executed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
