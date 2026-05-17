#!/usr/bin/env python3
"""SEA-AD metadata-only extraction.

This script reads only H5AD metadata groups (`obs` and `var`) via h5py.
It must never touch the expression matrix root (`X`) or any `X/...` arrays.
It has no Scanpy dependency, creates no in-memory analysis object, and does not run
normalization, PCA, UMAP, clustering, or model training.
"""

from __future__ import annotations

import argparse
import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import h5py


FORBIDDEN_ROOT_KEY = "X"

OBS_FIELDS = [
    "Donor ID",
    "Subclass",
    "Class",
    "Supertype",
    "Sex",
    "Gender",
    "Age at Death",
    "Braak",
    "CERAD score",
    "Overall AD neuropathological Change",
    "Cognitive Status",
    "PMI",
    "RIN",
    "Fraction mitochondrial UMIs",
    "Genes detected",
    "Number of UMIs",
    "Brain Region",
    "APOE Genotype",
    "Highest Lewy Body Disease",
    "LATE",
    "Thal",
    "Overall CAA Score",
    "Used in analysis",
    "Neurotypical reference",
]

OBS_OUTPUT = "sea_ad_obs_metadata_minimal.tsv"
VAR_OUTPUT = "sea_ad_var_genes.tsv"
SUMMARY_OUTPUT = "sea_ad_metadata_summary.tsv"
TABLE1_OUTPUT = "table1_sea_ad_cohort_cell_summary.tsv"


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


def forbid_expression_path(path: str) -> None:
    root = path.strip("/").split("/", 1)[0]
    if root == FORBIDDEN_ROOT_KEY:
        raise RuntimeError(
            "Forbidden expression-matrix access requested. "
            "Phase 2C may read obs/var metadata only."
        )


class SafeH5Reader:
    """Small path guard around h5py that refuses expression-matrix access."""

    def __init__(self, handle: h5py.File):
        self._handle = handle

    def get(self, path: str) -> Any:
        forbid_expression_path(path)
        node: Any = self._handle
        for part in [item for item in path.strip("/").split("/") if item]:
            node = node[part]
        return node

    def exists(self, path: str) -> bool:
        forbid_expression_path(path)
        node: Any = self._handle
        for part in [item for item in path.strip("/").split("/") if item]:
            if part not in node:
                return False
            node = node[part]
        return True


def decode_scalar(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item"):
        try:
            value = value.item()
        except ValueError:
            pass
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if value is None:
        return ""
    text = str(value)
    return "" if text == "nan" else text


def read_dataset_slice(dataset: h5py.Dataset, start: int | None, stop: int | None) -> list[str]:
    raw = dataset[()] if start is None or stop is None else dataset[start:stop]
    if getattr(raw, "shape", ()) == ():
        return [decode_scalar(raw)]
    return [decode_scalar(value) for value in raw]


def read_encoded_values(node: Any, start: int | None = None, stop: int | None = None) -> list[str]:
    """Read h5ad-encoded metadata values from a Dataset or categorical Group."""
    if isinstance(node, h5py.Dataset):
        return read_dataset_slice(node, start, stop)

    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        codes = node["codes"][()] if start is None or stop is None else node["codes"][start:stop]
        categories = read_encoded_values(node["categories"], None, None)
        values: list[str] = []
        for code in codes:
            code_int = int(code)
            values.append("" if code_int < 0 else categories[code_int])
        return values

    if isinstance(node, h5py.Group) and "values" in node:
        return read_encoded_values(node["values"], start, stop)

    raise TypeError(f"Unsupported metadata encoding at node {getattr(node, 'name', '<unknown>')}")


def obs_index_name(obs_group: h5py.Group) -> str:
    raw = obs_group.attrs.get("_index", "_index")
    return decode_scalar(raw)


def node_length(node: Any) -> int:
    if isinstance(node, h5py.Dataset):
        return int(node.shape[0])
    if isinstance(node, h5py.Group) and "codes" in node:
        return int(node["codes"].shape[0])
    if isinstance(node, h5py.Group) and "values" in node:
        return node_length(node["values"])
    raise TypeError(f"Cannot determine length for node {getattr(node, 'name', '<unknown>')}")


def selected_obs_fields(reader: SafeH5Reader) -> list[str]:
    return [field for field in OBS_FIELDS if reader.exists(f"obs/{field}")]


def write_obs_metadata(
    reader: SafeH5Reader,
    out_path: Path,
    max_obs_preview: int | None,
    chunk_size: int,
) -> dict[str, Any]:
    obs_group = reader.get("obs")
    index_name = obs_index_name(obs_group)
    index_node = reader.get(f"obs/{index_name}")
    available_fields = selected_obs_fields(reader)
    total_rows = node_length(index_node)
    rows_to_write = min(total_rows, max_obs_preview) if max_obs_preview is not None else total_rows

    stats: dict[str, Any] = {
        "cell_count": 0,
        "donors": set(),
        "subclasses": set(),
        "subclass_counts": Counter(),
        "ad_pathology_counts": Counter(),
        "cognitive_status_counts": Counter(),
        "class_counts": Counter(),
        "brain_region_counts": Counter(),
        "donors_by_category": defaultdict(set),
        "available_fields": available_fields,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["obs_index", *available_fields], delimiter="\t")
        writer.writeheader()

        for start in range(0, rows_to_write, chunk_size):
            stop = min(start + chunk_size, rows_to_write)
            chunk: dict[str, list[str]] = {
                "obs_index": read_encoded_values(index_node, start, stop)
            }
            for field in available_fields:
                chunk[field] = read_encoded_values(reader.get(f"obs/{field}"), start, stop)

            for offset in range(stop - start):
                row = {field: chunk[field][offset] for field in chunk}
                writer.writerow(row)
                update_stats(stats, row)

            logging.info("Wrote obs metadata rows %d-%d of %d", start + 1, stop, rows_to_write)

    stats["cell_count"] = rows_to_write
    stats["total_obs_rows_in_file"] = total_rows
    return stats


def update_stats(stats: dict[str, Any], row: dict[str, str]) -> None:
    donor_id = row.get("Donor ID", "")
    subclass = row.get("Subclass", "")
    ad_pathology = row.get("Overall AD neuropathological Change", "")
    cognitive_status = row.get("Cognitive Status", "")
    class_label = row.get("Class", "")
    brain_region = row.get("Brain Region", "")

    if donor_id:
        stats["donors"].add(donor_id)
    if subclass:
        stats["subclasses"].add(subclass)
        stats["subclass_counts"][subclass] += 1
        if donor_id:
            stats["donors_by_category"][("subclass", subclass)].add(donor_id)
    if ad_pathology:
        stats["ad_pathology_counts"][ad_pathology] += 1
        if donor_id:
            stats["donors_by_category"][("ad_pathology", ad_pathology)].add(donor_id)
    if cognitive_status:
        stats["cognitive_status_counts"][cognitive_status] += 1
        if donor_id:
            stats["donors_by_category"][("cognitive_status", cognitive_status)].add(donor_id)
    if class_label:
        stats["class_counts"][class_label] += 1
        if donor_id:
            stats["donors_by_category"][("class", class_label)].add(donor_id)
    if brain_region:
        stats["brain_region_counts"][brain_region] += 1
        if donor_id:
            stats["donors_by_category"][("brain_region", brain_region)].add(donor_id)


def write_var_genes(reader: SafeH5Reader, out_path: Path) -> int:
    var_index_path = "var/_index" if reader.exists("var/_index") else None
    gene_ids_path = "var/gene_ids" if reader.exists("var/gene_ids") else None
    if var_index_path is None and gene_ids_path is None:
        raise RuntimeError("No supported var gene index fields found.")

    index_values = read_encoded_values(reader.get(var_index_path), None, None) if var_index_path else []
    gene_ids = read_encoded_values(reader.get(gene_ids_path), None, None) if gene_ids_path else []
    row_count = max(len(index_values), len(gene_ids))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_symbol", "gene_id"], delimiter="\t")
        writer.writeheader()
        for index in range(row_count):
            writer.writerow(
                {
                    "gene_symbol": index_values[index] if index < len(index_values) else "",
                    "gene_id": gene_ids[index] if index < len(gene_ids) else "",
                }
            )
    return row_count


def write_summary_tables(stats: dict[str, Any], gene_count: int, tables_dir: Path) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)

    summary_path = tables_dir / SUMMARY_OUTPUT
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["cell_count", stats["cell_count"]])
        writer.writerow(["total_obs_rows_in_file", stats["total_obs_rows_in_file"]])
        writer.writerow(["donor_count", len(stats["donors"])])
        writer.writerow(["subclass_count", len(stats["subclasses"])])
        writer.writerow(["gene_count", gene_count])
        writer.writerow(["obs_fields_extracted", ",".join(stats["available_fields"])])
        for label, count in sorted(stats["ad_pathology_counts"].items()):
            writer.writerow([f"overall_ad_neuropathological_change:{label}", count])
        for label, count in sorted(stats["cognitive_status_counts"].items()):
            writer.writerow([f"cognitive_status:{label}", count])

    table1_path = tables_dir / TABLE1_OUTPUT
    with table1_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["summary_group", "label", "cell_count", "donor_count"],
            delimiter="\t",
        )
        writer.writeheader()
        for group_name, counter in [
            ("Class", stats["class_counts"]),
            ("Subclass", stats["subclass_counts"]),
            ("Brain Region", stats["brain_region_counts"]),
            ("Overall AD neuropathological Change", stats["ad_pathology_counts"]),
            ("Cognitive Status", stats["cognitive_status_counts"]),
        ]:
            category_key = group_name.lower().replace(" ", "_")
            if group_name == "Overall AD neuropathological Change":
                category_key = "ad_pathology"
            for label, cell_count in sorted(counter.items()):
                writer.writerow(
                    {
                        "summary_group": group_name,
                        "label": label,
                        "cell_count": cell_count,
                        "donor_count": len(stats["donors_by_category"].get((category_key, label), set())),
                    }
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract SEA-AD obs/var metadata only.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("data/interim/sea_ad"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/11_extract_sea_ad_metadata_only.log"),
    )
    parser.add_argument(
        "--max-obs-preview",
        type=int,
        default=None,
        help="Optional row limit for previewing obs metadata. Default extracts all selected obs rows.",
    )
    parser.add_argument("--chunk-size", type=int, default=100_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting SEA-AD metadata-only extraction.")
    logging.info("Input h5ad: %s", args.input)
    logging.info("Output metadata directory: %s", args.outdir)
    logging.info("Output tables directory: %s", args.tables_dir)
    logging.info("Max obs preview: %s", args.max_obs_preview)
    logging.info("Expression-matrix root '%s' is forbidden.", FORBIDDEN_ROOT_KEY)

    if not args.input.exists():
        logging.error("Input file not found: %s", args.input)
        return 2

    with h5py.File(args.input, "r") as raw_handle:
        reader = SafeH5Reader(raw_handle)
        obs_path = args.outdir / OBS_OUTPUT
        var_path = args.outdir / VAR_OUTPUT
        stats = write_obs_metadata(reader, obs_path, args.max_obs_preview, args.chunk_size)
        gene_count = write_var_genes(reader, var_path)
        write_summary_tables(stats, gene_count, args.tables_dir)

    logging.info("Donor count: %d", len(stats["donors"]))
    logging.info("Cell count: %d", stats["cell_count"])
    logging.info("Subclass count: %d", len(stats["subclasses"]))
    logging.info("Disease category counts: %s", dict(stats["ad_pathology_counts"]))
    logging.info("Cognitive status counts: %s", dict(stats["cognitive_status_counts"]))
    logging.info("Wrote %s", args.outdir / OBS_OUTPUT)
    logging.info("Wrote %s", args.outdir / VAR_OUTPUT)
    logging.info("Wrote %s", args.tables_dir / SUMMARY_OUTPUT)
    logging.info("Wrote %s", args.tables_dir / TABLE1_OUTPUT)
    logging.info("No X matrix arrays were read.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
