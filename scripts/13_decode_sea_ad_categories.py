#!/usr/bin/env python3
"""Decode SEA-AD metadata categories and prepare publication Table 1.

This script reads labels from `obs/__categories` in the local H5AD file and
decodes the metadata-only TSV produced by Phase 2C. It never reads expression
matrix arrays and never runs analysis workflows.
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

DECODED_METADATA_OUTPUT = "sea_ad_obs_metadata_decoded.tsv"
CATEGORY_MAPPING_OUTPUT = "sea_ad_category_mapping.tsv"
UNSUPPORTED_CATEGORY_NODES_OUTPUT = "sea_ad_unsupported_category_nodes.tsv"
PUBLICATION_TABLE_OUTPUT = "table1_sea_ad_publication_ready.tsv"

TABLE1_FIELDS = [
    "Class",
    "Subclass",
    "Supertype",
    "Sex",
    "Gender",
    "Age at Death",
    "Braak",
    "CERAD score",
    "Overall AD neuropathological Change",
    "Cognitive Status",
    "Brain Region",
    "APOE Genotype",
    "Highest Lewy Body Disease",
    "LATE",
    "Thal",
    "Overall CAA Score",
    "Used in analysis",
    "Neurotypical reference",
]


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
            "Phase 2D may read obs metadata and obs categories only."
        )


class SafeH5Reader:
    """Path guard around h5py that refuses expression-matrix access."""

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


def read_label_values(node: Any) -> list[str]:
    if isinstance(node, h5py.Dataset):
        raw = node[()]
        if getattr(raw, "shape", ()) == ():
            return [decode_scalar(raw)]
        return [decode_scalar(value) for value in raw]
    if isinstance(node, h5py.Group) and "categories" in node:
        return read_label_values(node["categories"])
    if isinstance(node, h5py.Group) and "values" in node:
        return read_label_values(node["values"])
    if isinstance(node, h5py.Group):
        keys = ", ".join(sorted(node.keys())) or "no child keys"
        raise TypeError(f"Unsupported category group encoding; child keys: {keys}")
    raise TypeError(f"Unsupported category node type: {type(node).__name__}")


def node_type_name(node: Any) -> str:
    if isinstance(node, h5py.Dataset):
        return "Dataset"
    if isinstance(node, h5py.Group):
        return "Group"
    return type(node).__name__


def load_category_mappings(
    reader: SafeH5Reader,
) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    category_root = "obs/__categories"
    if not reader.exists(category_root):
        logging.warning("No obs/__categories group found; decoded output will match input labels.")
        return {}, []

    root = reader.get(category_root)
    mappings: dict[str, dict[str, str]] = {}
    unsupported_nodes: list[dict[str, str]] = []
    for column in sorted(root.keys()):
        node = root[column]
        try:
            labels = read_label_values(node)
        except (TypeError, KeyError, ValueError, OSError) as exc:
            reason = str(exc)
            unsupported_nodes.append(
                {
                    "category_column": column,
                    "hdf5_path": getattr(node, "name", f"{category_root}/{column}"),
                    "node_type": node_type_name(node),
                    "reason": reason,
                }
            )
            logging.warning(
                "Skipping unsupported category node %s at %s: %s",
                column,
                getattr(node, "name", f"{category_root}/{column}"),
                reason,
            )
            continue
        mappings[column] = {str(index): label for index, label in enumerate(labels)}
    return mappings, unsupported_nodes


def decode_value(value: str, mapping: dict[str, str] | None) -> str:
    if mapping is None:
        return value
    stripped = value.strip()
    if stripped in {"", "-1", "nan", "None"}:
        return ""
    if stripped in mapping:
        return mapping[stripped]
    try:
        numeric_code = str(int(float(stripped)))
    except ValueError:
        return value
    return mapping.get(numeric_code, value)


def write_category_mapping(mappings: dict[str, dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["column", "code", "label"], delimiter="\t")
        writer.writeheader()
        for column, mapping in sorted(mappings.items()):
            for code, label in sorted(mapping.items(), key=lambda item: int(item[0])):
                writer.writerow({"column": column, "code": code, "label": label})


def write_unsupported_category_nodes(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["category_column", "hdf5_path", "node_type", "reason"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def decode_metadata_table(
    metadata_path: Path,
    decoded_path: Path,
    mappings: dict[str, dict[str, str]],
) -> dict[str, Any]:
    stats: dict[str, Any] = {
        "row_count": 0,
        "donors": set(),
        "cell_counts": Counter(),
        "donor_sets": defaultdict(set),
        "decoded_columns": [],
    }

    decoded_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("r", encoding="utf-8", newline="") as source, decoded_path.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Metadata TSV has no header: {metadata_path}")
        fieldnames = reader.fieldnames
        stats["decoded_columns"] = [field for field in fieldnames if field in mappings]

        writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            decoded = {
                field: decode_value(value, mappings.get(field))
                for field, value in row.items()
            }
            writer.writerow(decoded)
            update_table1_stats(stats, decoded)

    return stats


def update_table1_stats(stats: dict[str, Any], row: dict[str, str]) -> None:
    stats["row_count"] += 1
    donor_id = row.get("Donor ID", "")
    if donor_id:
        stats["donors"].add(donor_id)

    for field in TABLE1_FIELDS:
        if field not in row:
            continue
        label = row.get(field, "") or "missing"
        stats["cell_counts"][(field, label)] += 1
        if donor_id:
            stats["donor_sets"][(field, label)].add(donor_id)


def percent(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "0.00"
    return f"{(100 * numerator / denominator):.2f}"


def write_publication_table(stats: dict[str, Any], out_path: Path) -> None:
    total_cells = int(stats["row_count"])
    total_donors = len(stats["donors"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "table_section",
                "variable",
                "label",
                "cell_count",
                "cell_percent",
                "donor_count",
                "donor_percent",
                "notes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "table_section": "cohort",
                "variable": "Total",
                "label": "All nuclei/cells",
                "cell_count": total_cells,
                "cell_percent": "100.00",
                "donor_count": total_donors,
                "donor_percent": "100.00",
                "notes": "Metadata-only count; expression matrix not read.",
            }
        )
        for field in TABLE1_FIELDS:
            values = [
                ((variable, label), count)
                for (variable, label), count in stats["cell_counts"].items()
                if variable == field
            ]
            for (_, label), cell_count in sorted(values, key=lambda item: (-item[1], item[0][1])):
                donor_count = len(stats["donor_sets"].get((field, label), set()))
                writer.writerow(
                    {
                        "table_section": "metadata",
                        "variable": field,
                        "label": label,
                        "cell_count": cell_count,
                        "cell_percent": percent(cell_count, total_cells),
                        "donor_count": donor_count,
                        "donor_percent": percent(donor_count, total_donors),
                        "notes": "Decoded from obs metadata/category labels.",
                    }
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Decode SEA-AD obs category labels.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("data/interim/sea_ad/sea_ad_obs_metadata_minimal.tsv"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("data/interim/sea_ad"))
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/13_decode_sea_ad_categories.log"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)

    logging.info("Starting SEA-AD metadata category decoding.")
    logging.info("Input H5AD: %s", args.input)
    logging.info("Metadata TSV: %s", args.metadata)
    logging.info("Output metadata directory: %s", args.outdir)
    logging.info("Output tables directory: %s", args.tables_dir)

    if not args.input.exists():
        logging.error("Input H5AD file not found: %s", args.input)
        return 2
    if not args.metadata.exists():
        logging.error("Metadata TSV not found: %s", args.metadata)
        return 2

    with h5py.File(args.input, "r") as raw_handle:
        reader = SafeH5Reader(raw_handle)
        mappings, unsupported_nodes = load_category_mappings(reader)

    decoded_path = args.outdir / DECODED_METADATA_OUTPUT
    mapping_path = args.tables_dir / CATEGORY_MAPPING_OUTPUT
    unsupported_path = args.tables_dir / UNSUPPORTED_CATEGORY_NODES_OUTPUT
    table1_path = args.tables_dir / PUBLICATION_TABLE_OUTPUT

    write_category_mapping(mappings, mapping_path)
    write_unsupported_category_nodes(unsupported_nodes, unsupported_path)
    stats = decode_metadata_table(args.metadata, decoded_path, mappings)
    write_publication_table(stats, table1_path)

    logging.info("Decoded columns: %s", ", ".join(stats["decoded_columns"]))
    logging.info("Category mapping columns: %d", len(mappings))
    logging.info("Unsupported category nodes: %d", len(unsupported_nodes))
    logging.info("Rows decoded: %d", stats["row_count"])
    logging.info("Donor count: %d", len(stats["donors"]))
    logging.info("Wrote %s", decoded_path)
    logging.info("Wrote %s", mapping_path)
    logging.info("Wrote %s", unsupported_path)
    logging.info("Wrote %s", table1_path)
    logging.info("No expression matrix arrays were read.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
