#!/usr/bin/env python3
"""Metadata-only inspection for external H5AD cohorts.

This script is designed for Mathys 2019 onboarding. It reads metadata groups
only, refuses expression-matrix paths, and writes lightweight summaries.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any

import h5py


FORBIDDEN_ROOT_KEY = "X"
DEFAULT_OBS_KEYS = ["Donor ID", "donor", "individual", "subject", "diagnosis", "cell_type", "celltype", "braak", "cerad"]
DEFAULT_VAR_KEYS = ["_index", "gene_ids", "gene_symbols", "gene_symbol", "features", "gene_name"]


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
        raise RuntimeError("Forbidden expression-matrix access requested during metadata-only inspection.")


class SafeH5Reader:
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
    return "" if value is None else str(value)


def node_type(node: Any) -> str:
    if isinstance(node, h5py.Dataset):
        return "Dataset"
    if isinstance(node, h5py.Group):
        return "Group"
    return type(node).__name__


def node_shape(node: Any) -> str:
    shape = getattr(node, "shape", None)
    if shape is None:
        return ""
    return "x".join(str(item) for item in shape)


def node_dtype(node: Any) -> str:
    dtype = getattr(node, "dtype", None)
    return "" if dtype is None else str(dtype)


def encoding_type(node: Any) -> str:
    raw = getattr(node, "attrs", {}).get("encoding-type", "")
    return decode_scalar(raw)


def obs_index_name(obs_group: h5py.Group) -> str:
    return decode_scalar(obs_group.attrs.get("_index", "_index"))


def var_index_name(var_group: h5py.Group) -> str:
    return decode_scalar(var_group.attrs.get("_index", "_index"))


def node_length(node: Any) -> int:
    if isinstance(node, h5py.Dataset):
        return int(node.shape[0])
    if isinstance(node, h5py.Group) and "codes" in node:
        return int(node["codes"].shape[0])
    if isinstance(node, h5py.Group) and "values" in node:
        return node_length(node["values"])
    return 0


def read_dataset_values(dataset: h5py.Dataset) -> list[str]:
    raw = dataset[()]
    if getattr(raw, "shape", ()) == ():
        return [decode_scalar(raw)]
    return [decode_scalar(value) for value in raw]


def read_encoded_values(node: Any) -> list[str]:
    if isinstance(node, h5py.Dataset):
        return read_dataset_values(node)
    if isinstance(node, h5py.Group) and "categories" in node:
        return read_encoded_values(node["categories"])
    if isinstance(node, h5py.Group) and "values" in node:
        return read_encoded_values(node["values"])
    return []


def summarize_group(reader: SafeH5Reader, group_path: str, component: str) -> list[dict[str, str]]:
    if not reader.exists(group_path):
        return []
    group = reader.get(group_path)
    rows: list[dict[str, str]] = []
    for key in sorted(group.keys()):
        path = f"{group_path}/{key}"
        node = reader.get(path)
        rows.append(
            {
                "component": component,
                "key": key,
                "hdf5_path": path,
                "node_type": node_type(node),
                "shape": node_shape(node),
                "dtype": node_dtype(node),
                "encoding_type": encoding_type(node),
                "notes": "metadata key only; expression matrix not accessed",
            }
        )
    return rows


def shape_rows(reader: SafeH5Reader) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if reader.exists("obs"):
        obs_group = reader.get("obs")
        index_key = obs_index_name(obs_group)
        if reader.exists(f"obs/{index_key}"):
            rows.append(
                {
                    "component": "shape",
                    "key": "n_obs",
                    "hdf5_path": f"obs/{index_key}",
                    "node_type": "metadata_shape",
                    "shape": str(node_length(reader.get(f"obs/{index_key}"))),
                    "dtype": "",
                    "encoding_type": "",
                    "notes": "number of observation rows inferred from obs index only",
                }
            )
    if reader.exists("var"):
        var_group = reader.get("var")
        index_key = var_index_name(var_group)
        if reader.exists(f"var/{index_key}"):
            rows.append(
                {
                    "component": "shape",
                    "key": "n_vars",
                    "hdf5_path": f"var/{index_key}",
                    "node_type": "metadata_shape",
                    "shape": str(node_length(reader.get(f"var/{index_key}"))),
                    "dtype": "",
                    "encoding_type": "",
                    "notes": "number of variable rows inferred from var index only",
                }
            )
    return rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["component", "key", "hdf5_path", "node_type", "shape", "dtype", "encoding_type", "notes"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    logging.info("Wrote metadata summary: %s", path)


def choose_var_symbol_key(reader: SafeH5Reader) -> str | None:
    for key in DEFAULT_VAR_KEYS:
        if reader.exists(f"var/{key}"):
            return key
    return None


def write_var_genes(reader: SafeH5Reader, path: Path) -> int:
    if not reader.exists("var"):
        raise RuntimeError("No var group found in external H5AD metadata.")
    symbol_key = choose_var_symbol_key(reader)
    if symbol_key is None:
        raise RuntimeError("No supported var gene-symbol field found.")
    var_group = reader.get("var")
    index_key = var_index_name(var_group)
    symbol_values = read_encoded_values(reader.get(f"var/{symbol_key}"))
    index_values = read_encoded_values(reader.get(f"var/{index_key}")) if reader.exists(f"var/{index_key}") else []
    gene_id_values = read_encoded_values(reader.get("var/gene_ids")) if reader.exists("var/gene_ids") else []
    row_count = max(len(symbol_values), len(index_values), len(gene_id_values))

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_symbol", "gene_id", "var_index"], delimiter="\t")
        writer.writeheader()
        for index in range(row_count):
            symbol = symbol_values[index] if index < len(symbol_values) else ""
            gene_id = gene_id_values[index] if index < len(gene_id_values) else ""
            if not gene_id and index < len(index_values) and symbol_key != index_key:
                gene_id = index_values[index]
            writer.writerow({"gene_symbol": symbol, "gene_id": gene_id, "var_index": index})
    logging.info("Wrote lightweight var gene metadata rows: %d", row_count)
    return row_count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect external H5AD metadata without expression access.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/tables/mathys2019_metadata_summary.tsv"),
    )
    parser.add_argument(
        "--var-output",
        type=Path,
        default=Path("data/interim/external/mathys_2019/mathys_var_genes.tsv"),
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=Path("results/logs/36_inspect_mathys2019.log"),
    )
    parser.add_argument("--skip-var-output", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    logging.info("Starting metadata-only external H5AD inspection.")
    logging.info("Input: %s", args.input)
    with h5py.File(args.input, "r") as handle:
        reader = SafeH5Reader(handle)
        rows: list[dict[str, str]] = []
        rows.extend(shape_rows(reader))
        rows.extend(summarize_group(reader, "obs", "obs_key"))
        rows.extend(summarize_group(reader, "var", "var_key"))
        rows.extend(summarize_group(reader, "obs/__categories", "obs_category"))
        write_summary(args.output, rows)
        if not args.skip_var_output:
            write_var_genes(reader, args.var_output)

    logging.info("Metadata inspection complete. Expression matrix was not accessed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
