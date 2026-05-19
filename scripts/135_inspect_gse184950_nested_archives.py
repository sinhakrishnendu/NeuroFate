#!/usr/bin/env python3
"""Inspect nested GSE184950 per-sample archives without extracting files."""

from __future__ import annotations

import argparse
import csv
import logging
import tarfile
from pathlib import Path


OUTPUT_COLUMNS = [
    "outer_member_path",
    "gsm_accession",
    "sample_id",
    "nested_member_path",
    "nested_size_bytes",
    "extension",
    "likely_role",
    "complete_processed_matrix_set",
]


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


def extension_for(name: str) -> str:
    lower = name.lower()
    for suffix in [".matrix.mtx.gz", ".barcodes.tsv.gz", ".features.tsv.gz", ".genes.tsv.gz", ".fastq.gz", ".mtx.gz", ".tsv.gz", ".tar.gz"]:
        if lower.endswith(suffix):
            return suffix
    return Path(lower).suffix


def likely_role(name: str) -> str:
    lower = name.lower()
    if "matrix.mtx" in lower:
        return "tenx_matrix"
    if "barcodes.tsv" in lower:
        return "tenx_barcodes"
    if "features.tsv" in lower or "genes.tsv" in lower:
        return "tenx_features"
    if "fastq" in lower or lower.endswith((".fq", ".fq.gz")):
        return "raw_sequence_do_not_process"
    if "filtered_feature_bc_matrix" in lower:
        return "count_matrix_folder"
    return "other"


def parse_outer_name(name: str) -> tuple[str, str]:
    base = Path(name).name.replace(".tar.gz", "").replace(".tgz", "")
    parts = base.split("_", 1)
    if len(parts) == 2 and parts[0].startswith("GSM"):
        return parts[0], parts[1]
    return "", base


def expected_outer_members(manifest: list[dict[str, str]]) -> set[str]:
    expected = set()
    for row in manifest:
        for column in ("processed_tar_name", "expected_archive_member"):
            value = row.get(column, "")
            if value:
                expected.add(Path(value).name)
    return expected


def inspect_nested(raw_tar: Path, manifest: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    expected = expected_outer_members(manifest)
    with tarfile.open(raw_tar, "r:*") as outer:
        for outer_member in outer.getmembers():
            outer_base = Path(outer_member.name).name
            if not outer_member.isfile() or not outer_base.endswith((".tar.gz", ".tgz")):
                continue
            if expected and outer_base not in expected:
                continue
            gsm, sample = parse_outer_name(outer_base)
            nested_handle = outer.extractfile(outer_member)
            if nested_handle is None:
                continue
            with tarfile.open(fileobj=nested_handle, mode="r:*") as nested:
                nested_members = nested.getmembers()
                roles = {likely_role(member.name) for member in nested_members}
                complete = {"tenx_matrix", "tenx_barcodes", "tenx_features"}.issubset(roles)
                for member in nested_members:
                    rows.append(
                        {
                            "outer_member_path": outer_member.name,
                            "gsm_accession": gsm,
                            "sample_id": sample,
                            "nested_member_path": member.name,
                            "nested_size_bytes": str(member.size),
                            "extension": extension_for(member.name),
                            "likely_role": likely_role(member.name),
                            "complete_processed_matrix_set": str(complete).lower(),
                        }
                    )
    return rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    samples = sorted({row["sample_id"] for row in rows if row.get("sample_id")})
    complete = sorted({row["sample_id"] for row in rows if row.get("complete_processed_matrix_set") == "true"})
    role_counts: dict[str, int] = {}
    for row in rows:
        role = row.get("likely_role", "other")
        role_counts[role] = role_counts.get(role, 0) + 1
    lines = [
        "# Phase 26 GSE184950 Nested Archive Summary",
        "",
        "Nested per-sample archives were listed only; no files were extracted.",
        "",
        f"- Samples with nested members listed: {len(samples)}",
        f"- Samples with complete processed 10x matrix sets: {len(complete)}",
        "",
        "Role counts:",
    ]
    lines.extend(f"- {role}: {count}" for role, count in sorted(role_counts.items()))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect GSE184950 nested per-sample archives without extraction.")
    parser.add_argument("--raw-tar", type=Path, default=Path("data/raw/external/gse184950_pd_sn/GSE184950_RAW.tar"))
    parser.add_argument("--series-manifest", type=Path, default=Path("results/tables/phase25_gse184950_series_processed_file_manifest.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/tables/phase26_gse184950_nested_archive_inventory.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/reports/phase26_gse184950_nested_archive_summary.md"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/135_inspect_gse184950_nested_archives.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = inspect_nested(args.raw_tar, read_tsv(args.series_manifest))
    write_tsv(args.output, rows, OUTPUT_COLUMNS)
    write_summary(args.summary_output, rows)
    logging.info("Listed nested GSE184950 archive members rows=%d without extraction", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
