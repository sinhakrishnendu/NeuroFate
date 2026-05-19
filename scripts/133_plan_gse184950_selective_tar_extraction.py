#!/usr/bin/env python3
"""Plan manual selective GSE184950 processed-matrix extraction from RAW archive listings."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path


PLAN_COLUMNS = [
    "sample_name",
    "processed_tar_name",
    "archive_member_path",
    "manual_action",
    "output_directory",
    "fastq_handling",
    "status",
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


def basename(path: str) -> str:
    return Path(path).name


def load_axis_gene_count(axis_registry: Path) -> int:
    rows = read_tsv(axis_registry)
    genes: set[str] = set()
    for row in rows:
        genes.update(gene.strip().upper() for gene in row.get("gene_members", "").replace(",", ";").split(";") if gene.strip())
    return len(genes)


def build_plan(inventory: list[dict[str, str]], metadata: list[dict[str, str]], axis_registry: Path) -> list[dict[str, str]]:
    load_axis_gene_count(axis_registry)
    members = {basename(row.get("member_path", "")): row.get("member_path", "") for row in inventory}
    plan: list[dict[str, str]] = []
    for row in metadata:
        processed_tar = row.get("processed_tar_name", "")
        sample = row.get("sample_name", processed_tar.replace(".tar.gz", ""))
        member_path = members.get(processed_tar, "")
        status = "ready_for_manual_selective_processed_tar_extraction" if member_path else "expected_processed_tar_not_in_inventory"
        plan.append(
            {
                "sample_name": sample,
                "processed_tar_name": processed_tar,
                "archive_member_path": member_path,
                "manual_action": "extract_processed_tar_only_then_run_127" if member_path else "inspect_archive_inventory_or_download_raw_archive",
                "output_directory": f"data/interim/external/gse184950_pd_sn/processed_matrices/{sample}",
                "fastq_handling": "skip_fastq_no_raw_preprocessing",
                "status": status,
            }
        )
    return plan


def write_manual_script(path: Path, plan_rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    member_lines = "\n".join(
        f'# tar -xf "${{RAW_ARCHIVE}}" -C "${{OUT_ROOT}}" "{row["archive_member_path"]}"'
        for row in plan_rows
        if row.get("archive_member_path")
    )
    path.write_text(
        f"""#!/usr/bin/env bash
set -euo pipefail

echo "MANUAL_HEAVY: GSE184950 selective processed-matrix extraction."
echo "Do not run from Codex. Manual user execution only."
echo "This script extracts reviewed processed tar.gz members only; FASTQ/SRA processing is not used."

RUN_MANUAL_GSE184950_EXTRACTION="${{RUN_MANUAL_GSE184950_EXTRACTION:-NO}}"
if [[ "${{RUN_MANUAL_GSE184950_EXTRACTION}}" != "YES" ]]; then
  echo "Set RUN_MANUAL_GSE184950_EXTRACTION=YES only after reviewing phase25_gse184950_selective_extraction_plan.tsv."
  exit 1
fi

RAW_ARCHIVE="data/raw/external/gse184950_pd_sn/GSE184950_RAW.tar"
OUT_ROOT="data/interim/external/gse184950_pd_sn/processed_matrices"
mkdir -p "${{OUT_ROOT}}"

# Review and uncomment specific processed tar members only. Do not extract FASTQ files.
{member_lines}

# After selective extraction and per-sample unpacking, run the guarded axis-gene extractor manually:
# python scripts/127_extract_gse184950_axis_genes_from_10x.py --run-manual-extraction YES
""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan selective GSE184950 processed tar extraction.")
    parser.add_argument("--archive-inventory", type=Path, default=Path("results/tables/phase24_gse184950_raw_archive_inventory.tsv"))
    parser.add_argument("--series-metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--axis-registry", type=Path, default=Path("metadata/neurofate_axis_registry.tsv"))
    parser.add_argument("--output-plan", type=Path, default=Path("results/tables/phase25_gse184950_selective_extraction_plan.tsv"))
    parser.add_argument("--manual-script-output", type=Path, default=Path("results/logs/manual_phase25_gse184950_selective_extraction.sh"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/133_plan_gse184950_selective_tar_extraction.log"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.log_file)
    rows = build_plan(read_tsv(args.archive_inventory), read_tsv(args.series_metadata), args.axis_registry)
    write_tsv(args.output_plan, rows, PLAN_COLUMNS)
    write_manual_script(args.manual_script_output, rows)
    logging.info("Wrote GSE184950 selective extraction plan rows=%d without extracting archives", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
