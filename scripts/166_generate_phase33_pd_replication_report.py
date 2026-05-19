#!/usr/bin/env python3
"""Generate a conservative Phase 33 PD replication expansion report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_labels(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = row.get(column, "") or "missing"
        counts[key] = counts.get(key, 0) + 1
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Phase 33 PD replication expansion report.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase33_pd_replication_registry.tsv"))
    parser.add_argument("--inventory", type=Path, default=Path("results/reports/phase33_pd_replication_file_inventory.tsv"))
    parser.add_argument("--results-glob", default="results/tables/phase33_*_pd_axis_replication_statistics.tsv")
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase33_pd_replication_expansion_report.md"))
    args = parser.parse_args()

    registry = read_tsv(args.registry)
    inventory = read_tsv(args.inventory)
    result_paths = sorted(Path(".").glob(args.results_glob))
    result_rows = [(path, read_tsv(path)) for path in result_paths]

    lines = [
        "# Phase 33 PD Replication Expansion Report",
        "",
        "## 1. Why Phase 33 Was Needed",
        "AD replication is now nominally supported, while independent PD axis replication remains the main bottleneck for stronger cross-disease claims.",
        "",
        "## 2. Candidate PD Replication Cohorts",
    ]
    for row in registry:
        lines.append(f"- `{row['cohort_id']}` ({row['geo_accession']}): priority {row['priority']}; {row['tissue_or_cell_type']}; strategy: {row['first_use_strategy']}.")
    lines.extend(["", "## 3. Local Acquisition Status"])
    if inventory:
        for row in inventory:
            if row.get("likely_role") == "missing":
                lines.append(f"- `{row['cohort_id']}`: local files missing.")
    else:
        lines.append("- File inventory has not been generated yet.")
    lines.extend(["", "## 4. Expression/Metadata Availability"])
    lines.append("Phase 33 uses GEO series matrices and processed sample-level expression tables when present. FASTQ/SRA processing is out of scope.")
    lines.extend(["", "## 5. Axis Coverage"])
    coverage_paths = sorted(Path("results/tables").glob("phase33_*_axis_feature_coverage.tsv"))
    if coverage_paths:
        for path in coverage_paths:
            rows = read_tsv(path)
            ok = sum(1 for row in rows if row.get("status") == "ok")
            lines.append(f"- `{path.name}`: {ok}/{len(rows)} axes have at least one mapped gene.")
    else:
        lines.append("- No Phase 33 axis coverage tables are available yet.")
    lines.extend(["", "## 6. PD Replication Results If Available"])
    if result_rows:
        for path, rows in result_rows:
            counts = count_labels(rows, "evidence_label")
            lines.append(f"- `{path.name}`: {counts}.")
    else:
        lines.append("- No Phase 33 PD replication statistics are available yet.")
    lines.extend(["", "## 7. Whether Neuronal Vulnerability Replicates In PD"])
    neuronal = []
    for path, rows in result_rows:
        neuronal.extend((path.name, row) for row in rows if row.get("axis_id") == "neuronal_vulnerability_axis")
    if neuronal:
        for name, row in neuronal:
            lines.append(f"- `{name}`: effect={row.get('effect_size', '')}; p={row.get('pvalue', '')}; FDR={row.get('fdr', '')}; label={row.get('evidence_label', '')}.")
    else:
        lines.append("- Not yet tested in a Phase 33 cohort.")
    lines.extend(
        [
            "",
            "## 8. PNAS Implication",
            "The immediate implication is practical: prioritize GSE20141 first, because laser-captured substantia nigra pars compacta neurons are the most direct donor/sample-level test of the neuronal vulnerability axis.",
            "",
            "## 9. Next Acquisition Priority",
            "Acquire `GSE20141_series_matrix.txt.gz` and any processed supplementary expression table, then run Phase 33 axis extraction. Direction-only evidence must remain preliminary and must not be used for confirmed cross-disease biology.",
            "",
        ]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
