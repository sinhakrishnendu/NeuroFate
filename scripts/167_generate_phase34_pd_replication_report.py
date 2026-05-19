#!/usr/bin/env python3
"""Generate the Phase 34 PD microarray/bulk replication expansion report."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def counts(rows: list[dict[str, str]], key: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for row in rows:
        label = row.get(key, "") or "missing"
        out[label] = out.get(label, 0) + 1
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Phase 34 PD replication report.")
    parser.add_argument("--registry", type=Path, default=Path("metadata/phase34_pd_replication_registry.tsv"))
    parser.add_argument("--inventory", type=Path, default=Path("results/reports/phase34_pd_geo_file_inventory.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/reports/phase34_pd_replication_expansion_report.md"))
    args = parser.parse_args()

    registry = read_tsv(args.registry)
    inventory = read_tsv(args.inventory)
    stats_paths = [
        *sorted(Path("results/tables").glob("phase34_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase35_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase36_*_pd_axis_replication_statistics.tsv")),
        *sorted(Path("results/tables").glob("phase37_*_pd_axis_replication_statistics.tsv")),
    ]
    coverage_paths = [
        *sorted(Path("results/tables").glob("phase34_*_axis_feature_coverage.tsv")),
        *sorted(Path("results/tables").glob("phase35_*_axis_feature_coverage.tsv")),
        *sorted(Path("results/tables").glob("phase36_*_axis_feature_coverage.tsv")),
        *sorted(Path("results/tables").glob("phase37_*_axis_feature_coverage.tsv")),
    ]

    lines = [
        "# Phase 34 PD Replication Expansion Report",
        "",
        "## 1. Why Phase 34 Was Needed",
        "Independent AD replication is nominally supported, but independent PD axis replication remains the limiting bottleneck for stronger cross-disease claims.",
        "",
        "## 2. Candidate Datasets",
    ]
    for row in registry:
        lines.append(f"- `{row['cohort_id']}` ({row['geo_accession']}): priority {row['priority']}; {row['tissue_or_cell_type']}; expected {row['expected_case_count']} cases and {row['expected_control_count']} controls.")
    lines.extend(["", "## 3. Local File Status"])
    if inventory:
        for row in inventory:
            if row.get("likely_role") == "missing":
                lines.append(f"- `{row['cohort_id']}`: manual acquisition required.")
    else:
        lines.append("- File inventory has not been generated.")
    lines.extend(["", "## 4. Metadata Endpoint Availability"])
    metadata_paths = [
        *sorted(Path("results/tables").glob("phase34_*_sample_metadata.tsv")),
        *sorted(Path("results/tables").glob("phase35_*_sample_metadata.tsv")),
        *sorted(Path("results/tables").glob("phase36_*_sample_metadata.tsv")),
        *sorted(Path("results/tables").glob("phase37_*_sample_metadata.tsv")),
    ]
    for path in metadata_paths:
        rows = read_tsv(path)
        lines.append(f"- `{path.name}`: endpoint counts {counts(rows, 'label__pd_vs_control')}.")
    if not metadata_paths:
        lines.append("- No Phase 34 metadata tables are available yet.")
    lines.extend(["", "## 5. Probe/Gene Mapping Status"])
    for path in sorted(Path("results/tables").glob("phase34_*_axis_probe_mapping.tsv")):
        rows = read_tsv(path)
        lines.append(f"- `{path.name}`: {len(rows)} NeuroFate probe-to-gene mappings.")
    lines.extend(["", "## 6. Axis-Score Availability"])
    if coverage_paths:
        for path in coverage_paths:
            rows = read_tsv(path)
            lines.append(f"- `{path.name}`: {sum(row.get('status') == 'ok' for row in rows)}/{len(rows)} axes have mapped genes.")
    else:
        lines.append("- No Phase 34 axis-score coverage is available yet.")
    lines.extend(["", "## 7. PD Replication Statistics If Available"])
    if stats_paths:
        for path in stats_paths:
            rows = read_tsv(path)
            lines.append(f"- `{path.name}`: {counts(rows, 'evidence_label')}.")
    else:
        lines.append("- No Phase 34 replication statistics are available yet.")
    lines.extend(["", "## 8. Whether Neuronal Vulnerability Axis Replicates"])
    neuronal = []
    for path in stats_paths:
        neuronal.extend((path.name, row) for row in read_tsv(path) if row.get("axis_id") == "neuronal_vulnerability_axis")
    if neuronal:
        for name, row in neuronal:
            lines.append(f"- `{name}`: effect={row.get('effect_size', '')}; p={row.get('pvalue', '')}; FDR={row.get('fdr', '')}; label={row.get('evidence_label', '')}.")
    else:
        lines.append("- Not yet tested in Phase 34.")
    lines.extend(["", "## 9. PNAS Implication", "Phase 34 is designed to close the PD replication gap. Direction-only support remains preliminary and does not establish confirmed cross-disease biology.", "", "## 10. Next Manual Command", "`RUN_MANUAL_DOWNLOAD=YES scripts/manual_downloads/download_gse20141_manual.sh`", ""])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
