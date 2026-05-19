#!/usr/bin/env python3
"""Generate Phase 26 GSE184950 replication figures from lightweight TSV outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def to_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def save_bar(labels: list[str], values: list[float], title: str, ylabel: str, output: Path) -> None:
    import matplotlib.pyplot as plt

    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, values, color="#4C78A8")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def disease_counts(metadata: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    counts: dict[str, int] = {}
    for row in metadata:
        label = row.get("disease_state", "") or "missing"
        counts[label] = counts.get(label, 0) + 1
    return list(counts), [float(value) for value in counts.values()]


def coverage_counts(audit: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    labels = [row.get("sample_id", "") for row in audit[:20]]
    values = [to_float(row.get("genes_found")) for row in audit[:20]]
    return labels, values


def effect_counts(stats: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    labels = [row.get("axis_id", "").replace("_axis", "") for row in stats]
    values = [to_float(row.get("effect_size")) for row in stats]
    return labels, values


def readiness_counts(readiness: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    counts: dict[str, int] = {}
    for row in readiness:
        status = row.get("status", "") or "missing"
        counts[status] = counts.get(status, 0) + 1
    return list(counts), [float(value) for value in counts.values()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 26 GSE184950 figures.")
    parser.add_argument("--metadata", type=Path, default=Path("results/tables/phase25_gse184950_series_sample_metadata.tsv"))
    parser.add_argument("--axis-audit", type=Path, default=Path("results/tables/phase26_gse184950_axis_gene_extraction_audit.tsv"))
    parser.add_argument("--replication-stats", type=Path, default=Path("results/tables/phase25_gse184950_axis_replication_statistics.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("results/reports/phase23_pnas_readiness_matrix.tsv"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    save_bar(*disease_counts(read_tsv(args.metadata)), "GSE184950 Sample Summary", "Samples", args.figures_dir / "figure60_gse184950_sample_summary.png")
    save_bar(*coverage_counts(read_tsv(args.axis_audit)), "GSE184950 Axis Gene Coverage", "Genes found", args.figures_dir / "figure61_gse184950_axis_coverage.png")
    save_bar(*effect_counts(read_tsv(args.replication_stats)), "GSE184950 Axis Replication Effects", "Rank-biserial effect", args.figures_dir / "figure62_gse184950_axis_replication_effects.png")
    save_bar(*readiness_counts(read_tsv(args.readiness)), "PNAS Replication Status", "Criteria", args.figures_dir / "figure63_pnas_replication_status.png")
    print(f"Wrote Phase 26 figures to {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
