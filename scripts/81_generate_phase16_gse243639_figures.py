#!/usr/bin/env python3
"""Generate Phase 16 GSE243639 planning/validation figures with matplotlib only."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_bar(path: Path, labels: list[str], values: list[float], title: str, ylabel: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels or ["unavailable"], values or [0.0], color="#5876a3")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def figure_sample_summary(label_summary: Path, output: Path) -> None:
    rows = [row for row in read_tsv(label_summary) if row.get("label_field") == "diagnosis"]
    save_bar(output, [row["label"] for row in rows], [float(row["sample_count"]) for row in rows], "GSE243639 sample diagnosis summary", "sample count")


def figure_gene_overlap(audit: Path, output: Path) -> None:
    rows = read_tsv(audit)
    if rows:
        row = rows[0]
        labels = ["extracted", "missing"]
        values = [float(row.get("extracted_target_genes", 0) or 0), float(row.get("missing_target_genes", 0) or 0)]
    else:
        labels, values = ["unavailable"], [0.0]
    save_bar(output, labels, values, "GSE243639 NeuroFate target-gene overlap", "gene count")


def figure_validation(metrics: Path, output: Path) -> None:
    rows = read_tsv(metrics)
    labels = []
    values = []
    for row in rows:
        try:
            labels.append(row.get("validation_mode", "unknown"))
            values.append(float(row.get("auroc", "nan")))
        except ValueError:
            continue
    save_bar(output, labels, values, "GSE243639 PD validation AUROC", "AUROC")


def figure_feature_space(schema: Path, output: Path) -> None:
    rows = read_tsv(schema)
    counts = {"shared": 0, "schema_specific": 0}
    for row in rows:
        status = row.get("status", "schema_specific")
        counts[status] = counts.get(status, 0) + 1
    save_bar(output, list(counts), [float(value) for value in counts.values()], "GSE243639 feature schema alignment", "feature count")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 16 GSE243639 figures.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    figure_sample_summary(args.tables_dir / "phase16_gse243639_label_summary.tsv", args.figures_dir / "figure35_gse243639_pd_sample_summary.png")
    figure_gene_overlap(args.tables_dir / "phase16_gse243639_gene_extraction_audit.tsv", args.figures_dir / "figure36_gse243639_target_gene_overlap.png")
    figure_validation(args.tables_dir / "phase16_gse243639_external_validation_metrics.tsv", args.figures_dir / "figure37_gse243639_pd_validation.png")
    figure_feature_space(args.tables_dir / "phase16_gse243639_feature_schema_alignment.tsv", args.figures_dir / "figure38_gse243639_feature_space_summary.png")
    print(f"Wrote Phase 16 figures under {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
