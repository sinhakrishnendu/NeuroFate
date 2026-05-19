#!/usr/bin/env python3
"""Generate lightweight Phase 27 GSE184950 clean replication figures."""

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
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        output.with_suffix(output.suffix + ".skipped.txt").write_text(
            f"{title}\nMatplotlib is not installed in this environment; figure generation was skipped.\n",
            encoding="utf-8",
        )
        return
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(labels, values, color="#4C78A8")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output, dpi=200)
    plt.close(fig)


def sample_integrity_values(rows: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    row = rows[0] if rows else {}
    return ["expected", "observed", "valid", "invalid"], [
        to_float(row.get("expected_samples")),
        to_float(row.get("observed_axis_score_samples")),
        to_float(row.get("valid_axis_score_samples")),
        to_float(row.get("invalid_axis_score_samples")),
    ]


def replication_effects(rows: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    return [row.get("axis_id", "").replace("_axis", "") for row in rows], [to_float(row.get("effect_size")) for row in rows]


def evidence_counts(rows: list[dict[str, str]]) -> tuple[list[str], list[float]]:
    counts: dict[str, int] = {}
    for row in rows:
        label = row.get("evidence_label", "") or "missing"
        counts[label] = counts.get(label, 0) + 1
    return list(counts), [float(value) for value in counts.values()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Phase 27 clean GSE184950 figures.")
    parser.add_argument("--sample-integrity", type=Path, default=Path("results/tables/phase27_gse184950_sample_integrity_audit.tsv"))
    parser.add_argument("--replication-stats", type=Path, default=Path("results/tables/phase27_gse184950_axis_replication_statistics_clean.tsv"))
    parser.add_argument("--figures-dir", type=Path, default=Path("results/figures"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    save_bar(*sample_integrity_values(read_tsv(args.sample_integrity)), "GSE184950 Clean Sample Integrity", "Samples", args.figures_dir / "figure64_gse184950_clean_sample_integrity.png")
    save_bar(*replication_effects(read_tsv(args.replication_stats)), "GSE184950 Clean Axis Replication Effects", "Rank-biserial effect", args.figures_dir / "figure65_gse184950_clean_axis_replication_effects.png")
    save_bar(*evidence_counts(read_tsv(args.replication_stats)), "GSE184950 Replication Evidence Category", "Axes", args.figures_dir / "figure66_gse184950_replication_evidence_category.png")
    print(f"Wrote Phase 27 figures to {args.figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
