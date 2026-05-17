#!/usr/bin/env python3
"""Generate end-user NeuroFate report artifacts from existing results."""

from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


SECTIONS = [
    ("Dataset Summary", ["sea_ad_metadata_summary.tsv", "table1_sea_ad_publication_ready.tsv"]),
    ("Gene Panel", ["target_gene_panel_presence.tsv", "mathys_gene_overlap.tsv"]),
    ("SEA-AD Metadata", ["sea_ad_donor_summary.tsv", "sea_ad_celltype_by_ad_pathology.tsv"]),
    ("Sparse Extraction", ["target_gene_panel_presence.tsv"]),
    ("Phase 3 Biological Summaries", ["gene_by_celltype_summary.tsv", "microglial_activation_signature.tsv"]),
    ("Phase 4 Statistical Findings", ["phase4_gene_statistics.tsv", "phase4_composite_indices.tsv"]),
    ("Phase 5 Classical ML Models", ["phase5_model_metrics.tsv", "phase5_feature_importance.tsv"]),
    ("Phase 6 MPS Neural Model", ["phase6_mps_model_metrics.tsv", "phase6_mps_training_log.tsv"]),
    ("Phase 9 Mathys External Feasibility", ["phase9_mathys_external_validation_metrics.tsv", "mathys_2019_label_summary.tsv"]),
    ("Limitations", ["phase10_validation_warning_flags.tsv"]),
    ("Reproducibility Metadata", ["output_validation_report.tsv", "no_overclaiming_audit.tsv"]),
    ("Warnings and No-Overclaiming Notes", ["phase10_validation_warning_flags.tsv", "no_overclaiming_audit.tsv"]),
]


def read_tsv_preview(path: Path, max_rows: int) -> tuple[list[str], list[list[str]]]:
    if not path.exists():
        return [], []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        rows = []
        for index, row in enumerate(reader):
            if index >= max_rows:
                break
            rows.append(row)
    return header, rows


def markdown_table(header: list[str], rows: list[list[str]]) -> list[str]:
    if not header:
        return ["_Not available._"]
    lines = ["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |"]
    for row in rows:
        padded = row + [""] * (len(header) - len(row))
        lines.append("| " + " | ".join(padded[: len(header)]) + " |")
    return lines


def build_markdown(tables_dir: Path, max_rows: int) -> str:
    lines = [
        "# NeuroFate Analysis Report",
        "",
        "This report is generated from existing NeuroFate outputs. It does not run analysis.",
        "",
    ]
    for title, filenames in SECTIONS:
        lines.extend([f"## {title}", ""])
        for filename in filenames:
            path = tables_dir / filename
            lines.append(f"### `{filename}`")
            header, rows = read_tsv_preview(path, max_rows)
            lines.extend(markdown_table(header, rows))
            lines.append("")
    return "\n".join(lines)


def markdown_to_html(markdown: str) -> str:
    body = []
    for line in markdown.splitlines():
        if line.startswith("# "):
            body.append(f"<h1>{html.escape(line[2:])}</h1>")
        elif line.startswith("## "):
            body.append(f"<h2>{html.escape(line[3:])}</h2>")
        elif line.startswith("### "):
            body.append(f"<h3>{html.escape(line[4:])}</h3>")
        elif line.startswith("| "):
            body.append(f"<pre>{html.escape(line)}</pre>")
        elif not line:
            body.append("")
        else:
            body.append(f"<p>{html.escape(line)}</p>")
    return "<!doctype html><html><head><meta charset='utf-8'><title>NeuroFate Report</title></head><body>\n" + "\n".join(body) + "\n</body></html>\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate end-user NeuroFate analysis report.")
    parser.add_argument("--tables-dir", type=Path, default=Path("results/tables"))
    parser.add_argument("--reports-dir", type=Path, default=Path("results/reports"))
    parser.add_argument("--max-preview-rows", type=int, default=12)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.reports_dir.mkdir(parents=True, exist_ok=True)
    markdown = build_markdown(args.tables_dir, args.max_preview_rows)
    md_path = args.reports_dir / "neurofate_analysis_report.md"
    html_path = args.reports_dir / "neurofate_analysis_report.html"
    md_path.write_text(markdown, encoding="utf-8")
    html_path.write_text(markdown_to_html(markdown), encoding="utf-8")
    print(f"Wrote {md_path}")
    print(f"Wrote {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
