#!/usr/bin/env python3
"""Run a tiny, self-contained NeuroFate demo without external data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from neurofate.demo import build_demo_outputs


DEFAULT_METADATA = Path("examples/tiny_demo/tiny_metadata.tsv")
DEFAULT_PANEL = Path("examples/tiny_demo/tiny_gene_panel.tsv")
DEFAULT_EXPRESSION = Path("examples/tiny_demo/tiny_sparse_expression.tsv.gz")
DEFAULT_OUTDIR = Path("results/demo")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the tiny bundled NeuroFate demo.")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--gene-panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--expression", type=Path, default=DEFAULT_EXPRESSION)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        build_demo_outputs(args.metadata, args.gene_panel, args.expression, args.outdir)
    except FileNotFoundError as exc:
        print(f"NeuroFate demo could not start: {exc}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
