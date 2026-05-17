#!/usr/bin/env python3
"""Run configurable repeated donor-level baseline benchmarks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

from neurofate.benchmarking import load_config, load_donor_table, summarize_metric_rows, train_evaluate_split


METRICS = ["auroc", "auprc", "balanced_accuracy", "brier"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run repeated Phase 12 donor-level benchmarks.")
    parser.add_argument("--features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--config", type=Path, default=Path("configs/benchmark_config.yaml"))
    parser.add_argument("--per-seed-output", type=Path, default=Path("results/tables/phase12_repeated_benchmark_metrics.tsv"))
    parser.add_argument("--summary-output", type=Path, default=Path("results/tables/phase12_repeated_benchmark_summary.tsv"))
    parser.add_argument("--dry-run", action="store_true", help="Print planned benchmark size and exit.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.dry_run:
        print("Phase 12 repeated benchmark plan:")
        print(f"  features: {args.features}")
        print(f"  seeds: {config.get('seed_list', [])}")
        print(f"  models: {config.get('models', [])}")
        print(f"  tasks: {config.get('tasks', [])}")
        return 0
    frame = load_donor_table(args.features)
    rows = []
    for task in config["tasks"]:
        for model in config["models"]:
            for seed in config["seed_list"]:
                rows.append(
                    train_evaluate_split(
                        frame,
                        task=task,
                        model_name=model,
                        seed=int(seed),
                        test_size=float(config["test_size"]),
                        min_class_count=int(config["min_class_count"]),
                    )
                )
    args.per_seed_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.per_seed_output, sep="\t", index=False)
    pd.DataFrame(summarize_metric_rows(rows, METRICS)).to_csv(args.summary_output, sep="\t", index=False)
    print(f"Wrote {args.per_seed_output}")
    print(f"Wrote {args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
