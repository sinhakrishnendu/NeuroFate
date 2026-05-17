#!/usr/bin/env python3
"""Run donor-level feature-group ablation benchmarks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

from neurofate.benchmarking import (
    feature_groups,
    load_config,
    load_donor_table,
    select_feature_columns,
    summarize_metric_rows,
    train_evaluate_split,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 12 feature-group ablation.")
    parser.add_argument("--features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--config", type=Path, default=Path("configs/benchmark_config.yaml"))
    parser.add_argument("--metrics-output", type=Path, default=Path("results/tables/phase12_feature_ablation_metrics.tsv"))
    parser.add_argument("--importance-output", type=Path, default=Path("results/tables/phase12_feature_group_importance.tsv"))
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.dry_run:
        print("Phase 12 feature-ablation plan:")
        print(f"  features: {args.features}")
        print(f"  feature groups: gene, cell_fraction, celltype_index, inflammatory, astrocyte, neuronal, mitochondrial")
        return 0
    frame = load_donor_table(args.features)
    all_features = select_feature_columns(frame)
    groups = feature_groups(all_features)
    seed_list = [int(seed) for seed in config["seed_list"]]
    rows = []
    for task in config["tasks"]:
        for seed in seed_list:
            rows.append(
                {
                    "feature_group": "all_features",
                    **train_evaluate_split(
                        frame,
                        task=task,
                        model_name=args.model,
                        seed=seed,
                        test_size=float(config["test_size"]),
                        min_class_count=int(config["min_class_count"]),
                        feature_columns=all_features,
                    ),
                }
            )
            for group_name, group_columns in groups.items():
                retained = [column for column in all_features if column not in set(group_columns)]
                rows.append(
                    {
                        "feature_group": f"without_{group_name}",
                        **train_evaluate_split(
                            frame,
                            task=task,
                            model_name=args.model,
                            seed=seed,
                            test_size=float(config["test_size"]),
                            min_class_count=int(config["min_class_count"]),
                            feature_columns=retained,
                        ),
                    }
                )
    args.metrics_output.parent.mkdir(parents=True, exist_ok=True)
    metrics = pd.DataFrame(rows)
    metrics.to_csv(args.metrics_output, sep="\t", index=False)
    summary = pd.DataFrame(summarize_metric_rows(rows, ["auroc", "balanced_accuracy"]))
    all_summary = summary[summary["task"].notna()].copy()
    importance_rows = []
    if not all_summary.empty:
        for task in config["tasks"]:
            baseline = metrics[(metrics["task"] == task) & (metrics["feature_group"] == "all_features")]
            baseline_mean = pd.to_numeric(baseline["auroc"], errors="coerce").mean()
            for group_name in groups:
                ablated = metrics[
                    (metrics["task"] == task) & (metrics["feature_group"] == f"without_{group_name}")
                ]
                ablated_mean = pd.to_numeric(ablated["auroc"], errors="coerce").mean()
                importance_rows.append(
                    {
                        "task": task,
                        "feature_group": group_name,
                        "baseline_auroc_mean": baseline_mean,
                        "ablated_auroc_mean": ablated_mean,
                        "delta_auroc_when_removed": baseline_mean - ablated_mean,
                    }
                )
    pd.DataFrame(importance_rows).to_csv(args.importance_output, sep="\t", index=False)
    print(f"Wrote {args.metrics_output}")
    print(f"Wrote {args.importance_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
