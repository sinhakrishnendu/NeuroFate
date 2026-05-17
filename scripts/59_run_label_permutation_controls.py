#!/usr/bin/env python3
"""Run label-permutation controls for donor-level NeuroFate benchmarks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd

from neurofate.benchmarking import load_config, load_donor_table, task_labels, train_evaluate_split


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 12 label permutation controls.")
    parser.add_argument("--features", type=Path, default=Path("results/tables/phase5_donor_feature_table.tsv"))
    parser.add_argument("--config", type=Path, default=Path("configs/benchmark_config.yaml"))
    parser.add_argument("--metrics-output", type=Path, default=Path("results/tables/phase12_permutation_metrics.tsv"))
    parser.add_argument("--pvalues-output", type=Path, default=Path("results/tables/phase12_empirical_pvalues.tsv"))
    parser.add_argument("--model", default="logistic_regression")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    if args.dry_run:
        print("Phase 12 permutation-control plan:")
        print(f"  features: {args.features}")
        print(f"  permutations: {config['n_permutations']}")
        print(f"  model: {args.model}")
        return 0
    frame = load_donor_table(args.features)
    seed = int(config["seed_list"][0])
    rows = []
    pvalue_rows = []
    for task in config["tasks"]:
        observed = train_evaluate_split(
            frame,
            task=task,
            model_name=args.model,
            seed=seed,
            test_size=float(config["test_size"]),
            min_class_count=int(config["min_class_count"]),
        )
        observed_auroc = observed.get("auroc", np.nan)
        labels = task_labels(frame, task)
        valid = labels.notna()
        if observed.get("status") != "ok":
            pvalue_rows.append(
                {
                    "task": task,
                    "model": args.model,
                    "observed_auroc": observed_auroc,
                    "empirical_pvalue": np.nan,
                    "status": observed.get("status", "skipped"),
                }
            )
            continue
        rng = np.random.default_rng(seed)
        permuted_aurocs = []
        for permutation_id in range(int(config["n_permutations"])):
            permuted = frame.copy()
            shuffled = labels[valid].to_numpy().copy()
            rng.shuffle(shuffled)
            permuted.loc[valid, "__permuted_label__"] = shuffled
            original_task_labels = task_labels
            # Keep the model path identical by swapping the task label into Cognitive Status bins.
            if task == "dementia_vs_reference":
                permuted.loc[valid, "label__Cognitive_Status"] = np.where(
                    shuffled == 1, "Dementia", "Reference"
                )
            elif task == "high_vs_low_ad_neuropathology":
                permuted.loc[valid, "label__Overall_AD_neuropathological_Change"] = np.where(
                    shuffled == 1, "High", "Low"
                )
            elif task == "apoe_risk_prediction":
                permuted.loc[valid, "label__APOE_Genotype"] = np.where(shuffled == 1, "3/4", "3/3")
            else:
                permuted.loc[valid, "label__LATE"] = np.where(shuffled == 1, "LATE Stage 1", "Reference")
            metric = train_evaluate_split(
                permuted,
                task=task,
                model_name=args.model,
                seed=seed + permutation_id + 1,
                test_size=float(config["test_size"]),
                min_class_count=int(config["min_class_count"]),
            )
            auroc = metric.get("auroc", np.nan)
            permuted_aurocs.append(auroc)
            rows.append(
                {
                    "task": task,
                    "model": args.model,
                    "permutation_id": permutation_id,
                    "observed_auroc": observed_auroc,
                    "permuted_auroc": auroc,
                    "status": metric.get("status", "ok"),
                }
            )
        valid_perm = np.asarray([value for value in permuted_aurocs if not pd.isna(value)])
        empirical = (float(np.sum(valid_perm >= observed_auroc)) + 1.0) / (len(valid_perm) + 1.0)
        pvalue_rows.append(
            {
                "task": task,
                "model": args.model,
                "observed_auroc": observed_auroc,
                "empirical_pvalue": empirical,
                "status": "ok",
            }
        )
    args.metrics_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.metrics_output, sep="\t", index=False)
    pd.DataFrame(pvalue_rows).to_csv(args.pvalues_output, sep="\t", index=False)
    print(f"Wrote {args.metrics_output}")
    print(f"Wrote {args.pvalues_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
