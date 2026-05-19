#!/usr/bin/env python3
"""Run Phase 20 sample-level GSE243639 cell-type-aware PD validation."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


def load_phase91_module():
    script = Path(__file__).resolve().parent / "91_run_gse243639_repaired_celltype_pd_validation.py"
    spec = importlib.util.spec_from_file_location("phase91_pd_validation", script)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load script 91 validation module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 20 GSE243639 sample-level cell-type-aware PD validation.")
    parser.add_argument("--features", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_table.tsv"))
    parser.add_argument("--annotation-summary", type=Path, default=Path("results/tables/phase20_gse243639_feature_group_counts.tsv"))
    parser.add_argument("--metrics-output", type=Path, default=Path("results/tables/phase20_gse243639_celltype_validation_metrics.tsv"))
    parser.add_argument("--predictions-output", type=Path, default=Path("results/tables/phase20_gse243639_celltype_predictions.tsv"))
    parser.add_argument("--importance-output", type=Path, default=Path("results/tables/phase20_gse243639_celltype_feature_importance.tsv"))
    parser.add_argument("--log-file", type=Path, default=Path("results/logs/102_run_gse243639_phase20_celltype_pd_validation.log"))
    parser.add_argument("--n-permutations", type=int, default=100)
    parser.add_argument("--n-bootstrap", type=int, default=200)
    parser.add_argument("--n-repeats", type=int, default=10)
    parser.add_argument("--seed", type=int, default=201)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    phase91 = load_phase91_module()
    phase91.configure_logging(args.log_file)
    rows, y = phase91.filtered_table(phase91.read_tsv(args.features))
    counts = phase91.Counter(y.tolist())
    features = phase91.feature_columns(rows)
    match_rate = phase91.read_match_rate(args.annotation_summary)
    metric_fields = [
        "model",
        "validation_mode",
        "task_id",
        "auroc",
        "auroc_sd",
        "auroc_ci_low",
        "auroc_ci_high",
        "auprc",
        "balanced_accuracy",
        "brier_score",
        "empirical_permutation_pvalue",
        "null_auroc_mean",
        "n_samples",
        "positive_count",
        "negative_count",
        "feature_count",
        "annotation_match_rate",
        "reliability_flag",
        "notes",
    ]
    if len(rows) < 8 or len(counts) < 2 or min(counts.values()) < 2 or len(features) < phase91.MIN_REPAIRED_FEATURES or match_rate < phase91.MIN_MATCH_RATE:
        reason = "Phase 20 safe-map feature table did not meet sample, class, feature, or annotation match-rate requirements"
        phase91.write_tsv(args.metrics_output, phase91.failure_metric(reason, len(rows), counts.get(1, 0), counts.get(0, 0), len(features), match_rate), metric_fields)
        phase91.write_tsv(args.predictions_output, [], ["model", "validation_mode", "split_id", "row_index", "true_label", "predicted_probability", "predicted_label"])
        phase91.write_tsv(args.importance_output, [], ["rank", "feature", "importance", "model"])
        return 0
    X = phase91.matrix(rows, features)
    logistic_metrics, split_predictions = phase91.repeated_split_predictions(X, y, "logistic_regression", args.seed, args.n_repeats)
    rf_metrics, rf_predictions = phase91.repeated_split_predictions(X, y, "random_forest_baseline", args.seed + 500, args.n_repeats)
    loo_metrics, loo_predictions = phase91.leave_one_out_predictions(X, y, args.seed + 900)
    observed_auroc = float(phase91.np.mean([row["auroc"] for row in logistic_metrics]))
    empirical_p, null_auroc = phase91.permutation_control(X, y, args.seed + 1200, observed_auroc, args.n_permutations)
    loo_probabilities = phase91.np.asarray([float(row["predicted_probability"]) for row in loo_predictions], dtype=float)
    ci_low, ci_high = phase91.bootstrap_ci(y, loo_probabilities, args.seed + 1400, args.n_bootstrap)
    metric_rows = [
        phase91.summarize_metrics("logistic_regression", "repeated_stratified_split", logistic_metrics, len(rows), counts, len(features), match_rate, empirical_p, null_auroc, ci_low, ci_high),
        phase91.summarize_metrics("logistic_regression", "leave_one_out", [loo_metrics], len(rows), counts, len(features), match_rate, empirical_p, null_auroc, ci_low, ci_high),
        phase91.summarize_metrics("random_forest_baseline", "repeated_stratified_split", rf_metrics, len(rows), counts, len(features), match_rate),
    ]
    phase91.write_tsv(args.metrics_output, metric_rows, metric_fields)
    phase91.write_tsv(args.predictions_output, [*split_predictions, *rf_predictions, *loo_predictions], ["model", "validation_mode", "split_id", "row_index", "true_label", "predicted_probability", "predicted_label"])
    phase91.write_tsv(args.importance_output, phase91.feature_importance(features, X, y, args.seed + 1600), ["rank", "feature", "importance", "model"])
    print(f"Wrote {args.metrics_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
