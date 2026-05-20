# NeuroFate Tiny Demo

This directory contains bundled synthetic toy data for CLI smoke testing.

Files:

- `tiny_metadata.tsv`: pseudo-cell metadata with donor labels.
- `tiny_gene_panel.tsv`: small gene panel used by the demo.
- `tiny_sparse_expression.tsv` and `.gz`: toy long-form expression rows.
- `tiny_expected_output_summary.tsv`: expected output artifacts.

Run:

```bash
neurofate run-demo
```

Expected outputs are written under `results/demo/`:

- `demo_donor_feature_table.tsv`
- `demo_model_metrics.tsv`
- `axis_scores.tsv`
- `neurofate_risk_scores.tsv`
- `risk_score_report.md`
- `demo_report.md`

The tiny demo is not biological evidence. It only verifies installation, packaged resources, and CLI workflow plumbing.

