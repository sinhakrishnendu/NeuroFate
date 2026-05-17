# Input Formats

NeuroFate accepts sparse or tabular inputs:

- SEA-AD metadata TSVs.
- Sparse long gene-panel expression TSVs.
- Donor-level feature tables.
- Mathys 2019 CSV count and covariate files.

The preferred sparse gene-panel format is:

```text
cell_id or row_index    gene_symbol    expression_value
```

Full single-cell matrices are never required for model training.
