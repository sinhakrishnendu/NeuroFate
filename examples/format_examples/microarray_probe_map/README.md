# Microarray Probe-Map Example

Run:

```bash
neurofate run \
  --expression examples/format_examples/microarray_probe_map/expression.tsv \
  --metadata examples/format_examples/microarray_probe_map/metadata.tsv \
  --gene-map examples/format_examples/microarray_probe_map/probe_map.tsv \
  --outdir results/examples/microarray_probe_map
```

This example demonstrates compact probe-to-gene mapping without writing a
genome-wide expression table.
