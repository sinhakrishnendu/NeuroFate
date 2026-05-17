# Troubleshooting

Common checks:

- Run `neurofate doctor`.
- Confirm files are under the expected `data/` and `results/` paths.
- Inspect logs under `results/logs/`.
- Use `scripts/53_validate_neurofate_outputs.py` to find missing outputs.

If a command would touch large data, review memory settings and run it manually.
