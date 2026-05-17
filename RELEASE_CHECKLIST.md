# Release Checklist

- [ ] Confirm `python -m pip install -e .` works in a clean Python 3.11 environment.
- [ ] Run `python -m py_compile scripts/*.py neurofate/*.py`.
- [ ] Run `python -m pytest tests/`.
- [ ] Run `neurofate check-system`.
- [ ] Run `neurofate doctor`.
- [ ] Run `neurofate run-demo`.
- [ ] Confirm no bundled external datasets are present.
- [ ] Regenerate reports and reproducibility manifest.
- [ ] Run the no-overclaiming audit before manuscript or release text is published.
- [ ] Verify citation, CodeMeta, changelog, and data-use notes.
- [ ] Build PyPI artifacts with `python -m build`; `dist/` should contain only `.whl` and `.tar.gz` package files.
- [ ] Build reviewer/source ZIP artifacts with `python scripts/66_create_source_release_package.py` and `python scripts/67_create_results_review_package.py`; these write to `release_artifacts/`.
