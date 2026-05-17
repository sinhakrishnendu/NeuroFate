# PyPI Release Guide

NeuroFate is prepared for PyPI distribution as the `neurofate` package.

## Package Entry Points

After installation, users should be able to run:

```bash
neurofate check-system
neurofate doctor
neurofate run-demo
```

The tiny demo is bundled as package data under `neurofate.resources.tiny_demo`, so it works from an installed wheel and does not require the GitHub repository layout.

## Manual Release Commands

```bash
python -m pip install -e ".[dev]"
python -m py_compile scripts/*.py neurofate/*.py
python -m pytest tests/
python -m build
python -m twine check dist/*
```

Upload to TestPyPI first:

```bash
python -m twine upload --repository testpypi dist/*
```

Then upload to PyPI only after TestPyPI installation succeeds:

```bash
python -m twine upload dist/*
```

## GitHub Release Pairing

For each PyPI release:

1. Update `pyproject.toml`, `CHANGELOG.md`, `CITATION.cff`, and `codemeta.json` to the same version.
2. Commit the release changes.
3. Create a Git tag such as `v0.1.0`.
4. Push the tag to GitHub.
5. Upload to PyPI.

## Data Policy

The PyPI package must not include real SEA-AD, Mathys, ROSMAP, STRING, HMDB, KEGG, Synapse, or other external datasets. Only synthetic tiny demo data are bundled.
