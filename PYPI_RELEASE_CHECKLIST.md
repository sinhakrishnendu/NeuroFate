# PyPI Release Checklist

Use these commands manually from a clean checkout. Do not upload until the license, repository URL, version, and release notes are final.

## 1. Clean Local Checks

```bash
python -m pip install -e ".[dev]"
python -m py_compile scripts/*.py neurofate/*.py
python -m pytest tests/
python -m neurofate run-demo
python -m neurofate doctor
```

## 2. Build Source And Wheel Distributions

```bash
python -m build
```

Expected outputs:

- `dist/neurofate-0.1.0.tar.gz`
- `dist/neurofate-0.1.0-py3-none-any.whl`

## 3. Validate Distribution Metadata

```bash
python -m twine check dist/*
```

## 4. TestPyPI Dry Run

```bash
python -m twine upload --repository testpypi dist/*
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple neurofate
neurofate run-demo
```

## 5. PyPI Upload

```bash
python -m twine upload dist/*
```

## 6. Post-Release Verification

```bash
python -m pip install neurofate
neurofate check-system
neurofate run-demo
```

## Notes

- Do not bundle SEA-AD, Mathys, ROSMAP, STRING, or other external datasets.
- Ensure `CHANGELOG.md`, `CITATION.cff`, `codemeta.json`, and README release notes match the version.
- Ensure GitHub tags and PyPI versions match exactly.
