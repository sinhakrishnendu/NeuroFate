# Contributing

Thank you for helping improve NeuroFate.

## Development Setup

```bash
python -m pip install -e ".[dev]"
python -m py_compile scripts/*.py neurofate/*.py
python -m pytest tests/
neurofate run-demo
```

## Safety Rules

- Do not commit external datasets or generated large files.
- Do not add code paths that download data automatically.
- Keep H5AD expression-matrix access guarded, sparse, and explicit.
- Keep external validation claims proportional to cohort size and feature overlap.

## Pull Requests

Include a short summary, tests run, and any known limitations. For methods changes, update documentation and relevant registry templates.
