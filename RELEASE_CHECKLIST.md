# NeuroFate 0.3.0 Release Checklist

Use this checklist before publishing NeuroFate to TestPyPI, PyPI, GitHub Releases, or Zenodo.

## Version and Metadata

- [x] `pyproject.toml` version is `0.3.0`.
- [x] `neurofate/__init__.py` version is `0.3.0`.
- [x] `CHANGELOG.md` includes `0.3.0`.
- [x] `CITATION.cff` version is `0.3.0`.
- [x] `codemeta.json` version is `0.3.0`.
- [x] README states release-candidate version `0.3.0`.
- [x] Bioinformatics manuscript states version `0.3.0`.

## Documentation

- [x] README is a full user/reviewer manual.
- [x] `docs/cli_reference.md` documents stable public commands.
- [x] `docs/input_output_schema.md` documents supported inputs and outputs.
- [x] `docs/real_world_geo_smoke_test_gse20141.md` documents the real GEO smoke test.
- [x] `docs/reproducible_commands.md` lists install, demo, smoke-test, build, test, and manuscript commands.
- [x] Research-use-only language is present in README, reports, and manuscript.

## Public CLI

- [x] `neurofate check-system`
- [x] `neurofate doctor`
- [x] `neurofate run-demo`
- [x] `neurofate ingest --help`
- [x] `neurofate build-axis-scores --help`
- [x] `neurofate score-risk --help`
- [x] `neurofate run --help`
- [x] `neurofate adapt-endpoint --help`

## Tests

- [x] Python compilation passes.
- [x] Public CLI tests pass.
- [x] Ingestion tests pass.
- [x] Endpoint adapter tests pass.
- [x] Research-use-only report tests pass.
- [x] Bioinformatics manuscript claim tests pass.
- [x] Final no-overclaiming audit has zero high-severity flags.

## Real-World Smoke Test

- [x] GSE20141 series matrix downloaded manually for local smoke test.
- [x] GPL570 annotation downloaded manually for local smoke test.
- [x] `neurofate run` completed on the public GEO series matrix.
- [x] 18/18 samples matched.
- [x] 10 PD and 8 control labels retained.
- [x] 29/30 NeuroFate genes retained.
- [x] 10/10 axes scored.
- [x] Research-use risk scores generated.

## Package Build

- [x] `python -m build --outdir dist_final`
- [x] `python -m twine check dist_final/*`
- [x] `dist/` should contain only `.whl` and `.tar.gz` release artifacts when using the default build path.
- [ ] TestPyPI upload dry run.
- [ ] Test install from TestPyPI in a clean environment.
- [ ] PyPI release.

## GitHub and Archive

- [ ] Confirm GitHub repository is public.
- [ ] Create release tag `v0.3.0`.
- [ ] Publish GitHub release notes.
- [ ] Archive release on Zenodo.
- [ ] Add Zenodo DOI to README and `CITATION.cff`.

## Manuscript

- [x] Full Bioinformatics methods/software paper exists.
- [x] Author order and correspondence are updated.
- [x] Cover letter author/correspondence details are updated.
- [x] Figures compile.
- [x] Tables compile.
- [x] Bibliography compiles.
- [x] Manuscript PDF compiles.
- [ ] Final funding statement.
- [ ] Final author contributions.
- [ ] Final journal formatting after editorial decision.

## Data and Code Availability

- [x] GitHub URL is listed.
- [x] Public data accessions are listed in manuscript.
- [x] Large external data are not bundled.
- [x] Raw protected data are not redistributed.
- [x] Research-use-only scope is explicit.
