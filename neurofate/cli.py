"""Command-line interface for the NeuroFate platform."""

from __future__ import annotations

import argparse
import importlib.util
import platform
import shlex
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

WORKFLOWS = {
    "inspect-sea-ad": {
        "help": "Inspect SEA-AD metadata only; never reads expression matrix arrays.",
        "command": "python scripts/11_extract_sea_ad_metadata_only.py --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad --outdir data/interim/sea_ad --tables-dir results/tables --log-file results/logs/11_extract_sea_ad_metadata_only.log",
        "expected": ["data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"],
        "hint": "Place the SEA-AD H5AD at data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad, then rerun with --run.",
    },
    "extract-metadata": {
        "help": "Extract metadata-only SEA-AD obs/var tables with H5AD X guards.",
        "command": "python scripts/11_extract_sea_ad_metadata_only.py --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad --outdir data/interim/sea_ad --tables-dir results/tables --log-file results/logs/11_extract_sea_ad_metadata_only.log",
        "expected": ["data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad"],
        "hint": "Run the guarded SEA-AD manual acquisition first; no download is performed by NeuroFate.",
    },
    "extract-panel": {
        "help": "Extract only the configured target gene panel with chunked sparse safeguards.",
        "command": "python scripts/15_sparse_gene_extraction_safe.py --input data/raw/sea_ad/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad --gene-panel metadata/target_gene_panel_v1.tsv --out data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz --dry-run",
        "expected": ["metadata/target_gene_panel_v1.tsv"],
        "hint": "Use --run only after the local SEA-AD file is present and the gene panel plan has passed.",
    },
    "summarize": {
        "help": "Summarize sparse gene-panel expression and decoded metadata tables.",
        "command": "python scripts/16_compute_sparse_expression_statistics.py --expression data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz --metadata data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv --tables-dir results/tables --log-file results/logs/16_compute_sparse_expression_statistics.log",
        "expected": [
            "data/interim/sea_ad/sparse_gene_panel_expression.tsv.gz",
            "data/interim/sea_ad/sea_ad_obs_metadata_decoded.tsv",
        ],
        "hint": "Run metadata decoding and sparse panel extraction before summary generation.",
    },
    "train-baseline": {
        "help": "Run donor-level classical baseline models from Phase 5 tables.",
        "command": "python scripts/23_run_phase5_models.py --features results/tables/phase5_donor_feature_table.tsv --tables-dir results/tables --figures-dir results/figures",
        "expected": ["results/tables/phase5_donor_feature_table.tsv"],
        "hint": "Build the donor feature table first with scripts/22_build_donor_feature_table.py.",
    },
    "train-mps": {
        "help": "Run the small donor-level PyTorch MPS model when explicitly requested.",
        "command": "python scripts/27_train_neurofate_mps_model.py --features results/tables/phase5_donor_feature_table.tsv --config configs/neurofate_mps_model_config.yaml --task dementia_vs_reference",
        "expected": ["results/tables/phase5_donor_feature_table.tsv", "configs/neurofate_mps_model_config.yaml"],
        "hint": "This is a model-training command; run only deliberately on donor-level features.",
    },
    "validate-external": {
        "help": "Run Mathys external feasibility validation from prepared donor-level tables.",
        "command": "python scripts/43_run_mathys_external_validation.py --sea-ad-features results/tables/phase5_donor_feature_table.tsv --mathys-features results/tables/mathys_2019_phase5_donor_feature_table.tsv",
        "expected": [
            "results/tables/phase5_donor_feature_table.tsv",
            "results/tables/mathys_2019_phase5_donor_feature_table.tsv",
        ],
        "hint": "Prepare Mathys CSV feature tables first; external data are not bundled.",
    },
    "make-report": {
        "help": "Generate Markdown and HTML reports from existing results only.",
        "command": "python scripts/51_generate_end_user_report.py --tables-dir results/tables --reports-dir results/reports",
        "expected": [],
        "hint": "Reports can be generated any time; missing result tables are shown as unavailable.",
    },
    "audit-leakage": {
        "help": "Audit donor-level feature tables for label and identifier leakage.",
        "command": "python scripts/57_audit_feature_leakage.py --input results/tables/phase5_donor_feature_table.tsv --output results/reports/feature_leakage_audit.tsv",
        "expected": ["results/tables/phase5_donor_feature_table.tsv", "configs/benchmark_config.yaml"],
        "hint": "Run Phase 5 donor feature table generation first, then rerun with --run.",
    },
    "benchmark": {
        "help": "Run configurable repeated donor-level baseline benchmarks.",
        "command": "python scripts/58_run_repeated_baseline_benchmarks.py --features results/tables/phase5_donor_feature_table.tsv --config configs/benchmark_config.yaml",
        "expected": ["results/tables/phase5_donor_feature_table.tsv", "configs/benchmark_config.yaml"],
        "hint": "Use --run only when ready for the configured repeated benchmark loops.",
    },
    "benchmark-report": {
        "help": "Generate Phase 12 benchmark uncertainty and evidence reports from existing tables.",
        "command": "python scripts/61_generate_benchmark_uncertainty_report.py --summary results/tables/phase12_repeated_benchmark_summary.tsv --pvalues results/tables/phase12_empirical_pvalues.tsv --ablation results/tables/phase12_feature_group_importance.tsv --evidence results/reports/evidence_strength_matrix.tsv --output results/reports/phase12_benchmark_uncertainty_report.md",
        "expected": ["configs/benchmark_config.yaml"],
        "hint": "Run benchmark, permutation, ablation, and evidence-classification commands first for complete reports.",
    },
}


def package_status(name: str) -> str:
    return "available" if importlib.util.find_spec(name) is not None else "missing"


def missing_paths(paths: list[str]) -> list[str]:
    return [path for path in paths if not (PROJECT_ROOT / path).exists()]


def run_command(command: str, run: bool, expected: list[str] | None = None, hint: str = "") -> int:
    if not run:
        print("DRY RUN command template:")
        print(command)
        print("Pass --run to execute this lightweight wrapper deliberately.")
        if hint:
            print(f"Next step: {hint}")
        return 0
    expected_paths = list(expected or [])
    command_parts = shlex.split(command)
    if len(command_parts) > 1 and command_parts[1].startswith("scripts/"):
        expected_paths.append(command_parts[1])
    missing = missing_paths(expected_paths)
    if missing:
        print("NeuroFate cannot run this command because expected input is missing:")
        for path in missing:
            print(f"- {path}")
        if hint:
            print(f"Next step: {hint}")
        if any(path.startswith("scripts/") for path in missing):
            print("This workflow wrapper requires the GitHub repository checkout layout.")
        return 2
    try:
        completed = subprocess.run(shlex.split(command), cwd=PROJECT_ROOT, check=False)
    except FileNotFoundError as exc:
        print(f"NeuroFate could not start the command: {exc}")
        print("Check that Python is available and that you are running from an installed NeuroFate environment.")
        return 127
    return int(completed.returncode)


def check_system(_: argparse.Namespace) -> int:
    print(f"NeuroFate root: {PROJECT_ROOT}")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {platform.platform()}")
    for package in ["numpy", "scipy", "sklearn", "h5py", "yaml", "matplotlib", "torch"]:
        print(f"{package}: {package_status(package)}")
    return 0


def doctor(_: argparse.Namespace) -> int:
    from neurofate.demo import demo_resource

    package_missing = [
        name
        for name in ["tiny_metadata.tsv", "tiny_gene_panel.tsv", "tiny_sparse_expression.tsv"]
        if not demo_resource(name).is_file()
    ]
    if package_missing:
        print("NeuroFate doctor found missing packaged resources:")
        for name in package_missing:
            print(f"- neurofate.resources.tiny_demo/{name}")
        return 1

    repo_layout = (PROJECT_ROOT / "pyproject.toml").exists() and (PROJECT_ROOT / "scripts").exists()
    if not repo_layout:
        print("NeuroFate doctor: installed package resources are present.")
        print("Repository workflow scripts are available only in the GitHub checkout.")
        return 0

    required = [
        "README.md",
        "pyproject.toml",
        "LICENSE",
        "configs/project_config.yaml",
        "configs/benchmark_config.yaml",
        "configs/templates/sea_ad_minimal.yaml",
        "docs/quickstart.md",
        "docs/pypi_release.md",
        "docs/benchmarking_and_validation.md",
        "examples/tiny_demo/tiny_sparse_expression.tsv.gz",
        "CITATION.cff",
        "codemeta.json",
        "CHANGELOG.md",
        "PYPI_RELEASE_CHECKLIST.md",
        "RELEASE_CHECKLIST.md",
        "CONTRIBUTING.md",
        "CODE_OF_CONDUCT.md",
        "metadata/target_gene_panel_v1.tsv",
        "scripts/51_generate_end_user_report.py",
        "scripts/52_generate_reproducibility_manifest.py",
        "scripts/53_validate_neurofate_outputs.py",
        "scripts/54_no_overclaiming_audit.py",
        "scripts/55_run_tiny_demo.py",
        "scripts/56_inventory_outputs.py",
        "scripts/57_audit_feature_leakage.py",
        "scripts/58_run_repeated_baseline_benchmarks.py",
        "scripts/59_run_label_permutation_controls.py",
        "scripts/60_run_feature_ablation.py",
        "scripts/61_generate_benchmark_uncertainty_report.py",
        "scripts/62_generate_phase12_benchmark_figures.py",
        "scripts/63_classify_evidence_strength.py",
    ]
    missing = [path for path in required if not (PROJECT_ROOT / path).exists()]
    if missing:
        print("NeuroFate doctor found missing files:")
        for path in missing:
            print(f"- {path}")
        return 1
    print("NeuroFate doctor: core platform files are present.")
    return 0


def template_subcommand(name: str):
    def handler(args: argparse.Namespace) -> int:
        workflow = WORKFLOWS[name]
        return run_command(
            workflow["command"],
            args.run,
            expected=workflow["expected"],
            hint=workflow["hint"],
        )

    return handler


def run_demo(_: argparse.Namespace) -> int:
    from neurofate.demo import build_demo_outputs

    try:
        build_demo_outputs(outdir=Path.cwd() / "results/demo")
    except FileNotFoundError as exc:
        print(f"NeuroFate demo could not start: {exc}")
        print("Reinstall NeuroFate with package data, then rerun `neurofate run-demo`.")
        return 2
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="neurofate",
        description="NeuroFate donor-level neurodegeneration systems biology platform.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check-system", help="Report local Python/platform readiness.")
    check.set_defaults(func=check_system)

    for name, workflow in WORKFLOWS.items():
        subcommand = subparsers.add_parser(name, help=workflow["help"])
        subcommand.add_argument("--run", action="store_true", help="Execute the wrapper command.")
        subcommand.set_defaults(func=template_subcommand(name))

    demo_parser = subparsers.add_parser("run-demo", help="Run the bundled tiny synthetic demo.")
    demo_parser.set_defaults(func=run_demo)

    doctor_parser = subparsers.add_parser("doctor", help="Validate core platform file layout.")
    doctor_parser.set_defaults(func=doctor)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
