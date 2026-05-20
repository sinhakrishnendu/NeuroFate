"""Command-line interface for the NeuroFate research software package."""

from __future__ import annotations

import argparse
import importlib.util
import platform
import shlex
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_AXIS_REGISTRY = PACKAGE_ROOT / "resources" / "neurofate_axis_registry.tsv"
DEFAULT_ALIAS_TABLE = PACKAGE_ROOT / "resources" / "neurofate_axis_gene_aliases.tsv"

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
    "external-triage": {
        "help": "Triage Phase 15 external validation candidates without internet access.",
        "command": "python scripts/69_triage_external_validation_candidates.py --registry metadata/phase15_external_validation_candidates.tsv --output results/reports/phase15_external_dataset_triage.tsv --summary results/reports/phase15_external_dataset_priority_summary.md",
        "expected": ["metadata/phase15_external_validation_candidates.tsv"],
        "hint": "This is a lightweight registry-only step; it does not download or inspect dataset contents.",
    },
    "inspect-external": {
        "help": "Inventory a manually acquired external cohort directory by filename and size only.",
        "command": "python scripts/70_inspect_external_dataset_files.py --dataset-id gse243639_pd_snpc --input-dir data/raw/external/gse243639_pd_snpc --output-summary results/reports/phase15_gse243639_pd_snpc_file_inventory.tsv --format-output results/reports/phase15_gse243639_pd_snpc_format_recommendation.tsv --log-file results/logs/70_inspect_external_dataset_files.log",
        "expected": ["data/raw/external/gse243639_pd_snpc"],
        "hint": "Manually place external files under data/raw/external/<dataset_id>/, then run the underlying script with the desired --dataset-id and --input-dir.",
    },
    "plan-external-extraction": {
        "help": "Create a manual external sparse-extraction plan; no extraction is executed.",
        "command": "python scripts/73_prepare_external_sparse_extraction_plan.py --dataset-id gse243639_pd_snpc --format h5ad_csr --input-matrix data/raw/external/gse243639_pd_snpc/COUNTS_OR_CONTAINER_FILE --metadata-file data/raw/external/gse243639_pd_snpc/METADATA_FILE --feature-file data/raw/external/gse243639_pd_snpc/FEATURE_FILE --panel metadata/target_gene_panel_v1.tsv --output-plan results/tables/phase15_gse243639_pd_snpc_sparse_extraction_plan.tsv --manual-script-output results/logs/manual_phase15_gse243639_pd_snpc_extraction_template.sh --log-file results/logs/73_prepare_external_sparse_extraction_plan.log",
        "expected": ["metadata/target_gene_panel_v1.tsv", "data/raw/external/gse243639_pd_snpc/COUNTS_OR_CONTAINER_FILE"],
        "hint": "Use the underlying script with real local file paths after manual acquisition and format inspection.",
    },
    "validate-multi-external": {
        "help": "Show or run donor-level multi-external validation from prepared feature tables only.",
        "command": "python scripts/75_run_multi_external_validation.py --sea-ad-features results/tables/phase5_donor_feature_table.tsv --tables-dir results/tables --reports-dir results/reports --log-file results/logs/75_run_multi_external_validation.log --dry-run",
        "expected": [],
        "hint": "Provide --external-feature-table dataset_id=path directly to scripts/75_run_multi_external_validation.py only after donor/sample-level external feature tables exist.",
    },
    "external-report": {
        "help": "Generate the Phase 15 external validation expansion report from existing planning outputs.",
        "command": "python scripts/77_generate_phase15_external_validation_report.py --output results/reports/phase15_external_validation_report.md",
        "expected": ["metadata/phase15_external_validation_candidates.tsv"],
        "hint": "Run `neurofate external-triage --run` first for a complete readiness section.",
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
    for package in ["numpy", "pandas", "scipy", "sklearn", "yaml", "matplotlib", "torch"]:
        print(f"{package}: {package_status(package)}")
    print("research_use_only: true")
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
        "metadata/phase15_external_validation_candidates.tsv",
        "docs/external_validation_expansion.md",
        "scripts/69_triage_external_validation_candidates.py",
        "scripts/70_inspect_external_dataset_files.py",
        "scripts/71_inspect_external_metadata_safe.py",
        "scripts/72_plan_external_target_gene_overlap.py",
        "scripts/73_prepare_external_sparse_extraction_plan.py",
        "scripts/74_build_external_feature_table_generic.py",
        "scripts/75_run_multi_external_validation.py",
        "scripts/76_generate_phase15_external_validation_figures.py",
        "scripts/77_generate_phase15_external_validation_report.py",
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


def build_axis_scores(args: argparse.Namespace) -> int:
    from neurofate.axis import build_axis_score_tables

    try:
        outputs = build_axis_score_tables(
            expression=Path(args.expression),
            metadata=Path(args.metadata),
            axis_registry=Path(args.axis_registry),
            outdir=Path(args.outdir),
            sample_id_column=args.sample_id_column,
            endpoint_column=args.endpoint_column,
            positive_class=args.positive_class,
            negative_class=args.negative_class,
            orientation=args.orientation,
            gene_column=args.gene_column,
        )
    except Exception as exc:  # noqa: BLE001 - CLI should print clean user-facing failures.
        print(f"NeuroFate axis scoring failed: {exc}")
        return 2
    for label, path in outputs.items():
        print(f"Wrote {label}: {path}")
    return 0


def score_risk(args: argparse.Namespace) -> int:
    from neurofate.axis import RESEARCH_USE_NOTICE, score_research_risk

    try:
        outputs = score_research_risk(Path(args.axis_scores), Path(args.outdir))
    except Exception as exc:  # noqa: BLE001 - CLI should print clean user-facing failures.
        print(f"NeuroFate research risk scoring failed: {exc}")
        return 2
    print(RESEARCH_USE_NOTICE)
    for label, path in outputs.items():
        print(f"Wrote {label}: {path}")
    return 0


def ingest_data(args: argparse.Namespace) -> int:
    from neurofate.ingest import IngestConfig, run_ingest

    try:
        result = run_ingest(
            IngestConfig(
                expression=Path(args.expression),
                metadata=Path(args.metadata),
                outdir=Path(args.outdir),
                sample_id_column=args.sample_id_column,
                endpoint_column=args.endpoint_column,
                positive_class=args.positive_class,
                negative_class=args.negative_class,
                gene_id_column=args.gene_id_column,
                orientation=args.orientation,
                gene_map=Path(args.gene_map) if args.gene_map else None,
                axis_registry=Path(args.axis_registry),
                alias_table=Path(args.alias_table) if args.alias_table else None,
                assist=args.assist,
                min_axis_genes=args.min_axis_genes,
            )
        )
    except Exception as exc:  # noqa: BLE001 - CLI should print clean user-facing failures.
        print(f"NeuroFate ingest failed: {exc}")
        return 2
    for label, path in result.__dict__.items():
        print(f"Wrote {label}: {path}")
    return 0


def run_full_workflow(args: argparse.Namespace) -> int:
    from neurofate.ingest import IngestConfig, run_complete_workflow

    try:
        outputs = run_complete_workflow(
            IngestConfig(
                expression=Path(args.expression),
                metadata=Path(args.metadata),
                outdir=Path(args.outdir),
                sample_id_column=args.sample_id_column,
                endpoint_column=args.endpoint_column,
                positive_class=args.positive_class,
                negative_class=args.negative_class,
                gene_id_column=args.gene_id_column,
                orientation=args.orientation,
                gene_map=Path(args.gene_map) if args.gene_map else None,
                axis_registry=Path(args.axis_registry),
                alias_table=Path(args.alias_table) if args.alias_table else None,
                assist=args.assist,
                min_axis_genes=args.min_axis_genes,
            )
        )
    except Exception as exc:  # noqa: BLE001 - CLI should print clean user-facing failures.
        print(f"NeuroFate run failed: {exc}")
        return 2
    for label, path in outputs.items():
        print(f"Wrote {label}: {path}")
    return 0


def adapt_endpoint(args: argparse.Namespace) -> int:
    from neurofate.adapters import adapt_endpoint_metadata

    try:
        outputs = adapt_endpoint_metadata(
            metadata=Path(args.metadata),
            outdir=Path(args.outdir),
            task=args.task,
            endpoint_column=args.endpoint_column,
        )
    except Exception as exc:  # noqa: BLE001 - CLI should print clean user-facing failures.
        print(f"NeuroFate endpoint adaptation failed: {exc}")
        return 2
    for label, path in outputs.__dict__.items():
        print(f"Wrote {label}: {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="neurofate",
        description=(
            "NeuroFate command-line research software for endpoint-locked "
            "transcriptomic neurodegeneration-axis scoring."
        ),
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

    axis_parser = subparsers.add_parser(
        "build-axis-scores",
        help="Build donor/sample-level NeuroFate axis scores from a compact expression matrix.",
    )
    axis_parser.add_argument("--expression", required=True, help="Sample-level or gene-row TSV/TSV.GZ.")
    axis_parser.add_argument("--metadata", required=True, help="Sample metadata TSV.")
    axis_parser.add_argument("--axis-registry", required=True, help="NeuroFate axis registry TSV.")
    axis_parser.add_argument("--sample-id-column", required=True, help="Sample identifier column.")
    axis_parser.add_argument("--endpoint-column", required=True, help="Endpoint column in metadata.")
    axis_parser.add_argument("--positive-class", required=True, help="Positive endpoint class.")
    axis_parser.add_argument("--negative-class", required=True, help="Negative endpoint class.")
    axis_parser.add_argument("--outdir", required=True, help="Output directory.")
    axis_parser.add_argument(
        "--orientation",
        choices=["auto", "genes_rows", "samples_rows"],
        default="auto",
        help="Expression table orientation.",
    )
    axis_parser.add_argument(
        "--gene-column",
        default="gene_symbol",
        help="Gene identifier column for genes_rows orientation.",
    )
    axis_parser.set_defaults(func=build_axis_scores)

    risk_parser = subparsers.add_parser(
        "score-risk",
        help="Compute a research-use NeuroFate risk score from axis scores.",
    )
    risk_parser.add_argument("--axis-scores", required=True, help="Axis score TSV.")
    risk_parser.add_argument("--outdir", required=True, help="Output directory.")
    risk_parser.set_defaults(func=score_risk)

    ingest_parser = subparsers.add_parser(
        "ingest",
        help="Inspect and standardize user transcriptomic tables for NeuroFate scoring.",
    )
    ingest_parser.add_argument("--expression", required=True, help="User expression table: CSV/TSV/TXT/GZ.")
    ingest_parser.add_argument("--metadata", required=True, help="User sample metadata table: CSV/TSV/TXT/GZ.")
    ingest_parser.add_argument("--outdir", required=True, help="Output directory for standardized files.")
    ingest_parser.add_argument("--sample-id-column", default="auto", help="Metadata sample ID column or auto.")
    ingest_parser.add_argument("--endpoint-column", default="auto", help="Endpoint/group column or auto.")
    ingest_parser.add_argument("--positive-class", default="auto", help="Positive endpoint class or auto.")
    ingest_parser.add_argument("--negative-class", default="auto", help="Negative endpoint class or auto.")
    ingest_parser.add_argument("--gene-id-column", default="auto", help="Gene/probe ID column or auto.")
    ingest_parser.add_argument(
        "--orientation",
        choices=["auto", "genes_rows", "samples_rows", "long"],
        default="auto",
        help="Expression orientation.",
    )
    ingest_parser.add_argument("--gene-map", default="", help="Optional probe/gene mapping table.")
    ingest_parser.add_argument(
        "--axis-registry",
        default=str(DEFAULT_AXIS_REGISTRY),
        help="NeuroFate axis registry TSV.",
    )
    ingest_parser.add_argument(
        "--alias-table",
        default=str(DEFAULT_ALIAS_TABLE),
        help="Optional gene symbol/Ensembl alias table.",
    )
    ingest_parser.add_argument("--assist", action="store_true", help="Record assisted-mode intent for ambiguous runs.")
    ingest_parser.add_argument("--min-axis-genes", type=int, default=3, help="Minimum retained axis genes.")
    ingest_parser.set_defaults(func=ingest_data)

    run_parser = subparsers.add_parser(
        "run",
        help="Run ingest, axis scoring, research-use risk scoring, and report generation.",
    )
    run_parser.add_argument("--expression", required=True, help="User expression table: CSV/TSV/TXT/GZ.")
    run_parser.add_argument("--metadata", required=True, help="User sample metadata table: CSV/TSV/TXT/GZ.")
    run_parser.add_argument("--outdir", required=True, help="Output directory for the complete run.")
    run_parser.add_argument("--sample-id-column", default="auto", help="Metadata sample ID column or auto.")
    run_parser.add_argument("--endpoint-column", default="auto", help="Endpoint/group column or auto.")
    run_parser.add_argument("--positive-class", default="auto", help="Positive endpoint class or auto.")
    run_parser.add_argument("--negative-class", default="auto", help="Negative endpoint class or auto.")
    run_parser.add_argument("--gene-id-column", default="auto", help="Gene/probe ID column or auto.")
    run_parser.add_argument(
        "--orientation",
        choices=["auto", "genes_rows", "samples_rows", "long"],
        default="auto",
        help="Expression orientation.",
    )
    run_parser.add_argument("--gene-map", default="", help="Optional probe/gene mapping table.")
    run_parser.add_argument(
        "--axis-registry",
        default=str(DEFAULT_AXIS_REGISTRY),
        help="NeuroFate axis registry TSV.",
    )
    run_parser.add_argument(
        "--alias-table",
        default=str(DEFAULT_ALIAS_TABLE),
        help="Optional gene symbol/Ensembl alias table.",
    )
    run_parser.add_argument("--assist", action="store_true", help="Record assisted-mode intent for ambiguous runs.")
    run_parser.add_argument("--min-axis-genes", type=int, default=3, help="Minimum retained axis genes.")
    run_parser.set_defaults(func=run_full_workflow)

    adapter_parser = subparsers.add_parser(
        "adapt-endpoint",
        help="Create explicit endpoint-label aliases for validation scripts.",
    )
    adapter_parser.add_argument("--metadata", required=True, help="Standardized metadata TSV.")
    adapter_parser.add_argument("--outdir", required=True, help="Output directory for adapted metadata.")
    adapter_parser.add_argument(
        "--endpoint-column",
        default="auto",
        help="Binary endpoint label column, usually label__endpoint.",
    )
    adapter_parser.add_argument(
        "--task",
        choices=["generic", "pd_vs_control", "ad_vs_control"],
        default="generic",
        help="Task-specific alias set to create.",
    )
    adapter_parser.set_defaults(func=adapt_endpoint)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
