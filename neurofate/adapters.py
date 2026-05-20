"""Adapters between public NeuroFate CLI outputs and legacy validation tables.

The public ingestion workflow writes a deliberately generic endpoint label,
``label__endpoint``.  Older research-validation scripts often expect cohort-
specific names such as ``label__pd_vs_control``.  This module creates explicit
aliases without changing the underlying 0/1 label semantics.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from neurofate.axis import RESEARCH_USE_NOTICE


TASK_ALIASES = {
    "pd_vs_control": ["label__pd_vs_control"],
    "ad_vs_control": ["label__ad_vs_control"],
    "generic": [],
}
GENERIC_ALIASES = ["endpoint_label", "label"]
DEFAULT_LABEL_COLUMNS = (
    "label__endpoint",
    "endpoint_label",
    "label",
    "label__pd_vs_control",
    "label__ad_vs_control",
)


@dataclass
class EndpointAdapterResult:
    adapted_metadata: Path
    endpoint_adapter_report: Path
    endpoint_aliases: Path


def _read_table(input_table: Path | pd.DataFrame) -> pd.DataFrame:
    if isinstance(input_table, pd.DataFrame):
        return input_table.copy()
    return pd.read_csv(input_table, sep="\t", dtype=str)


def _label_key(value: object) -> str:
    return str(value).strip().casefold()


def _normalize_binary_label(value: object) -> str | None:
    key = _label_key(value)
    if key in {"1", "1.0", "true", "case", "positive"}:
        return "1"
    if key in {"0", "0.0", "false", "control", "negative"}:
        return "0"
    return None


def normalize_endpoint_column(
    input_table: Path | pd.DataFrame,
    endpoint_column: str | None = None,
) -> tuple[pd.DataFrame, str]:
    """Return a copy with a validated binary ``label__endpoint`` column.

    The function never infers biological meaning from free-text disease labels.
    It only accepts an explicit or already standardized binary label column.
    """

    table = _read_table(input_table)
    selected = endpoint_column
    if selected in {"", "auto"}:
        selected = None
    if selected is None:
        for candidate in DEFAULT_LABEL_COLUMNS:
            if candidate in table.columns:
                selected = candidate
                break
    if selected is None or selected not in table.columns:
        raise ValueError(
            "Could not find a binary endpoint label column. Pass --endpoint-column explicitly; "
            f"available columns are {list(table.columns)}."
        )

    normalized = [_normalize_binary_label(value) for value in table[selected]]
    invalid = [str(value) for value, label in zip(table[selected], normalized, strict=False) if label is None]
    if invalid:
        examples = "; ".join(invalid[:8])
        raise ValueError(
            f"Endpoint column {selected!r} contains values that are not unambiguous binary labels: {examples}"
        )
    table = table.copy()
    table["label__endpoint"] = normalized
    return table, selected


def map_public_label_to_internal(label_column: str, task: str = "generic") -> list[str]:
    """Return explicit alias columns for a task-specific validation script."""

    if task not in TASK_ALIASES:
        raise ValueError(f"Unsupported task {task!r}; expected one of {sorted(TASK_ALIASES)}")
    aliases = ["label__endpoint", *TASK_ALIASES[task], *GENERIC_ALIASES]
    return list(dict.fromkeys(alias for alias in aliases if alias != label_column or alias == "label__endpoint"))


def ensure_label_endpoint_aliases(
    table: Path | pd.DataFrame,
    task: str = "generic",
    endpoint_column: str | None = None,
) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    """Create explicit endpoint-label aliases and document every mapping."""

    adapted, source_column = normalize_endpoint_column(table, endpoint_column=endpoint_column)
    alias_rows: list[dict[str, str]] = []
    for alias in map_public_label_to_internal("label__endpoint", task=task):
        adapted[alias] = adapted["label__endpoint"]
        alias_rows.append(
            {
                "source_column": source_column,
                "alias_column": alias,
                "task": task,
                "mapping_rule": "copied_binary_0_1_without_semantic_reinterpretation",
                "unique_values": ";".join(sorted(set(adapted[alias].astype(str)))),
            }
        )
    return adapted, alias_rows


def _write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_endpoint_adapter_report(
    path: Path,
    metadata_path: Path,
    source_column: str,
    task: str,
    alias_rows: list[dict[str, str]],
    sample_count: int,
) -> None:
    lines = [
        "# NeuroFate Endpoint Adapter Report",
        "",
        RESEARCH_USE_NOTICE,
        "",
        f"- Input metadata: `{metadata_path}`",
        f"- Source endpoint column: `{source_column}`",
        f"- Task: `{task}`",
        f"- Samples adapted: {sample_count}",
        "",
        "## Aliases Created",
    ]
    for row in alias_rows:
        lines.append(
            f"- `{row['alias_column']}` copied from `{row['source_column']}` "
            f"({row['mapping_rule']}; values={row['unique_values']})"
        )
    lines.extend(
        [
            "",
            "No biological label direction was changed by this adapter.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def adapt_endpoint_metadata(
    metadata: Path,
    outdir: Path,
    task: str = "generic",
    endpoint_column: str | None = None,
) -> EndpointAdapterResult:
    outdir.mkdir(parents=True, exist_ok=True)
    adapted, alias_rows = ensure_label_endpoint_aliases(
        metadata,
        task=task,
        endpoint_column=endpoint_column,
    )
    source_column = alias_rows[0]["source_column"] if alias_rows else endpoint_column or "label__endpoint"
    adapted_path = outdir / "adapted_metadata.tsv"
    aliases_path = outdir / "endpoint_aliases.tsv"
    report_path = outdir / "endpoint_adapter_report.md"
    adapted.to_csv(adapted_path, sep="\t", index=False)
    _write_rows(
        aliases_path,
        alias_rows,
        ["source_column", "alias_column", "task", "mapping_rule", "unique_values"],
    )
    write_endpoint_adapter_report(
        report_path,
        metadata_path=metadata,
        source_column=source_column,
        task=task,
        alias_rows=alias_rows,
        sample_count=len(adapted),
    )
    return EndpointAdapterResult(
        adapted_metadata=adapted_path,
        endpoint_adapter_report=report_path,
        endpoint_aliases=aliases_path,
    )
