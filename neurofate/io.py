"""Input/output helpers with safety-first defaults."""

from __future__ import annotations

from pathlib import Path


HEAVY_SUFFIXES = {".h5", ".hdf5", ".h5ad", ".zarr", ".loom", ".mtx", ".parquet"}


def is_heavy_data_path(path: str | Path) -> bool:
    """Return True when a path looks like a large biological data artifact."""
    return Path(path).suffix.lower() in HEAVY_SUFFIXES


def assert_lightweight_path(path: str | Path) -> Path:
    """Refuse paths that should only be handled by explicit heavy workflows."""
    resolved = Path(path)
    if is_heavy_data_path(resolved):
        raise ValueError(f"Refusing heavy data path in lightweight context: {resolved}")
    return resolved
