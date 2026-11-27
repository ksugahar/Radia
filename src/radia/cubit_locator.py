"""Locate a Coreform Cubit ``bin`` directory for the Radia Cubit panels.

Guards: ``radia`` stays independent of ``radia_mcp``; the Cubit panels need
only Cubit's own Python package, not the MCP server or installer.
Assumes: ``CUBIT_BIN_DIR`` or ``CUBIT_INSTALL_DIR`` when set, else the
standard install locations; the newest matching install wins.
"""
from __future__ import annotations

import glob
import os
from pathlib import Path

_DEFAULT_PATTERNS = (
    "C:/Program Files/Coreform Cubit */bin",
    "/opt/Coreform-Cubit-*/bin",
)


def _has_cubit_python(bin_dir: Path) -> bool:
    return (bin_dir / "cubit.py").is_file() or (bin_dir / "cubit").is_dir()


def find_cubit_bin() -> str | None:
    """Return the Cubit ``bin`` directory as a string, or None when none is found.

    An environment variable that names a directory without Cubit's Python
    package raises instead of falling through to another install.
    """
    for name, suffix in (("CUBIT_BIN_DIR", ""), ("CUBIT_INSTALL_DIR", "bin")):
        value = os.environ.get(name)
        if value:
            path = Path(value) / suffix if suffix else Path(value)
            if not _has_cubit_python(path):
                raise FileNotFoundError(f"{name}={value}: no Cubit Python package under {path}")
            return str(path.resolve())
    for pattern in _DEFAULT_PATTERNS:
        for found in sorted(glob.glob(pattern), reverse=True):
            if _has_cubit_python(Path(found)):
                return str(Path(found).resolve())
    return None
