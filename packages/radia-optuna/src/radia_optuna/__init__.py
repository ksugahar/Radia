"""Locate and verify the separately packaged MATLAB Optuna component."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path


__version__ = "0.1.0"
ORACLE_VERSION = "4.9.0"


def matlab_path() -> Path:
    """Return the installed directory that MATLAB must add to its path."""
    return Path(str(files(__package__).joinpath("matlab")))


def mex_path() -> Path:
    """Return the installed Windows x64 Optuna MEX path."""
    return matlab_path() / "optuna_mex.mexw64"


def matlab_addpath_command() -> str:
    """Return a copy-pasteable MATLAB addpath command."""
    escaped = str(matlab_path()).replace("'", "''")
    return f"addpath('{escaped}')"


__all__ = [
    "ORACLE_VERSION",
    "__version__",
    "matlab_addpath_command",
    "matlab_path",
    "mex_path",
]
