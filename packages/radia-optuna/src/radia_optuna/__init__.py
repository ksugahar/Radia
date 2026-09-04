"""Locate and verify the separately packaged MATLAB Optuna component."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

__version__ = "0.1.6"
ORACLE_VERSION = "4.9.0"


def _staged_matlab_path() -> Path:
    return Path(str(files(__package__).joinpath("matlab")))


def _checkout_matlab_path() -> Path | None:
    candidate = Path(__file__).resolve().parents[4] / "matlab"
    if (candidate / "+radia" / "+optuna").is_dir():
        return candidate
    return None


def layout() -> str:
    """Report whether the MATLAB tree came from a wheel or checkout."""
    if _staged_matlab_path().is_dir():
        return "wheel"
    if _checkout_matlab_path() is not None:
        return "checkout"
    return "missing"


def matlab_path() -> Path:
    """Return the staged wheel tree or the editable monorepo MATLAB tree."""
    staged = _staged_matlab_path()
    if staged.is_dir():
        return staged
    checkout = _checkout_matlab_path()
    return checkout if checkout is not None else staged


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
    "layout",
    "matlab_addpath_command",
    "matlab_path",
    "mex_path",
]
