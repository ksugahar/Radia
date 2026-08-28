"""Locate and verify the separately packaged MATLAB Optuna component."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path


__version__ = "0.1.0"
ORACLE_VERSION = "4.9.0"


def _staged_matlab_path() -> Path:
    return Path(str(files(__package__).joinpath("matlab")))


def _checkout_matlab_path() -> Path | None:
    """Locate the monorepo tree a wheel build would have staged.

    ``setup.py`` copies ``<repo>/matlab`` into the package only while
    building a wheel, so an editable install has no staged tree.
    """
    candidate = Path(__file__).resolve().parents[4] / "matlab"
    if (candidate / "+radia" / "+optuna").is_dir():
        return candidate
    return None


def layout() -> str:
    """Report which MATLAB tree ``matlab_path`` resolved: never guess."""
    if _staged_matlab_path().is_dir():
        return "wheel"
    if _checkout_matlab_path() is not None:
        return "checkout"
    return "missing"


def matlab_path() -> Path:
    """Return the directory that MATLAB must add to its path.

    A published wheel stages the canonical monorepo sources as
    ``radia_optuna/matlab``. An editable install has no staged tree, so the
    checked-out ``<repo>/matlab`` directory is used instead; ``layout()``
    always says which one was resolved.
    """
    staged = _staged_matlab_path()
    if staged.is_dir():
        return staged
    checkout = _checkout_matlab_path()
    if checkout is not None:
        return checkout
    return staged


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
