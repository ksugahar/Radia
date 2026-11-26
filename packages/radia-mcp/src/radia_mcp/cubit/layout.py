"""Installed Cubit layout for the package-based 2026.8 runtime.

``daemon.py`` loads this file by path inside Cubit's bundled Python, so it
must stay standard-library only and import nothing from ``radia_mcp``.
"""
from pathlib import Path
import os


SUPPORTED_CUBIT_RELEASE = "2026.8"


def build_info_path(root: Path) -> Path:
    """The release manifest of a Cubit install root."""
    return Path(root) / "BUILD_INFO"


def cubit_release(root: Path) -> str | None:
    """The `short` field of BUILD_INFO under a Cubit install root, or None.

    This is the only reader of BUILD_INFO: the session and the daemon both go
    through it, so they cannot disagree on what a release is.  The directory
    name is never evidence.
    """
    try:
        text = build_info_path(root).read_text(encoding="utf-8")
    except FileNotFoundError:
        return None
    for line in text.splitlines():
        key, sep, value = line.partition("=")
        if sep and not line.lstrip().startswith("#") and key.strip() == "short":
            return value.strip() or None
    return None


def is_supported_cubit_root(root: Path) -> bool:
    """True only when BUILD_INFO records the supported release."""
    return cubit_release(root) == SUPPORTED_CUBIT_RELEASE


def has_python_binding(bin_dir: Path) -> bool:
    return any((bin_dir / name).is_file() for name in (
        'cubit.py', '_cubit3.pyd', '_cubit3.so')) or (
            (bin_dir / 'cubit/__init__.py').is_file()
            and (bin_dir / 'cubit/_cubit3.pyd').is_file())


def runtime_environment(bin_dir: Path) -> dict[str, str]:
    """Process environment for Cubit's bundled Python (ACIS DLLs on PATH)."""
    env = os.environ.copy()
    acis = Path(bin_dir).parent / 'acis/code/bin'
    if acis.is_dir():
        env['PATH'] = str(acis) + os.pathsep + env.get('PATH', '')
    return env
