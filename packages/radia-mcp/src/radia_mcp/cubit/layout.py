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


def install_root(path: Path) -> Path:
    """The Cubit install root for either the root itself or its ``bin``."""
    root = Path(path)
    return root.parent if root.name.lower() == 'bin' else root


def plugin_directory(install_dir: Path) -> Path:
    """The directory Cubit scans for command plugins."""
    root = install_root(install_dir)
    if (root / 'bin/cubitx.exe').is_file():
        return root / 'plugins'
    return root / 'bin/plugins'


def console_executable(bin_dir: Path) -> Path | None:
    """Cubit's headless console executable under ``bin_dir``, or None."""
    for name in ('cubitx.exe', 'coreform_cubit.com', 'coreform_cubit', 'cubit'):
        candidate = Path(bin_dir) / name
        if candidate.is_file():
            return candidate
    return None


COMMAND_PLUGIN = 'cubit_mesh_export.dll'


def command_plugin_path(install_dir: Path) -> Path:
    """Where the installer deploys the ``export netgen`` command plugin."""
    return plugin_directory(install_dir) / COMMAND_PLUGIN


def runtime_environment(bin_dir: Path) -> dict[str, str]:
    """Process environment for Cubit's bundled Python (ACIS DLLs on PATH)."""
    env = os.environ.copy()
    acis = Path(bin_dir).parent / 'acis/code/bin'
    if acis.is_dir():
        env['PATH'] = str(acis) + os.pathsep + env.get('PATH', '')
    return env
