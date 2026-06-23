"""Helpers for Coreform Cubit 2025.12 export tests.

These tests exercise the current public Cubit path:

    cubit.cmd('export netgen "model.vol" order N overwrite')

The old top-level ``extract_curved_mesh`` Python helper is retired from the
public API.  ``cubit_mesh_curver`` remains a low-level pybind module used by
the .ccm plugin.
"""

from __future__ import annotations

import os
import re
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
RADIA_SRC = REPO_ROOT / "src" / "radia"
CME_PKG_DIR = REPO_ROOT / "packages" / "cubit-mesh-export" / "src" / "cubit_mesh_export"
_CUBIT_MODULE = None
_DLL_HANDLES = []


def find_cubit_bin() -> Path:
    """Return the supported Coreform Cubit 2025.12+ bin directory."""
    if str(RADIA_SRC) not in sys.path:
        sys.path.insert(0, str(RADIA_SRC))
    from install_panels import find_cubit_bin as _find_cubit_bin

    cubit_bin = _find_cubit_bin()
    if not cubit_bin:
        raise RuntimeError("Coreform Cubit 2025.12+ not found; set CUBIT_PATH")
    return Path(cubit_bin)


def add_cubit_to_path() -> Path:
    """Make Cubit's Python module and DLL directory visible."""
    cubit_bin = find_cubit_bin()
    if str(cubit_bin) not in sys.path:
        sys.path.insert(0, str(cubit_bin))
    os.environ["PATH"] = str(cubit_bin) + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        for dll_dir in (cubit_bin, cubit_bin / "python3"):
            if dll_dir.is_dir():
                _DLL_HANDLES.append(os.add_dll_directory(str(dll_dir)))
    return cubit_bin


def _required_python_minor(cubit_bin: Path) -> tuple[int, int] | None:
    """Infer Cubit's embedded Python minor version from its DLL name."""
    python_dir = cubit_bin / "python3"
    for dll in python_dir.glob("python3*.dll"):
        match = re.fullmatch(r"python(\d)(\d{2})\.dll", dll.name)
        if match:
            return int(match.group(1)), int(match.group(2))
    return None


def init_cubit():
    """Import and initialize Cubit in batch mode."""
    global _CUBIT_MODULE
    if _CUBIT_MODULE is not None:
        return _CUBIT_MODULE

    cubit_bin = add_cubit_to_path()
    required = _required_python_minor(cubit_bin)
    if required and sys.version_info[:2] != required:
        import pytest

        pytest.skip(
            "Coreform Cubit 2025.12 embeds Python "
            f"{required[0]}.{required[1]}; real Cubit import tests must run "
            "with Cubit's bundled python3/python.exe"
        )

    import cubit

    plugin_dir = cubit_bin / "plugins"
    cubit.init([
        "cubit",
        "-nojournal",
        "-batch",
        "-nographics",
        "-commandplugindir",
        str(plugin_dir),
    ])
    _CUBIT_MODULE = cubit
    return cubit


def add_cubit_mesh_curver_to_path() -> None:
    """Expose the low-level pybind module for module-shape tests only."""
    if str(CME_PKG_DIR) not in sys.path:
        sys.path.insert(0, str(CME_PKG_DIR))


def export_netgen(cubit, stem: str, order: int = 1) -> Path:
    """Run the public 2025.12 Netgen export command and return the .vol path."""
    out_dir = Path(tempfile.mkdtemp(prefix="radia_cubit_test_"))
    out = out_dir / f"{stem}_o{order}.vol"
    cmd_path = str(out).replace("\\", "/")
    cubit.cmd(f'export netgen "{cmd_path}" order {order} overwrite')
    if not out.is_file() or out.stat().st_size == 0:
        raise AssertionError(f"export netgen did not create {out}")
    return out


def load_ngsolve_mesh(path: Path):
    """Load a .vol file with NGSolve."""
    from ngsolve import Mesh

    return Mesh(str(path))
