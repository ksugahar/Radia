"""
pytest configuration and shared fixtures for Radia tests

Usage:
  pytest tests/                     # Run all tests
  pytest tests/ -m basic            # Run only basic tests
  pytest tests/ -m "not slow"       # Skip slow tests
  pytest tests/ -m "not slow and not golden"  # Fast CI/simple gate
  pytest validation_test/           # Heavy validation / golden / GUI / Cubit checks
"""

import sys
import os
from pathlib import Path
import pytest


def setup_radia_path():
    """Setup Python path to import radia from src/radia/."""
    current = Path(__file__).resolve().parent
    while current.parent != current:
        if (current / 'CMakeLists.txt').exists():
            project_root = current
            break
        current = current.parent
    else:
        project_root = Path(__file__).resolve().parent.parent

    src_path = project_root / 'src'
    if src_path.exists():
        sys.path.insert(0, str(src_path))

    # Note: `sys.path.insert(0, src/radia)` used to be added here so bare
    # imports like `from coil_from_cad import ...` worked.  Removed
    # because bare imports + canonical `radia.coil_from_cad` simul-
    # taneously load peec_matrices.pyd under two module keys, triggering
    # pybind11's "type X is already registered" error.  All tests now
    # use `from radia.X import Y` exclusively.

    mcp_src = project_root / 'packages' / 'radia-mcp' / 'src'
    if mcp_src.exists():
        sys.path.insert(0, str(mcp_src))

    # Add MKL DLL directory for peec_matrices.pyd etc.
    mkl_bin = os.path.join(sys.prefix, 'Library', 'bin')
    if os.path.isdir(mkl_bin) and hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(mkl_bin)

    return project_root


# ---------------------------------------------------------------
# Qt6 DLL preload (POLICY 2026-05-30): must happen BEFORE
# setup_radia_path() adds MKL bin to the Windows DLL search path.
#
# Why: setup_radia_path() calls os.add_dll_directory(<sys.prefix>/Library/bin)
# so radia's C++ extension (peec_matrices.pyd etc.) finds MKL.  But
# MKL's bin ships an older vcruntime / msvcp that the Windows DLL
# loader picks BEFORE PySide6's site-packages copy once MKL bin is on
# the search path -- and Qt6Core's required symbols differ across
# vcruntime versions, so Qt6Core.dll fails to load with 0xc0000139
# ("specified procedure not found") the first time any test module
# imports PySide6.QtCore.  Historically this made all panel pytest
# tests crash on LAB ("LAB caveat" in panel-qt-test SKILL.md).
#
# Fix: load Qt6Core BEFORE we touch the DLL search path.  Once
# Qt6Core is mapped into the process, the MKL vcruntime cannot shadow
# it (Windows does not reload an already-loaded DLL).  Verified
# 2026-05-30: panel pytest goes from "DLL load failed while importing
# QtCore" to 11 / 11 passing in 2.97 s on LAB.
#
# pytest-qt's qapp fixture also preloads QtCore, but pytest collection
# (where test modules' top-level imports execute) happens BEFORE
# fixtures fire -- so the preload MUST be at conftest top, not in a
# fixture.  No PySide6 -> try/except handles it; panel tests are
# auto-skipped via collect_ignore further down.
try:
    from PySide6 import QtCore as _qt6_preload  # noqa: F401
except ImportError:
    pass


PROJECT_ROOT = setup_radia_path()


def _load_ci_slow_nodeids():
    """Load measured-slow tests so fast CI can keep the runner responsive."""
    path = PROJECT_ROOT / "tests" / "ci_slow_nodeids.txt"
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return set()

    nodeids = set()
    for line in lines:
        line = line.split("#", 1)[0].strip()
        if line:
            nodeids.add(line.replace("\\", "/"))
    return nodeids


CI_SLOW_NODEIDS = _load_ci_slow_nodeids()


# ---------------------------------------------------------------
# collect_ignore: Skip test files that import unavailable packages.
# This runs BEFORE pytest tries to import test modules, preventing
# DLL load failures and access violations from crashing the process.
# ---------------------------------------------------------------
def _check_module(name):
    """Check if a module can be imported safely (with src/radia in path)."""
    try:
        __import__(name)
        return True
    except (ImportError, OSError, Exception):
        return False

# Auto-detect Cubit via install_panels.find_cubit_bin() (single source of truth)
try:
    from radia.install_panels import find_cubit_bin
    _cubit_path = find_cubit_bin()
except ImportError:
    _cubit_path = None
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

_OPTIONAL_DEPS = {
    "ngsolve": _check_module("ngsolve"),
    "netgen": _check_module("netgen"),
    # radia_ngsolve is no longer a separate module; RadiaField is in radia
    "magpylib": _check_module("magpylib"),
    "cubit": _check_module("cubit"),
}

# Build the exclusion list by scanning test files for top-level imports
_tests_dir = Path(__file__).parent
collect_ignore = []

for _tf in sorted(_tests_dir.glob("test_*.py")):
    try:
        _content = _tf.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        continue

    _skip = False
    for _line in _content.split("\n"):
        _stripped = _line.lstrip()
        # Stop scanning at first function/class definition
        if _stripped.startswith("def ") or _stripped.startswith("class "):
            break
        # Check for imports of unavailable modules at any top-level indent
        # (including inside try/except, because DLL load can segfault)
        for _mod, _avail in _OPTIONAL_DEPS.items():
            # `from {_mod}.` catches submodule imports (e.g. `from netgen.occ
            # import ...`) that the trailing-space form `from {_mod} ` misses.
            if not _avail and (f"from {_mod} " in _stripped or
                               f"from {_mod}." in _stripped or
                               f"import {_mod}" in _stripped):
                _skip = True
                break
        if _skip:
            break
        # Check for sys.exit() at module level (crashes pytest collection)
        if "sys.exit(" in _stripped:
            _skip = True
            break

    if _skip:
        collect_ignore.append(str(_tf))


# ---------------------------------------------------------------
# Markers
# ---------------------------------------------------------------
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line("markers", "basic: Basic functionality tests (fast)")
    config.addinivalue_line("markers", "comprehensive: Comprehensive test suite")
    config.addinivalue_line("markers", "advanced: Advanced features and edge cases")
    config.addinivalue_line("markers", "performance: Performance and scaling tests")
    config.addinivalue_line("markers", "slow: Tests that take more than 10 seconds")
    config.addinivalue_line("markers", "golden: Golden/reference tests separated from the simple CI gate")
    config.addinivalue_line("markers", "benchmark: Performance benchmarks")
    config.addinivalue_line("markers", "ngsolve: Tests requiring NGSolve")
    config.addinivalue_line("markers", "validation: Heavy validation_test suite")


# ---------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------
@pytest.fixture(scope="session")
def radia_module():
    """Provides the radia module with clean state."""
    import radia as rad
    rad.UtiDelAll()
    yield rad
    rad.UtiDelAll()


@pytest.fixture
def radia_clean():
    """Provides a clean radia state for each test."""
    import radia as rad
    rad.UtiDelAll()
    yield rad
    rad.UtiDelAll()


@pytest.fixture(scope="session")
def project_root():
    """Provides the project root path."""
    return PROJECT_ROOT


# ---------------------------------------------------------------
# Auto-markers
# ---------------------------------------------------------------
def pytest_collection_modifyitems(config, items):
    """Add markers based on test location."""
    for item in items:
        nodeid = item.nodeid.replace("\\", "/")
        if nodeid in CI_SLOW_NODEIDS:
            item.add_marker(pytest.mark.slow)
        if "golden" in nodeid.lower():
            item.add_marker(pytest.mark.golden)
        if "benchmarks" in str(item.fspath):
            item.add_marker(pytest.mark.benchmark)
            item.add_marker(pytest.mark.slow)
        if "ngsolve" in item.name.lower() or "ngsolve" in str(item.fspath).lower():
            item.add_marker(pytest.mark.ngsolve)
