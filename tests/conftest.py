"""
pytest configuration and shared fixtures for Radia tests

Usage:
  pytest tests/                     # Run all tests
  pytest tests/ -m basic            # Run only basic tests
  pytest tests/ -m "not slow"       # Skip slow tests
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

    # Add MKL DLL directory for peec_matrices.pyd etc.
    mkl_bin = os.path.join(sys.prefix, 'Library', 'bin')
    if os.path.isdir(mkl_bin) and hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(mkl_bin)

    return project_root


PROJECT_ROOT = setup_radia_path()


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
            if not _avail and (f"from {_mod} " in _stripped or
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
    config.addinivalue_line("markers", "benchmark: Performance benchmarks")
    config.addinivalue_line("markers", "ngsolve: Tests requiring NGSolve")


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
        if "benchmarks" in str(item.fspath):
            item.add_marker(pytest.mark.benchmark)
            item.add_marker(pytest.mark.slow)
        if "ngsolve" in item.name.lower() or "ngsolve" in str(item.fspath).lower():
            item.add_marker(pytest.mark.ngsolve)
