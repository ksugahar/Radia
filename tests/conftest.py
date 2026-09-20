"""
pytest configuration and shared fixtures for Radia tests

Usage:
  pytest tests/                     # Run all tests
  pytest tests/ -m basic            # Run only basic tests
  pytest tests/                     # Short CI/debug suite
  pytest validation_test/           # Heavy validation / golden / GUI / Cubit checks
"""

import os
import sys
from pathlib import Path

import pytest


def setup_radia_path():
    """Add the repository's package roots without exposing package internals."""
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


PROJECT_ROOT = setup_radia_path()


# ---------------------------------------------------------------
# collect_ignore: Skip test files that import unavailable packages.
# This runs BEFORE pytest tries to import test modules, preventing
# DLL load failures and access violations from crashing the process.
# ---------------------------------------------------------------
def _check_module(name):
    """Import a module, returning None on success or the reason it failed."""
    try:
        __import__(name)
        return None
    except BaseException as exc:  # a DLL load failure is not an ImportError
        return "%s: %s" % (type(exc).__name__, exc)

# Auto-detect Cubit through its owning distribution.  This probe is for an
# optional extra, so nothing it raises may end collection -- a broken netgen
# reaches this import first and used to kill the session with a raw traceback,
# ahead of the check below that can actually say what to fix.
try:
    from cubit_mesh_export.toolbar_install import find_cubit_bin
    _cubit_path = find_cubit_bin()
except BaseException:
    _cubit_path = None
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

# `ngsolve` and `netgen` are pinned, non-optional dependencies in
# pyproject.toml, so a failed import there means a broken environment, not an
# absent extra.  Treating them as optional silently shrank a 2026-09-20 run by
# 529 of 5313 tests -- ~10 % of the suite vanished and the summary was still
# green, which is indistinguishable from a full pass.  A partial run now has
# to be asked for by name.
_REQUIRED_MODULES = ("ngsolve", "netgen")
# radia_ngsolve is no longer a separate module; RadiaField is in radia
_OPTIONAL_MODULES = ("magpylib", "cubit")
ALLOW_PARTIAL_ENV = "RADIA_TESTS_ALLOW_PARTIAL"

_IMPORT_FAILURES = {
    name: _check_module(name)
    for name in _REQUIRED_MODULES + _OPTIONAL_MODULES
}
_OPTIONAL_DEPS = {
    name: reason is None for name, reason in _IMPORT_FAILURES.items()
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


_MISSING_REQUIRED = {
    name: reason
    for name in _REQUIRED_MODULES
    if (reason := _IMPORT_FAILURES[name]) is not None
}
_ALLOW_PARTIAL = os.environ.get(ALLOW_PARTIAL_ENV, "") not in ("", "0")

if _MISSING_REQUIRED and not _ALLOW_PARTIAL:
    raise pytest.UsageError(
        "pinned dependencies failed to import, so %d of the %d test files "
        "in %s would be skipped without the run looking any different from "
        "a full pass:%s%s%sFix the environment (both are pinned in "
        "pyproject.toml), or ask for the partial run by name with %s=1."
        % (len(collect_ignore), len(list(_tests_dir.glob("test_*.py"))),
           _tests_dir, os.linesep,
           os.linesep.join("  %s -> %s" % (name, reason)
                           for name, reason in sorted(_MISSING_REQUIRED.items())),
           os.linesep, ALLOW_PARTIAL_ENV))


def pytest_report_header(config):
    """Say out loud how much of the suite this environment can collect."""
    del config
    lines = []
    if _MISSING_REQUIRED:
        lines.append(
            "PARTIAL RUN (%s=1): pinned %s unavailable"
            % (ALLOW_PARTIAL_ENV, ", ".join(sorted(_MISSING_REQUIRED))))
    absent = sorted(name for name in _OPTIONAL_MODULES
                    if _IMPORT_FAILURES[name] is not None)
    if absent:
        lines.append("optional modules unavailable: %s" % ", ".join(absent))
    if collect_ignore:
        lines.append("test files not collected: %d" % len(collect_ignore))
    return lines


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
        if "golden" in nodeid.lower():
            item.add_marker(pytest.mark.golden)
        if "benchmarks" in str(item.fspath):
            item.add_marker(pytest.mark.benchmark)
            item.add_marker(pytest.mark.slow)
        if "ngsolve" in item.name.lower() or "ngsolve" in str(item.fspath).lower():
            item.add_marker(pytest.mark.ngsolve)
