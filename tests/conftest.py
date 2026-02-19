"""
pytest configuration and shared fixtures for Radia tests

This file provides:
- Automatic path setup for importing radia module
- Shared fixtures for common test resources
- pytest hooks for test discovery

Usage:
  pytest tests/                     # Run all tests
  pytest tests/ -m basic            # Run only basic tests
  pytest tests/ -m "not slow"       # Skip slow tests
  pytest tests/ -k "test_import"    # Run tests matching pattern
"""

import sys
import os
from pathlib import Path
import pytest

def setup_radia_path():
    """
    Setup Python path to import radia module from build directory.

    This function works regardless of where the test is run from:
    - Project root: python tests/test_simple.py
    - Tests directory: python test_simple.py
    - Benchmarks directory: python benchmark_openmp.py

    Priority order for finding _radia.pyd:
    1. src/radia/ (package directory - Build.ps1 copies here)
    2. build-msvc/ (MSVC build output)
    3. build-intel/ (Intel build output)
    4. build/lib/Release/ (legacy build location)
    5. build/Release/ (alternative build location)
    """
    # Find project root by looking for CMakeLists.txt
    current = Path(__file__).resolve().parent

    # Go up from tests/ directory to find project root
    while current.parent != current:
        if (current / 'CMakeLists.txt').exists():
            project_root = current
            break
        current = current.parent
    else:
        # Fallback: assume we're in tests/ and go up one level
        project_root = Path(__file__).resolve().parent.parent

    # Add src/ directory so that 'import radia' finds src/radia/__init__.py
    src_path = project_root / 'src'
    if src_path.exists():
        sys.path.insert(0, str(src_path))

    return project_root

# Setup path when this module is imported
PROJECT_ROOT = setup_radia_path()

# pytest configuration
def pytest_configure(config):
    """Configure pytest with custom markers"""
    # Markers are also defined in pyproject.toml, but we add them here too for completeness
    config.addinivalue_line("markers", "basic: Basic functionality tests (fast)")
    config.addinivalue_line("markers", "comprehensive: Comprehensive test suite")
    config.addinivalue_line("markers", "advanced: Advanced features and edge cases")
    config.addinivalue_line("markers", "performance: Performance and scaling tests")
    config.addinivalue_line("markers", "slow: Tests that take more than 10 seconds")
    config.addinivalue_line("markers", "benchmark: Performance benchmarks (not run by default)")
    config.addinivalue_line("markers", "ngsolve: Tests requiring NGSolve integration")


# Shared fixtures
@pytest.fixture(scope="session")
def radia_module():
    """
    Fixture that provides the radia module.
    Ensures proper import and cleanup.
    """
    import radia as rad
    rad.UtiDelAll()  # Clean state
    yield rad
    rad.UtiDelAll()  # Cleanup after all tests


@pytest.fixture
def radia_clean():
    """
    Fixture that provides a clean radia state for each test.
    """
    import radia as rad
    rad.UtiDelAll()  # Clean before test
    yield rad
    rad.UtiDelAll()  # Cleanup after test


@pytest.fixture(scope="session")
def project_root():
    """Fixture providing the project root path."""
    return PROJECT_ROOT


def pytest_collection_modifyitems(config, items):
    """
    Modify test collection to add markers based on test location.
    """
    for item in items:
        # Add benchmark marker to tests in benchmarks/ directory
        if "benchmarks" in str(item.fspath):
            item.add_marker(pytest.mark.benchmark)
            item.add_marker(pytest.mark.slow)

        # Add ngsolve marker to tests with ngsolve in name
        if "ngsolve" in item.name.lower() or "ngsolve" in str(item.fspath).lower():
            item.add_marker(pytest.mark.ngsolve)
