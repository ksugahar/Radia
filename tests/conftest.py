"""
pytest configuration and shared fixtures for Radia tests

This file provides:
- Automatic path setup for importing radia module
- Shared fixtures for common test resources
- pytest hooks for test discovery
"""

import sys
import os
from pathlib import Path

def setup_radia_path():
    """
    Setup Python path to import radia module from build directory.

    This function works regardless of where the test is run from:
    - Project root: python tests/test_simple.py
    - Tests directory: python test_simple.py
    - Benchmarks directory: python benchmark_openmp.py
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

    # Add build output directory to Python path
    build_dir = project_root / 'build' / 'lib' / 'Release'
    if build_dir.exists():
        sys.path.insert(0, str(build_dir))

    # Also try dist directory
    dist_dir = project_root / 'dist'
    if dist_dir.exists():
        sys.path.insert(0, str(dist_dir))

    return project_root

# Setup path when this module is imported
PROJECT_ROOT = setup_radia_path()

# pytest configuration
def pytest_configure(config):
    """Configure pytest"""
    config.addinivalue_line("markers", "basic: basic functionality tests")
    config.addinivalue_line("markers", "comprehensive: comprehensive test suite")
    config.addinivalue_line("markers", "slow: tests that take more than 10 seconds")
    config.addinivalue_line("markers", "benchmark: performance benchmarks")
