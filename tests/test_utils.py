"""
Utility functions for Radia tests

Provides common utilities for:
- Path setup
- Module importing
- Test helpers
"""

import sys
import os
from pathlib import Path

def get_project_root():
	"""
	Get the project root directory.

	Returns:
	    Path: Absolute path to project root directory
	"""
	# Find project root by looking for CMakeLists.txt
	current = Path(__file__).resolve().parent

	# Go up from tests/ directory to find project root
	while current.parent != current:
	    if (current / 'CMakeLists.txt').exists():
	        return current
	    current = current.parent

	# Fallback: assume we're in tests/ and go up one level
	return Path(__file__).resolve().parent.parent

def setup_radia_import():
	"""
	Setup Python path to import radia module.

	This works from anywhere:
	- Project root
	- tests/ directory
	- tests/benchmarks/ directory

	Returns:
	    Path: Path to project root
	"""
	project_root = get_project_root()

	# Add build output directory to Python path
	build_dir = project_root / 'build' / 'lib' / 'Release'
	if build_dir.exists() and str(build_dir) not in sys.path:
	    sys.path.insert(0, str(build_dir))

	# Also try dist directory
	dist_dir = project_root / 'dist'
	if dist_dir.exists() and str(dist_dir) not in sys.path:
	    sys.path.insert(0, str(dist_dir))

	return project_root

def import_radia():
	"""
	Import radia module with automatic path setup.

	Returns:
	    module: The radia module

	Raises:
	    ImportError: If radia module cannot be imported
	"""
	setup_radia_import()

	try:
	    import radia
	    return radia
	except ImportError as e:
	    project_root = get_project_root()
	    build_dir = project_root / 'build' / 'lib' / 'Release'
	    raise ImportError(
	        f"Cannot import radia module. "
	        f"Please build the module first.\n"
	        f"Expected location: {build_dir}\n"
	        f"Original error: {e}"
	    )
