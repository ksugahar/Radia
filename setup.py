#!/usr/bin/env python
"""
Setup script for Radia-NGSolve package

This setup.py is used for building binary distributions (wheels) that include
the pre-compiled C++ extension module (_radia_pybind.pyd).
radia_ngsolve is now pure Python (no C++ build needed).

For development builds, use the CMake build scripts:
- Build.ps1 for _radia_pybind.pyd (main module)
- radia_ngsolve.py is pure Python (no build needed)

Module naming: The C++ extension is named '_radia_pybind' (with underscore) so that
'import radia' uses the Python package with __init__.py, which then imports
from _radia_pybind. This follows NGSolve's pattern.

pybind11 Migration Complete (2026-01):
All bindings now use pybind11 exclusively via _radia_pybind.pyd.
"""

from setuptools import setup, find_packages
from pathlib import Path
import os
import shutil
import sys

# Version is read from pyproject.toml by setuptools
# This fallback is only used if pyproject.toml is not available
version = "4.3.0"

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

def prepare_package_data():
	"""
	Prepare package data by copying built extension modules to the package directory

	Note: Package directory is src/radia so that
	'import radia' works correctly after pip install.

	Build path: build-msvc/ (MSVC + Intel MKL build)
	"""
	package_dir = Path(__file__).parent / "src" / "radia"
	package_dir.mkdir(parents=True, exist_ok=True)

	build_dir = Path(__file__).parent / "build-msvc"

	# Copy _radia_pybind.pyd (main C++ extension module)
	radia_pyd_found = False
	# Try versioned name first
	radia_pyd = build_dir / "_radia_pybind.cp312-win_amd64.pyd"
	if radia_pyd.exists():
		shutil.copy2(radia_pyd, package_dir / "_radia_pybind.pyd")
		print(f"Copied {radia_pyd} to {package_dir}")
		radia_pyd_found = True
	else:
		# Try simple name
		radia_pyd = build_dir / "_radia_pybind.pyd"
		if radia_pyd.exists():
			shutil.copy2(radia_pyd, package_dir / "_radia_pybind.pyd")
			print(f"Copied {radia_pyd} to {package_dir}")
			radia_pyd_found = True

	if not radia_pyd_found:
		print(f"Warning: _radia_pybind.pyd not found. Run Build.ps1 first.")

	# radia_ngsolve is now pure Python (radia_ngsolve.py), no .pyd needed

	# Copy Cubit plugin files (Build.ps1 normally handles this)
	cubit_files = [
		("src/cubit_plugin/build-pyd/radia_cubit_mesh.cp312-win_amd64.pyd",
		 "radia_cubit_mesh.pyd"),
		("src/cubit_plugin/build-ccm/radia_cubit.ccm",
		 "radia_cubit.ccm"),
	]
	root = Path(__file__).parent
	for src_rel, dst_name in cubit_files:
		src = root / src_rel
		dst = package_dir / dst_name
		if src.exists() and not dst.exists():
			shutil.copy2(src, dst)
			print(f"Copied {src.name} to {package_dir}")
		elif not dst.exists():
			print(f"Warning: {dst_name} not found (Cubit features disabled)")

	return package_dir

# Prepare package data before setup
# Skip for editable installs (pyd already in src/radia/) and sdist
# pip runs setup.py via build_meta without "editable" in sys.argv,
# so we check if the pyd is already present and not locked.
_skip_copy = "sdist" in sys.argv or os.environ.get("RADIA_SKIP_PYD_COPY")
if not _skip_copy:
	try:
		package_dir = prepare_package_data()
	except PermissionError:
		print("Warning: .pyd files locked (editable install?). Skipping copy.")
		package_dir = Path(__file__).parent / "src" / "radia"

setup(
	name="radia",
	version=version,
	description="Radia 3D Magnetostatics with NGSolve Integration and TaskManager Parallelization",
	long_description=long_description,
	long_description_content_type="text/markdown",
	author="Oleg Chubar, Pascal Elleaume",
	author_email="chubar@bnl.gov",
	maintainer="Radia Development Team",
	url="https://github.com/ksugahar/Radia",
	project_urls={
		"Homepage": "https://github.com/ksugahar/Radia",
		"Documentation": "https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia",
		"Repository": "https://github.com/ksugahar/Radia",
		"Issues": "https://github.com/ksugahar/Radia/issues",
	},
	license="BSD-style AND MIT",
	classifiers=[
		"Development Status :: 5 - Production/Stable",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: BSD License",
		"License :: OSI Approved :: MIT License",
		"Programming Language :: Python :: 3.12",
		"Programming Language :: Python :: 3 :: Only",
		"Programming Language :: C++",
		"Topic :: Scientific/Engineering :: Physics",
		"Operating System :: Microsoft :: Windows",
	],
	keywords=["magnetostatics", "magnetic field", "radia", "synchrotron", "ngsolve", "fem"],
	python_requires=">=3.12",
	packages=find_packages(where="src"),
	package_dir={"": "src"},
	package_data={
		"radia": [
			"*.pyd",  # C++ extensions (_radia_pybind.pyd, radia_cubit_mesh.pyd, etc.)
			"*.ccm",  # Cubit plugin (radia_cubit.ccm)
			"*.py",   # Python utility modules
			# MKL DLLs are NOT bundled - installed via pip dependency (mkl>=2024.2.0)
			"panels/**/*.py",   # Cubit toolbar panel scripts
			"panels/**/*.ui",   # Cubit toolbar panel UI files
			"panels/**/*.txt",  # Cubit toolbar panel marker files
			"resources/*.png",  # App icons
		],
	},
	entry_points={
		"console_scripts": [
			"cubit-install-panels=radia.install_panels:main",
			"radia-setup=radia.setup_cubit:main",
			"radia-app=radia.radia_app:main",
			# MCP server entry-points moved to packages/radia-mcp (public)
			# and s:\mcp-server\mcp-server-{ih,electromagnet,peec} (LAB-private)
			# as of 2026-04-21. radia wheel no longer ships MCP servers.
		],
	},
	include_package_data=True,
	install_requires=[
		"numpy>=1.20",
		"mkl>=2024.2.0",
		"intel-cmplr-lib-rt>=2024.2.0",
		"ngsolve==6.2.2603",
		"netgen-mesher==6.2.2603",
	],
	extras_require={
		"viz": [
			"pyvista>=0.40",
			"matplotlib>=3.5",
		],
		"test": [
			"pytest>=7.0",
			"pytest-cov>=4.0",
		],
		"dev": [
			"pytest>=7.0",
			"pytest-cov>=4.0",
			"pyvista>=0.40",
			"matplotlib>=3.5",
			"build>=0.10",
			"twine>=4.0",
		],
	},
	zip_safe=False,  # Don't zip the package (needed for .pyd files)
)
