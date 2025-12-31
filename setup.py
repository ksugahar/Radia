#!/usr/bin/env python
"""
Setup script for Radia-NGSolve package

This setup.py is used for building binary distributions (wheels) that include
the pre-compiled C++ extension modules (radia.pyd, rad_ngsolve.pyd).

For development builds, use the CMake build scripts:
- Build.ps1 for radia.pyd
- Build_NGSolve.ps1 for rad_ngsolve.pyd
"""

from setuptools import setup, find_packages
from pathlib import Path
import shutil
import sys

# Read version from pyproject.toml
version = "1.3.16"

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

def prepare_package_data():
	"""
	Prepare package data by copying built extension modules to the package directory

	Note: Package directory is now src/radia (not src/python) so that
	'import radia' works correctly after pip install.

	Build paths searched (in order):
	1. build-msvc/ (MSVC + Intel MKL build - PREFERRED)
	2. build/Release/ (legacy CMake build)
	"""
	package_dir = Path(__file__).parent / "src" / "radia"
	package_dir.mkdir(parents=True, exist_ok=True)

	# Try build-msvc first (MSVC + Intel MKL), then fall back to build/Release
	build_dirs = [
		Path(__file__).parent / "build-msvc",  # MSVC build (preferred)
		Path(__file__).parent / "build" / "Release",  # Legacy CMake build
	]

	# Copy radia.pyd
	radia_pyd_found = False
	for build_dir in build_dirs:
		# Try versioned name first
		radia_pyd = build_dir / "radia.cp312-win_amd64.pyd"
		if radia_pyd.exists():
			shutil.copy2(radia_pyd, package_dir / "radia.pyd")
			print(f"Copied {radia_pyd} to {package_dir}")
			radia_pyd_found = True
			break
		# Try simple name
		radia_pyd = build_dir / "radia.pyd"
		if radia_pyd.exists():
			shutil.copy2(radia_pyd, package_dir / "radia.pyd")
			print(f"Copied {radia_pyd} to {package_dir}")
			radia_pyd_found = True
			break

	if not radia_pyd_found:
		print(f"Warning: radia.pyd not found. Run BuildMSVC.ps1 first.")

	# Copy radia_ngsolve.pyd
	for build_dir in build_dirs:
		radia_ngsolve_pyd = build_dir / "radia_ngsolve.pyd"
		if radia_ngsolve_pyd.exists():
			shutil.copy2(radia_ngsolve_pyd, package_dir / "radia_ngsolve.pyd")
			print(f"Copied {radia_ngsolve_pyd} to {package_dir}")
			break
	else:
		print(f"Info: radia_ngsolve.pyd not found. This is optional.")

	return package_dir

# Prepare package data before setup
if "sdist" not in sys.argv:
	# Only copy files for binary distributions, not source distributions
	package_dir = prepare_package_data()

setup(
	name="radia",
	version=version,
	description="Radia 3D Magnetostatics with NGSolve Integration and OpenMP Parallelization",
	long_description=long_description,
	long_description_content_type="text/markdown",
	author="Oleg Chubar, Pascal Elleaume",
	author_email="chubar@bnl.gov",
	maintainer="Radia Development Team",
	url="https://github.com/ksugahar/Radia_NGSolve",
	project_urls={
		"Homepage": "https://github.com/ksugahar/Radia_NGSolve",
		"Documentation": "https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia",
		"Repository": "https://github.com/ksugahar/Radia_NGSolve",
		"Issues": "https://github.com/ksugahar/Radia_NGSolve/issues",
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
			"*.pyd",  # Include all .pyd files (radia.pyd, radia_ngsolve.pyd)
			"*.dll",  # Include Intel MKL DLLs (mkl_rt.2.dll, mkl_core.2.dll, libiomp5md.dll, etc.)
			"*.py",   # Include all Python utility modules
		],
	},
	include_package_data=True,
	install_requires=[
		"numpy>=1.20",
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
