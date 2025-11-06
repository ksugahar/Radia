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
version = "1.0.5"

# Read the README file
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

def prepare_package_data():
	"""
	Prepare package data by copying built extension modules to the package directory
	"""
	package_dir = Path(__file__).parent / "src" / "python"
	package_dir.mkdir(parents=True, exist_ok=True)

	# Copy radia.pyd from dist/ if it exists
	radia_pyd = Path(__file__).parent / "dist" / "radia.pyd"
	if radia_pyd.exists():
		shutil.copy2(radia_pyd, package_dir / "radia.pyd")
		print(f"Copied {radia_pyd} to {package_dir}")
	else:
		print(f"Warning: {radia_pyd} not found. Run Build.ps1 first.")

	# Copy rad_ngsolve.pyd from build/Release/ if it exists
	rad_ngsolve_pyd = Path(__file__).parent / "build" / "Release" / "rad_ngsolve.pyd"
	if rad_ngsolve_pyd.exists():
		shutil.copy2(rad_ngsolve_pyd, package_dir / "rad_ngsolve.pyd")
		print(f"Copied {rad_ngsolve_pyd} to {package_dir}")
	else:
		print(f"Info: {rad_ngsolve_pyd} not found. This is optional.")

	return package_dir

# Prepare package data before setup
if "sdist" not in sys.argv:
	# Only copy files for binary distributions, not source distributions
	package_dir = prepare_package_data()

setup(
	name="radia-ngsolve",
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
	license="LGPL-2.1",
	classifiers=[
		"Development Status :: 5 - Production/Stable",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)",
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
		"python": [
			"*.pyd",  # Include all .pyd files (radia.pyd, rad_ngsolve.pyd)
			"radia_pyvista_viewer.py",
			"radia_vtk_export.py",
			"nastran_reader.py",
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
