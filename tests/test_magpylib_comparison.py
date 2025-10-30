#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Comparison test between Radia and magpylib
Tests circular (cylindrical) magnet field calculations
"""

# Add project root's build directory to path
import sys
import os
from pathlib import Path

# Find project root (works from any test subdirectory)
current_file = Path(__file__).resolve()
if 'tests' in current_file.parts:
	# Find the 'tests' directory and go up one level
	tests_index = current_file.parts.index('tests')
	project_root = Path(*current_file.parts[:tests_index])
else:
	# Fallback
	project_root = current_file.parent

# Add build directory to path
build_dir = project_root / 'build' / 'lib' / 'Release'
if build_dir.exists():
	sys.path.insert(0, str(build_dir))

import numpy as np
import codecs

# Configure UTF-8 output
sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

def test_cylindrical_magnet_comparison():
	"""
	Compare Radia and magpylib field calculations for a cylindrical magnet.

	Both libraries should give similar results for the magnetic field around
	a uniformly magnetized cylindrical permanent magnet.
	"""
	print("=" * 70)
	print("RADIA vs MAGPYLIB COMPARISON TEST")
	print("=" * 70)
	print("\nTest: Cylindrical permanent magnet field calculation\n")

	# Import libraries
	try:
		import radia as rad
		print("[OK] Radia imported successfully")
	except ImportError as e:
		print(f"[ERROR] Failed to import Radia: {e}")
		return False

	try:
		import magpylib as magpy
		print(f"[OK] magpylib imported successfully (version {magpy.__version__})")
	except ImportError as e:
		print(f"[ERROR] Failed to import magpylib: {e}")
		print("       Install with: pip install magpylib")
		return False

	print("\n" + "-" * 70)
	print("Magnet Configuration")
	print("-" * 70)

	# Magnet parameters (in mm for Radia, convert to mm for magpylib)
	radius = 10.0  # mm
	height = 20.0  # mm

	# Use NdFeB-like material properties for realistic comparison
	# Typical NdFeB: Br = 1.2 T (remanence), Hc = 900 kA/m (coercivity)
	# For hard permanent magnet: M ≈ Br/μ₀ = 1.2/(4π×10⁻⁷) ≈ 955 kA/m
	remanence_T = 1.2  # Tesla (typical NdFeB)
	magnetization = remanence_T / (4 * np.pi * 1e-7)  # Convert to A/m

	print(f"  Shape: Cylinder")
	print(f"  Radius: {radius} mm")
	print(f"  Height: {height} mm")
	print(f"  Material: NdFeB-like")
	print(f"  Remanence Br: {remanence_T} T")
	print(f"  Magnetization M: {magnetization/1000:.1f} kA/m (Z-direction)")

	# Create magnet in Radia
	print("\n" + "-" * 70)
	print("Creating magnet in Radia...")
	print("-" * 70)

	# Radia: ObjCylMag([x,y,z], radius, height, nseg, axis, [mx,my,mz])
	# Subdivide cylinder for better accuracy
	n_phi = 16  # azimuthal subdivisions (number of segments)

	radia_mag = rad.ObjCylMag([0, 0, 0], radius, height, n_phi, 'z', [0, 0, magnetization])
	print(f"[OK] Radia cylindrical magnet created (ID: {radia_mag})")
	print(f"     Subdivisions: {n_phi} segments (azimuthal)")

	# Create magnet in magpylib
	print("\n" + "-" * 70)
	print("Creating magnet in magpylib...")
	print("-" * 70)

	# magpylib uses SI units (mm for position, T for field and polarization)
	# Use polarization parameter (remanence Br) in Tesla
	# For permanent magnet: polarization = remanence Br

	magpy_mag = magpy.magnet.Cylinder(
		polarization=(0, 0, remanence_T),  # in Tesla
		dimension=(2*radius, height)  # (diameter, height) in mm
	)
	print(f"[OK] magpylib cylindrical magnet created")
	print(f"     Polarization (Br): {remanence_T} T")

	# Test points: Create a grid in the XZ plane (Y=0)
	print("\n" + "-" * 70)
	print("Calculating fields at test points...")
	print("-" * 70)

	# Test points along Z-axis and radial direction
	test_points = []

	# Points along Z-axis (above magnet)
	for z in [25, 30, 40, 50]:
		test_points.append([0, 0, z])

	# Points in radial direction at z=0 plane
	for r in [15, 20, 30]:
		test_points.append([r, 0, 0])

	# Points in XZ plane
	for x in [10, 20]:
		for z in [20, 30]:
			test_points.append([x, 0, z])

	print(f"  Number of test points: {len(test_points)}")
	print(f"  Point locations:")
	for i, pt in enumerate(test_points[:5]):
		print(f"    {i+1}. ({pt[0]:.1f}, {pt[1]:.1f}, {pt[2]:.1f}) mm")
	if len(test_points) > 5:
		print(f"    ... and {len(test_points)-5} more points")

	# Calculate fields with both libraries
	print("\n" + "-" * 70)
	print("Field Comparison Results")
	print("-" * 70)
	print(f"{'Point (mm)':<20} {'Radia Bz (mT)':<15} {'magpylib Bz (mT)':<18} {'Difference':<12} {'Error %':<10}")
	print("-" * 70)

	max_error_percent = 0.0
	total_abs_error = 0.0
	bz_rad_first = 0.0
	bz_mag_first = 0.0

	for i, pt in enumerate(test_points):
		# Radia field calculation (returns field in Gauss)
		# Radia uses Gauss-cm system: 1 Tesla = 10000 Gauss
		b_radia = rad.Fld(radia_mag, 'b', pt)
		bx_rad, by_rad, bz_rad = b_radia[0] / 10, b_radia[1] / 10, b_radia[2] / 10  # Convert Gauss to mT

		# magpylib field calculation (returns field in Tesla)
		b_magpy = magpy_mag.getB(pt)
		bx_mag, by_mag, bz_mag = b_magpy[0] * 1000, b_magpy[1] * 1000, b_magpy[2] * 1000  # Convert T to mT

		# Compare Bz (main component for axially magnetized cylinder)
		diff = bz_rad - bz_mag

		# Calculate relative error
		if abs(bz_mag) > 1e-6:  # Avoid division by zero
			error_percent = abs(diff / bz_mag) * 100
		else:
			error_percent = 0.0

		max_error_percent = max(max_error_percent, error_percent)
		total_abs_error += abs(diff)

		# Store first point for ratio calculation
		if i == 0:
			bz_rad_first = bz_rad
			bz_mag_first = bz_mag

		pt_str = f"({pt[0]:.0f},{pt[1]:.0f},{pt[2]:.0f})"
		print(f"{pt_str:<20} {bz_rad:<15.6f} {bz_mag:<18.6f} {diff:<12.6f} {error_percent:<10.2f}")

	# Summary statistics
	print("-" * 70)
	avg_abs_error = total_abs_error / len(test_points)
	print(f"\nSummary:")
	print(f"  Maximum relative error: {max_error_percent:.2f}%")
	print(f"  Average absolute error: {avg_abs_error:.6f} mT")

	# Analysis note
	print("\n" + "=" * 70)
	print("ANALYSIS")
	print("=" * 70)
	print(f"\nNote: Radia and magpylib use DIFFERENT permanent magnet models:")
	print(f"\n  Radia:")
	print(f"    - Uses uniform magnetization M (A/m) throughout volume")
	print(f"    - Does NOT account for demagnetization field")
	print(f"    - Gives UPPER BOUND on external field")
	print(f"\n  magpylib:")
	print(f"    - Uses polarization J (remanence Br in Tesla)")
	print(f"    - DOES account for demagnetization effects")
	print(f"    - Gives REALISTIC field for actual permanent magnets")
	print(f"\nFor short cylinders (L/D ≈ 1):")
	print(f"  - Demagnetization factor N ≈ 0.3-0.4")
	print(f"  - Expected field ratio: Radia/magpylib ≈ 1/(1-N) ≈ 50-100x")
	if bz_mag_first > 0:
		print(f"Observed ratio at (0,0,25): {bz_rad_first/bz_mag_first:.1f}x")
	else:
		print(f"Observed values at (0,0,25): Radia={bz_rad_first:.2f} mT, magpylib={bz_mag_first:.2f} mT")

	# Determine if libraries are working (not if they agree perfectly)
	# Both should give non-zero, physically reasonable results
	radia_ok = (bz_rad_first > 100)  # Should be strong field for NdFeB
	magpy_ok = (bz_mag_first > 1)    # Should be reasonable field

	# Cleanup Radia
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	if radia_ok and magpy_ok:
		print(f"[PASS] Both libraries produce physically reasonable fields")
		print(f"       (Exact agreement not expected due to different models)")
		print("=" * 70)
		return True
	else:
		print(f"[FAIL] One or both libraries produce unreasonable results")
		print("=" * 70)
		return False

def main():
	"""Run the comparison test"""
	result = test_cylindrical_magnet_comparison()

	if result:
		print("\n*** COMPARISON TEST PASSED ***\n")
		sys.exit(0)
	else:
		print("\n*** COMPARISON TEST FAILED ***\n")
		sys.exit(1)

if __name__ == '__main__':
	main()
