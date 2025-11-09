#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Square coil magnetic field analytical validation test
Compares Radia with analytical Biot-Savart solution
"""

import sys
import os
from pathlib import Path
import numpy as np
import codecs

# Find project root (works from any test subdirectory)
current_file = Path(__file__).resolve()
if 'tests' in current_file.parts:
	tests_index = current_file.parts.index('tests')
	project_root = Path(*current_file.parts[:tests_index])
else:
	project_root = current_file.parent

# Add build directory to path
build_dir = project_root / 'build' / 'lib' / 'Release'
if build_dir.exists():
	sys.path.insert(0, str(build_dir))

# Configure UTF-8 output
sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

def analytical_square_coil_field(I, a, z):
	"""
	Calculate magnetic field on axis of square coil using analytical formula.

	User's formula: H = 8*I/(4*pi*sqrt(a^2+y^2)*a/sqrt(2*a^2+y^2))
	where y is NORMALIZED: y = z/a

	Args:
		I: Current in Amperes
		a: Half-side length in meters
		z: Axial distance in meters (on axis)

	Returns:
		Bz: Magnetic field in Tesla on axis
	"""
	mu0 = 4 * np.pi * 1e-7  # T·m/A

	# Normalized distance
	y = z / a

	# Analytical formula from user (with y = z/a)
	# H = 8*I/(4*pi) * 1/(sqrt(a^2+y^2)*a/sqrt(2*a^2+y^2))
	# With y = z/a:
	# H = 8*I/(4*pi) * 1/(sqrt(a^2+(z/a)^2)*a/sqrt(2*a^2+(z/a)^2))
	# Simplify: sqrt(a^2+(z/a)^2) = sqrt(a^2+z^2/a^2) = sqrt((a^4+z^2)/a^2) = sqrt(a^4+z^2)/a
	# Similarly: sqrt(2*a^2+(z/a)^2) = sqrt(2*a^4+z^2)/a
	# So: H = 8*I/(4*pi) * 1/((sqrt(a^4+z^2)/a)*a/(sqrt(2*a^4+z^2)/a))
	#       = 8*I/(4*pi) * sqrt(2*a^4+z^2)/(sqrt(a^4+z^2)*a^2)

	# Actually, let me use the direct form with y = z/a:
	# Hy = 1./sqrt(1+y.^2)./sqrt(2+y.^2)  (normalized form)
	# This suggests: H/H0 = 1/sqrt(1+y^2)/sqrt(2+y^2)
	# where H0 = 2*I/(pi*a) (characteristic field scale)

	# Try standard textbook formula for square coil:
	# Bz = (2*mu0*I*a^2) / (pi*(a^2+z^2)*sqrt(2*a^2+z^2))
	Bz_standard = (2 * mu0 * I * a**2) / (np.pi * (a**2 + z**2) * np.sqrt(2*a**2 + z**2))

	# Also try user's formula with different interpretations:
	# H = (2*I/(pi*a)) / (sqrt(1+y^2)*sqrt(2+y^2)) where y=z/a
	H = (2 * I / (np.pi * a)) / (np.sqrt(1 + y**2) * np.sqrt(2 + y**2))
	Bz_user = mu0 * H

	# Compare: are they equivalent?
	# Bz_standard = (2*mu0*I*a^2) / (pi*(a^2+z^2)*sqrt(2*a^2+z^2))
	#             = (2*mu0*I*a^2) / (pi*a^2*(1+(z/a)^2)*a*sqrt(2+( z/a)^2))
	#             = (2*mu0*I) / (pi*a*sqrt(1+y^2)*sqrt(2+y^2))
	#
	# Bz_user = mu0 * (2*I/(pi*a)) / (sqrt(1+y^2)*sqrt(2+y^2))
	#         = (2*mu0*I) / (pi*a*sqrt(1+y^2)*sqrt(2+y^2))
	#
	# They are THE SAME! Use standard formula.

	return Bz_standard

def test_square_coil_analytical():
	"""
	Test Radia square coil against analytical Biot-Savart solution.
	"""
	print("=" * 70)
	print("SQUARE COIL ANALYTICAL VALIDATION TEST")
	print("=" * 70)
	print("\nComparing Radia with Biot-Savart analytical solution\n")

	# Import Radia
	try:
		import radia as rad
		print("[OK] Radia imported successfully")
	except ImportError as e:
		print(f"[ERROR] Failed to import Radia: {e}")
		return False

	print("\n" + "-" * 70)
	print("Coil Configuration")
	print("-" * 70)

	# Square coil parameters
	a = 25.0  # mm (half side length, so total side = 50mm)
	I = 1.0   # A (total current)

	# Wire dimensions (thin wire approximation)
	wx = 0.1  # mm (wire thickness in x for vertical segments)
	wy = 0.1  # mm (wire thickness in y for horizontal segments)
	h = 0.1   # mm (wire height in z)

	print(f"  Shape: Square coil (4 straight segments)")
	print(f"  Half-side length a: {a} mm")
	print(f"  Total side length: {2*a} mm")
	print(f"  Total current I: {I} A")
	print(f"  Wire dimensions: {wx}×{wy}×{h} mm")

	# Calculate current density J = I / Area
	# For thin segments, area ≈ wy * h (horizontal) or wx * h (vertical)
	J_horizontal = I / (wy * h)  # A/mm²
	J_vertical = I / (wx * h)    # A/mm²

	print(f"  Current density J (horizontal): {J_horizontal:.3f} A/mm²")
	print(f"  Current density J (vertical): {J_vertical:.3f} A/mm²")

	print("\n" + "-" * 70)
	print("Creating square coil in Radia...")
	print("-" * 70)

	# Create 4 rectangular current segments forming a square
	# Counter-clockwise current (viewed from +z) produces +z field at center
	# Square corners at (±a, ±a, 0)

	# Bottom segment: from (-a,-a) to (+a,-a), current in +x direction
	seg_bottom = rad.ObjRecCur([0, -a, 0], [2*a-wx, wy, h], [J_horizontal, 0, 0])

	# Right segment: from (+a,-a) to (+a,+a), current in +y direction
	seg_right = rad.ObjRecCur([a, 0, 0], [wx, 2*a-wy, h], [0, J_vertical, 0])

	# Top segment: from (+a,+a) to (-a,+a), current in -x direction
	seg_top = rad.ObjRecCur([0, a, 0], [2*a-wx, wy, h], [-J_horizontal, 0, 0])

	# Left segment: from (-a,+a) to (-a,-a), current in -y direction
	seg_left = rad.ObjRecCur([-a, 0, 0], [wx, 2*a-wy, h], [0, -J_vertical, 0])

	# Create container
	coil = rad.ObjCnt([seg_bottom, seg_right, seg_top, seg_left])

	print(f"[OK] Square coil created (ID: {coil})")
	print(f"     4 segments in counter-clockwise configuration")

	# Test points along Z-axis
	print("\n" + "-" * 70)
	print("Calculating fields along Z-axis...")
	print("-" * 70)

	z_values_mm = [1, 2, 5, 10, 20, 30, 50, 75, 100]

	print(f"{'Z (mm)':<10} {'Radia Bz (mT)':<16} {'Analytical Bz (mT)':<20} {'Difference':<12} {'Error %':<10}")
	print("-" * 70)

	max_error_percent = 0.0
	total_abs_error = 0.0

	for z_mm in z_values_mm:
		# Radia calculation
		pt = [0, 0, z_mm]
		b_radia = rad.Fld(coil, 'b', pt)
		bz_radia_mT = b_radia[2] * 1000  # Convert T to mT

		# Analytical calculation
		a_m = a / 1000  # Convert mm to m
		z_m = z_mm / 1000  # Convert mm to m
		bz_analytical_T = analytical_square_coil_field(I, a_m, z_m)
		bz_analytical_mT = bz_analytical_T * 1000  # Convert T to mT

		# Compare
		diff = bz_radia_mT - bz_analytical_mT
		if abs(bz_analytical_mT) > 1e-6:
			error_percent = abs(diff / bz_analytical_mT) * 100
		else:
			error_percent = 0.0

		max_error_percent = max(max_error_percent, error_percent)
		total_abs_error += abs(diff)

		print(f"{z_mm:<10.1f} {bz_radia_mT:<16.6f} {bz_analytical_mT:<20.6f} {diff:<12.6f} {error_percent:<10.2f}")

	# Summary
	print("-" * 70)
	avg_abs_error = total_abs_error / len(z_values_mm)
	print(f"\nSummary:")
	print(f"  Maximum relative error: {max_error_percent:.2f}%")
	print(f"  Average absolute error: {avg_abs_error:.6f} mT")

	# Analysis
	print("\n" + "=" * 70)
	print("ANALYSIS")
	print("=" * 70)
	print(f"\nAnalytical formula: H = (2*I)/(pi*a) * 1/[sqrt(1+y²)*sqrt(2+y²)]")
	print(f"where y = z/a (normalized axial distance)")
	print(f"\nThis represents the Biot-Savart integral for a square current loop.")
	print(f"The formula is exact for an infinitely thin wire on the axis.")
	print(f"\nRadia uses finite wire thickness ({wx}×{wy}×{h} mm), so small")
	print(f"differences are expected, especially at small z.")

	# Determine pass/fail
	tolerance_percent = 10.0  # Allow 10% for finite wire effects
	passed = max_error_percent < tolerance_percent

	# Cleanup
	rad.UtiDelAll()

	print("\n" + "=" * 70)
	if passed:
		print(f"[PASS] Radia matches analytical solution within {tolerance_percent}% tolerance")
		print(f"       Maximum error: {max_error_percent:.2f}%")
		print("=" * 70)
	else:
		print(f"[FAIL] Radia differs from analytical by {max_error_percent:.2f}%")
		print(f"       Exceeds {tolerance_percent}% tolerance")
		print("=" * 70)

	# Use assertion for pytest
	assert passed, f"Radia differs from analytical by {max_error_percent:.2f}%, exceeds {tolerance_percent}% tolerance"

def main():
	"""Run the analytical validation test"""
	result = test_square_coil_analytical()

	if result:
		print("\n*** ANALYTICAL VALIDATION TEST PASSED ***\n")
		sys.exit(0)
	else:
		print("\n*** ANALYTICAL VALIDATION TEST FAILED ***\n")
		sys.exit(1)

if __name__ == '__main__':
	main()
