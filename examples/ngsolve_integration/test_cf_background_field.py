"""
Test CoefficientFunction-based background field

Tests for radia.ObjBckgCF() functionality
"""

import sys
sys.path.insert(0, 'S:/radia/01_GitHub/src/python')

import radia as rad
from radia_ngsolve_field import (
	create_uniform_field,
	create_analytical_field,
	dipole_field,
	solenoid_field
)
import numpy as np


def test_uniform_field():
	"""Test 1: Uniform field equivalence"""
	print("\n" + "="*70)
	print("Test 1: Uniform Field Equivalence")
	print("="*70)

	# Create using standard method
	B0 = [0, 0, 1.0]  # 1 Tesla in z
	field1 = rad.ObjBckg(B0)

	# Create using new CF method
	field2 = create_uniform_field(B0[0], B0[1], B0[2])

	# Test points
	test_points = [
		[0, 0, 0],
		[10, 0, 0],
		[0, 10, 0],
		[0, 0, 10],
		[5, 5, 5]
	]

	print(f"Uniform field: Bz = {B0[2]} T")
	print(f"\nTesting at various points:")

	max_error = 0
	for pt in test_points:
		B1 = rad.Fld(field1, 'b', pt)
		B2 = rad.Fld(field2, 'b', pt)

		error = np.linalg.norm(np.array(B1) - np.array(B2))
		max_error = max(max_error, error)

		print(f"Point {pt}:")
		print(f"  ObjBckg:   {B1}")
		print(f"  ObjBckgCF: {B2}")
		print(f"  Error: {error:.2e} T")

	print(f"\nMaximum error: {max_error:.2e} T")

	if max_error < 1e-10:
		print("PASS: Uniform field equivalence")
		return True
	else:
		print("FAIL: Uniform field equivalence")
		return False


def test_linear_gradient():
	"""Test 2: Linear gradient field"""
	print("\n" + "="*70)
	print("Test 2: Linear Gradient Field")
	print("="*70)

	# Define linear gradient: Bz = 1.0 + 0.01 * z (in Tesla, z in mm)
	def gradient_field(x, y, z):
		B0 = 1.0  # T
		grad = 0.001  # T/mm
		return [0, 0, B0 + grad * z]

	field = create_analytical_field(gradient_field, unit='mm')

	# Test at various z positions
	test_z = [0, 10, 20, 50, 100]  # mm

	print("Field: Bz(z) = 1.0 + 0.001*z [T], z in mm")
	print("\nField values:")

	all_pass = True
	for z in test_z:
		B = rad.Fld(field, 'b', [0, 0, z])
		expected_Bz = 1.0 + 0.001 * z
		error = abs(B[2] - expected_Bz)

		print(f"z = {z:3.0f} mm: Bz = {B[2]:.6f} T, expected = {expected_Bz:.6f} T, error = {error:.2e}")

		if error > 1e-10:
			all_pass = False

	if all_pass:
		print("PASS: Linear gradient field")
		return True
	else:
		print("FAIL: Linear gradient field")
		return False


def test_quadrupole():
	"""Test 3: Quadrupole field"""
	print("\n" + "="*70)
	print("Test 3: Quadrupole Field")
	print("="*70)

	# Define quadrupole: Bx = g*y, By = g*x (g in T/m, x,y in m)
	gradient = 10.0  # T/m

	def quadrupole(x, y, z):
		# x, y, z in meters
		Bx = gradient * y
		By = gradient * x
		Bz = 0.0
		return [Bx, By, Bz]

	field = create_analytical_field(quadrupole, unit='m')

	# Test points (mm)
	test_points = [
		[0, 0, 0],	 # Origin: should be zero
		[10, 0, 0],	# x-axis: Bx=0, By=g*x
		[0, 10, 0],	# y-axis: Bx=g*y, By=0
		[10, 10, 0],   # Diagonal
	]

	print(f"Quadrupole field: Bx = {gradient}*y, By = {gradient}*x [T/m]")
	print("\nField values:")

	all_pass = True
	for pt in test_points:
		B = rad.Fld(field, 'b', pt)
		x_m = pt[0] / 1000.0
		y_m = pt[1] / 1000.0

		expected_Bx = gradient * y_m
		expected_By = gradient * x_m
		expected_Bz = 0.0

		error_x = abs(B[0] - expected_Bx)
		error_y = abs(B[1] - expected_By)
		error_z = abs(B[2] - expected_Bz)

		print(f"Point ({pt[0]}, {pt[1]}, {pt[2]}) mm:")
		print(f"  B = [{B[0]:.6f}, {B[1]:.6f}, {B[2]:.6f}] T")
		print(f"  Expected = [{expected_Bx:.6f}, {expected_By:.6f}, {expected_Bz:.6f}] T")
		print(f"  Error = [{error_x:.2e}, {error_y:.2e}, {error_z:.2e}] T")

		if max(error_x, error_y, error_z) > 1e-10:
			all_pass = False

	if all_pass:
		print("PASS: Quadrupole field")
		return True
	else:
		print("FAIL: Quadrupole field")
		return False


def test_with_magnet():
	"""Test 4: Background field + permanent magnet"""
	print("\n" + "="*70)
	print("Test 4: Background Field + Permanent Magnet")
	print("="*70)

	# Create uniform background field
	B_background = [0, 0, 0.5]  # 0.5 T in z
	field = create_uniform_field(*B_background)

	# Create permanent magnet
	magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.0])

	# Create containers
	container_bg = rad.ObjCnt([field])
	container_mag = rad.ObjCnt([magnet])
	container_both = rad.ObjCnt([magnet, field])

	# Evaluate at test point
	test_point = [0, 0, 20]  # 20 mm above magnet

	B_bg = rad.Fld(container_bg, 'b', test_point)
	B_mag = rad.Fld(container_mag, 'b', test_point)
	B_both = rad.Fld(container_both, 'b', test_point)

	print(f"Background field: {B_background} T")
	print(f"Test point: {test_point} mm")
	print(f"\nField values:")
	print(f"  Background only: {B_bg}")
	print(f"  Magnet only:	 {B_mag}")
	print(f"  Combined:		{B_both}")

	# Check superposition
	B_expected = [B_bg[i] + B_mag[i] for i in range(3)]
	error = np.linalg.norm(np.array(B_both) - np.array(B_expected))

	print(f"  Expected (sum):  {B_expected}")
	print(f"  Error: {error:.2e} T")

	if error < 1e-6:  # Relaxed tolerance for magnet field
		print("PASS: Superposition with magnet")
		return True
	else:
		print("FAIL: Superposition with magnet")
		return False


def test_dipole_field():
	"""Test 5: Dipole field"""
	print("\n" + "="*70)
	print("Test 5: Magnetic Dipole Field")
	print("="*70)

	# Create dipole at origin
	moment = [0, 0, 1.0]  # 1 A*m^2 in z
	field = dipole_field(moment, position=(0, 0, 0), unit='m')

	# Test along z-axis
	# For dipole along z at origin: Bz(0,0,z) = (μ₀/2π) * m / z³
	mu0_2pi = 2e-7  # T*m/A

	test_z = [0.1, 0.2, 0.5, 1.0]  # meters

	print("Dipole moment: [0, 0, 1.0] A*m^2")
	print("Testing along z-axis:")

	all_pass = True
	for z in test_z:
		B = rad.Fld(field, 'b', [0, 0, z*1000])  # Convert to mm
		expected_Bz = mu0_2pi * moment[2] / (z**3)
		error = abs(B[2] - expected_Bz) / expected_Bz * 100  # percent

		print(f"z = {z:.1f} m: Bz = {B[2]:.6e} T, expected = {expected_Bz:.6e} T, error = {error:.2f}%")

		if error > 5.0:  # 5% tolerance (numerical approximation)
			all_pass = False

	if all_pass:
		print("PASS: Dipole field")
		return True
	else:
		print("FAIL: Dipole field")
		return False


def main():
	"""Run all tests"""
	print("\n" + "="*70)
	print("RADIA CF BACKGROUND FIELD TEST SUITE")
	print("="*70)

	tests = [
		("Uniform field equivalence", test_uniform_field),
		("Linear gradient field", test_linear_gradient),
		("Quadrupole field", test_quadrupole),
		("Background + magnet", test_with_magnet),
		("Dipole field", test_dipole_field),
	]

	results = []
	for name, test_func in tests:
		try:
			passed = test_func()
			results.append((name, passed))
		except Exception as e:
			print(f"\nFAIL EXCEPTION in {name}: {e}")
			import traceback
			traceback.print_exc()
			results.append((name, False))

	# Summary
	print("\n" + "="*70)
	print("TEST SUMMARY")
	print("="*70)

	passed_count = sum(1 for _, passed in results if passed)
	total_count = len(results)

	for name, passed in results:
		status = "PASS" if passed else "FAIL"
		print(f"{status}: {name}")

	print(f"\nPassed: {passed_count}/{total_count}")

	if passed_count == total_count:
		print("\n*** ALL TESTS PASSED!")
		return 0
	else:
		print(f"\nFAIL {total_count - passed_count} TEST(S) FAILED")
		return 1


if __name__ == "__main__":
	exit(main())
