"""
Advanced Radia Test
Tests more complex magnetic field calculations
"""

import sys
import os
import math

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

# Configure UTF-8 output
import codecs
sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

def test_dipole_magnet():
	"""Create and test a simple dipole magnet"""
	print("\n" + "=" * 60)
	print("TEST: Dipole Magnet")
	print("=" * 60)

	import radia as rad

	# Create two blocks with opposite magnetization (dipole)
	print("\n1. Creating dipole magnet...")
	gap = 20  # mm
	block_size = [50, 50, 20]

	# Upper pole (magnetized downward)
	upper = rad.ObjRecMag([0, 0, gap/2 + block_size[2]/2], block_size)
	rad.ObjSetM(upper, [0, 0, -1000])
	print(f"   Upper pole created: ID={upper}")

	# Lower pole (magnetized upward)
	lower = rad.ObjRecMag([0, 0, -(gap/2 + block_size[2]/2)], block_size)
	rad.ObjSetM(lower, [0, 0, 1000])
	print(f"   Lower pole created: ID={lower}")

	# Create container
	container = rad.ObjCnt([upper, lower])
	print(f"   Container created: ID={container}")

	# Calculate field in the gap
	print("\n2. Calculating field along Z-axis in the gap...")
	z_points = [-10, -5, 0, 5, 10]
	print(f"\n   {'Z (mm)':>10} {'Bx (T)':>12} {'By (T)':>12} {'Bz (T)':>12}")
	print(f"   {'-'*10} {'-'*12} {'-'*12} {'-'*12}")

	for z in z_points:
		field = rad.Fld(container, 'b', [0, 0, z])
		print(f"   {z:10.1f} {field[0]:12.6f} {field[1]:12.6f} {field[2]:12.6f}")

	# Check field uniformity
	center_field = rad.Fld(container, 'b', [0, 0, 0])
	print(f"\n   ✓ Center field (Bz): {center_field[2]:.6f} T")

	# Cleanup
	rad.UtiDelAll()
	print(f"\n   ✓ Test completed")

def test_quadrupole():
	"""Create and test a quadrupole magnet"""
	print("\n" + "=" * 60)
	print("TEST: Quadrupole Magnet")
	print("=" * 60)

	import radia as rad

	print("\n1. Creating quadrupole magnet (4 poles)...")
	size = [20, 20, 100]
	offset = 25  # mm from center

	# Create 4 blocks arranged in quadrupole configuration
	# Magnetization pattern creates quadrupole field
	poles = []

	# Pole 1: +X position, magnetized in +Y
	p1 = rad.ObjRecMag([offset, 0, 0], size)
	rad.ObjSetM(p1, [0, 1000, 0])
	poles.append(p1)

	# Pole 2: +Y position, magnetized in -X
	p2 = rad.ObjRecMag([0, offset, 0], size)
	rad.ObjSetM(p2, [-1000, 0, 0])
	poles.append(p2)

	# Pole 3: -X position, magnetized in -Y
	p3 = rad.ObjRecMag([-offset, 0, 0], size)
	rad.ObjSetM(p3, [0, -1000, 0])
	poles.append(p3)

	# Pole 4: -Y position, magnetized in +X
	p4 = rad.ObjRecMag([0, -offset, 0], size)
	rad.ObjSetM(p4, [1000, 0, 0])
	poles.append(p4)

	container = rad.ObjCnt(poles)
	print(f"   Quadrupole created with 4 poles")

	# Calculate field along X-axis
	print("\n2. Calculating field along X-axis...")
	x_points = [0, 2, 4, 6, 8, 10]
	print(f"\n   {'X (mm)':>10} {'Bx (T)':>12} {'By (T)':>12} {'Gradient':>12}")
	print(f"   {'-'*10} {'-'*12} {'-'*12} {'-'*12}")

	prev_bx = None
	for x in x_points:
		field = rad.Fld(container, 'b', [x, 0, 0])
		if prev_bx is not None and x > 0:
			gradient = (field[0] - prev_bx) / 2.0  # T/mm
		else:
			gradient = 0
		print(f"   {x:10.1f} {field[0]:12.6f} {field[1]:12.6f} {gradient:12.6f}")
		prev_bx = field[0]

	# Cleanup
	rad.UtiDelAll()
	print(f"\n   ✓ Test completed")

def test_iron_core():
	"""Test with iron core (nonlinear material)"""
	print("\n" + "=" * 60)
	print("TEST: Iron Core Magnet")
	print("=" * 60)

	import radia as rad

	print("\n1. Creating iron core with coil...")

	# Create iron core
	core = rad.ObjRecMag([0, 0, 0], [40, 40, 100])
	print(f"   Core created: ID={core}")

	# Create iron material
	mat_iron = rad.MatStd('Steel37', 2000)
	rad.MatApl(core, mat_iron)
	print(f"   Iron material applied")

	# Set external field (simulating coil)
	rad.ObjSetM(core, [0, 0, 500])
	print(f"   External field applied")

	print("\n2. Solving with relaxation...")
	precision = 0.001
	max_iter = 1000
	result = rad.Solve(core, precision, max_iter)
	print(f"   Solver result: {result}")

	print("\n3. Calculating field...")
	points = [[0, 0, 0], [0, 0, 60]]
	for pt in points:
		field = rad.Fld(core, 'b', pt)
		print(f"   Point {pt}: Bz = {field[2]:.6f} T")

	# Cleanup
	rad.UtiDelAll()
	print(f"\n   ✓ Test completed")

def test_field_integral():
	"""Test field integral calculation"""
	print("\n" + "=" * 60)
	print("TEST: Field Integral")
	print("=" * 60)

	import radia as rad

	print("\n1. Creating simple magnet...")
	magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 50])
	rad.ObjSetM(magnet, [0, 0, 1000])

	print("\n2. Calculating field integrals...")

	# Field at several points
	z_start = -30
	z_end = 30
	n_points = 21

	print(f"\n   Integrating field from Z={z_start} to Z={z_end} mm")

	integral_bz = 0.0
	dz = (z_end - z_start) / (n_points - 1)

	for i in range(n_points):
		z = z_start + i * dz
		field = rad.Fld(magnet, 'b', [0, 0, z])
		integral_bz += field[2] * dz

	print(f"   ∫Bz·dz = {integral_bz:.6f} T·mm = {integral_bz/1000:.6f} T·m")

	# Cleanup
	rad.UtiDelAll()
	print(f"\n   ✓ Test completed")

def main():
	"""Run all advanced tests"""
	print("\n")
	print("╔" + "=" * 58 + "╗")
	print("║" + " " * 12 + "RADIA ADVANCED TEST SUITE" + " " * 19 + "║")
	print("╚" + "=" * 58 + "╝")

	try:
		import radia as rad
		print(f"\nRadia version: {rad.UtiVer()}")

		tests = [
			test_dipole_magnet,
			test_quadrupole,
			test_iron_core,
			test_field_integral,
		]

		for test in tests:
			try:
				test()
			except Exception as e:
				print(f"\n✗ Test failed: {e}")
				import traceback
				traceback.print_exc()

		print("\n" + "=" * 60)
		print("✓ ADVANCED TESTS COMPLETED")
		print("=" * 60)

	except ImportError as e:
		print(f"\n✗ ERROR: Cannot import radia module")
		print(f"  {e}")
		sys.exit(1)

if __name__ == "__main__":
	main()
