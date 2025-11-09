"""
Radia Python Module Test Suite
Tests basic functionality of the radia.pyd module
"""

import sys
import os

# Set UTF-8 encoding for output
if sys.platform == 'win32':
	import codecs

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
build_dir = project_root / 'build' / 'Release'
if build_dir.exists():
	sys.path.insert(0, str(build_dir))

sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

def test_import():
	"""Test 1: Module import"""
	print("=" * 60)
	print("Test 1: Module Import")
	print("=" * 60)
	try:
		import radia as rad
		print("[OK] SUCCESS: radia module imported successfully")
		print(f"  Module location: {rad.__file__}")
		return rad
	except ImportError as e:
		print(f"[FAIL] FAILED: Cannot import radia module")
		print(f"  Error: {e}")
		print(f"  sys.path: {sys.path}")
		return None

def test_version(rad):
	"""Test 2: Version information"""
	print("\n" + "=" * 60)
	print("Test 2: Version Information")
	print("=" * 60)
	try:
		version = rad.UtiVer()
		print(f"[OK] SUCCESS: Radia version: {version}")
		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Cannot get version")
		print(f"  Error: {e}")
		return False

def test_basic_geometry(rad):
	"""Test 3: Basic geometry creation"""
	print("\n" + "=" * 60)
	print("Test 3: Basic Geometry Creation")
	print("=" * 60)
	try:
		# Create a simple rectangular block with magnetization
		# Parameters: center point, dimensions
		block = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
		print(f"[OK] SUCCESS: Created rectangular block")
		print(f"  Block ID: {block}")

		# Set magnetization
		magnetization = [0, 0, 1000]  # 1000 A/m in z-direction
		rad.ObjSetM(block, magnetization)
		print(f"[OK] SUCCESS: Set magnetization to {magnetization}")

		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Cannot create basic geometry")
		print(f"  Error: {e}")
		return False

def test_material(rad):
	"""Test 4: Material definition"""
	print("\n" + "=" * 60)
	print("Test 4: Material Definition")
	print("=" * 60)
	try:
		# Create a material (Steel37 equivalent using MatSatIsoFrm)
		# Steel37: C<0.13% steel with nonlinear magnetic properties
		mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
		print(f"[OK] SUCCESS: Created material")
		print(f"  Material ID: {mat}")

		# Create object and apply material
		block = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
		rad.MatApl(block, mat)
		print(f"[OK] SUCCESS: Applied material to object")

		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Cannot create/apply material")
		print(f"  Error: {e}")
		return False

def test_field_calculation(rad):
	"""Test 5: Magnetic field calculation"""
	print("\n" + "=" * 60)
	print("Test 5: Magnetic Field Calculation")
	print("=" * 60)
	try:
		# Create a simple magnet
		block = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
		rad.ObjSetM(block, [0, 0, 1000])

		# Calculate field at a point
		point = [0, 0, 20]  # 20mm above the magnet
		field = rad.Fld(block, 'b', point)
		print(f"[OK] SUCCESS: Calculated magnetic field")
		print(f"  Point: {point} mm")
		print(f"  Field: Bx={field[0]:.6f}, By={field[1]:.6f}, Bz={field[2]:.6f} T")

		# Verify field is reasonable (should be positive in z-direction)
		if field[2] > 0:
			print(f"[OK] Field direction is correct (Bz > 0)")
		else:
			print(f"[WARN] Warning: Field direction unexpected")

		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Cannot calculate field")
		print(f"  Error: {e}")
		return False

def test_solve(rad):
	"""Test 6: Relaxation/Solve"""
	print("\n" + "=" * 60)
	print("Test 6: Relaxation Solver")
	print("=" * 60)
	try:
		# Create a simple magnetic system
		block = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
		mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])  # Steel37 equivalent
		rad.MatApl(block, mat)

		# Solve the system
		precision = 0.001
		max_iter = 1000
		result = rad.Solve(block, precision, max_iter)
		print(f"[OK] SUCCESS: Solver completed")
		print(f"  Result: {result}")

		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Solver error")
		print(f"  Error: {e}")
		return False

def test_transformation(rad):
	"""Test 7: Geometric transformation"""
	print("\n" + "=" * 60)
	print("Test 7: Geometric Transformation")
	print("=" * 60)
	try:
		# Create object
		block = rad.ObjRecMag([0, 0, 0], [10, 10, 10])

		# Create translation transformation
		trans = rad.TrfTrsl([10, 0, 0])
		print(f"[OK] SUCCESS: Created translation transformation")

		# Apply transformation
		rad.TrfOrnt(block, trans)
		print(f"[OK] SUCCESS: Applied transformation to object")

		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Transformation error")
		print(f"  Error: {e}")
		return False

def test_cleanup(rad):
	"""Test 8: Memory cleanup"""
	print("\n" + "=" * 60)
	print("Test 8: Memory Cleanup")
	print("=" * 60)
	try:
		# Delete all objects
		rad.UtiDelAll()
		print(f"[OK] SUCCESS: All objects deleted")
		return True
	except Exception as e:
		print(f"[FAIL] FAILED: Cleanup error")
		print(f"  Error: {e}")
		return False

def main():
	"""Run all tests"""
	print("\n")
	print("╔" + "=" * 58 + "╗")
	print("║" + " " * 14 + "RADIA MODULE TEST SUITE" + " " * 21 + "║")
	print("╚" + "=" * 58 + "╝")
	print()

	# Test 1: Import
	rad = test_import()
	if rad is None:
		print("\n[FAIL] CRITICAL: Cannot proceed without module import")
		return False

	# Run tests
	tests = [
		("Version", test_version),
		("Basic Geometry", test_basic_geometry),
		("Material", test_material),
		("Field Calculation", test_field_calculation),
		("Solver", test_solve),
		("Transformation", test_transformation),
		("Cleanup", test_cleanup),
	]

	results = []
	for name, test_func in tests:
		try:
			result = test_func(rad)
			results.append((name, result))
		except Exception as e:
			print(f"\n[FAIL] EXCEPTION in {name}: {e}")
			results.append((name, False))

	# Summary
	print("\n" + "=" * 60)
	print("TEST SUMMARY")
	print("=" * 60)

	passed = sum(1 for _, result in results if result)
	total = len(results)

	for name, result in results:
		status = "[PASS]" if result else "[FAIL]"
		print(f"  {status}: {name}")

	print("-" * 60)
	print(f"  Total: {passed}/{total} tests passed ({passed/total*100:.1f}%)")
	print("=" * 60)

	if passed == total:
		print("\n*** ALL TESTS PASSED! ***")
		return True
	else:
		print(f"\n[WARNING] {total - passed} test(s) failed")
		return False

if __name__ == "__main__":
	success = main()
	sys.exit(0 if success else 1)
