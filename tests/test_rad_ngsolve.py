"""
rad_ngsolve Module Test
Tests the NGSolve integration with Radia (RadiaField is now in the main radia module)
"""

import sys
import os
from pathlib import Path
import pytest

# Set UTF-8 encoding for output
if sys.platform == 'win32':
	import codecs
	sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
	sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

# Find project root and add src to path
current_file = Path(__file__).resolve()
if 'tests' in current_file.parts:
	tests_index = current_file.parts.index('tests')
	project_root = Path(*current_file.parts[:tests_index])
else:
	project_root = current_file.parent

src_dir = project_root / 'src'
if src_dir.exists():
	sys.path.insert(0, str(src_dir))


def check_ngsolve_available():
	"""Check if NGSolve is installed"""
	try:
		import ngsolve
		return True
	except ImportError:
		return False


@pytest.mark.skipif(not check_ngsolve_available(),
	               reason="NGSolve not installed")
class TestRadNGSolve:
	"""Test suite for RadiaField (integrated into radia module)"""

	def test_import(self):
		"""Test 1: Module import"""
		print("\n[Test 1] Importing RadiaField...")

		from radia import RadiaField
		print(f"  [OK] RadiaField imported")

		# Check it has the right attributes
		assert hasattr(RadiaField, '__call__')
		print(f"  [OK] RadiaField has expected methods")

	def test_radiafield_api(self):
		"""Test 2: RadiaField API check"""
		print("\n[Test 2] Checking RadiaField API...")

		from radia import RadiaField
		import radia as rad

		rad.UtiDelAll()
		magnet = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])

		# Create RadiaField
		bf = RadiaField(magnet, 'b')
		assert callable(bf), "RadiaField must be callable"
		assert hasattr(bf, 'dim'), "RadiaField must have dim attribute"
		assert bf.dim == 3
		print(f"  [OK] RadiaField('b') created, dim={bf.dim}")

		# Test other field types
		hf = RadiaField(magnet, 'h')
		assert callable(hf) and hf.dim == 3
		print(f"  [OK] RadiaField('h') works")

		af = RadiaField(magnet, 'a')
		assert callable(af) and af.dim == 3
		print(f"  [OK] RadiaField('a') works")

		mf = RadiaField(magnet, 'm')
		assert callable(mf) and mf.dim == 3
		print(f"  [OK] RadiaField('m') works")

		pf = RadiaField(magnet, 'phi')
		assert callable(pf) and pf.dim == 1
		print(f"  [OK] RadiaField('phi') works, dim={pf.dim}")

		rad.UtiDelAll()

	def test_integration_with_radia(self):
		"""Test 3: Integration with Radia magnetic field"""
		print("\n[Test 3] Testing integration with Radia...")

		from radia import RadiaField
		import radia as rad

		rad.UtiDelAll()

		# Create a simple permanent magnet (M in A/m, Br=1.2T -> M=954930 A/m)
		magnet = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])
		print(f"  [OK] Radia magnet created: ID={magnet}")

		# Create RadiaField
		B_cf = RadiaField(magnet, 'b')
		print(f"  [OK] RadiaField created")

		# Verify callable
		assert callable(B_cf)
		print(f"  [OK] RadiaField is callable")

		# Verify Radia field values at external point
		B_ext = rad.Fld(magnet, 'b', [0.05, 0, 0])
		assert abs(B_ext[2]) > 1e-5, f"Expected non-zero Bz, got {B_ext[2]}"
		print(f"  [OK] Field at (0.05,0,0): Bz = {B_ext[2]:.6e} T")

		# Cleanup
		rad.UtiDelAll()
		print(f"  [OK] Radia objects cleaned up")

	def test_all_field_types(self):
		"""Test 4: All field types (b, h, a, m, phi)"""
		print("\n[Test 4] Testing all field types...")

		from radia import RadiaField
		import radia as rad

		rad.UtiDelAll()

		magnet = rad.magnet_box([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])

		# Test all field types
		for ftype in ['b', 'h', 'a', 'm', 'phi']:
			field = RadiaField(magnet, ftype)
			assert callable(field)
			assert field.field_type == ftype
			expected_dim = 1 if ftype == 'phi' else 3
			assert field.dim == expected_dim
			print(f"  [OK] RadiaField('{ftype}') works, dim={field.dim}")

		rad.UtiDelAll()


# Standalone test function for non-pytest execution
def run_standalone_test():
	"""Run standalone test without pytest"""
	print("=" * 70)
	print("rad_ngsolve Module Test (Pure Python)")
	print("=" * 70)

	if not check_ngsolve_available():
		print("\n[SKIP] NGSolve not installed")
		print("Install with: pip install ngsolve")
		return 1

	print("\n[OK] Prerequisites satisfied")

	# Run tests
	test = TestRadNGSolve()

	try:
		test.test_import()
		test.test_radiafield_api()
		test.test_integration_with_radia()
		test.test_all_field_types()

		print("\n" + "=" * 70)
		print("[OK] ALL TESTS PASSED!")
		print("=" * 70)
		return 0

	except Exception as e:
		print(f"\n[FAIL] ERROR: {e}")
		import traceback
		traceback.print_exc()
		return 1


if __name__ == '__main__':
	sys.exit(run_standalone_test())
