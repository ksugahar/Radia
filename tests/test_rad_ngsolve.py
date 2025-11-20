"""
rad_ngsolve Module Test
Tests the NGSolve CoefficientFunction integration with Radia
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

# Find project root and add build directories to path
current_file = Path(__file__).resolve()
if 'tests' in current_file.parts:
	tests_index = current_file.parts.index('tests')
	project_root = Path(*current_file.parts[:tests_index])
else:
	project_root = current_file.parent

# Add build directories
build_dir = project_root / 'build' / 'lib' / 'Release'
if build_dir.exists():
	sys.path.insert(0, str(build_dir))

dist_dir = project_root / 'dist'
if dist_dir.exists():
	sys.path.insert(0, str(dist_dir))

# Add NGSolve module build directory
ngsolve_build_dir = project_root / 'build' / 'Release'
if ngsolve_build_dir.exists():
	sys.path.insert(0, str(ngsolve_build_dir))


def check_ngsolve_available():
	"""Check if NGSolve is installed"""
	try:
	    import ngsolve
	    return True
	except ImportError:
	    return False


def check_rad_ngsolve_available():
	"""Check if rad_ngsolve module is built"""
	try:
	    import ngsolve  # Must import first
	    import rad_ngsolve
	    return True
	except ImportError:
	    return False


@pytest.mark.skipif(not check_ngsolve_available(),
	               reason="NGSolve not installed")
@pytest.mark.skipif(not check_rad_ngsolve_available(),
	               reason="rad_ngsolve module not built")
class TestRadNGSolve:
	"""Test suite for rad_ngsolve module"""

	def test_import(self):
	    """Test 1: Module import"""
	    print("\n[Test 1] Importing rad_ngsolve...")

	    import ngsolve
	    print("  [OK] ngsolve imported")

	    import rad_ngsolve
	    print(f"  [OK] rad_ngsolve imported from {rad_ngsolve.__file__}")

	    # Check available functions
	    funcs = [name for name in dir(rad_ngsolve) if not name.startswith('_')]
	    assert 'RadiaField' in funcs, "RadiaField not found"
	    print(f"  [OK] Available functions: {funcs}")

	def test_coefficient_function_type(self):
	    """Test 2: CoefficientFunction type check"""
	    print("\n[Test 2] Checking CoefficientFunction type...")

	    import ngsolve
	    from ngsolve import CoefficientFunction
	    import rad_ngsolve

	    # Create RadiaField with dummy Radia object ID
	    bf = rad_ngsolve.RadiaField(1, 'b')

	    assert isinstance(bf, CoefficientFunction), "RadiaField is not a CoefficientFunction"
	    print(f"  [OK] RadiaField is CoefficientFunction: {type(bf)}")

	    # Test other field types
	    hf = rad_ngsolve.RadiaField(1, 'h')
	    assert isinstance(hf, CoefficientFunction), "RadiaField('h') is not a CoefficientFunction"
	    print(f"  [OK] RadiaField('h') is CoefficientFunction")

	    af = rad_ngsolve.RadiaField(1, 'a')
	    assert isinstance(af, CoefficientFunction), "RadiaField('a') is not a CoefficientFunction"
	    print(f"  [OK] RadiaField('a') is CoefficientFunction")

	    mf = rad_ngsolve.RadiaField(1, 'm')
	    assert isinstance(mf, CoefficientFunction), "RadiaField('m') is not a CoefficientFunction"
	    print(f"  [OK] RadiaField('m') is CoefficientFunction")

	def test_integration_with_radia(self):
	    """Test 3: Integration with Radia magnetic field"""
	    print("\n[Test 3] Testing integration with Radia...")

	    import ngsolve
	    from ngsolve import CoefficientFunction
	    import rad_ngsolve
	    import radia as rad

	    # Set Radia to use meters (required for NGSolve integration)
	    rad.FldUnits('m')

	    # Create a simple Radia magnet with permanent magnet material
	    magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 0])
	    # NdFeB: Br=1.2T, Hc=900kA/m, magnetization axis in z-direction
	    rad.MatApl(magnet, rad.MatPM(1.2, 900000, [0, 0, 1]))
	    rad.Solve(magnet, 0.0001, 10000)
	    print(f"  [OK] Radia magnet created: ID={magnet}")

	    # Create CoefficientFunction
	    B_cf = rad_ngsolve.RadiaField(magnet, 'b')
	    print(f"  [OK] RadiaField CoefficientFunction created")

	    # Verify it's a CoefficientFunction
	    assert isinstance(B_cf, CoefficientFunction)
	    print(f"  [OK] Type verified: {type(B_cf)}")

	    # Verify Radia field values
	    B_center = rad.Fld(magnet, 'b', [0, 0, 0])
	    assert B_center[2] > 0.5, f"Expected Bz > 0.5T, got {B_center[2]}"
	    print(f"  [OK] Field at center: Bz = {B_center[2]:.4f} T")

	    # Cleanup
	    rad.UtiDelAll()
	    print(f"  [OK] Radia objects cleaned up")

	def test_all_field_types(self):
	    """Test 4: All field types (b, h, a, m)"""
	    print("\n[Test 4] Testing all field types...")

	    import ngsolve
	    from ngsolve import CoefficientFunction
	    import rad_ngsolve
	    import radia as rad

	    # Set Radia to use meters (required for NGSolve integration)
	    rad.FldUnits('m')

	    # Create magnet with permanent magnet material
	    magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 0])
	    # NdFeB: Br=1.2T, Hc=900kA/m
	    rad.MatApl(magnet, rad.MatPM(1.2, 900000, [0, 0, 1]))
	    rad.Solve(magnet, 0.0001, 10000)

	    # Test all field types
	    field_types = ['b', 'h', 'a', 'm']
	    for ftype in field_types:
	        field = rad_ngsolve.RadiaField(magnet, ftype)
	        assert isinstance(field, CoefficientFunction)
	        assert field.field_type == ftype
	        print(f"  [OK] RadiaField('{ftype}') works")

	    rad.UtiDelAll()


# Standalone test function for non-pytest execution
def run_standalone_test():
	"""Run standalone test without pytest"""
	print("=" * 70)
	print("rad_ngsolve Module Test")
	print("=" * 70)

	if not check_ngsolve_available():
	    print("\n[SKIP] NGSolve not installed")
	    print("Install with: pip install ngsolve")
	    return 1

	if not check_rad_ngsolve_available():
	    print("\n[SKIP] rad_ngsolve module not built")
	    print("Build with: cmake --build build --target rad_ngsolve")
	    return 1

	print("\n[OK] Prerequisites satisfied")

	# Run tests
	test = TestRadNGSolve()

	try:
	    test.test_import()
	    test.test_coefficient_function_type()
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
