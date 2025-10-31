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
	    assert 'RadBfield' in funcs, "RadBfield not found"
	    assert 'RadHfield' in funcs, "RadHfield not found"
	    assert 'RadAfield' in funcs, "RadAfield not found"
	    print(f"  [OK] Available functions: {funcs}")

	def test_coefficient_function_type(self):
	    """Test 2: CoefficientFunction type check"""
	    print("\n[Test 2] Checking CoefficientFunction type...")

	    import ngsolve
	    from ngsolve import CoefficientFunction
	    import rad_ngsolve

	    # Create RadBfield with dummy Radia object ID
	    bf = rad_ngsolve.RadBfield(1)

	    assert isinstance(bf, CoefficientFunction), "RadBfield is not a CoefficientFunction"
	    print(f"  [OK] RadBfield is CoefficientFunction: {type(bf)}")

	    # Test other field types
	    hf = rad_ngsolve.RadHfield(1)
	    assert isinstance(hf, CoefficientFunction), "RadHfield is not a CoefficientFunction"
	    print(f"  [OK] RadHfield is CoefficientFunction: {type(hf)}")

	    af = rad_ngsolve.RadAfield(1)
	    assert isinstance(af, CoefficientFunction), "RadAfield is not a CoefficientFunction"
	    print(f"  [OK] RadAfield is CoefficientFunction: {type(af)}")

	def test_integration_with_radia(self):
	    """Test 3: Integration with Radia magnetic field"""
	    print("\n[Test 3] Testing integration with Radia...")

	    import ngsolve
	    from ngsolve import CoefficientFunction
	    import rad_ngsolve
	    import radia as rad

	    # Create a simple Radia magnet
	    magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.2])
	    rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
	    rad.Solve(magnet, 0.0001, 10000)
	    print(f"  [OK] Radia magnet created: ID={magnet}")

	    # Create CoefficientFunction
	    B_cf = rad_ngsolve.RadBfield(magnet)
	    print(f"  [OK] RadBfield CoefficientFunction created")

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

	def test_mesh_evaluation(self):
	    """Test 4: Field evaluation on NGSolve mesh"""
	    print("\n[Test 4] Testing field evaluation on mesh...")

	    import ngsolve
	    from ngsolve import Mesh, H1, GridFunction
	    from netgen.geom2d import unit_square
	    import rad_ngsolve
	    import radia as rad

	    # Create Radia magnet
	    magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1.2])
	    rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
	    rad.Solve(magnet, 0.0001, 10000)
	    print(f"  [OK] Radia magnet created")

	    # Create RadBfield CoefficientFunction
	    B_cf = rad_ngsolve.RadBfield(magnet)
	    print(f"  [OK] RadBfield created")

	    # Create simple 2D mesh
	    mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))
	    print(f"  [OK] Mesh created: {mesh.ne} elements")

	    # Create finite element space
	    fes = H1(mesh, order=2, dim=3)
	    print(f"  [OK] FE space created: {fes.ndof} DOFs")

	    # Evaluate field on mesh
	    gf = GridFunction(fes)
	    gf.Set(B_cf)
	    print(f"  [OK] Field evaluated on mesh")

	    # Get field data
	    import numpy as np
	    vec_data = gf.vec.FV().NumPy()
	    print(f"  [OK] Field range: [{vec_data.min():.6f}, {vec_data.max():.6f}]")

	    # Cleanup
	    rad.UtiDelAll()
	    print(f"  [OK] Test completed")

	def test_field_components(self):
	    """Test 5: Access individual field components"""
	    print("\n[Test 5] Testing field component access...")

	    import ngsolve
	    from ngsolve import Integrate, dx
	    from ngsolve import Mesh
	    from netgen.geom2d import unit_square
	    import rad_ngsolve
	    import radia as rad

	    # Create Radia magnet (with material to avoid solve error)
	    magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 20], [0, 0, 1.0])
	    rad.MatApl(magnet, rad.MatStd('NdFeB', 1.0))
	    rad.Solve(magnet, 0.0001, 10000)
	    print(f"  [OK] Radia magnet created")

	    # Create RadBfield
	    B_cf = rad_ngsolve.RadBfield(magnet)
	    print(f"  [OK] RadBfield created")

	    # Create mesh
	    mesh = Mesh(unit_square.GenerateMesh(maxh=0.5))
	    print(f"  [OK] Mesh created")

	    # Access individual components
	    Bx = B_cf[0]
	    By = B_cf[1]
	    Bz = B_cf[2]
	    print(f"  [OK] Component access: Bx, By, Bz")

	    # Compute field magnitude squared
	    B_squared = Bx**2 + By**2 + Bz**2
	    print(f"  [OK] Computed |B|^2")

	    # Integrate Bz over domain
	    integral_Bz = Integrate(Bz, mesh)
	    print(f"  [OK] Integral of Bz = {integral_Bz:.6e}")

	    # Cleanup
	    rad.UtiDelAll()
	    print(f"  [OK] Test completed")


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
	    print("Build with: .\\build_ngsolve.ps1")
	    return 1

	print("\n[OK] Prerequisites satisfied")

	# Run tests
	test = TestRadNGSolve()

	try:
	    test.test_import()
	    test.test_coefficient_function_type()
	    test.test_integration_with_radia()
	    test.test_mesh_evaluation()
	    test.test_field_components()

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
