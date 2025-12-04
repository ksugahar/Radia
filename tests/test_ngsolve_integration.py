"""
NGSolve Integration Test Suite

Tests the NGSolve integration with Radia according to CLAUDE.md best practices:
- rad.FldUnits('m') is REQUIRED for NGSolve integration
- HDiv(mesh, order=2) for best accuracy
- Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- Use CoefficientFunction directly for maximum accuracy near boundaries

This test suite validates:
1. Module import and CoefficientFunction creation
2. Field types (b, h, a, m)
3. HDiv function space integration
4. Field accuracy at various distances from magnet
5. Comparison between GridFunction and direct Radia evaluation
"""

import sys
import os
from pathlib import Path
import pytest
import numpy as np

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

# Add build directories (multiple locations for .pyd files)
build_release = project_root / 'build' / 'Release'
if build_release.exists():
    sys.path.insert(0, str(build_release))

build_lib_radia = project_root / 'build' / 'lib' / 'radia'
if build_lib_radia.exists():
    sys.path.insert(0, str(build_lib_radia))

src_radia = project_root / 'src' / 'radia'
if src_radia.exists():
    sys.path.insert(0, str(src_radia))

src_python = project_root / 'src' / 'python'
if src_python.exists():
    sys.path.insert(0, str(src_python))


def check_ngsolve_available():
    """Check if NGSolve is installed"""
    try:
        import ngsolve
        return True
    except ImportError:
        return False


def check_radia_ngsolve_available():
    """Check if radia_ngsolve module is built"""
    try:
        import ngsolve  # Must import first
        import radia_ngsolve
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not check_ngsolve_available(),
                   reason="NGSolve not installed")
@pytest.mark.skipif(not check_radia_ngsolve_available(),
                   reason="radia_ngsolve module not built")
class TestNGSolveIntegration:
    """Test suite for NGSolve integration following CLAUDE.md best practices"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup: import modules and create test magnet"""
        import radia as rad
        import radia_ngsolve
        from ngsolve import Mesh, HDiv, GridFunction, CoefficientFunction
        from netgen.csg import CSGeometry, OrthoBrick, Pnt

        self.rad = rad
        self.radia_ngsolve = radia_ngsolve
        self.Mesh = Mesh
        self.HDiv = HDiv
        self.GridFunction = GridFunction
        self.CoefficientFunction = CoefficientFunction
        self.CSGeometry = CSGeometry
        self.OrthoBrick = OrthoBrick
        self.Pnt = Pnt

        # CRITICAL: Set units to meters (REQUIRED for NGSolve integration)
        rad.UtiDelAll()
        rad.FldUnits('m')

        # Create test magnet (permanent magnet)
        self.magnet_center = [0, 0, 0]  # meters
        self.magnet_size = [0.020, 0.020, 0.030]  # 20mm x 20mm x 30mm
        self.magnet = rad.ObjRecMag(
            self.magnet_center,
            self.magnet_size,
            [0, 0, 1.2]  # Magnetization 1.2 T in z-direction
        )
        # Apply NdFeB material
        rad.MatApl(self.magnet, rad.MatPM(1.2, 900000, [0, 0, 1]))
        rad.Solve(self.magnet, 0.0001, 10000)

        yield

        # Cleanup
        rad.UtiDelAll()

    def test_units_are_meters(self):
        """Test 1: Verify FldUnits is set to meters"""
        print("\n[Test 1] Verifying FldUnits('m') is set")

        units_str = self.rad.FldUnits()
        assert 'Length:  m' in units_str, f"Expected meters, got: {units_str}"
        print("  [OK] Units set to meters")

    def test_radiafield_is_coefficientfunction(self):
        """Test 2: RadiaField returns NGSolve CoefficientFunction"""
        print("\n[Test 2] Checking RadiaField type")

        B_cf = self.radia_ngsolve.RadiaField(self.magnet, 'b')
        assert isinstance(B_cf, self.CoefficientFunction), \
            f"Expected CoefficientFunction, got {type(B_cf)}"
        print(f"  [OK] RadiaField is CoefficientFunction: {type(B_cf)}")

    def test_all_field_types(self):
        """Test 3: All field types (b, h, a, m) work correctly"""
        print("\n[Test 3] Testing all field types")

        field_types = ['b', 'h', 'a', 'm']
        for ftype in field_types:
            field = self.radia_ngsolve.RadiaField(self.magnet, ftype)
            assert isinstance(field, self.CoefficientFunction)
            assert field.field_type == ftype
            print(f"  [OK] RadiaField('{ftype}') works")

    def test_hdiv_gridfunction_projection(self):
        """Test 4: HDiv GridFunction projection (CLAUDE.md best practice)"""
        print("\n[Test 4] HDiv GridFunction projection (order=2)")

        # Create mesh outside magnet region (evaluate in air only)
        geo = self.CSGeometry()
        # Mesh region offset from magnet
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.03, -0.03),
            self.Pnt(0.08, 0.03, 0.03)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.01))

        print(f"  Mesh: {mesh.ne} elements, {mesh.nv} vertices")

        # CLAUDE.md recommends HDiv with order=2
        fes = self.HDiv(mesh, order=2)
        B_gf = self.GridFunction(fes)

        # Create CoefficientFunction
        B_cf = self.radia_ngsolve.RadiaField(self.magnet, 'b')

        # Project to GridFunction
        B_gf.Set(B_cf)

        print(f"  FES DOFs: {fes.ndof}")
        print("  [OK] HDiv GridFunction projection successful")

    def test_field_accuracy_far_from_magnet(self):
        """Test 5: Field accuracy at distance > 1 mesh cell from magnet"""
        print("\n[Test 5] Field accuracy at distance from magnet surface")

        # Create mesh far from magnet
        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.04, -0.02, -0.02),
            self.Pnt(0.08, 0.02, 0.02)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.008))

        # Use HDiv order=2 as recommended
        fes = self.HDiv(mesh, order=2)
        B_gf = self.GridFunction(fes)
        B_cf = self.radia_ngsolve.RadiaField(self.magnet, 'b')
        B_gf.Set(B_cf)

        # Test points far from magnet (> 1 mesh cell = 8mm from surface)
        # Magnet surface is at x = 0.01m, so test at x >= 0.05m
        test_points = [
            (0.05, 0.0, 0.0),   # 40mm from magnet center, 30mm from surface
            (0.06, 0.0, 0.0),
            (0.07, 0.0, 0.0),
        ]

        max_rel_error = 0.0
        print(f"  {'Point':<25s} {'Radia Bx':>12s} {'NGSolve Bx':>12s} {'Error %':>10s}")
        print("  " + "-" * 65)

        for pt in test_points:
            # Direct Radia evaluation
            B_radia = self.rad.Fld(self.magnet, 'b', list(pt))

            # NGSolve GridFunction evaluation
            B_ngsolve = B_gf(mesh(*pt))

            # Calculate relative error for Bx component
            if abs(B_radia[0]) > 1e-6:
                rel_error = abs(B_radia[0] - B_ngsolve[0]) / abs(B_radia[0]) * 100
            else:
                rel_error = abs(B_radia[0] - B_ngsolve[0]) * 100

            max_rel_error = max(max_rel_error, rel_error)
            print(f"  {str(pt):<25s} {B_radia[0]:>12.6f} {B_ngsolve[0]:>12.6f} {rel_error:>9.2f}%")

        # At > 1 mesh cell distance, error should be < 5%
        assert max_rel_error < 5.0, f"Max relative error {max_rel_error:.2f}% exceeds 5%"
        print(f"  [OK] Max relative error: {max_rel_error:.2f}% (< 5%)")

    def test_coefficient_function_direct_evaluation(self):
        """Test 6: Direct CoefficientFunction evaluation (most accurate)"""
        print("\n[Test 6] Direct CoefficientFunction evaluation")

        # CLAUDE.md states: Use CoefficientFunction directly for maximum
        # accuracy near boundaries

        B_cf = self.radia_ngsolve.RadiaField(self.magnet, 'b')

        # Create a simple mesh to get mesh object
        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.02, -0.02),
            self.Pnt(0.06, 0.02, 0.02)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.01))

        # Test point closer to magnet surface (within 1 mesh cell)
        # This is where direct CoefficientFunction should be used
        test_point = (0.025, 0.0, 0.0)  # 15mm from magnet center

        # Direct Radia evaluation
        B_radia = self.rad.Fld(self.magnet, 'b', list(test_point))

        # Direct CoefficientFunction evaluation
        B_cf_val = B_cf(mesh(*test_point))

        print(f"  Test point: {test_point}")
        print(f"  Radia B: [{B_radia[0]:.6f}, {B_radia[1]:.6f}, {B_radia[2]:.6f}]")
        print(f"  CF B:    [{B_cf_val[0]:.6f}, {B_cf_val[1]:.6f}, {B_cf_val[2]:.6f}]")

        # Direct CoefficientFunction should match Radia exactly
        for i in range(3):
            diff = abs(B_radia[i] - B_cf_val[i])
            assert diff < 1e-6, f"Component {i}: diff={diff}"

        print("  [OK] CoefficientFunction matches Radia directly")

    def test_field_type_attribute(self):
        """Test 7: RadiaField has field_type attribute"""
        print("\n[Test 7] field_type attribute")

        for ftype in ['b', 'h', 'a', 'm']:
            field = self.radia_ngsolve.RadiaField(self.magnet, ftype)
            assert hasattr(field, 'field_type')
            assert field.field_type == ftype
            print(f"  [OK] RadiaField('{ftype}').field_type = '{field.field_type}'")


@pytest.mark.skipif(not check_ngsolve_available(),
                   reason="NGSolve not installed")
@pytest.mark.skipif(not check_radia_ngsolve_available(),
                   reason="radia_ngsolve module not built")
class TestNGSolveFunctionSpaces:
    """Test different NGSolve function spaces"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup for function space tests"""
        import radia as rad
        import radia_ngsolve
        from ngsolve import Mesh, HDiv, HCurl, VectorH1, GridFunction
        from netgen.csg import CSGeometry, OrthoBrick, Pnt

        self.rad = rad
        self.radia_ngsolve = radia_ngsolve
        self.Mesh = Mesh
        self.HDiv = HDiv
        self.HCurl = HCurl
        self.VectorH1 = VectorH1
        self.GridFunction = GridFunction
        self.CSGeometry = CSGeometry
        self.OrthoBrick = OrthoBrick
        self.Pnt = Pnt

        rad.UtiDelAll()
        rad.FldUnits('m')

        self.magnet = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.03], [0, 0, 1.2])
        rad.MatApl(self.magnet, rad.MatPM(1.2, 900000, [0, 0, 1]))
        rad.Solve(self.magnet, 0.0001, 10000)

        yield
        rad.UtiDelAll()

    def test_hdiv_space(self):
        """Test HDiv function space (CLAUDE.md recommended)"""
        print("\n[Test] HDiv function space")

        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.02, -0.02),
            self.Pnt(0.06, 0.02, 0.02)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.01))

        fes = self.HDiv(mesh, order=2)
        gf = self.GridFunction(fes)
        B_cf = self.radia_ngsolve.RadiaField(self.magnet, 'b')
        gf.Set(B_cf)

        print(f"  HDiv DOFs: {fes.ndof}")
        print("  [OK] HDiv projection successful")

    def test_hcurl_space(self):
        """Test HCurl function space (for vector potential A)"""
        print("\n[Test] HCurl function space")

        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.02, -0.02),
            self.Pnt(0.06, 0.02, 0.02)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.01))

        fes = self.HCurl(mesh, order=2)
        gf = self.GridFunction(fes)
        A_cf = self.radia_ngsolve.RadiaField(self.magnet, 'a')
        gf.Set(A_cf)

        print(f"  HCurl DOFs: {fes.ndof}")
        print("  [OK] HCurl projection successful")

    def test_vectorh1_space(self):
        """Test VectorH1 function space (continuous vector field)"""
        print("\n[Test] VectorH1 function space")

        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.02, -0.02),
            self.Pnt(0.06, 0.02, 0.02)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.01))

        fes = self.VectorH1(mesh, order=2)
        gf = self.GridFunction(fes)
        B_cf = self.radia_ngsolve.RadiaField(self.magnet, 'b')
        gf.Set(B_cf)

        print(f"  VectorH1 DOFs: {fes.ndof}")
        print("  [OK] VectorH1 projection successful")


# Standalone test function
def run_standalone_test():
    """Run standalone test without pytest"""
    print("=" * 70)
    print("NGSolve Integration Test Suite")
    print("Following CLAUDE.md Best Practices")
    print("=" * 70)

    if not check_ngsolve_available():
        print("\n[SKIP] NGSolve not installed")
        print("Install with: pip install ngsolve")
        return 1

    if not check_radia_ngsolve_available():
        print("\n[SKIP] radia_ngsolve module not built")
        print("Build with: cmake --build build --target radia_ngsolve")
        return 1

    print("\n[OK] Prerequisites satisfied")

    try:
        # Run integration tests
        test = TestNGSolveIntegration()
        test.setup()
        test.test_units_are_meters()
        test.test_radiafield_is_coefficientfunction()
        test.test_all_field_types()
        test.test_hdiv_gridfunction_projection()
        test.test_field_accuracy_far_from_magnet()
        test.test_coefficient_function_direct_evaluation()
        test.test_field_type_attribute()

        # Run function space tests
        test2 = TestNGSolveFunctionSpaces()
        test2.setup()
        test2.test_hdiv_space()
        test2.test_hcurl_space()
        test2.test_vectorh1_space()

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
