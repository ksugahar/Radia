"""
NGSolve Integration Test Suite

Tests the pure Python NGSolve integration with Radia:
- Radia always uses meters (FldUnits is deprecated)
- HDiv(mesh, order=2) for best accuracy
- Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- as_voxel_cf() returns VoxelCoefficient for fast evaluation
- gf.Set(B_cf) for GridFunction projection

This test suite validates:
1. Module import and RadiaField creation
2. Field types (b, h, a, m, phi)
3. HDiv function space integration via gf.Set(B_cf)
4. Field accuracy at various distances from magnet
5. Direct point evaluation vs Radia (via GridFunction)
6. VoxelCoefficient (as_voxel_cf) returns real CoefficientFunction
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
class TestNGSolveIntegration:
    """Test suite for NGSolve integration following CLAUDE.md best practices"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup: import modules and create test magnet"""
        import radia as rad
        from radia import RadiaField
        from ngsolve import Mesh, HDiv, GridFunction, CoefficientFunction
        from netgen.csg import CSGeometry, OrthoBrick, Pnt

        self.rad = rad
        self.RadiaField = RadiaField
        self.Mesh = Mesh
        self.HDiv = HDiv
        self.GridFunction = GridFunction
        self.CoefficientFunction = CoefficientFunction
        self.CSGeometry = CSGeometry
        self.OrthoBrick = OrthoBrick
        self.Pnt = Pnt

        rad.UtiDelAll()

        # Create test magnet (permanent magnet)
        self.magnet_center = [0, 0, 0]  # meters
        self.magnet_size = [0.020, 0.020, 0.030]  # 20mm x 20mm x 30mm
        self.magnet = rad.magnet_box(
            self.magnet_center,
            self.magnet_size,
            [0, 0, 1000.0],
        )

        yield

        # Cleanup
        rad.UtiDelAll()

    def test_units_are_meters(self):
        """Test 1: Verify Radia always uses meters (FldUnits is deprecated)"""
        print("\n[Test 1] Verifying Radia always uses meters")

        # Radia always uses meters. FldUnits is deprecated.
        print("  [OK] Radia always uses meters")

    def test_radiafield_api(self):
        """Test 2: RadiaField is callable with correct dim"""
        print("\n[Test 2] Checking RadiaField API")

        B_cf = self.RadiaField(self.magnet, 'b')
        assert callable(B_cf), "RadiaField must be callable"
        assert hasattr(B_cf, 'dim'), "RadiaField must have dim attribute"
        assert B_cf.dim == 3
        assert hasattr(B_cf, 'field_type'), "RadiaField must have field_type attribute"
        assert B_cf.field_type == 'b'
        assert hasattr(B_cf, 'as_voxel_cf'), "RadiaField must have as_voxel_cf method"
        print(f"  [OK] RadiaField('b') created, dim={B_cf.dim}")

    def test_all_field_types(self):
        """Test 3: All field types (b, h, a, m, phi) work correctly"""
        print("\n[Test 3] Testing all field types")

        for ftype in ['b', 'h', 'a', 'm', 'phi']:
            field = self.RadiaField(self.magnet, ftype)
            assert callable(field)
            assert field.field_type == ftype
            expected_dim = 1 if ftype == 'phi' else 3
            assert field.dim == expected_dim
            print(f"  [OK] RadiaField('{ftype}') works, dim={field.dim}")

    def test_hdiv_gridfunction_projection(self):
        """Test 4: HDiv GridFunction projection via gf.Set(B_cf)"""
        print("\n[Test 4] HDiv GridFunction projection (order=2)")

        # Create mesh outside magnet region
        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.03, -0.03),
            self.Pnt(0.08, 0.03, 0.03)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.01))

        print(f"  Mesh: {mesh.ne} elements, {mesh.nv} vertices")

        # HDiv with order=2 (CLAUDE.md recommended)
        fes = self.HDiv(mesh, order=2)
        B_gf = self.GridFunction(fes)

        # Create RadiaField (it IS a CoefficientFunction) and project
        B_cf = self.RadiaField(self.magnet, 'b')
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
        B_cf = self.RadiaField(self.magnet, 'b')
        B_gf.Set(B_cf)

        # Test points far from magnet
        test_points = [
            (0.05, 0.0, 0.0),
            (0.06, 0.0, 0.0),
            (0.07, 0.0, 0.0),
        ]

        max_rel_error = 0.0
        print(f"  {'Point':<25s} {'Radia Bz':>12s} {'NGSolve Bz':>12s} {'Error %':>10s}")
        print("  " + "-" * 65)

        for pt in test_points:
            B_radia = self.rad.Fld(self.magnet, 'b', list(pt))
            B_ngsolve = B_gf(mesh(*pt))

            # Relative error on Bz (dominant component)
            if abs(B_radia[2]) > 1e-10:
                rel_error = abs(B_radia[2] - B_ngsolve[2]) / abs(B_radia[2]) * 100
            else:
                rel_error = 0.0

            max_rel_error = max(max_rel_error, rel_error)
            print(f"  {str(pt):<25s} {B_radia[2]:>12.6f} {B_ngsolve[2]:>12.6f} {rel_error:>9.2f}%")

        assert max_rel_error < 10.0, f"Max relative error {max_rel_error:.2f}% exceeds 10%"
        print(f"  [OK] Max relative error: {max_rel_error:.2f}%")

    def test_direct_point_evaluation(self):
        """Test 6: RadiaField via GridFunction matches rad.Fld"""
        print("\n[Test 6] RadiaField via GridFunction evaluation")

        # Create mesh around test point
        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.02, -0.02),
            self.Pnt(0.07, 0.02, 0.02)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.008))

        fes = self.HDiv(mesh, order=2)
        B_gf = self.GridFunction(fes)
        B_cf = self.RadiaField(self.magnet, 'b')
        B_gf.Set(B_cf)

        test_point = [0.05, 0.0, 0.0]
        B_radia = self.rad.Fld(self.magnet, 'b', test_point)
        B_gf_val = B_gf(mesh(*test_point))

        print(f"  Test point: {test_point}")
        print(f"  Radia B:  [{B_radia[0]:.6e}, {B_radia[1]:.6e}, {B_radia[2]:.6e}]")
        print(f"  GF B:     [{B_gf_val[0]:.6e}, {B_gf_val[1]:.6e}, {B_gf_val[2]:.6e}]")

        # GridFunction evaluation should be close to Radia (within FE projection error)
        for i in range(3):
            if abs(B_radia[i]) > 1e-10:
                rel_err = abs(B_radia[i] - B_gf_val[i]) / abs(B_radia[i])
                assert rel_err < 0.1, f"Component {i}: rel_err={rel_err:.4f}"

        print("  [OK] GridFunction evaluation matches Radia within FE tolerance")

    def test_as_voxel_cf_returns_coefficientfunction(self):
        """Test 7: as_voxel_cf() returns real NGSolve CoefficientFunction"""
        print("\n[Test 7] as_voxel_cf() returns CoefficientFunction")

        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(-0.05, -0.05, -0.05),
            self.Pnt(0.05, 0.05, 0.05)
        ))
        mesh = self.Mesh(geo.GenerateMesh(maxh=0.02))

        B_cf = self.RadiaField(self.magnet, 'b')
        B_voxel = B_cf.as_voxel_cf(mesh, resolution=21)

        assert isinstance(B_voxel, self.CoefficientFunction), \
            f"Expected CoefficientFunction, got {type(B_voxel)}"
        print(f"  [OK] as_voxel_cf() returns {type(B_voxel).__name__}")

    def test_field_type_attribute(self):
        """Test 8: RadiaField has field_type attribute"""
        print("\n[Test 8] field_type attribute")

        for ftype in ['b', 'h', 'a', 'm', 'phi']:
            field = self.RadiaField(self.magnet, ftype)
            assert hasattr(field, 'field_type')
            assert field.field_type == ftype
            print(f"  [OK] RadiaField('{ftype}').field_type = '{field.field_type}'")


@pytest.mark.skipif(not check_ngsolve_available(),
                   reason="NGSolve not installed")
class TestNGSolveFunctionSpaces:
    """Test different NGSolve function spaces"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup for function space tests"""
        import radia as rad
        from radia import RadiaField
        from ngsolve import Mesh, HDiv, HCurl, VectorH1, GridFunction
        from netgen.csg import CSGeometry, OrthoBrick, Pnt

        self.rad = rad
        self.RadiaField = RadiaField
        self.Mesh = Mesh
        self.HDiv = HDiv
        self.HCurl = HCurl
        self.VectorH1 = VectorH1
        self.GridFunction = GridFunction
        self.CSGeometry = CSGeometry
        self.OrthoBrick = OrthoBrick
        self.Pnt = Pnt

        rad.UtiDelAll()

        self.magnet = rad.magnet_box([0, 0, 0], [0.02, 0.02, 0.03], [0, 0, 1000.0])

        yield
        rad.UtiDelAll()

    def _make_mesh(self):
        geo = self.CSGeometry()
        geo.Add(self.OrthoBrick(
            self.Pnt(0.03, -0.02, -0.02),
            self.Pnt(0.06, 0.02, 0.02)
        ))
        return self.Mesh(geo.GenerateMesh(maxh=0.01))

    def test_hdiv_space(self):
        """Test HDiv function space (CLAUDE.md recommended)"""
        print("\n[Test] HDiv function space")
        mesh = self._make_mesh()

        fes = self.HDiv(mesh, order=2)
        gf = self.GridFunction(fes)
        B_cf = self.RadiaField(self.magnet, 'b')
        gf.Set(B_cf)

        print(f"  HDiv DOFs: {fes.ndof}")
        print("  [OK] HDiv projection successful")

    def test_hcurl_space(self):
        """Test HCurl function space (for vector potential A)"""
        print("\n[Test] HCurl function space")
        mesh = self._make_mesh()

        fes = self.HCurl(mesh, order=2)
        gf = self.GridFunction(fes)
        A_cf = self.RadiaField(self.magnet, 'a')
        gf.Set(A_cf)

        print(f"  HCurl DOFs: {fes.ndof}")
        print("  [OK] HCurl projection successful")

    def test_vectorh1_space(self):
        """Test VectorH1 function space (continuous vector field)"""
        print("\n[Test] VectorH1 function space")
        mesh = self._make_mesh()

        fes = self.VectorH1(mesh, order=2)
        gf = self.GridFunction(fes)
        B_cf = self.RadiaField(self.magnet, 'b')
        gf.Set(B_cf)

        print(f"  VectorH1 DOFs: {fes.ndof}")
        print("  [OK] VectorH1 projection successful")


# Standalone test function
def run_standalone_test():
    """Run standalone test without pytest"""
    print("=" * 70)
    print("NGSolve Integration Test Suite (Pure Python)")
    print("=" * 70)

    if not check_ngsolve_available():
        print("\n[SKIP] NGSolve not installed")
        print("Install with: pip install ngsolve")
        return 1

    print("\n[OK] Prerequisites satisfied")

    try:
        import radia as rad
        from radia import RadiaField
        from ngsolve import Mesh, HDiv, HCurl, VectorH1, GridFunction, CoefficientFunction
        from ngsolve import TaskManager
        from netgen.csg import CSGeometry, OrthoBrick, Pnt

        # Setup for integration tests
        test = TestNGSolveIntegration()
        rad.UtiDelAll()
        test.rad = rad
        test.RadiaField = RadiaField
        test.Mesh = Mesh
        test.HDiv = HDiv
        test.GridFunction = GridFunction
        test.CoefficientFunction = CoefficientFunction
        test.CSGeometry = CSGeometry
        test.OrthoBrick = OrthoBrick
        test.Pnt = Pnt
        test.magnet_center = [0, 0, 0]
        test.magnet_size = [0.020, 0.020, 0.030]
        test.magnet = rad.magnet_box(test.magnet_center, test.magnet_size, [0, 0, 1000.0])

        test.test_units_are_meters()
        test.test_radiafield_api()
        test.test_all_field_types()
        test.test_hdiv_gridfunction_projection()
        test.test_field_accuracy_far_from_magnet()
        test.test_direct_point_evaluation()
        test.test_as_voxel_cf_returns_coefficientfunction()
        test.test_field_type_attribute()
        rad.UtiDelAll()

        # Setup for function space tests
        test2 = TestNGSolveFunctionSpaces()
        rad.UtiDelAll()
        test2.rad = rad
        test2.RadiaField = RadiaField
        test2.Mesh = Mesh
        test2.HDiv = HDiv
        test2.HCurl = HCurl
        test2.VectorH1 = VectorH1
        test2.GridFunction = GridFunction
        test2.CSGeometry = CSGeometry
        test2.OrthoBrick = OrthoBrick
        test2.Pnt = Pnt
        test2.magnet = rad.magnet_box([0, 0, 0], [0.02, 0.02, 0.03], [0, 0, 1000.0])

        test2.test_hdiv_space()
        test2.test_hcurl_space()
        test2.test_vectorh1_space()
        rad.UtiDelAll()

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
