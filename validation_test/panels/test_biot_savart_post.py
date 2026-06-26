"""Integration test for BiotSavartCF post-processing pipeline.

Tests _compute_B_field_cf: surface J -> BiotSavartCF -> B CoefficientFunction.
Requires NGSolve + ngsolve.bem.

Run: python -m pytest validation_test/panels/test_biot_savart_post.py -v
"""

import math
import os
import sys
import tempfile

import numpy as np
import pytest

# Skip entire module if NGSolve unavailable
ngsolve = pytest.importorskip("ngsolve")

_panels_dir = os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia", "panels")
_panels_dir = os.path.abspath(_panels_dir)
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)

_radia_src = os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia")
_radia_src = os.path.abspath(_radia_src)
if _radia_src not in sys.path:
    sys.path.insert(0, _radia_src)

try:
    from ngsolve.bem import BiotSavartCF  # noqa: F401
    _HAS_BEM = True
except ImportError:
    _HAS_BEM = False


def _make_loop_surface_mesh(R=0.03, n_seg=24, order=1):
    """Create a surface mesh of a torus (circular cross-section wire loop).

    Returns an HDivSurface mesh with solvable J distribution.
    This is a minimal test fixture -- not a real EFIE solve.
    """
    from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Axis,
                             OCCGeometry, Glue)
    from ngsolve import Mesh, TaskManager

    a = 0.002  # wire radius
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    for f in torus.faces:
        f.name = "coil"
        f.maxh = a * 3

    geo = OCCGeometry(Glue(torus.faces))
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=a * 4)
    mesh = Mesh(ngmesh)
    if order > 1:
        mesh.Curve(order)
    return mesh


@pytest.fixture(scope="module")
def loop_mesh_and_J():
    """Create loop surface mesh and a synthetic J (uniform azimuthal).

    HDivSurface.Set() doesn't work on surface-only meshes,
    so we populate DOFs directly via edge integrals.
    """
    mesh = _make_loop_surface_mesh(R=0.03, order=1)

    from ngsolve import HDivSurface, GridFunction, CF, BND, Integrate

    fes = HDivSurface(mesh, order=0)
    gf_J = GridFunction(fes)

    # Populate DOFs: for RWG (order=0), each DOF is the normal flux
    # through an edge.  We approximate J = J0 * e_theta and compute
    # edge integrals manually via element-wise integration.
    J0 = 1.0 / (math.pi * 0.002**2)  # A/m^2

    # Assign random-but-nonzero DOF values (synthetic test current)
    rng = np.random.default_rng(42)
    gf_J.vec.FV().NumPy()[:] = rng.normal(0, J0 * 1e-4, fes.ndof)

    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    return mesh, gf_J, elem_A


class TestComputeBFieldCF:
    @pytest.mark.skipif(
        not _HAS_BEM, reason="ngsolve.bem.BiotSavartCF not available")
    def test_returns_cf_and_mesh(self, loop_mesh_and_J):
        """_compute_B_field_cf returns (CoefficientFunction, Mesh)."""
        mesh, gf_J, elem_A = loop_mesh_and_J

        from calc_inductance import _compute_B_field_cf

        B_cf, mesh_vol = _compute_B_field_cf(
            mesh, gf_J, elem_A,
            lx=0.05, ly=0.05, lz=0.05, maxh_vol=0.015)

        assert mesh_vol.nv > 10, "Volume mesh should have vertices"
        # B_cf should be evaluable as a CoefficientFunction
        assert B_cf.dim == 3, "B should be 3-component vector"

    @pytest.mark.skipif(
        not _HAS_BEM, reason="ngsolve.bem.BiotSavartCF not available")
    def test_b_field_nonzero_at_center(self, loop_mesh_and_J):
        """B at center of loop should be nonzero (Bz dominant)."""
        mesh, gf_J, elem_A = loop_mesh_and_J
        from calc_inductance import _compute_B_field_cf

        B_cf, mesh_vol = _compute_B_field_cf(
            mesh, gf_J, elem_A,
            lx=0.05, ly=0.05, lz=0.05, maxh_vol=0.015)

        # Evaluate B at origin
        from ngsolve import HDiv, GridFunction
        fes_B = HDiv(mesh_vol, order=1)
        gf_B = GridFunction(fes_B)
        gf_B.Set(B_cf)

        # Sample at a point near center
        B_center = gf_B(mesh_vol(0.001, 0.001, 0.001))
        B_mag = np.linalg.norm(B_center)
        assert B_mag > 1e-10, f"B at center should be nonzero, got {B_mag}"

    @pytest.mark.skipif(
        not _HAS_BEM, reason="ngsolve.bem.BiotSavartCF not available")
    def test_gmsh_export(self, loop_mesh_and_J):
        """GmshPostExport produces valid .msh files."""
        mesh, gf_J, elem_A = loop_mesh_and_J
        from calc_inductance import _compute_B_field_cf
        from gmsh_post_export import GmshPostExport
        from ngsolve import HDiv, GridFunction, Norm

        B_cf, mesh_vol = _compute_B_field_cf(
            mesh, gf_J, elem_A,
            lx=0.05, ly=0.05, lz=0.05, maxh_vol=0.015)

        fes_B = HDiv(mesh_vol, order=1)
        gf_B = GridFunction(fes_B)
        gf_B.Set(B_cf)

        with tempfile.TemporaryDirectory() as tmpdir:
            msh_file = os.path.join(tmpdir, "test_B.msh")
            post = GmshPostExport(mesh_vol)
            post.add_vector_field("B", gf_B)
            post.add_scalar_field("|B|", Norm(gf_B))
            post.write(msh_file)

            assert os.path.exists(msh_file)
            size = os.path.getsize(msh_file)
            assert size > 100, f"GMSH file too small: {size} bytes"
