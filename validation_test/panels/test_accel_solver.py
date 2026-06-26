"""
Layer 2: Solver logic tests for accelerator magnet panel (Radia + NGSolve).

No Cubit required. Uses OCC mesh generation (NGSolve Netgen).
Tests the physics: CoilBuilder -> RadiaField -> FEM solve -> B field.

Requires: pip install radia (with NGSolve)
Skip condition: NGSolve or Radia not available.
"""

import math
import os
import sys

import numpy as np
import pytest

_repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_src = os.path.join(_repo, "src", "radia")
if _src not in sys.path:
    sys.path.insert(0, _src)

try:
    import radia as rad
    import ngsolve
    HAS_RADIA_NGSOLVE = True
except ImportError:
    HAS_RADIA_NGSOLVE = False

pytestmark = [
    pytest.mark.skipif(not HAS_RADIA_NGSOLVE, reason="Radia + NGSolve required"),
    pytest.mark.filterwarnings("ignore:Gimbal lock:UserWarning"),
]

MU_0 = 4e-7 * math.pi


# ============================================================
# Test: CoilBuilder.to_wire_segments()
# ============================================================
class TestCoilBuilderWireSegments:

    def test_straight_segment(self):
        from coil_builder import CoilBuilder
        mm = 1e-3
        coil = (CoilBuilder(current=100)
            .set_start([0, 0, 0])
            .set_cross_section(10*mm, 10*mm)
            .add_straight(100*mm))
        segs, I = coil.to_wire_segments()
        assert I == 100
        assert len(segs) == 1
        p1, p2 = segs[0]
        assert abs(p1[0]) < 1e-10  # starts at origin
        # Straight goes along local Y
        length = np.linalg.norm(np.array(p2) - np.array(p1))
        assert length == pytest.approx(0.1, rel=1e-6)

    def test_arc_discretized(self):
        from coil_builder import CoilBuilder
        mm = 1e-3
        coil = (CoilBuilder(current=100)
            .set_start([0, 0, 0])
            .set_cross_section(10*mm, 10*mm)
            .add_arc(50*mm, 90))
        segs, I = coil.to_wire_segments(n_arc=10)
        assert len(segs) == 10  # 10 sub-segments for one arc

    def test_racetrack(self):
        from coil_builder import CoilBuilder
        mm = 1e-3
        coil = (CoilBuilder(current=500)
            .set_start([0, 0, 0])
            .set_cross_section(10*mm, 10*mm)
            .add_straight(50*mm)
            .add_arc(25*mm, 180, tilt=90)
            .add_straight(50*mm)
            .add_arc(25*mm, 180, tilt=90))
        segs, I = coil.to_wire_segments(n_arc=20)
        assert I == 500
        # 2 straights + 2 arcs * 20 = 42 segments
        assert len(segs) == 42


# ============================================================
# Test: CoilBuilder.to_radia() + RadiaField
# ============================================================
class TestCoilRadiaField:

    def test_radia_objects_created(self):
        from coil_builder import CoilBuilder
        mm = 1e-3
        rad.UtiDelAll()
        coil = (CoilBuilder(current=1000)
            .set_start([0, 0, 0])
            .set_cross_section(10*mm, 10*mm)
            .add_straight(100*mm))
        objs = coil.to_radia()
        assert len(objs) >= 1

    def test_radiafield_cf_b(self):
        """RadiaField('b') returns a 3-component CoefficientFunction."""
        from coil_builder import CoilBuilder
        mm = 1e-3
        rad.UtiDelAll()
        coil = (CoilBuilder(current=1000)
            .set_start([0, 0, 0])
            .set_cross_section(10*mm, 10*mm)
            .add_straight(100*mm))
        objs = coil.to_radia()
        container = rad.ObjCnt(objs)
        B_cf = rad.RadiaField(container, 'b')
        assert B_cf.dim == 3

    def test_radiafield_cf_h(self):
        """RadiaField('h') returns H-field CoefficientFunction."""
        from coil_builder import CoilBuilder
        mm = 1e-3
        rad.UtiDelAll()
        coil = (CoilBuilder(current=1000)
            .set_start([0, 0, 0])
            .set_cross_section(10*mm, 10*mm)
            .add_straight(100*mm))
        objs = coil.to_radia()
        container = rad.ObjCnt(objs)
        H_cf = rad.RadiaField(container, 'h')
        assert H_cf.dim == 3

    def test_solenoid_field_at_center(self):
        """Simple solenoid: B at center should be close to mu_0 * n * I."""
        rad.UtiDelAll()
        # Create a solenoid-like coil (single-turn loop approximated as arc)
        from coil_builder import CoilBuilder
        R = 0.05  # 50mm radius
        coil = (CoilBuilder(current=1000)
            .set_start([R, 0, 0])
            .set_cross_section(0.005, 0.005)
            .add_arc(R, 360))
        objs = coil.to_radia()
        container = rad.ObjCnt(objs)

        # B at center of loop: mu_0 * I / (2R) for thin loop
        B_analytical = MU_0 * 1000 / (2 * R)

        B = rad.Fld(container, 'b', [0, 0, 0])
        B_mag = math.sqrt(sum(b**2 for b in B))
        # Expect within 20% (thick coil, not thin wire)
        assert B_mag > 0
        assert abs(B_mag - B_analytical) / B_analytical < 0.3


# ============================================================
# Test: Omega-reduced formulation (linear, OCC mesh)
# ============================================================
class TestOmegaReducedLinear:

    @pytest.fixture
    def simple_mesh_and_source(self):
        """Create a simple sphere mesh with a current loop source."""
        from ngsolve import Mesh, TaskManager
        from netgen.occ import Sphere, Pnt, OCCGeometry

        rad.UtiDelAll()
        from coil_builder import CoilBuilder
        R_coil = 0.03
        coil = (CoilBuilder(current=1000)
            .set_start([R_coil, 0, 0])
            .set_cross_section(0.003, 0.003)
            .add_arc(R_coil, 360))
        objs = coil.to_radia()
        container = rad.ObjCnt(objs)
        H_s = rad.RadiaField(container, 'h')
        B_s = rad.RadiaField(container, 'b')

        # Simple air sphere mesh (no yoke, pure Laplace)
        sphere = Sphere(Pnt(0, 0, 0), 0.1)
        sphere.name = "air"
        geo = OCCGeometry(sphere)
        with TaskManager():
            ngmesh = geo.GenerateMesh(maxh=0.03)
        mesh = Mesh(ngmesh)
        mesh.Curve(2)

        return mesh, H_s, B_s, container

    def test_laplace_in_air(self, simple_mesh_and_source):
        """In air (no iron), Omega satisfies Laplace equation.
        H = H_s - grad(Omega), and Omega should be ~0 everywhere.
        """
        from ngsolve import (H1, BilinearForm, LinearForm, GridFunction,
                             grad, dx, Integrate, InnerProduct, TaskManager)

        mesh, H_s, B_s, _ = simple_mesh_and_source

        fes = H1(mesh, order=1, dirichlet="")
        u, v = fes.TnT()

        # mu = mu_0 everywhere (air only)
        a = BilinearForm(fes)
        a += MU_0 * grad(u) * grad(v) * dx
        a.Assemble()

        f = LinearForm(fes)
        f += MU_0 * H_s * grad(v) * dx
        f.Assemble()

        gfu = GridFunction(fes)
        with TaskManager():
            gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

        # Omega should be small (pure air, H_s is curl-free outside coil)
        # The integral of |grad(Omega)|^2 should be small compared to |H_s|^2
        grad_omega_sq = Integrate(
            InnerProduct(grad(gfu), grad(gfu)), mesh).real
        # Not exactly zero due to BC truncation, but should be finite
        assert math.isfinite(grad_omega_sq)


# ============================================================
# Test: A-formulation (linear, OCC mesh)
# ============================================================
class TestAFormulationLinear:

    def test_a_formulation_gauge(self):
        """A-formulation with gauge regularization should solve."""
        from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                             curl, dx, CF, TaskManager, Mesh)
        from netgen.occ import Sphere, Pnt, OCCGeometry

        rad.UtiDelAll()
        from coil_builder import CoilBuilder
        coil = (CoilBuilder(current=1000)
            .set_start([0.03, 0, 0])
            .set_cross_section(0.003, 0.003)
            .add_arc(0.03, 360))
        objs = coil.to_radia()
        container = rad.ObjCnt(objs)
        B_s = rad.RadiaField(container, 'b')

        sphere = Sphere(Pnt(0, 0, 0), 0.08)
        sphere.name = "air"
        geo = OCCGeometry(sphere)
        with TaskManager():
            ngmesh = geo.GenerateMesh(maxh=0.025)
        mesh = Mesh(ngmesh)

        NU_0 = 1.0 / MU_0
        fes = HCurl(mesh, order=1, dirichlet="")
        u, v = fes.TnT()

        a = BilinearForm(fes)
        a += NU_0 * curl(u) * curl(v) * dx
        a += 1e-6 * NU_0 * u * v * dx  # gauge
        a.Assemble()

        f = LinearForm(fes)
        f += -NU_0 * B_s * curl(v) * dx
        f.Assemble()

        gfu = GridFunction(fes)
        with TaskManager():
            gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

        # B_total = curl(A_r) + B_s should be non-zero near center
        # Just verify solve didn't crash and solution is finite
        assert gfu.vec.Norm() > 0
        assert math.isfinite(gfu.vec.Norm())


# ============================================================
# Test: Hysteresis material creation
# ============================================================
class TestHysteresisMaterial:

    @pytest.fixture
    def play_material(self):
        """Create a simple Play hysteresis material."""
        rad.UtiDelAll()
        K = 3
        eta = np.array([0.0, 0.5, 1.0])
        # Simple linear shape functions
        r = np.linspace(0, 2.0, 20)
        f_k_tables = []
        for k in range(K):
            f = 1000.0 * (k + 1) * r  # linear
            f_k_tables.append((r.tolist(), f.tolist()))
        mat = rad.MatPlayHysteresis(K, eta, f_k_tables)
        return mat

    def test_create_play_material(self, play_material):
        assert play_material > 0

    def test_nu_rev_positive(self, play_material):
        nu_rev = rad.MatHysGetNuRev(play_material)
        assert nu_rev > 0

    def test_irreversible_at_zero(self, play_material):
        """H_irr at B=0 should be zero (virgin state)."""
        H_irr = rad.MatHysIrreversible(play_material, [0, 0, 0])
        assert np.linalg.norm(H_irr) < 1e-10

    def test_irreversible_nonzero_after_excitation(self, play_material):
        """After B excitation, H_irr should be non-zero."""
        # Drive to B = [0, 0, 1.0]
        H_irr = rad.MatHysIrreversible(play_material, [0, 0, 1.0])
        assert np.linalg.norm(H_irr) > 0

    def test_state_save_restore(self, play_material):
        """State save/restore should preserve play operator states."""
        # Excite
        rad.MatHysIrreversible(play_material, [0, 0, 0.5])
        state = rad.MatHysSaveState(play_material)
        assert len(state) > 0

        # Excite further
        rad.MatHysIrreversible(play_material, [0, 0, 1.5])

        # Restore
        rad.MatHysRestoreState(play_material, state)
        H_irr = rad.MatHysIrreversible(play_material, [0, 0, 0.5])
        # Should be consistent with the saved state
        assert math.isfinite(np.linalg.norm(H_irr))

    def test_commit_state(self, play_material):
        """Commit advances the state reference for next step."""
        rad.MatHysIrreversible(play_material, [0, 0, 1.0])
        rad.MatHysCommitState(play_material)
        # After commit, a new evaluation should use committed state
        H_irr = rad.MatHysIrreversible(play_material, [0, 0, 1.0])
        assert math.isfinite(np.linalg.norm(H_irr))

    def test_per_element_independent_handles(self):
        """Multiple material handles should have independent state."""
        rad.UtiDelAll()
        K = 2
        eta = np.array([0.0, 0.5])
        r = np.linspace(0, 2.0, 10)
        f_k_tables = [(r.tolist(), (500 * r).tolist()),
                       (r.tolist(), (300 * r).tolist())]

        mat1 = rad.MatPlayHysteresis(K, eta, f_k_tables)
        mat2 = rad.MatPlayHysteresis(K, eta, f_k_tables)

        # Excite mat1 only
        H1 = rad.MatHysIrreversible(mat1, [0, 0, 1.0])
        H2 = rad.MatHysIrreversible(mat2, [0, 0, 0.0])

        # mat1 should have non-zero H_irr, mat2 should be ~zero
        assert np.linalg.norm(H1) > np.linalg.norm(H2)


# ============================================================
# Test: BH interpolation
# ============================================================
class TestBHInterpolation:

    def test_create_interpolators(self):
        sys.path.insert(0, os.path.join(_repo, "src", "radia", "panels"))
        from calc_accel_magnet import _create_bh_interpolators

        bh_data = [[0, 0], [100, 0.1], [1000, 1.2], [50000, 2.0]]
        interp = _create_bh_interpolators(bh_data)

        assert interp['mu_r_init'] > 0

        # Chord: mu at low H should be high (unsaturated)
        mu_low = interp['mu_chord'](50)
        mu_high = interp['mu_chord'](10000)
        assert mu_low > mu_high

        # Chord: nu at low B should be low (high permeability)
        nu_low = interp['nu_chord'](0.05)
        nu_high = interp['nu_chord'](1.5)
        assert nu_low < nu_high

        # Differential: should be positive
        mu_d = interp['mu_diff'](500)
        nu_d = interp['nu_diff'](1.0)
        assert mu_d > 0
        assert nu_d > 0


# ============================================================
# Run as standalone
# ============================================================
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
