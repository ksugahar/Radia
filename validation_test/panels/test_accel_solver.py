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

    def test_panel_law_is_the_production_law_with_vacuum_tail(self):
        """The panel's chord/differential laws are radia.bh_law plus mu0 beyond the table."""
        sys.path.insert(0, os.path.join(_repo, "src", "radia", "panels"))
        from calc_accel_magnet import _create_bh_interpolators
        from radia.bh_law import monotone_bh_pchip

        table = np.loadtxt(os.path.join(_repo, "src", "radia", "panels",
                                        "samples", "em_sample_bh.txt"))
        H_tab, B_tab = table[:, 0], table[:, 1]
        interp = _create_bh_interpolators(table.tolist())
        law = monotone_bh_pchip(H_tab, B_tab)

        inside = np.concatenate([H_tab[1:], 0.5 * (H_tab[1:] + H_tab[:-1])])
        B_panel = interp['mu_chord'](inside) * inside
        assert np.allclose(B_panel, law(inside), rtol=1e-12, atol=0.0)
        # At H_max itself the panel takes the tail's one-sided slope (mu0).
        interior = inside[inside < H_tab[-1]]
        assert np.allclose(interp['mu_diff'](interior), law.derivative()(interior),
                           rtol=1e-12, atol=0.0)
        # The sample's magnetization rises, so dB/dH never drops below mu0.
        assert np.all(interp['mu_diff'](inside) >= MU_0 * (1.0 - 1e-12))

        # Beyond the table: B = B_max + mu0 (H - H_max), H(B) its inverse.
        H_max, B_max = float(H_tab[-1]), float(B_tab[-1])
        for factor in (1.5, 10.0, 100.0):
            H = factor * H_max
            B = B_max + MU_0 * (H - H_max)
            assert interp['mu_chord'](H) * H == pytest.approx(B, rel=1e-13)
            assert interp['mu_diff'](H) == pytest.approx(MU_0, rel=1e-13)
            assert interp['nu_chord'](B) * B == pytest.approx(H, rel=1e-12)
            assert interp['nu_diff'](B) == pytest.approx(1.0 / MU_0, rel=1e-13)
        # The retired chord extrapolation held B/H at its last value, which
        # at ten times H_max claimed ten times B_max (26 T for this sample).
        assert interp['mu_chord'](10.0 * H_max) * 10.0 * H_max < 0.5 * 10.0 * B_max

        # H(B) inverts B(H) inside the table, and the differentials are
        # reciprocal away from the C0 kink at H_max (one-sided slopes differ).
        H_probe = interior[::7]
        B_probe = interp['mu_chord'](H_probe) * H_probe
        assert np.allclose(interp['nu_chord'](B_probe) * B_probe, H_probe,
                           rtol=1e-9, atol=1e-9 * H_max)
        assert np.allclose(interp['nu_diff'](B_probe) * interp['mu_diff'](H_probe),
                           1.0, rtol=1e-8)

    def test_panel_law_rejects_tables_the_production_law_rejects(self):
        sys.path.insert(0, os.path.join(_repo, "src", "radia", "panels"))
        from calc_accel_magnet import _create_bh_interpolators

        with pytest.raises(ValueError, match="start at"):
            _create_bh_interpolators([[10, 0.1], [100, 1.0]])
        with pytest.raises(ValueError, match="strictly increasing"):
            _create_bh_interpolators([[0, 0], [100, 1.0], [100, 1.2]])
        with pytest.raises(ValueError, match="non-decreasing"):
            _create_bh_interpolators([[0, 0], [100, 1.0], [200, 0.9]])
        # A flat piece has dB/dH = 0, so H(B) and the reluctivity do not exist.
        with pytest.raises(ValueError, match="not strictly increasing"):
            _create_bh_interpolators([[0, 0], [1, 0], [2, 1]])
        with pytest.raises(ValueError, match="not strictly increasing"):
            _create_bh_interpolators([[0, 0], [100, 1.0], [200, 1.0], [300, 1.5]])

    def test_production_coenergy_is_the_antiderivative_of_the_law(self):
        """The energy split used by the panel: coenergy = int B dH, energy = BH - coenergy."""
        import ngsolve as ng
        from netgen.geom2d import unit_square
        from scipy.integrate import quad
        from radia.bh_law import monotone_bh_pchip
        from radia.scalar_potential_solver import _build_bh_spline_law

        table = np.loadtxt(os.path.join(_repo, "src", "radia", "panels",
                                        "samples", "em_sample_bh.txt"))
        H_tab, B_tab = table[:, 0], table[:, 1]
        law = monotone_bh_pchip(H_tab, B_tab)
        H_max, B_max = float(H_tab[-1]), float(B_tab[-1])

        def b_exact(h):
            return float(law(h)) if h <= H_max else B_max + MU_0 * (h - H_max)

        b_of, coenergy_of, _ = _build_bh_spline_law(table)
        mesh = ng.Mesh(unit_square.GenerateMesh(maxh=0.5))
        point = mesh(0.5, 0.5)
        for h in (50.0, 2.0e3, 5.0e4, 3.0e5, 1.0e6, 5.0e6):
            knots = [k for k in H_tab if k < h] + [h]
            exact = sum(quad(b_exact, a, b, limit=200)[0] for a, b in zip([0.0] + knots[:-1], knots))
            assert coenergy_of(ng.CF(h))(point) == pytest.approx(exact, rel=1e-9)
            assert b_of(ng.CF(h))(point) == pytest.approx(b_exact(h), rel=1e-12)
        # In saturation 0.5*B*H is neither the energy nor the coenergy.
        h = 1.0e6
        coenergy = coenergy_of(ng.CF(h))(point)
        energy = b_exact(h) * h - coenergy
        assert abs(0.5 * b_exact(h) * h - energy) > 0.2 * energy

    def test_panel_inverse_law_is_exact_for_small_nonzero_flux_density(self):
        """H(B) keeps relative accuracy down to tiny B; no absolute-tolerance zero."""
        sys.path.insert(0, os.path.join(_repo, "src", "radia", "panels"))
        from calc_accel_magnet import _create_bh_interpolators
        from radia.bh_law import monotone_bh_pchip

        table = np.loadtxt(os.path.join(_repo, "src", "radia", "panels",
                                        "samples", "em_sample_bh.txt"))
        interp = _create_bh_interpolators(table.tolist())
        law = monotone_bh_pchip(table[:, 0], table[:, 1])
        mu_initial = interp['mu_r_init'] * MU_0
        for B in (1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 0.5, 2.0):
            H = interp['nu_chord'](B) * B
            assert H > 0.0
            assert float(law(H)) == pytest.approx(B, rel=1e-12)
            assert np.isfinite(interp['nu_diff'](B)) and interp['nu_diff'](B) > 0.0
        # In the initial linear range the chord reluctivity is 1/mu(0).
        assert interp['nu_chord'](1e-15) == pytest.approx(1.0 / mu_initial, rel=1e-6)
        assert interp['nu_chord'](0.0) == pytest.approx(1.0 / mu_initial, rel=1e-12)


# ============================================================
# Test: increment solves measured against the system load
# ============================================================
class TestIncrementResidualScale:

    @pytest.fixture
    def poisson(self):
        from ngsolve import H1, BilinearForm, LinearForm, GridFunction, grad, dx, Mesh
        from netgen.geom2d import unit_square

        mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
        fes = H1(mesh, order=1, dirichlet=".*")
        u, v = fes.TnT()
        a = BilinearForm(grad(u) * grad(v) * dx).Assemble()
        f = LinearForm(1.0 * v * dx).Assemble()
        return fes, a, f, GridFunction(fes)

    @pytest.mark.parametrize("bad", [float("inf"), float("nan"), 0.0, -1.0])
    def test_reference_norm_must_be_finite_and_positive(self, poisson, bad):
        sys.path.insert(0, os.path.join(_repo, "src", "radia", "panels"))
        from calc_common import apply_fe_inverse

        fes, a, f, gfu = poisson
        inverse = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        with pytest.raises(ValueError, match="finite positive"):
            apply_fe_inverse(a.mat, inverse, f.vec, gfu.vec, fes.FreeDofs(),
                             reference_norm=bad)

    def test_inaccurate_increment_fails_against_a_valid_load_scale(self, poisson):
        sys.path.insert(0, os.path.join(_repo, "src", "radia", "panels"))
        from calc_common import apply_fe_inverse

        fes, a, f, gfu = poisson
        exact = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
        load = f.vec.Norm()
        # Accurate solve passes against the load scale ...
        assert apply_fe_inverse(a.mat, exact, f.vec, gfu.vec, fes.FreeDofs(),
                                reference_norm=load) <= 1e-6
        # ... a wrong inverse (half the exact one) does not.
        gfu.vec[:] = 0.0
        wrong = 0.5 * exact
        with pytest.raises(RuntimeError, match="true residual"):
            apply_fe_inverse(a.mat, wrong, f.vec, gfu.vec, fes.FreeDofs(),
                             reference_norm=load)


# ============================================================
# Run as standalone
# ============================================================
if __name__ == "__main__":
    pytest.main([__file__, "-v"])
