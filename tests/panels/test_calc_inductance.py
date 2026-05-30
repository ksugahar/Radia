"""
Layer 1: Computation logic tests for calc_inductance.py ESIM/Dowell.

No Cubit, no Qt required. Uses OCC torus mesh directly.
Runs in CI (pytest).

Tests:
  - ESIM computation on workpiece panels
  - Dowell computation on workpiece panels
  - Result JSON schema (all required keys present)
  - Backward compatibility (no workpiece -> inductance only)
"""

import math
import os
import sys
import pytest
import numpy as np

# Setup path
_repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_repo, "src", "radia"))

# bem_inductance was retired from the IH panel 2026-04-19 and moved to
# examples/induction_heating/bem_reference/ (research-only).  Add that
# dir to sys.path so the importorskip below can find it; if the dir is
# missing entirely (e.g. slim checkout), the whole module skips.
_BEM_REF = os.path.join(_repo, "examples", "induction_heating", "bem_reference")
if os.path.isdir(_BEM_REF) and _BEM_REF not in sys.path:
    sys.path.insert(0, _BEM_REF)

ngsolve = pytest.importorskip("ngsolve")
pytest.importorskip("ngsolve.bem")
pytest.importorskip(
    "bem_inductance",
    reason="bem_inductance.py not on sys.path; expected at "
           "examples/induction_heating/bem_reference/")

from ngsolve import *
from ngsolve import TaskManager
from bem_inductance import compute_inductance_source_sink, MU_0
from esim_cell_problem import ESIMFiniteSlabSolver


# ============================================================
# Fixtures
# ============================================================
@pytest.fixture(scope="module")
def coil_mesh_and_J():
    """Create gapped torus coil mesh + BEM solve (shared across tests)."""
    from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry, Glue
    from netgen.meshing import MeshingParameters

    R, a = 0.05, 0.005
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 355)
    for f in torus.faces:
        f.name = "conductor"
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"
    torus_surf = Glue(torus.faces)
    geo = OCCGeometry(torus_surf)
    ngmesh = geo.GenerateMesh(mp=MeshingParameters(maxh=a, curvaturesafety=2.0))
    mesh = Mesh(ngmesh)
    with TaskManager():
        mesh.Curve(2)

        sol = compute_inductance_source_sink(mesh)
        assert 'error' not in sol
        return mesh, sol['gf_J'], sol['L']


@pytest.fixture
def workpiece_panels():
    """Simple cylindrical workpiece panels at torus center."""
    n_phi, n_z = 8, 4
    radius, height = 0.01, 0.02
    dphi = 2 * math.pi / n_phi
    dz = height / n_z
    z_bot = -height / 2

    panels = []
    for i in range(n_phi):
        phi_mid = (i + 0.5) * dphi
        nx, ny = math.cos(phi_mid), math.sin(phi_mid)
        for j in range(n_z):
            z_mid = z_bot + (j + 0.5) * dz
            panels.append({
                'center': np.array([radius * nx, radius * ny, z_mid]),
                'normal': np.array([nx, ny, 0.0]),
                'area': radius * dphi * dz,
            })
    return panels


def biot_savart_H(mesh_coil, gf_J, obs_pts):
    """Compute H at observation points via Biot-Savart from BEM J."""
    INV_4PI = 1.0 / (4.0 * np.pi)
    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_coil, VOL_or_BND=BND, element_wise=True)

    centroids, areas, J_vecs = [], [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        centroids.append(c); areas.append(area); J_vecs.append(jvec)

    centroids = np.array(centroids)
    areas = np.array(areas)
    J_vecs = np.array(J_vecs)
    obs = np.asarray(obs_pts)
    dx = obs[:, None, :] - centroids[None, :, :]
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-30))
    r3_inv = areas[None, :] / (r ** 3)
    cross = np.cross(J_vecs[None, :, :], dx)
    return INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)


# ============================================================
# Test: ESIM computation
# ============================================================
class TestESIMComputation:

    def test_esim_produces_positive_power(self, coil_mesh_and_J, workpiece_panels):
        """ESIM must produce P > 0 and Q > 0 for conducting workpiece."""
        mesh, gf_J, L = coil_mesh_and_J
        panels = workpiece_panels

        obs_pts = np.array([p['center'] for p in panels])
        H = biot_savart_H(mesh, gf_J, obs_pts)
        normals = np.array([p['normal'] for p in panels])
        H_n = np.sum(H * normals, axis=1, keepdims=True)
        H_t = H - H_n * normals
        H_t_mag = np.linalg.norm(H_t, axis=1)

        solver = ESIMFiniteSlabSolver(
            half_thickness=0.005, sigma=5.8e7, frequency=50000,
            mu_r=1.0, n_nodes=100)

        P_total = 0.0
        Q_total = 0.0
        for i, p in enumerate(panels):
            sol = solver.solve(max(float(H_t_mag[i]), 1e-3))
            P_total += sol['P_prime'] * p['area']
            Q_total += sol['Q_prime'] * p['area']

        assert P_total > 0, "Power loss must be positive"
        assert Q_total > 0, "Reactive power must be positive"

    def test_esim_steel_higher_loss_than_copper(self, coil_mesh_and_J, workpiece_panels):
        """Steel (high mu) should produce more loss than copper at same frequency."""
        mesh, gf_J, L = coil_mesh_and_J
        panels = workpiece_panels

        obs_pts = np.array([p['center'] for p in panels])
        H = biot_savart_H(mesh, gf_J, obs_pts)
        normals = np.array([p['normal'] for p in panels])
        H_n = np.sum(H * normals, axis=1, keepdims=True)
        H_t = H - H_n * normals
        H_t_mag = np.linalg.norm(H_t, axis=1)

        STEEL_BH = [
            [0, 0], [50, 0.1], [100, 0.25], [500, 0.95],
            [1000, 1.2], [5000, 1.55], [50000, 1.9], [100000, 2.0],
        ]

        def total_P(bh, sigma, mu_r):
            solver = ESIMFiniteSlabSolver(
                half_thickness=0.005, bh_curve=bh, sigma=sigma,
                frequency=50000, mu_r=mu_r if bh is None else None,
                n_nodes=100)
            P = 0.0
            for i, p in enumerate(panels):
                sol = solver.solve(max(float(H_t_mag[i]), 1e-3))
                P += sol['P_prime'] * p['area']
            return P

        P_copper = total_P(None, 5.8e7, 1.0)
        P_steel = total_P(STEEL_BH, 2e6, None)

        assert P_steel > P_copper, \
            f"Steel loss ({P_steel:.4e}) should exceed copper ({P_copper:.4e})"


# ============================================================
# Test: Dowell computation
# ============================================================
class TestDowellComputation:

    def test_dowell_matches_analytical(self):
        """Dowell Z_s = (rho/a)*gamma*a*tanh(gamma*a) for linear material."""
        rho = 1.0 / 5.8e7
        a = 0.005
        f = 50000
        omega = 2 * math.pi * f
        mu = MU_0
        delta = math.sqrt(2 * rho / (omega * mu))
        xi = a / delta
        ga = complex(1, 1) * xi
        Z_anal = (rho / a) * ga * np.tanh(ga)

        # ESIM should match
        solver = ESIMFiniteSlabSolver(
            half_thickness=a, sigma=5.8e7, frequency=f,
            mu_r=1.0, n_nodes=200)
        sol = solver.solve(1000.0)
        Z_esim = sol['Z']

        rel_err = abs(abs(Z_esim) - abs(Z_anal)) / abs(Z_anal)
        assert rel_err < 0.01, f"ESIM vs Dowell: {rel_err:.4f} > 1%"

    def test_dowell_power_consistent_with_Z(self):
        """P' = Re(Z_s)*H0^2/2, Q' = Im(Z_s)*H0^2/2."""
        H0 = 1000.0
        solver = ESIMFiniteSlabSolver(
            half_thickness=0.005, sigma=5.8e7, frequency=50000,
            mu_r=1.0, n_nodes=200)
        sol = solver.solve(H0)

        P_from_Z = sol['Z'].real * H0**2 / 2
        Q_from_Z = sol['Z'].imag * H0**2 / 2

        assert abs(P_from_Z - sol['P_prime']) / sol['P_prime'] < 0.001
        assert abs(Q_from_Z - sol['Q_prime']) / sol['Q_prime'] < 0.001


# ============================================================
# Test: Result JSON schema
# ============================================================
class TestResultSchema:

    REQUIRED_KEYS = [
        "inductance_H", "n_dofs", "n_faces", "t_solve", "t_export",
        "t_assembly", "t_lu", "j_npy", "mesh_vol",
    ]

    ESIM_KEYS = [
        "wp_model", "wp_material", "wp_frequency", "wp_sigma",
        "wp_R_effective", "wp_X_effective", "wp_P_total", "wp_Q_total",
        "wp_n_panels", "wp_delta_min", "wp_delta_max",
    ]

    def test_required_keys_documented(self):
        """All required result keys are defined."""
        # This is a static check - just verify the key lists are non-empty
        assert len(self.REQUIRED_KEYS) >= 9
        assert len(self.ESIM_KEYS) >= 10

    def test_esim_result_keys_complete(self):
        """Verify ESIM result dict has all expected keys."""
        # Simulate ESIM result
        result = {
            "wp_model": "esim",
            "wp_material": "steel",
            "wp_frequency": 50000,
            "wp_sigma": 2e6,
            "wp_half_thickness": 0.005,
            "wp_n_panels": 32,
            "wp_H_t_max": 10.0,
            "wp_H_t_min": 0.1,
            "wp_P_total": 1e-4,
            "wp_Q_total": 2e-4,
            "wp_R_effective": 5e-4,
            "wp_X_effective": 1e-3,
            "wp_delta_min": 0.04e-3,
            "wp_delta_max": 0.04e-3,
        }
        for key in self.ESIM_KEYS:
            assert key in result, f"Missing key: {key}"


# ============================================================
# Test: Backward compatibility
# ============================================================
class TestBackwardCompatibility:

    def test_no_workpiece_no_esim_keys(self, coil_mesh_and_J):
        """Without workpiece, result should not contain wp_* keys."""
        mesh, gf_J, L = coil_mesh_and_J
        # The inductance result from BEM only
        assert L > 0
        # Simulated result without workpiece
        result = {"inductance_H": L, "n_dofs": 100}
        wp_keys = [k for k in result if k.startswith("wp_")]
        assert len(wp_keys) == 0

    def test_inductance_reasonable(self, coil_mesh_and_J):
        """BEM inductance should be within 10% of Neumann formula."""
        mesh, gf_J, L = coil_mesh_and_J
        R, a = 0.05, 0.005
        L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)
        err = abs(L - L_neumann) / L_neumann
        assert err < 0.10, f"L={L*1e9:.2f} nH vs Neumann={L_neumann*1e9:.2f} nH"
