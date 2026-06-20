"""test_omega_reduced_omega.py — pytest tests for Omega-ReducedOmega
total+reduced split with Kelvin and coil source.

Tests ScalarPotentialSolver.solve_total_reduced_potential() against
Radia BEM reference for a mu_r=100 cylinder in a closed-loop coil field.
"""
import math
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src', 'radia'))

import radia as rad
from ngsolve import Mesh, TaskManager
from netgen.occ import (Cylinder, Sphere, Pnt, Z, Vertex, Glue,
                         OCCGeometry, IdentificationType)

from scalar_potential_solver import ScalarPotentialSolver
from netgen_mesh_import import netgen_mesh_to_radia

MU_0 = 4 * math.pi * 1e-7

# --- Geometry parameters ---
CYL_R = 0.03
CYL_H = 0.06
MU_R = 100
R_COIL = 0.08
A_COIL = 0.008
Z_COIL = 0.10
I_TOTAL = 1000.0
PHYS_R = 0.18
KELVIN_CENTER = (0.50, 0.0, 0.0)
MAXH_IRON = 0.008
MAXH_AIR = 0.025
MAXH_KELVIN = 0.03
FE_ORDER = 2


def _build_coil(gap_deg=0.0):
    """Build ObjArcCur coil (no UtiDelAll)."""
    arc = math.radians(360.0 - gap_deg)
    J0 = I_TOTAL / (2 * A_COIL) ** 2
    return rad.ObjArcCur(
        [0, 0, Z_COIL],
        [R_COIL - A_COIL, R_COIL + A_COIL],
        [0, arc], 2 * A_COIL, 200, 'man', 'z', J0)


def _build_kelvin_mesh():
    """Build iron + air_inner + air_outer (Kelvin) mesh."""
    half_h = CYL_H / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z, r=CYL_R, h=CYL_H)
    mag_cyl.mat("iron"); mag_cyl.maxh = MAXH_IRON
    for f in mag_cyl.faces:
        f.name = "iron_surf"
    inner = Sphere(Pnt(0, 0, 0), PHYS_R); inner.maxh = MAXH_AIR
    for f in inner.faces:
        f.name = "kelvin_int"
    inner_air = inner - mag_cyl; inner_air.mat("air_inner")
    outer = Sphere(Pnt(*KELVIN_CENTER), PHYS_R)
    outer.maxh = MAXH_KELVIN; outer.mat("air_outer")
    for f in outer.faces:
        f.name = "kelvin_ext"
    gnd = Vertex(Pnt(*KELVIN_CENTER)); gnd.name = "GND"
    geo = Glue([inner_air, mag_cyl, outer, gnd])
    geo.solids[0].name = "air_inner"
    geo.solids[1].name = "iron"
    geo.solids[2].name = "air_outer"
    k_int = k_ext = None
    for s in geo.solids:
        for f in s.faces:
            if f.name == "kelvin_int" and k_int is None:
                k_int = f
            elif f.name == "kelvin_ext" and k_ext is None:
                k_ext = f
    k_int.Identify(k_ext, "periodic", IdentificationType.PERIODIC)
    with TaskManager():
        ngmesh = OCCGeometry(geo).GenerateMesh(maxh=MAXH_AIR, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(FE_ORDER)
    return mesh


def _build_radia_reference(gap_deg=0.0):
    """Build Radia BEM reference (coil + iron cylinder, solved)."""
    rad.UtiDelAll()
    coil = _build_coil(gap_deg)
    half_h = CYL_H / 2
    mag_cyl = Cylinder(Pnt(0, 0, -half_h), Z, r=CYL_R, h=CYL_H)
    mag_cyl.mat("magnetic"); mag_cyl.maxh = 0.004
    geo = OCCGeometry(mag_cyl)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=0.004)
    mesh_rad = Mesh(ngmesh)
    cyl_obj = netgen_mesh_to_radia(mesh_rad,
                                    material={'magnetization': [0, 0, 0]},
                                    units='m', material_filter='magnetic',
                                    verbose=False)
    rad.MatApl(cyl_obj, rad.MatLin(MU_R))
    grp = rad.ObjCnt([coil, cyl_obj])
    rad.SolverConfig(bicgstab_tol=1e-8, relax_param=0.0)
    rad.Solve(grp, 1e-7, 2000, 1)
    return grp


def _sample_points():
    half_h = CYL_H / 2
    pts = ([(0, 0, zv) for zv in np.linspace(-half_h + 0.01, half_h - 0.01, 5)]
         + [(xv, 0, 0) for xv in np.linspace(0.005, CYL_R - 0.005, 4)]
         + [(0, 0, zv) for zv in np.linspace(half_h + 0.02, PHYS_R - 0.04, 4)])
    return pts


def _compare_B(solver, grp_ref, pts, mesh):
    B_cf = solver.get_B()
    rels = []
    for pt in pts:
        try:
            mip = mesh(*pt)
            B_fem = np.array([B_cf(mip)[i] for i in range(3)])
        except Exception:
            continue
        B_rad = np.array(rad.Fld(grp_ref, 'b', list(pt)))
        m_fem = float(np.linalg.norm(B_fem))
        m_rad = float(np.linalg.norm(B_rad))
        rels.append((m_fem - m_rad) / max(m_rad, 1e-12))
    return np.array(rels)


# 3 of the 4 B-accuracy tests below have been observed to fail intermittently
# in CI (tag-CI runs of v4.44.0; main-CI on the same commit passes; LAB-local
# always passes).  Cause not pinned but consistent with MKL/NGSolve threading
# non-determinism on the self-hosted runner.  pytest-rerunfailures retries up
# to 2 times to absorb the rare bad iteration without masking a real
# regression.  See cf576fa4 / eeb664c3 commit messages for the L_coil
# convergence study where the same runner-state effect showed up.
@pytest.mark.slow
@pytest.mark.flaky(reruns=2, reruns_delay=2)
class TestOmegaReducedOmegaKelvin:

    @pytest.fixture(scope="class")
    def radia_ref(self):
        return _build_radia_reference(gap_deg=0.0)

    @pytest.fixture(scope="class")
    def fem_solver(self):
        mesh = _build_kelvin_mesh()
        coil = _build_coil(gap_deg=0.0)
        solver = ScalarPotentialSolver(
            mesh, iron_domains='iron', mu_r=MU_R, order=FE_ORDER,
            kelvin_region='air_outer', kelvin_radius=PHYS_R,
            kelvin_center=list(KELVIN_CENTER))
        solver.set_source_from_radia(coil)
        solver.solve()
        return solver, mesh

    def test_auto_selects_total_reduced(self, fem_solver):
        solver, _ = fem_solver
        assert solver._phi_gf is not None

    def test_B_accuracy_inside_iron(self, fem_solver, radia_ref):
        solver, mesh = fem_solver
        half_h = CYL_H / 2
        pts = [(0, 0, zv) for zv in np.linspace(-half_h + 0.01,
                                                   half_h - 0.01, 5)]
        rels = _compare_B(solver, radia_ref, pts, mesh)
        assert np.max(np.abs(rels)) < 0.05, \
            f"Inside iron max|rel| = {np.max(np.abs(rels))*100:.2f}%"

    def test_B_accuracy_outside(self, fem_solver, radia_ref):
        solver, mesh = fem_solver
        half_h = CYL_H / 2
        pts = [(0, 0, zv) for zv in np.linspace(half_h + 0.02,
                                                   PHYS_R - 0.04, 4)]
        rels = _compare_B(solver, radia_ref, pts, mesh)
        assert np.max(np.abs(rels)) < 0.01, \
            f"Outside max|rel| = {np.max(np.abs(rels))*100:.2f}%"

    def test_overall_rms(self, fem_solver, radia_ref):
        solver, mesh = fem_solver
        rels = _compare_B(solver, radia_ref, _sample_points(), mesh)
        rms = np.sqrt(np.mean(rels**2))
        assert rms < 0.03, f"RMS = {rms*100:.2f}%"

    def test_phi_convention(self):
        """Verify H = -grad(Phi) for full-circle coil."""
        rad.UtiDelAll()
        coil = _build_coil(gap_deg=0.0)
        pt = np.array([0.0, 0.0, 0.0])
        eps = 1e-5
        H = np.array(rad.Fld(coil, 'h', pt.tolist()))
        grad_phi = np.zeros(3)
        for i in range(3):
            pp = pt.copy(); pp[i] += eps
            pm = pt.copy(); pm[i] -= eps
            grad_phi[i] = (rad.Fld(coil, 'phi', pp.tolist())
                          - rad.Fld(coil, 'phi', pm.tolist())) / (2*eps)
        # H = -grad(Phi) → H + grad(Phi) ≈ 0
        residual = np.linalg.norm(H + grad_phi) / np.linalg.norm(H)
        assert residual < 0.01, f"H + grad(Phi) residual = {residual:.4f}"
