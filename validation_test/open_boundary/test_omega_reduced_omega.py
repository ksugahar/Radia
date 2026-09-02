"""test_omega_reduced_omega.py — pytest tests for Omega-ReducedOmega
total+reduced split with Kelvin and coil source.

Checks automatic selection of the total/reduced Kelvin solve and the Radia
scalar-potential sign convention. Quantitative cross-formulation accuracy
belongs to the maintained HDiv-MMM / reduced-A / omega validation lane.
"""
import math

import numpy as np
import pytest

import radia as rad
from ngsolve import Mesh, TaskManager
from netgen.occ import (Cylinder, Sphere, Pnt, Z, Vertex, Glue,
                         OCCGeometry, IdentificationType)

from radia.scalar_potential_solver import ScalarPotentialSolver

pytestmark = pytest.mark.slow

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


class TestOmegaReducedOmegaKelvin:

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
