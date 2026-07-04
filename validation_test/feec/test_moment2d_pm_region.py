"""Golden lock for EMBEDDED permanent-magnet regions (design B) in the 2D planar MMMM.

A PM segment is a REGION of the SAME mesh as the soft iron (a real PM-motor rotor: magnets embedded
in the iron), NOT a separate body (design A, test_moment2d_magnet.py).  The PM is RIGID (fixed M): its
field (shared planar_charges.exterior_field) sources the soft-iron demag, then the iron is solved --
so design B is design A applied to a single partitioned mesh (fields superpose; no new C++).

Gates:
 1. RIGID PM: a uniformly magnetised disk's exterior field == the analytic 2D line dipole
    H_ext,x(r,0) = a^2 M / (2 r^2)  (M along x).
 2. design B == design A: same iron + PM geometry, embedded-region solve == separate-body magnets=.
 3. PM + iron == a MONOLITHIC magnetostatic FEM (equivalent-current A_z formulation) on the iron M.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2
import radia.planar_charges as pc

MU0 = 4e-7 * np.pi
A = 0.01
MREM = 8.0e5     # remanent magnetisation of the PM (A/m)


def _disk(cx, mat, a=A, maxh=0.1):
    geo = SplineGeometry(); geo.AddCircle((cx, 0.0), r=a, bc="e")
    geo.SetMaterial(1, mat)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def _two_region(cx_iron, cx_pm, a=A, maxh=0.1):
    """One mesh: an iron disk (domain 1) + a PM disk (domain 2) -- embedded design-B geometry."""
    geo = SplineGeometry()
    geo.AddCircle((cx_iron, 0.0), r=a, leftdomain=1, rightdomain=0, bc="iron_e")
    geo.AddCircle((cx_pm, 0.0), r=a, leftdomain=2, rightdomain=0, bc="pm_e")
    geo.SetMaterial(1, "iron"); geo.SetMaterial(2, "pm")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_rigid_pm_exterior_is_2d_dipole():
    """A uniformly magnetised disk (M along x) has exterior H_x(r,0) = a^2 M / (2 r^2)."""
    with ng.TaskManager():
        mag = _disk(0.0, "pm", maxh=A / 10)
        M = np.tile([MREM, 0.0], (mag.ne, 1))
        r = np.array([3 * A, 5 * A, 8 * A])
        P = np.stack([r, np.zeros_like(r)], axis=1)
        H = pc.exterior_field(mag, M, P)
    exact = A ** 2 * MREM / (2 * r ** 2)
    rel = np.abs(H[:, 0] - exact) / exact
    assert np.all(rel < 5e-3), (H[:, 0], exact, rel)
    assert np.all(np.abs(H[:, 1]) < 5e-3 * exact)          # on-axis: no transverse field


def test_design_b_equals_design_a():
    """Embedded-region PM (design B) == separate-body PM (design A) for the same geometry."""
    cx_iron, cx_pm = 1.6 * A, -1.6 * A
    with ng.TaskManager():
        # design B: one mesh, iron region + pm region
        mesh = _two_region(cx_iron, cx_pm, maxh=A / 7)
        rB = m2.solve_planar_demag(mesh, mu_r={"iron": 80.0}, H_ext=(0.0, 0.0),
                                   pm={"pm": [MREM, 0.0]})
        # design A: iron mesh alone + separate-body magnet
        iron = _disk(cx_iron, "iron", maxh=A / 7)
        mag = _disk(cx_pm, "pm", maxh=A / 7)
        Mm = np.tile([MREM, 0.0], (mag.ne, 1))
        rA = m2.solve_planar_demag(iron, mu_r=80.0, H_ext=(0.0, 0.0), magnets=[(mag, Mm)])
    # compare the iron magnetisation (design B M includes the PM elements; take the iron region)
    iron_ids = np.array([i for i, m in enumerate(m2._element_materials(mesh)) if m == "iron"], int)
    MB_iron_avg = rB["M"][iron_ids].mean(axis=0)
    MA_iron_avg = rA["M"].mean(axis=0)
    rel = np.linalg.norm(MB_iron_avg - MA_iron_avg) / np.linalg.norm(MA_iron_avg)
    assert rel < 1e-9, (MB_iron_avg, MA_iron_avg, rel)       # same mesh + same kernel -> identical
    assert rB["pm"] is True and rB["demag_factors"] is None


def _monolithic_pm(mesh, mu_r, Mpm, order=4):
    """Reference: magnetostatic A_z with a PM source (equivalent current):
        int nu grad(A).grad(v) = int_pm (Mx dv/dy - My dv/dx),  A=0 on outer.  Returns iron <M>."""
    nu = ng.CoefficientFunction([1.0 / (mu_r * MU0) if m == "iron" else 1.0 / MU0
                                 for m in mesh.GetMaterials()])
    fes = ng.H1(mesh, order=order, dirichlet="outer")
    u, v = fes.TnT()
    af = ng.BilinearForm(fes, symmetric=True)
    af += nu * ng.grad(u) * ng.grad(v) * ng.dx
    af.Assemble()
    lf = ng.LinearForm(fes)
    lf += (Mpm[0] * ng.grad(v)[1] - Mpm[1] * ng.grad(v)[0]) * ng.dx("pm")
    lf.Assemble()
    gfu = ng.GridFunction(fes)
    gfu.vec.data = af.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * lf.vec
    area = ng.Integrate(ng.CF(1), mesh, definedon=mesh.Materials("iron"))
    Bx = ng.Integrate(ng.grad(gfu)[1], mesh, definedon=mesh.Materials("iron")) / area
    By = ng.Integrate(-ng.grad(gfu)[0], mesh, definedon=mesh.Materials("iron")) / area
    return (mu_r - 1.0) / (mu_r * MU0) * np.array([Bx, By])


def _two_region_air(cx_iron, cx_pm, R, a=A, maxh_body=A / 8, maxh_air=None):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), r=R, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((cx_iron, 0.0), r=a, leftdomain=2, rightdomain=1, bc="iron_e")
    geo.AddCircle((cx_pm, 0.0), r=a, leftdomain=3, rightdomain=1, bc="pm_e")
    for i, m in enumerate(("air", "iron", "pm"), start=1):
        geo.SetMaterial(i, m)
    geo.SetDomainMaxH(2, maxh_body); geo.SetDomainMaxH(3, maxh_body)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh_air or R / 12))


def test_pm_iron_matches_monolithic_magnetostatic():
    """Design B PM+iron reproduces a monolithic magnetostatic FEM on the iron magnetisation."""
    mu_r = 60.0
    cx_iron, cx_pm = 1.6 * A, -1.6 * A
    R = 40 * A
    with ng.TaskManager():
        # MMMM design B (iron + pm regions, no air)
        mesh = _two_region(cx_iron, cx_pm, maxh=A / 8)
        r = m2.solve_planar_demag(mesh, mu_r={"iron": mu_r}, H_ext=(0.0, 0.0), pm={"pm": [MREM, 0.0]})
        iron_ids = np.array([i for i, m in enumerate(m2._element_materials(mesh)) if m == "iron"], int)
        M_mmmm = r["M"][iron_ids].mean(axis=0)
        # monolithic FEM (iron + pm + air)
        fmesh = _two_region_air(cx_iron, cx_pm, R)
        M_mono = _monolithic_pm(fmesh, mu_r, [MREM, 0.0])
    rel = abs(M_mmmm[0] - M_mono[0]) / abs(M_mono[0])
    assert rel < 2e-2, (M_mmmm, M_mono, rel)                # cross-method (charge-cloud vs FEM)
    assert M_mmmm[0] > 0 and abs(M_mmmm[1]) < 0.1 * abs(M_mmmm[0])
