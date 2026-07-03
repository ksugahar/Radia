"""Golden lock for the 2D planar collocation MMMM (radia.mmmm2d, C++ rad_moment2d).

The 2D twin of the 3D moment method: per-unit-length motor-cross-section soft-iron demag on a mesh of
triangles / quadrilaterals with the 2D Laplace kernel -ln(r)/(2 pi).  Each element carries one uniform
line-charge DOF per EDGE; M = chi H is imposed on the field MOMENTS about the centroid (1 monopole +
2 dipole + (nEdge-3) quadrupole).  Triangle = no quad row (2D simplex); quad = 1 quad row.

Locked physics (analytic 2D demag; the C++ core reproduces a numpy PoC to machine precision):
  * disk:         Dx = Dy = 1/2
  * ellipse a:b:  Dx = b/(a+b), Dy = a/(a+b),  Dx + Dy = 1   (the anisotropy discriminator)
  * linear resp:  <M> = chi H0 / (1 + chi D)  across a chi range
  * quad element == triangulated same geometry (shared-geometry cross-check)
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

import radia.mmmm2d as m2

MU0 = 4e-7 * np.pi


def _disk(maxh, quad=False):
    geo = SplineGeometry(); geo.AddCircle((0.0, 0.0), r=1.0, bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh, quad_dominated=quad))


def _ellipse(a, b, maxh, quad=False):
    geo = SplineGeometry()
    n = 128
    pts = [(a * np.cos(t), b * np.sin(t)) for t in np.linspace(0, 2 * np.pi, n, endpoint=False)]
    pid = [geo.AppendPoint(*p) for p in pts]
    for i in range(n):
        geo.Append(["line", pid[i], pid[(i + 1) % n]], bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh, quad_dominated=quad))


def _rect(hx, hy, maxh, quad):
    geo = SplineGeometry(); geo.AddRectangle((-hx, -hy), (hx, hy), bc="outer")
    return ng.Mesh(geo.GenerateMesh(maxh=maxh, quad_dominated=quad))


def test_disk_demag_half():
    """Circle: Dx = Dy = 1/2 (both axes), tiny cross-coupling."""
    with ng.TaskManager():
        Dx, Dy = m2.demag_factors(_disk(0.08), chi=3.0)
    assert abs(Dx - 0.5) < 2e-3, Dx
    assert abs(Dy - 0.5) < 2e-3, Dy


def test_ellipse_2to1_demag():
    """Ellipse a=2 (x), b=1 (y): Dx = 1/3, Dy = 2/3, sum = 1 -- the anisotropy discriminator."""
    with ng.TaskManager():
        Dx, Dy = m2.demag_factors(_ellipse(2.0, 1.0, 0.07), chi=3.0)
    assert abs(Dx - 1.0 / 3.0) < 3e-3, Dx
    assert abs(Dy - 2.0 / 3.0) < 3e-3, Dy
    assert abs((Dx + Dy) - 1.0) < 1e-3, Dx + Dy


def test_linear_chi_response_disk():
    """<M> = chi H0/(1 + chi/2) across a chi range (disk, unit H along x)."""
    with ng.TaskManager():
        mesh = _disk(0.08)
        for chi in (0.5, 3.0, 30.0, 300.0):
            r = m2.solve_planar_demag(mesh, mu_r=chi + 1.0, H_ext=(1.0, 0.0))
            Mx = r["M_avg"][0]
            ref = chi / (1.0 + chi * 0.5)
            assert abs(Mx / ref - 1.0) < 5e-3, (chi, Mx, ref)


def test_quad_element_matches_triangulated():
    """A quad-meshed square and a triangulated square give the SAME effective demag (the quad row
    is correct; a square's D_eff is NOT 1/2 -- non-ellipsoidal bodies have non-uniform interior M)."""
    with ng.TaskManager():
        rq = m2.solve_planar_demag(_rect(1.0, 1.0, 0.15, quad=True), mu_r=4.0, H_ext=(1.0, 0.0))
        rt = m2.solve_planar_demag(_rect(1.0, 1.0, 0.15, quad=False), mu_r=4.0, H_ext=(1.0, 0.0))
    Dq = 1.0 / rq["M_avg"][0] - 1.0 / 3.0
    Dt = 1.0 / rt["M_avg"][0] - 1.0 / 3.0
    assert rq["ndof"] > 0 and rt["ndof"] > 0
    assert abs(Dq - Dt) < 5e-3, (Dq, Dt)


def test_nonlinear_demag_limited():
    """A high-permeability soft-iron disk driven at H_ext is demag-limited to <M> ~ H_ext / D = 2 H_ext
    (D=1/2), NOT material-saturated -- the nonlinear Picard path converges to it."""
    H = np.array([0, 50, 200, 1000, 5000, 50000, 500000.0])
    B = np.array([0, 0.6, 1.2, 1.6, 1.85, 1.98, 2.05])
    with ng.TaskManager():
        r = m2.solve_planar_demag(_disk(0.1), bh_table=np.c_[H, B], H_ext=(2000.0, 0.0), nl_tol=1e-4)
    assert r["nonlinear"] and r["iters"] >= 1
    assert abs(r["M_avg"][0] - 4000.0) < 40.0, r["M_avg"]   # ~2*H_ext, demag-limited
