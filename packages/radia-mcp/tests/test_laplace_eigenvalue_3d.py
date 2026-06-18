"""3-D scalar Laplace-Dirichlet eigenvalues -- the modal decay spectrum of the heat equation
and the squared acoustic/quantum eigen-wavenumbers. Gated on closed forms: the unit CUBE
(lambda = pi^2(p^2+q^2+r^2), lowest 3 pi^2) and a BALL of radius R (lowest (pi/R)^2, next
(4.493409/R)^2 x3). The 3-D scalar companion of the 2-D waveguide cutoff eigensolver.
"""
import math
import os
import sys

import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.occ import unit_cube, Sphere, Pnt, OCCGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import laplace_dirichlet_eigenvalues


def test_unit_cube_lowest_eigenvalue():
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.18))
    lam = laplace_dirichlet_eigenvalues(mesh, 4, order=3)
    # lowest is 3 pi^2 (p=q=r=1); next is 6 pi^2 (one index =2, x3 degenerate)
    assert lam[0] == pytest.approx(3.0 * math.pi**2, rel=3e-3)
    assert lam[1] == pytest.approx(6.0 * math.pi**2, rel=1e-2)


def test_ball_lowest_eigenvalue():
    R = 1.0
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.3))
    mesh.Curve(4)                                   # curved tets for the round boundary
    lam = laplace_dirichlet_eigenvalues(mesh, 5, order=3)
    assert lam[0] == pytest.approx((math.pi / R)**2, rel=3e-3)        # radial s-mode
    # next eigenvalue is (4.493409/R)^2 (first l=1 mode), triply degenerate
    assert lam[1] == pytest.approx((4.493409 / R)**2, rel=1e-2)


def test_ball_eigenvalue_scales_with_radius():
    # lambda_1 = (pi/R)^2 scales as 1/R^2: a radius-2 ball has 1/4 the lowest eigenvalue
    def lam1(R):
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), R)).GenerateMesh(maxh=0.3 * R))
        mesh.Curve(4)
        return laplace_dirichlet_eigenvalues(mesh, 1, order=3)[0]
    l1 = lam1(1.0)
    l2 = lam1(2.0)
    assert l2 / l1 == pytest.approx(0.25, rel=2e-2)


def test_eigenvalues_ascending_and_positive():
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.2))
    lam = laplace_dirichlet_eigenvalues(mesh, 5, order=2)
    assert all(x > 0 for x in lam)
    assert all(b >= a for a, b in zip(lam, lam[1:]))
