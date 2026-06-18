"""Axisymmetric DC resistance (current-flow dual of the axi capacitance) via the power method.
Gated on closed forms: the uniform cylinder R = h/(sigma pi b^2) (exact, 1-D current) and the
disk-contact SPREADING / constriction resistance R = 1/(4 sigma a) (Holm) -- a partial-electrode
mixed-BC axisymmetric problem with a current-crowding rim singularity.
"""
import math
import os
import sys

import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.scalar_fem2d import solve_poisson_axi, resistance_axi


def test_uniform_cylinder_resistance_exact():
    # cylinder radius b, height h, axial current: V=1 at z=h, V=0 at z=0, sides insulated.
    # exact R = h/(sigma pi b^2) (uniform 1-D field) -> resistance_axi reproduces it exactly.
    b, h, sigma = 0.4, 1.3, 2.5
    geo = SplineGeometry()
    p = [geo.AppendPoint(*xy) for xy in [(0, 0), (b, 0), (b, h), (0, h)]]
    geo.Append(["line", p[0], p[1]], bc="bot", leftdomain=1, rightdomain=0)   # z=0  -> V=0
    geo.Append(["line", p[1], p[2]], bc="side", leftdomain=1, rightdomain=0)  # r=b  insulated
    geo.Append(["line", p[2], p[3]], bc="top", leftdomain=1, rightdomain=0)   # z=h  -> V=1
    geo.Append(["line", p[3], p[0]], bc="axis", leftdomain=1, rightdomain=0)  # r=0  natural
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.1))
    V = solve_poisson_axi(mesh, sigma, {"top": 1.0, "bot": 0.0}, order=2)
    R = resistance_axi(V, mesh, sigma, 1.0)
    R_exact = h / (sigma * math.pi * b * b)
    assert R == pytest.approx(R_exact, rel=1e-3)


def _spreading_R(a=1.0, sigma=1.0, L=120.0, order=3):
    geo = SplineGeometry()
    p = [geo.AppendPoint(*xy) for xy in [(0, 0), (a, 0), (L, 0), (L, L), (0, L)]]
    geo.Append(["line", p[0], p[1]], bc="contact", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[1], p[2]], bc="ins", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[2], p[3]], bc="far", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[3], p[4]], bc="far", leftdomain=1, rightdomain=0)
    geo.Append(["line", p[4], p[0]], bc="axis", leftdomain=1, rightdomain=0)
    mesh = ng.Mesh(geo.GenerateMesh(maxh=L / 15.0))
    for _ in range(3):                          # refine toward the current-crowding rim (a,0)
        for el in mesh.Elements(ng.VOL):
            cx = sum(mesh[v].point[0] for v in el.vertices) / 3.0
            cy = sum(mesh[v].point[1] for v in el.vertices) / 3.0
            mesh.SetRefinementFlag(el, ((cx - a) ** 2 + cy ** 2) < 9.0)
        mesh.Refine()
    V = solve_poisson_axi(mesh, sigma, {"contact": 1.0, "far": 0.0}, order=order)
    return resistance_axi(V, mesh, sigma, 1.0)


def test_disk_contact_spreading_resistance_holm():
    # Holm constriction resistance of a circular contact radius a on a half-space: R = 1/(4 sigma a)
    for a, sigma in ((1.0, 1.0), (0.5, 1.0)):
        R = _spreading_R(a=a, sigma=sigma)
        R_exact = 1.0 / (4.0 * sigma * a)
        assert R == pytest.approx(R_exact, rel=3e-2)      # finite domain + rim singularity ~1%
        assert R < R_exact                                 # finite far boundary -> slightly low


def test_resistance_scales_inverse_sigma():
    # R = 1/(4 sigma a): doubling conductivity halves the resistance
    r1 = _spreading_R(a=1.0, sigma=1.0)
    r2 = _spreading_R(a=1.0, sigma=2.0)
    assert r2 / r1 == pytest.approx(0.5, rel=1e-2)
