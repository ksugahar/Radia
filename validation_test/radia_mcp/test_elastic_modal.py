"""Elastic free-vibration (modal) eigenvalues -- the dynamic companion of the static
linear-elasticity solver. Gated on the 1-D axial-rod wave modes (nu=0, transverse motion
suppressed): fixed-free omega_n=(2n-1)pi/(2L)sqrt(E/rho); fixed-fixed omega_n=n pi/L sqrt(E/rho).
The structural-NVH modes that an electromagnetic force spectrum can excite.
"""
import math
import os
import sys

import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.elasticity import elastic_modal_2d

L, H = 1.0, 0.08
E, RHO = 1.0, 1.0
C = math.sqrt(E / RHO)


def _rod_mesh(maxh=0.05):
    geo = SplineGeometry()
    p = [geo.AppendPoint(*xy) for xy in [(0, 0), (L, 0), (L, H), (0, H)]]
    geo.Append(["line", p[0], p[1]], bc="bot", leftdomain=1)
    geo.Append(["line", p[1], p[2]], bc="right", leftdomain=1)
    geo.Append(["line", p[2], p[3]], bc="top", leftdomain=1)
    geo.Append(["line", p[3], p[0]], bc="left", leftdomain=1)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_fixed_free_axial_rod():
    mesh = _rod_mesh()
    w = elastic_modal_2d(mesh, E, 0.0, RHO, 4, dirichletx="left",
                         dirichlety="left|right|top|bot", order=3)
    for n, wn in enumerate(w, start=1):
        exact = (2 * n - 1) * math.pi / (2 * L) * C
        assert wn == pytest.approx(exact, rel=1e-4)


def test_fixed_fixed_axial_rod():
    mesh = _rod_mesh()
    w = elastic_modal_2d(mesh, E, 0.0, RHO, 4, dirichletx="left|right",
                         dirichlety="left|right|top|bot", order=3)
    for n, wn in enumerate(w, start=1):
        exact = n * math.pi / L * C
        assert wn == pytest.approx(exact, rel=1e-4)


def test_frequency_scales_with_wave_speed():
    # omega ~ sqrt(E/rho): 4x stiffer material -> 2x frequency
    mesh = _rod_mesh()
    w1 = elastic_modal_2d(mesh, 1.0, 0.0, 1.0, 1, dirichletx="left",
                          dirichlety="left|right|top|bot", order=3)[0]
    w2 = elastic_modal_2d(mesh, 4.0, 0.0, 1.0, 1, dirichletx="left",
                          dirichlety="left|right|top|bot", order=3)[0]
    assert w2 / w1 == pytest.approx(2.0, rel=1e-3)


def test_modes_ascending_positive():
    mesh = _rod_mesh()
    w = elastic_modal_2d(mesh, E, 0.0, RHO, 5, dirichletx="left",
                         dirichlety="left|right|top|bot", order=2)
    assert all(x > 0 for x in w)
    assert all(b >= a for a, b in zip(w, w[1:]))
