"""Frequency-domain (harmonic) heat conduction = thermal waves -- the complex-phasor sibling of
the transient heat solve. Gated on the exact 1-D semi-infinite damped wave
T(x) = T0 exp(-(1+i) x/delta), delta = sqrt(2 alpha/omega): |T| = exp(-x/delta), phase lag x/delta.
"""
import cmath
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

from radia_mcp.radia_ngsolve.multiphysics import solve_heat_harmonic, thermal_penetration_depth


def _strip(L=1.0, H=0.1, maxh=0.02):
    geo = SplineGeometry()
    p = [geo.AppendPoint(*xy) for xy in [(0, 0), (L, 0), (L, H), (0, H)]]
    geo.Append(["line", p[0], p[1]], bc="bot", leftdomain=1)
    geo.Append(["line", p[1], p[2]], bc="far", leftdomain=1)
    geo.Append(["line", p[2], p[3]], bc="top", leftdomain=1)
    geo.Append(["line", p[3], p[0]], bc="src", leftdomain=1)
    return ng.Mesh(geo.GenerateMesh(maxh=maxh))


def test_penetration_depth_formula():
    assert thermal_penetration_depth(1.0, 1.0, 200.0) == pytest.approx(math.sqrt(2.0 / 200.0))
    # delta scales as 1/sqrt(omega): 4x frequency -> half depth
    d1 = thermal_penetration_depth(2.0, 3.0, 50.0)
    d2 = thermal_penetration_depth(2.0, 3.0, 200.0)
    assert d2 / d1 == pytest.approx(0.5, rel=1e-12)


def test_thermal_wave_complex_profile():
    k = rho_cp = 1.0
    delta = 0.1
    omega = 2.0 * (k / rho_cp) / delta**2
    mesh = _strip()
    gfu = solve_heat_harmonic(mesh, k, rho_cp, omega, {"src": 1.0, "far": 0.0}, order=3)
    H = 0.1
    for xn in (1.0, 2.0, 3.0):                       # x in penetration depths
        val = gfu(mesh(xn * delta, H / 2))
        exact = cmath.exp(-(1 + 1j) * xn)            # exp(-(1+i) x/delta)
        assert abs(val - exact) < 5e-3, (xn, val, exact)
        assert abs(val) == pytest.approx(math.exp(-xn), rel=5e-3)   # |T| = exp(-x/delta)


def test_higher_frequency_decays_faster():
    k = rho_cp = 1.0
    mesh = _strip()
    g_lo = solve_heat_harmonic(mesh, k, rho_cp, 100.0, {"src": 1.0, "far": 0.0}, order=3)
    g_hi = solve_heat_harmonic(mesh, k, rho_cp, 400.0, {"src": 1.0, "far": 0.0}, order=3)
    xp = 0.15
    assert abs(g_hi(mesh(xp, 0.05))) < abs(g_lo(mesh(xp, 0.05)))   # higher omega -> shallower
