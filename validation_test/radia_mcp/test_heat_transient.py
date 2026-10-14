"""Transient (time-domain) heat conduction -- the theta-scheme complement of the steady
solve_heat_steady. Gated on the EXACT single Fourier-mode decay of the unit square
(T=0 on the boundary): T0=sin(pi x)sin(pi y) -> T = T0 exp(-alpha 2 pi^2 t), alpha=k/rho_cp.
Checks value agreement, Crank-Nicolson 2nd-order temporal convergence, and implicit-Euler stability.
"""
import math
import os
import sys

import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen")
from netgen.occ import unit_square

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import solve_heat_transient


def _setup():
    mesh = ng.Mesh(unit_square.GenerateMesh(maxh=0.05))
    T0 = ng.sin(math.pi * ng.x) * ng.sin(math.pi * ng.y)
    return mesh, T0


def _center(gfu, mesh):
    return gfu(mesh(0.5, 0.5))


def test_single_mode_decay_value():
    mesh, T0 = _setup()
    k = rho_cp = 1.0
    t_end = 0.05
    gfu = solve_heat_transient(mesh, T0, k, rho_cp, ".*", t_end, 100, order=3, theta=0.5)
    exact = 1.0 * math.exp(-(k / rho_cp) * 2.0 * math.pi**2 * t_end)   # sin(pi/2)^2 = 1
    assert _center(gfu, mesh) == pytest.approx(exact, rel=5e-5)


def test_crank_nicolson_second_order():
    mesh, T0 = _setup()
    k = rho_cp = 1.0
    t_end = 0.05
    exact = math.exp(-2.0 * math.pi**2 * t_end)
    errs = []
    for nst in (25, 50, 100):
        gfu = solve_heat_transient(mesh, T0, k, rho_cp, ".*", t_end, nst, order=3, theta=0.5)
        errs.append(abs(_center(gfu, mesh) - exact) / exact)
    # halving dt should quarter the error (2nd order); allow margin
    assert errs[0] / errs[1] > 3.0
    assert errs[1] / errs[2] > 3.0
    assert errs[-1] < 1e-4


def test_diffusivity_scaling():
    # doubling alpha = k/rho_cp doubles the decay rate -> T(t) with 2*alpha == T(2t) baseline
    mesh, T0 = _setup()
    t_end = 0.03
    g1 = solve_heat_transient(mesh, T0, 1.0, 1.0, ".*", t_end, 80, order=3)
    g2 = solve_heat_transient(mesh, T0, 2.0, 1.0, ".*", t_end, 80, order=3)   # alpha=2
    e1 = math.exp(-1.0 * 2 * math.pi**2 * t_end)
    e2 = math.exp(-2.0 * 2 * math.pi**2 * t_end)
    assert _center(g1, mesh) == pytest.approx(e1, rel=1e-3)
    assert _center(g2, mesh) == pytest.approx(e2, rel=1e-3)
    assert _center(g2, mesh) < _center(g1, mesh)        # faster decay


def test_implicit_euler_stable_and_decays():
    # theta=1 (backward Euler) is unconditionally stable: a coarse-dt run still decays monotonically
    mesh, T0 = _setup()
    gfu = solve_heat_transient(mesh, T0, 1.0, 1.0, ".*", 0.1, 10, order=2, theta=1.0)
    val = _center(gfu, mesh)
    assert 0.0 < val < 1.0                               # decayed, no blow-up (unconditionally stable)
    exact = math.exp(-2.0 * math.pi**2 * 0.1)
    # backward Euler is only 1st-order and this dt is coarse (rate*dt~0.2) -> ~20% high, but tracks
    assert val == pytest.approx(exact, rel=0.25)
