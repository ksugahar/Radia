"""Golden test: NONLINEAR HDiv-VIM damped Newton with a REAL tabulated BH curve (production #3 / real
material).

The earlier nonlinear tests used the analytic saturating curve M(H)=chi0 H/(1+chi0|H|/Msat).  Real soft
iron is given as a [[H,B]] TABLE (the same data Radia's MatSatIsoTab consumes).  solve_nonlinear_newton
now accepts bh_table=(Harr,Barr): it builds M(H) via a monotone PCHIP interpolation of B(H)
(M=B/mu0-H), with the physically-correct SATURATION beyond the table (|H|>Hmax: B extends with slope
mu0 => M -> M(Hmax) const, dM/dH -> 0 -- matching Radia's linear-B extension; raw PCHIP extrapolation
would blow up).

Locks that the table-driven operator Newton (routed through the production C++ analytic charge Gram)
reproduces the analytic uniform-sphere fixed point (with the SAME table M(H), demag D ~ 1/3) across
moderate -> deep saturation, and saturates to the table's Msat.
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
pytest.importorskip("scipy.interpolate")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

# a realistic soft-iron BH table [[H(A/m), B(T)], ...] (mu_r ~ 4000 at low H, saturating ~2.2 T)
_TABLE = [[0, 0], [50, 0.25], [100, 0.6], [150, 0.9], [200, 1.10], [300, 1.30], [500, 1.45],
          [800, 1.55], [1500, 1.65], [3000, 1.74], [6000, 1.82], [12000, 1.90], [30000, 1.98],
          [80000, 2.08], [2e5, 2.2]]
_HARR = [r[0] for r in _TABLE]
_BARR = [r[1] for r in _TABLE]


def _sphere(h=0.4):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_real_bh_table_newton_matches_analytic():
    mesh = _sphere()
    Mof, _, _, _, Mmax = nl._bh_table_funcs(_HARR, _BARR)
    with ng.TaskManager():
        for H0 in (1e4, 1e6, 3e6):
            Mh, nit, _ = nl.solve_nonlinear_newton(mesh, 0, 0, H0, bh_table=(_HARR, _BARR), maxit=150)
            Ma = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
            assert abs(Mh - Ma) < 1e-2 * abs(Ma), \
                f"table Newton {Mh:.1f} vs analytic {Ma:.1f} at H0={H0:.0e} (rel {abs(Mh-Ma)/abs(Ma):.1e})"
            # RT1 deep drive (H0>=3e6) runs the M-form limit cycle longer before settled-acceptance (~59 iters
            # at H0=3e6 vs <30 at RT0); the M is still the correct best-energy iterate (asserted above).
            assert nit < 80, f"table Newton not converging at H0={H0:.0e}: {nit} iters"


def test_real_bh_table_saturates_to_table_Msat():
    """Deep drive saturates M to the table's Msat = M(Hmax) (the beyond-table extension is physical,
    not a PCHIP blow-up)."""
    mesh = _sphere()
    _, _, _, _, Mmax = nl._bh_table_funcs(_HARR, _BARR)
    with ng.TaskManager():
        Mh, _, _ = nl.solve_nonlinear_newton(mesh, 0, 0, 5e6, bh_table=(_HARR, _BARR), maxit=150)
    assert 0.95 * Mmax < Mh < 1.02 * Mmax, f"deep-drive M {Mh:.1f} not saturated to table Msat {Mmax:.1f}"
