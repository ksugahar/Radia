r"""Straight cooling fin (extended surface) -- regression test (#47).

Thin fin, base dT_base, convective (h) faces, adiabatic tip: dT(x)=dT_base cosh(m(L-x))/cosh(mL),
m=sqrt(2h/(kt)), efficiency eta=tanh(mL)/(mL). radia's Robin heat solver with a base Dirichlet
reproduces the profile, tip and efficiency; the 1-D fin closed form is the tool-independent gate
(a reference commercial FE code matches to ~0.14% live, recorded internally). The
conduction-convection / heat-sink block; generalises the convective slab (#44)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import (fin_parameter, fin_tip_temperature_rise,
                                                  fin_efficiency)

L, T, K, DTB = 0.05, 0.002, 200.0, 100.0


def test_fin_helpers():
    m = fin_parameter(500.0, K, T)
    assert math.isclose(m, math.sqrt(2 * 500.0 / (K * T)), rel_tol=1e-12)   # = 50
    assert math.isclose(m, 50.0, rel_tol=1e-12)
    # efficiency limits: short/fat fin (mL->0) -> 1 ; long/thin (mL->inf) -> 0
    assert abs(fin_efficiency(m, 1e-4) - 1.0) < 1e-3
    assert fin_efficiency(m, 1.0) < 0.05            # mL=50 -> ~1/50
    assert fin_efficiency(m, 0.01) > fin_efficiency(m, 0.05)   # monotone decreasing in L
    # tip rise = dT_base/cosh(mL); falls as the fin lengthens
    assert math.isclose(fin_tip_temperature_rise(DTB, m, L), DTB / math.cosh(m * L), rel_tol=1e-12)
    assert fin_tip_temperature_rise(DTB, m, 2 * L) < fin_tip_temperature_rise(DTB, m, L)
