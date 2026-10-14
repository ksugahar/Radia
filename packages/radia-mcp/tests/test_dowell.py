r"""Dowell's equation -- AC resistance of an m-layer winding -- regression test (#60).

m foils in series; the leakage field builds layer by layer -> skin + proximity loss. Dowell:
F_R = Delta[(sinh2D+sin2D)/(cosh2D-cos2D) + (2/3)(m^2-1)(sinhD-sinD)/(coshD+cosD)], D=Delta=h/delta.
m=1 reduces to the deep-bar factor (#55); the proximity term grows ~m^2. Closed-form helper
(tool-independent) + an m-conductor series eddy FE."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (dowell_resistance_factor, deep_bar_resistance_factor,
                                           deep_bar_skin_ratio)

SIGMA = 5.8e7
W, H = 0.010, 0.002


def test_dowell_closed_form():
    # m=1 IS the single-conductor deep-bar factor (continuity with #55)
    for D in (0.5, 1.0, 2.0):
        assert math.isclose(dowell_resistance_factor(D, 1), deep_bar_resistance_factor(D), rel_tol=1e-12)
    # Delta -> 0: F_R -> 1 (DC) for any m
    assert math.isclose(dowell_resistance_factor(1e-4, 5), 1.0, rel_tol=1e-6)
    # proximity grows with layers (monotone in m) and with Delta
    assert dowell_resistance_factor(1.0, 1) < dowell_resistance_factor(1.0, 2) < dowell_resistance_factor(1.0, 4)
    assert dowell_resistance_factor(0.5, 4) < dowell_resistance_factor(1.0, 4) < dowell_resistance_factor(1.5, 4)
    # the proximity term is (2/3)(m^2-1)*Delta*prox -> F_R(m) - F_R(1) scales as (m^2-1)
    D = 1.0
    extra2 = dowell_resistance_factor(D, 2) - dowell_resistance_factor(D, 1)
    extra4 = dowell_resistance_factor(D, 4) - dowell_resistance_factor(D, 1)
    assert math.isclose(extra4 / extra2, (4**2 - 1) / (2**2 - 1), rel_tol=1e-9)   # 15/3 = 5
