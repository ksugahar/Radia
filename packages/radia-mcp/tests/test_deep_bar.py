r"""Deep-bar (slot-conductor) AC resistance & reactance -- regression test (#55).

A solid conductor filling an iron slot, at frequency, has its current pushed toward the slot
opening by the slot-leakage field (1-D skin effect over the slot depth). Field's closed form:
k_R = Rac/Rdc = xi(sinh2xi+sin2xi)/(cosh2xi-cos2xi), k_X = Lac/Ldc = (3/2xi)(sinh2xi-sin2xi)/
(cosh2xi-cos2xi), xi=h/delta. Closed-form helpers (tool-independent) + a current-driven eddy FE.
The AC sequel to the DC slot leakage (#54, mu0 h/3w = the xi->0 limit)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (deep_bar_skin_ratio, deep_bar_resistance_factor,
                                           deep_bar_reactance_factor, slot_leakage_permeance,
                                           MU0)

SIGMA = 5.8e7
W, H = 0.004, 0.020


def test_deep_bar_closed_form():
    # xi -> 0: DC limit, both factors -> 1
    assert math.isclose(deep_bar_resistance_factor(1e-4), 1.0, rel_tol=1e-9)
    assert math.isclose(deep_bar_reactance_factor(1e-4), 1.0, rel_tol=1e-9)
    # monotonic: k_R rises, k_X falls with xi
    xis = [0.5, 1.0, 2.0, 4.0, 8.0]
    kR = [deep_bar_resistance_factor(x) for x in xis]
    kX = [deep_bar_reactance_factor(x) for x in xis]
    assert all(b > a for a, b in zip(kR, kR[1:]))
    assert all(b < a for a, b in zip(kX, kX[1:]))
    # high-xi asymptotes: k_R -> xi, k_X -> 3/(2 xi)
    assert math.isclose(deep_bar_resistance_factor(10.0), 10.0, rel_tol=2e-2)
    assert math.isclose(deep_bar_reactance_factor(10.0), 3.0 / (2 * 10.0), rel_tol=2e-2)
    # skin ratio: xi = h sqrt(pi f mu0 sigma); doubling h doubles xi, 4x f doubles xi
    xi = deep_bar_skin_ratio(H, 50.0, SIGMA)
    assert math.isclose(xi, H * math.sqrt(math.pi * 50.0 * MU0 * SIGMA), rel_tol=1e-12)
    assert math.isclose(deep_bar_skin_ratio(2 * H, 50.0, SIGMA), 2 * xi, rel_tol=1e-12)
    assert math.isclose(deep_bar_skin_ratio(H, 200.0, SIGMA), 2 * xi, rel_tol=1e-12)
