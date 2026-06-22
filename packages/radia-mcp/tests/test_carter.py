r"""Carter coefficient (slotting -> effective air-gap increase) -- regression test (#57).

A slot opening adds reluctance to the air gap, so it behaves magnetically larger: g_eff = k_C g,
k_C = tau_s/(tau_s - gamma g), gamma = (b_o/g)^2/(5 + b_o/g). Closed-form helper (tool-independent)
+ a scalar-magnetic-potential gap-permeance FE (k_C = P_smooth/P_slot)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import carter_coefficient, effective_air_gap


def test_carter_closed_form():
    tau_s, g = 0.012, 0.001
    # no slot opening -> k_C = 1 (no correction)
    assert math.isclose(carter_coefficient(tau_s, g, 0.0), 1.0, rel_tol=1e-12)
    # k_C >= 1 and grows with the opening
    kcs = [carter_coefficient(tau_s, g, bo) for bo in (0.001, 0.002, 0.004, 0.006)]
    assert all(k >= 1.0 for k in kcs)
    assert all(b > a for a, b in zip(kcs, kcs[1:]))
    # effective gap = k_C * g
    bo = 0.004
    assert math.isclose(effective_air_gap(tau_s, g, bo), carter_coefficient(tau_s, g, bo) * g, rel_tol=1e-12)
    # textbook value (b_o/g = 4): gamma = 16/9, k_C = 12/(12 - 16/9) ~ 1.174
    assert math.isclose(carter_coefficient(0.012, 0.001, 0.004), 1.1739, abs_tol=2e-3)
