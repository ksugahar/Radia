r"""Salient synchronous-machine power-angle (delta) torque + pull-out -- #43.

T(d) = (3p/we)[V E/Xd sin d + (V^2/2)(1/Xq-1/Xd) sin 2d] (field + reluctance torque). The
closed-form pull-out (delta_max, T_max) is checked against an INDEPENDENT numeric angle
sweep, plus the limits: non-salient -> 90 deg, pure reluctance (E=0) -> 45 deg, and the
field+reluctance decomposition. Pure circuit theory -> tool-independent. The grid/voltage-
frame companion to the current-frame MTPA (#26) and the induction-machine torque-slip (#40)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (synchronous_power_angle_torque,
                                           synchronous_pullout)

V, E, XD, XQ = 230.0, 276.0, 2.0, 1.2
WE, P = 2 * math.pi * 50.0, 2


def _numeric_max(V_, E_, Xd_, Xq_, n=400000):
    return max(((i / n * math.pi, synchronous_power_angle_torque(V_, E_, Xd_, Xq_, WE, P, i / n * math.pi))
                for i in range(1, n)), key=lambda t: t[1])


def test_field_plus_reluctance_decomposition():
    # total torque is exactly the field (sin d) + reluctance (sin 2d) terms
    d = math.radians(40)
    Tf = (3 * P / WE) * V * E / XD * math.sin(d)
    Tr = (3 * P / WE) * 0.5 * V * V * (1 / XQ - 1 / XD) * math.sin(2 * d)
    assert math.isclose(synchronous_power_angle_torque(V, E, XD, XQ, WE, P, d), Tf + Tr, rel_tol=1e-12)
    # reluctance term vanishes for a non-salient machine
    assert math.isclose(synchronous_power_angle_torque(V, E, 2.0, 2.0, WE, P, d),
                        (3 * P / WE) * V * E / 2.0 * math.sin(d), rel_tol=1e-12)
    # T(0)=0, generating torque rises with small delta
    assert abs(synchronous_power_angle_torque(V, E, XD, XQ, WE, P, 0.0)) < 1e-9
    assert synchronous_power_angle_torque(V, E, XD, XQ, WE, P, 0.1) > 0
