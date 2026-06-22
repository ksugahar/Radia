r"""Induction-machine equivalent-circuit torque-slip + breakdown -- regression test (#40).

The single-cage Thevenin equivalent circuit T(s) and its closed-form breakdown
(s_max, T_max) are checked against an INDEPENDENT numeric slip sweep, plus the classic
invariant T_max independent of rotor resistance R2 (only s_max ∝ R2). Pure circuit theory
-> tool-independent. The induction-machine companion to the PM dq blocks (#26 MTPA, #37 FW)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (induction_machine_thevenin,
                                           induction_machine_torque,
                                           induction_machine_breakdown)

# representative 4-pole 50 Hz IM (SI, referred to stator)
V1, R1, X1, R2, X2, Xm = 230.0, 0.5, 1.0, 0.3, 1.0, 30.0
OMEGA_S = 2 * math.pi * 50.0 / 2     # synchronous mechanical speed (pole_pairs=2)


def _numeric_breakdown(R2_=R2, n=400000):
    best = max(((i / n, induction_machine_torque(V1, R1, X1, R2_, X2, Xm, OMEGA_S, i / n))
                for i in range(1, n + 1)), key=lambda t: t[1])
    return best


def test_thevenin_reduction():
    Vth, Rth, Xth = induction_machine_thevenin(V1, R1, X1, Xm)
    assert 0 < Vth < V1                       # source is attenuated by Xm/(X1+Xm)
    assert Rth > 0 and Xth > 0
    # large Xm -> Vth -> V1, Rth -> R1, Xth -> X1
    Vth2, Rth2, Xth2 = induction_machine_thevenin(V1, R1, X1, 1e7)
    assert abs(Vth2 - V1) / V1 < 1e-3 and abs(Rth2 - R1) < 1e-3 and abs(Xth2 - X1) < 1e-3
