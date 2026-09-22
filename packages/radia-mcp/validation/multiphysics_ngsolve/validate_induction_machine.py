r"""Induction-machine equivalent-circuit torque-slip + breakdown -- regression test (#40).

The single-cage Thevenin equivalent circuit T(s) and its closed-form breakdown
(s_max, T_max) are checked against an INDEPENDENT numeric slip sweep, plus the classic
invariant T_max independent of rotor resistance R2 (only s_max ∝ R2). Pure circuit theory
-> tool-independent. The induction-machine companion to the PM dq blocks (#26 MTPA, #37 FW)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
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


def validate_thevenin_reduction():
    Vth, Rth, Xth = induction_machine_thevenin(V1, R1, X1, Xm)
    assert 0 < Vth < V1                       # source is attenuated by Xm/(X1+Xm)
    assert Rth > 0 and Xth > 0
    # large Xm -> Vth -> V1, Rth -> R1, Xth -> X1
    Vth2, Rth2, Xth2 = induction_machine_thevenin(V1, R1, X1, 1e7)
    assert abs(Vth2 - V1) / V1 < 1e-3 and abs(Rth2 - R1) < 1e-3 and abs(Xth2 - X1) < 1e-3


def validate_breakdown_matches_numeric():
    s_max, T_max = induction_machine_breakdown(V1, R1, X1, R2, X2, Xm, OMEGA_S)
    s_num, T_num = _numeric_breakdown()
    assert abs(s_max - s_num) / s_num < 1e-3, f"s_max {s_max} vs {s_num}"
    assert abs(T_max - T_num) / T_num < 1e-4, f"T_max {T_max} vs {T_num}"
    # T(s_max) == T_max
    assert abs(induction_machine_torque(V1, R1, X1, R2, X2, Xm, OMEGA_S, s_max) - T_max) / T_max < 1e-9


def validate_breakdown_torque_is_R2_invariant():
    s0, T0 = induction_machine_breakdown(V1, R1, X1, R2, X2, Xm, OMEGA_S)
    for k in (2.0, 3.0, 5.0):
        sk, Tk = induction_machine_breakdown(V1, R1, X1, k * R2, X2, Xm, OMEGA_S)
        assert abs(Tk - T0) / T0 < 1e-12               # T_max independent of R2
        assert abs(sk - k * s0) / (k * s0) < 1e-12     # s_max proportional to R2
        # the numeric peak of the R2-scaled curve also equals T0
        _, Tk_num = _numeric_breakdown(k * R2)
        assert abs(Tk_num - T0) / T0 < 1e-3


def validate_torque_curve_shape():
    s_max, T_max = induction_machine_breakdown(V1, R1, X1, R2, X2, Xm, OMEGA_S)
    T = lambda s: induction_machine_torque(V1, R1, X1, R2, X2, Xm, OMEGA_S, s)
    assert T(1e-6) < 1e-2 * T_max                       # -> 0 at synchronism
    assert 0 < T(1.0) < T_max                           # starting torque (s=1), below breakdown
    # monotone rising below s_max, falling above
    assert T(0.5 * s_max) < T(s_max) and T(2 * s_max) < T(s_max)


if __name__ == "__main__":
    validate_thevenin_reduction()
    validate_breakdown_matches_numeric()
    validate_breakdown_torque_is_R2_invariant()
    validate_torque_curve_shape()
    print("[OK] induction-machine torque-slip + breakdown validated (closed form == numeric).")
