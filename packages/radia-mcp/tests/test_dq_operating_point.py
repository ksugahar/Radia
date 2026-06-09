r"""dq steady-state operating point (voltage, power, power factor) -- regression test (#46).

Pure dq circuit theory -> tool-independent. Checks: power balance P_in = P_cu + P_em (machine
precision), P_em = torque*omega_mech, torque matches dq_torque (#26), the lossless (R=0) limit
P_in=P_em with efficiency 1, the phasor voltage relations, and PF in [0,1]. The terminal-quantity
(V, I, P, PF) layer that MTPA (#26) / field weakening (#37) feed into."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (dq_voltages, dq_operating_point, dq_torque,
                                           mtpa_operating_point)

LM, LD, LQ, R, P = 0.1, 0.5e-3, 1.5e-3, 0.05, 4
WE = 2 * math.pi * 200.0


def test_power_balance_and_torque():
    _, id_, iq, _ = mtpa_operating_point(LM, LD, LQ, 15.0, P)
    op = dq_operating_point(R, LD, LQ, LM, id_, iq, WE, P)
    # input power = copper loss + air-gap power (exact)
    assert abs(op["P_in"] - (op["P_cu"] + op["P_em"])) < 1e-6 * abs(op["P_in"])
    # air-gap power = torque * mechanical speed
    assert abs(op["P_em"] - op["torque"] * op["omega_mech"]) < 1e-6 * abs(op["P_em"])
    # torque matches the dq torque (#26)
    assert math.isclose(op["torque"], dq_torque(LM, LD, LQ, id_, iq, P), rel_tol=1e-12)
    # power factor in (0, 1]
    assert 0.0 < op["power_factor"] <= 1.0 + 1e-12


def test_lossless_limit():
    _, id_, iq, _ = mtpa_operating_point(LM, LD, LQ, 15.0, P)
    op = dq_operating_point(0.0, LD, LQ, LM, id_, iq, WE, P)   # R=0
    assert abs(op["P_in"] - op["P_em"]) < 1e-6 * abs(op["P_em"])   # no copper loss
    assert abs(op["efficiency"] - 1.0) < 1e-9
    assert abs(op["P_cu"]) < 1e-12


def test_phasor_relations():
    # vd = R id - we Lq iq ; vq = R iq + we(Ld id + lm)
    id_, iq = -3.0, 12.0
    vd, vq = dq_voltages(R, LD, LQ, LM, id_, iq, WE)
    assert math.isclose(vd, R * id_ - WE * LQ * iq, rel_tol=1e-12)
    assert math.isclose(vq, R * iq + WE * (LD * id_ + LM), rel_tol=1e-12)
    # non-salient pure q-axis (id=0): vd = -we Lq iq, vq = R iq + we lm
    vd0, vq0 = dq_voltages(R, 1e-3, 1e-3, LM, 0.0, 10.0, WE)
    assert math.isclose(vd0, -WE * 1e-3 * 10.0, rel_tol=1e-12)
    assert math.isclose(vq0, R * 10.0 + WE * LM, rel_tol=1e-12)
    # |V| grows ~linearly with speed (back-EMF dominated)
    _, id_, iq, _ = mtpa_operating_point(LM, LD, LQ, 15.0, P)
    v1 = dq_operating_point(R, LD, LQ, LM, id_, iq, WE, P)["Vmag"]
    v2 = dq_operating_point(R, LD, LQ, LM, id_, iq, 2 * WE, P)["Vmag"]
    assert v2 > 1.8 * v1   # roughly doubles


if __name__ == "__main__":
    test_power_balance_and_torque()
    test_lossless_limit()
    test_phasor_relations()
    print("[OK] dq operating point: power balance, lossless limit, phasor relations validated.")
