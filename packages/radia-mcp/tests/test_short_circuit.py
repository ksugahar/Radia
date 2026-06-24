r"""PM machine 3-phase short-circuit current + braking torque + demag -- regression test (#49).

Pure dq circuit theory -> tool-independent. Checks: the shorted dq equations (vd=vq=0) are
satisfied, the high-speed limit |Isc| -> Ich = lambda_m/Ld (pure d-axis demagnetising), the
braking torque -> 0 at high speed and peaks at the critical speed (exactly omega_e = R/Ls for a
non-salient machine), and Ich = lambda_m/Ld. The fault/protection regime alongside the
operating-point #46."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (short_circuit_dq_currents, short_circuit_operating_point,
                                           characteristic_current,
                                           dq_torque)

LM, LD, LQ, R, P = 0.1, 0.5e-3, 1.5e-3, 0.05, 4


def test_characteristic_current():
    assert math.isclose(characteristic_current(LM, LD), LM / LD, rel_tol=1e-12)
    assert math.isclose(characteristic_current(LM, LD), 200.0, rel_tol=1e-12)


def test_short_circuit_equations_and_limit():
    we = 2 * math.pi * 100.0
    id_, iq = short_circuit_dq_currents(R, LD, LQ, LM, we)
    op = short_circuit_operating_point(R, LD, LQ, LM, we, P)
    assert op["id"] == pytest.approx(id_)
    assert op["iq"] == pytest.approx(iq)
    assert op["current_magnitude"] == pytest.approx(math.hypot(id_, iq))
    assert op["torque"] == pytest.approx(dq_torque(LM, LD, LQ, id_, iq, P))
    assert op["mechanical_power"] == pytest.approx(op["torque"] * we / P)
    assert abs(op["vd_residual"]) < 1e-9
    assert abs(op["vq_residual"]) < 1e-9
    # the shorted dq voltage equations hold: 0 = R id - we Lq iq ; 0 = R iq + we(Ld id + lm)
    assert abs(R * id_ - we * LQ * iq) < 1e-9
    assert abs(R * iq + we * (LD * id_ + LM)) < 1e-9
    # high-speed limit: id -> -Ich, iq -> 0, |Isc| -> Ich
    idh, iqh = short_circuit_dq_currents(R, LD, LQ, LM, 2 * math.pi * 1e6)
    Ich = characteristic_current(LM, LD)
    assert abs(idh + Ich) / Ich < 1e-3
    assert abs(iqh) < 1e-3 * Ich
    assert abs(math.hypot(idh, iqh) - Ich) / Ich < 1e-3
    oph = short_circuit_operating_point(R, LD, LQ, LM, 2 * math.pi * 1e6, P)
    assert abs(oph["current_ratio_to_characteristic"] - 1.0) < 1e-3
    assert abs(oph["d_axis_demag_fraction"] - 1.0) < 1e-3


def test_braking_torque_peak_and_decay():
    Tb = lambda we: dq_torque(LM, LD, LQ, *short_circuit_dq_currents(R, LD, LQ, LM, we), P)
    # braking torque (negative) -> 0 at high speed
    assert abs(Tb(2 * math.pi * 1e5)) < abs(Tb(2 * math.pi * 50))
    # it has an interior peak (rises from ~0 at low speed, falls at high speed)
    assert abs(Tb(2 * math.pi * 20)) > abs(Tb(2 * math.pi * 1))
    assert abs(Tb(2 * math.pi * 20)) > abs(Tb(2 * math.pi * 500))


def test_nonsalient_critical_speed():
    # non-salient (Ld=Lq=Ls): T_sc = (3/2)p lm iq, iq ~ we/(R^2+we^2 Ls^2) peaks at we = R/Ls
    Ls = 1.0e-3
    Tb = lambda we: abs(dq_torque(LM, Ls, Ls, *short_circuit_dq_currents(R, Ls, Ls, LM, we), P))
    we_crit = R / Ls
    # the analytic peak: torque just below/above we_crit is smaller
    assert Tb(we_crit) > Tb(0.7 * we_crit) and Tb(we_crit) > Tb(1.4 * we_crit)
    # |Isc| -> lm/Ls at high speed
    idh, iqh = short_circuit_dq_currents(R, Ls, Ls, LM, 2 * math.pi * 1e6)
    assert abs(math.hypot(idh, iqh) - LM / Ls) / (LM / Ls) < 1e-3


if __name__ == "__main__":
    test_characteristic_current()
    test_short_circuit_equations_and_limit()
    test_braking_torque_peak_and_decay()
    test_nonsalient_critical_speed()
    print("[OK] short-circuit currents -> Ich, braking-torque peak/decay, critical speed validated.")
