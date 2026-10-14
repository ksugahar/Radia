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
                                           field_weakening_speed_capability,
                                           dq_torque,
                                           pm_temperature_demag_sweep_summary)

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


def test_continuous_loop_slot45_pm_short_circuit_fault_table():
    R, Ld, Lq, lm, p, Imax = 0.05, 0.008, 0.016, 0.1, 4, 20.0
    cap = field_weakening_speed_capability(lm, Ld, Imax)
    assert cap["characteristic_current"] == pytest.approx(12.5)
    assert cap["infinite_speed_possible"] is True
    assert cap["mtpv_possible"] is True
    assert cap["current_margin"] == pytest.approx(7.5)
    assert cap["current_ratio"] == pytest.approx(1.6)

    omega_e_values = [1.0, 2.0, 5.0, 6.25, 10.0, 20.0, 50.0, 100.0, 500.0, 1000.0, 5000.0]
    rows = [short_circuit_operating_point(R, Ld, Lq, lm, omega_e, p) | {
        "omega_e": omega_e,
        "omega_mech": omega_e / p,
    } for omega_e in omega_e_values]
    peak = max(rows, key=lambda row: abs(row["torque"]))
    high = rows[-1]

    assert peak["omega_e"] == pytest.approx(6.25)
    assert peak["id"] == pytest.approx(-8.333333333333332)
    assert peak["iq"] == pytest.approx(-4.166666666666666)
    assert peak["current_ratio_to_characteristic"] == pytest.approx(0.7453559924999298)
    assert peak["d_axis_demag_fraction"] == pytest.approx(0.6666666666666665)
    assert peak["torque"] == pytest.approx(-4.166666666666666)

    assert high["omega_e"] == pytest.approx(5000.0)
    assert high["id"] == pytest.approx(-12.49999023438263)
    assert high["iq"] == pytest.approx(-0.0078124938964891436)
    assert high["current_ratio_to_characteristic"] == pytest.approx(0.9999994140629387)
    assert high["d_axis_demag_fraction"] == pytest.approx(0.9999992187506104)
    assert high["torque"] == pytest.approx(-0.009374989013683319)

    max_terminal_residual = max(max(abs(row["vd_residual"]), abs(row["vq_residual"])) for row in rows)
    max_speed_contract_error = max(abs(row["omega_e"] - p * row["omega_mech"]) for row in rows)
    assert max_terminal_residual <= 1.0e-12
    assert max_speed_contract_error <= 1.0e-12


def test_continuous_loop_slot54_fault_current_demag_screening_gate():
    R, Ld, Lq, lm, p, Imax = 0.05, 0.008, 0.016, 0.1, 4, 20.0
    cap = field_weakening_speed_capability(lm, Ld, Imax)
    high = short_circuit_operating_point(R, Ld, Lq, lm, 5000.0, p)
    sweep = pm_temperature_demag_sweep_summary(
        Br_20C=1.2,
        H_knee_20C=-9.0e5,
        temperature_C=120.0,
        magnet_len=0.004,
        gaps=(0.0005, 0.001, 0.002, 0.004, 0.008),
        iron_path=0.08,
        mu_r=1000.0,
        mu_rec=1.05,
    )

    assert cap["characteristic_current"] == pytest.approx(12.5)
    assert cap["current_margin"] == pytest.approx(7.5)
    assert high["id"] == pytest.approx(-12.49999023438263)
    assert high["d_axis_demag_fraction"] == pytest.approx(0.9999992187506104)
    assert high["current_ratio_to_characteristic"] == pytest.approx(0.9999994140629387)
    assert sweep["first_unsafe_gap_m"] == pytest.approx(0.008)
    assert sweep["safe_prefix_count"] == 4
    assert sweep["risk_label"] == "red"
    assert sweep["risk_summary"]["minimum_demag_margin_A_per_m"] == pytest.approx(-54988.18063331157)
    assert all(sweep["checks"].values())


if __name__ == "__main__":
    test_characteristic_current()
    test_short_circuit_equations_and_limit()
    test_braking_torque_peak_and_decay()
    test_nonsalient_critical_speed()
    test_continuous_loop_slot45_pm_short_circuit_fault_table()
    test_continuous_loop_slot54_fault_current_demag_screening_gate()
    print("[OK] short-circuit currents -> Ich, braking-torque peak/decay, critical speed validated.")
