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
                                           induction_machine_breakdown,
                                           induction_machine_slip_sweep_summary,
                                           induction_machine_noload_lockedrotor_summary)

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


def test_slip_sweep_summary_matches_femm_teaching_contract():
    summary = induction_machine_slip_sweep_summary(
        V1, R1, X1, R2, X2, Xm,
        line_frequency_hz=50.0,
        pole_pairs=2,
        slips=[0.005, 0.03, 0.1, 0.2, 1.0],
    )
    rows = {row["slip"]: row for row in summary["rows"]}

    assert summary["status"] == "ok"
    assert summary["policy"] == "induction_machine_slip_frequency_equivalent_circuit_gate"
    assert summary["synchronous_speed_rpm"] == pytest.approx(1500.0)
    assert rows[0.03]["slip_frequency_hz"] == pytest.approx(1.5)
    assert rows[1.0]["mechanical_speed_rpm"] == pytest.approx(0.0)
    assert summary["breakdown"]["slip"] == pytest.approx(0.147782600685332)
    assert summary["breakdown"]["torque_Nm"] == pytest.approx(189.32798248108452)
    assert summary["rotor_resistance_scaling"]["breakdown_slip"] == pytest.approx(
        2.0 * summary["breakdown"]["slip"]
    )
    assert summary["rotor_resistance_scaling"]["breakdown_torque_Nm"] == pytest.approx(
        summary["breakdown"]["torque_Nm"]
    )


def test_noload_lockedrotor_summary_extracts_equivalent_circuit_parameters():
    summary = induction_machine_noload_lockedrotor_summary(
        no_load_voltage_ll=400.0,
        no_load_current_line=5.0,
        no_load_power_total=900.0,
        locked_voltage_ll=90.0,
        locked_current_line=20.0,
        locked_power_total=1200.0,
        line_frequency_hz=50.0,
        pole_pairs=2,
        stator_resistance_ohm=0.4,
    )

    assert summary["status"] == "ok"
    assert summary["policy"] == "induction_machine_noload_lockedrotor_parameter_gate"
    assert summary["synchronous_speed_rpm"] == pytest.approx(1500.0)
    assert summary["no_load"]["phase_voltage_V"] == pytest.approx(400.0 / math.sqrt(3.0))
    assert summary["no_load"]["power_factor"] == pytest.approx(0.25980762113533157)
    assert summary["no_load"]["Rc_ohm"] == pytest.approx(177.7777777777778)
    assert summary["no_load"]["Xm_ohm"] == pytest.approx(47.830502044476084)
    assert summary["locked_rotor"]["Req_ohm"] == pytest.approx(1.0)
    assert summary["locked_rotor"]["Xeq_ohm"] == pytest.approx(2.3979157616563596)
    assert summary["locked_rotor"]["R2_referred_est_ohm"] == pytest.approx(0.6)
    assert summary["locked_rotor"]["equal_split_leakage_reactance_ohm"] == pytest.approx(1.1989578808281798)
