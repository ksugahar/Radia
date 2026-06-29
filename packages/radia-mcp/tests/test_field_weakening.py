r"""Field weakening operating regions for salient PM machines.

Pure dq theory: below base speed MTPA is voltage-feasible; above it the optimum
moves onto the current-circle / voltage-ellipse boundary, and a wide-CPSR machine
can enter MTPV when the characteristic current lies inside the current limit.
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")          # solve.py imports ngsolve at module load

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (
    base_speed_electrical,
    characteristic_current,
    dq_torque,
    field_weakening_operating_point,
    field_weakening_speed_capability,
    pm_drive_operating_point,
    pm_drive_efficiency_map_health,
    pm_drive_speed_sweep,
)


def _numeric_fw_argmax(lambda_m, Ld, Lq, Imax, Vmax, omega_e, pole_pairs):
    """Independent 1-D search over id; iq is set by the tighter active limit."""
    vlim = Vmax / omega_e

    def torque_at(id_):
        current_left = Imax * Imax - id_ * id_
        voltage_left = vlim * vlim - (Ld * id_ + lambda_m) ** 2
        if current_left < -1e-12 or voltage_left < -1e-12:
            return None
        iq_current = math.sqrt(max(0.0, current_left))
        iq_voltage = math.sqrt(max(0.0, voltage_left)) / Lq
        iq = min(iq_current, iq_voltage)
        torque = dq_torque(lambda_m, Ld, Lq, id_, iq, pole_pairs)
        if torque < 0.0:
            iq, torque = 0.0, 0.0
        return id_, iq, torque

    samples = 5000
    best = None
    best_i = None
    for i in range(samples + 1):
        id_ = -Imax + 2.0 * Imax * i / samples
        out = torque_at(id_)
        if out is not None and (best is None or out[2] > best[2]):
            best = out
            best_i = i
    if best is None:
        return None

    step = 2.0 * Imax / samples
    lo = max(-Imax, -Imax + 2.0 * Imax * max(0, best_i - 2) / samples)
    hi = min(Imax, -Imax + 2.0 * Imax * min(samples, best_i + 2) / samples)
    for _ in range(80):
        m1 = lo + (hi - lo) / 3.0
        m2 = hi - (hi - lo) / 3.0
        t1 = torque_at(m1)
        t2 = torque_at(m2)
        v1 = -math.inf if t1 is None else t1[2]
        v2 = -math.inf if t2 is None else t2[2]
        if v1 < v2:
            lo = m1
        else:
            hi = m2
    return torque_at(0.5 * (lo + hi))


def test_field_weakening_speed_capability_from_characteristic_current():
    wide = field_weakening_speed_capability(0.1, 8.0e-3, 20.0)
    assert wide["characteristic_current"] == pytest.approx(12.5)
    assert wide["infinite_speed_possible"] is True
    assert wide["finite_max_speed"] is False
    assert wide["mtpv_possible"] is True
    assert wide["current_margin"] == pytest.approx(7.5)
    assert wide["current_ratio"] == pytest.approx(1.6)

    finite = field_weakening_speed_capability(0.1, 0.8e-3, 20.0)
    assert finite["characteristic_current"] == pytest.approx(125.0)
    assert finite["infinite_speed_possible"] is False
    assert finite["finite_max_speed"] is True
    assert finite["mtpv_possible"] is False

    boundary = field_weakening_speed_capability(0.1, 5.0e-3, 20.0)
    assert boundary["characteristic_current"] == pytest.approx(20.0)
    assert boundary["infinite_speed_possible"] is True
    assert boundary["mtpv_possible"] is False
    assert boundary["current_margin"] == pytest.approx(0.0)


@pytest.mark.parametrize(
    "name,params,speed_multiple,region",
    [
        ("wide_mtpa", (0.1, 8.0e-3, 16.0e-3, 20.0, 120.0, 4), 0.5, "MTPA"),
        ("wide_fw", (0.1, 8.0e-3, 16.0e-3, 20.0, 120.0, 4), 2.0, "FW"),
        ("wide_mtpv", (0.1, 8.0e-3, 16.0e-3, 20.0, 120.0, 4), 10.0, "MTPV"),
        ("finite_fw", (0.1, 0.8e-3, 1.6e-3, 20.0, 120.0, 4), 1.2, "FW"),
    ],
)
def test_field_weakening_matches_independent_numeric_argmax(name, params, speed_multiple, region):
    lambda_m, Ld, Lq, Imax, Vmax, pole_pairs = params
    omega_base = base_speed_electrical(lambda_m, Ld, Lq, Imax, Vmax, pole_pairs)
    omega_e = speed_multiple * omega_base

    closed = field_weakening_operating_point(lambda_m, Ld, Lq, Imax, Vmax, omega_e, pole_pairs)
    numeric = _numeric_fw_argmax(lambda_m, Ld, Lq, Imax, Vmax, omega_e, pole_pairs)

    assert closed is not None, name
    assert numeric is not None, name
    id_, iq, torque, got_region = closed
    assert got_region == region
    assert torque == pytest.approx(numeric[2], rel=2e-4, abs=1e-9)
    assert id_ == pytest.approx(numeric[0], abs=2e-3)
    assert iq == pytest.approx(numeric[1], abs=2e-3)
    assert id_ * id_ + iq * iq <= Imax * Imax * (1.0 + 1e-9)
    assert (Ld * id_ + lambda_m) ** 2 + (Lq * iq) ** 2 <= (Vmax / omega_e) ** 2 * (1.0 + 1e-9)


def test_wide_cpsr_asymptote_and_finite_speed_cutoff():
    wide = (0.1, 8.0e-3, 16.0e-3, 20.0, 120.0, 4)
    lambda_m, Ld, Lq, Imax, Vmax, pole_pairs = wide
    omega_base = base_speed_electrical(lambda_m, Ld, Lq, Imax, Vmax, pole_pairs)
    id_, iq, torque, region = field_weakening_operating_point(
        lambda_m, Ld, Lq, Imax, Vmax, 1000.0 * omega_base, pole_pairs)
    assert region == "MTPV"
    assert id_ == pytest.approx(-characteristic_current(lambda_m, Ld), rel=5e-6)
    assert abs(iq) < 1e-3 * Imax
    assert torque < 0.02

    finite = (0.1, 0.8e-3, 1.6e-3, 20.0, 120.0, 4)
    lambda_m, Ld, Lq, Imax, Vmax, pole_pairs = finite
    omega_base = base_speed_electrical(lambda_m, Ld, Lq, Imax, Vmax, pole_pairs)
    assert field_weakening_operating_point(
        lambda_m, Ld, Lq, Imax, Vmax, 2.0 * omega_base, pole_pairs) is None
    assert _numeric_fw_argmax(lambda_m, Ld, Lq, Imax, Vmax, 2.0 * omega_base, pole_pairs) is None


def test_pm_drive_operating_point_rows_match_selector_and_terminal_quantities():
    lambda_m, Ld, Lq, Imax, Vmax, pole_pairs = (0.1, 8.0e-3, 16.0e-3, 20.0, 120.0, 4)
    omega_base = base_speed_electrical(lambda_m, Ld, Lq, Imax, Vmax, pole_pairs)

    mtpa = pm_drive_operating_point(lambda_m, Ld, Lq, Imax, Vmax, 0.5 * omega_base, pole_pairs)
    fw = pm_drive_operating_point(lambda_m, Ld, Lq, Imax, Vmax, 2.0 * omega_base, pole_pairs)
    mtpv = pm_drive_operating_point(lambda_m, Ld, Lq, Imax, Vmax, 10.0 * omega_base, pole_pairs)

    assert mtpa["region"] == "MTPA"
    assert fw["region"] == "FW"
    assert mtpv["region"] == "MTPV"
    for row in (mtpa, fw, mtpv):
        assert row["feasible"] is True
        assert row["selected_torque"] == pytest.approx(row["torque"], rel=1e-12)
        assert row["current_utilization"] <= 1.0 + 1e-12
        assert row["voltage_utilization_lossless"] <= 1.0 + 1e-12
        assert row["voltage_utilization"] == pytest.approx(row["voltage_utilization_lossless"])
        assert row["omega_mech"] == pytest.approx(row["omega_e"] / pole_pairs)
    assert fw["current_utilization"] == pytest.approx(1.0, rel=1e-9)
    assert fw["voltage_utilization"] == pytest.approx(1.0, rel=1e-9)
    assert mtpv["current_utilization"] < 1.0


def test_pm_drive_speed_sweep_records_infeasible_rows():
    finite = pm_drive_speed_sweep(
        0.1,
        0.8e-3,
        1.6e-3,
        20.0,
        120.0,
        4,
        speed_multiples=(0.5, 1.0, 1.2, 2.0),
    )
    assert finite["speed_capability"]["finite_max_speed"] is True
    assert finite["rows"][0]["region"] == "MTPA"
    assert finite["rows"][2]["region"] == "FW"
    assert finite["rows"][3]["region"] == "infeasible"
    assert finite["rows"][3]["feasible"] is False
    assert finite["rows"][3]["speed_multiple"] == pytest.approx(2.0)

    wide = pm_drive_speed_sweep(
        0.1,
        8.0e-3,
        16.0e-3,
        20.0,
        120.0,
        4,
        speed_multiples=(0.5, 2.0, 10.0),
        R=0.05,
    )
    assert [row["region"] for row in wide["rows"]] == ["MTPA", "FW", "MTPV"]
    assert wide["rows"][0]["P_cu"] > 0.0
    assert wide["rows"][0]["voltage_utilization"] > wide["rows"][0]["voltage_utilization_lossless"]


def test_jmag_pm_drive_speed_sweep_gate_records_region_and_voltage_contract():
    sweep = pm_drive_speed_sweep(
        0.1,
        8.0e-3,
        16.0e-3,
        20.0,
        120.0,
        4,
        speed_multiples=(0.5, 1.0, 2.0, 10.0),
        R=0.05,
    )
    rows = sweep["rows"]

    assert sweep["omega_base"] == pytest.approx(455.32785624862885)
    assert [row["region"] for row in rows] == ["MTPA", "FW", "FW", "MTPV"]
    assert rows[0]["torque"] == pytest.approx(18.851963084409263)
    assert rows[2]["torque"] == pytest.approx(11.410462301017828)
    assert rows[3]["id"] == pytest.approx(-12.91999276662178)
    assert rows[3]["iq"] == pytest.approx(1.633723744208235)
    assert rows[3]["current_utilization"] < 0.7
    assert rows[1]["voltage_utilization_lossless"] == pytest.approx(1.0)
    assert rows[1]["voltage_utilization"] > 1.0
    assert rows[2]["voltage_utilization"] > rows[2]["voltage_utilization_lossless"]


def test_jmag_efficiency_map_health_gate_summarizes_drive_rows():
    sweep = pm_drive_speed_sweep(
        0.1,
        8.0e-3,
        16.0e-3,
        20.0,
        120.0,
        4,
        speed_multiples=(0.5, 1.0, 2.0, 10.0),
        R=0.05,
    )
    health = pm_drive_efficiency_map_health(sweep)

    assert health["status"] == "ok"
    assert health["policy"] == "pm_drive_efficiency_map_health_gate"
    assert health["region_sequence"] == ["MTPA", "FW", "FW", "MTPV"]
    assert health["region_counts"] == {"MTPA": 1, "FW": 2, "MTPV": 1}
    assert health["max_efficiency"] == pytest.approx(0.994425732150756)
    assert health["max_efficiency_region"] == "MTPV"
    assert health["max_torque_Nm"] == pytest.approx(18.85196308440927)
    assert health["max_torque_region"] in {"MTPA", "FW"}
    assert health["max_output_power_W"] == pytest.approx(2597.750669164122)
    assert health["max_output_power_region"] == "FW"
    assert health["max_power_balance_rel_error"] < 1.0e-15
    assert health["max_speed_contract_rel_error"] == pytest.approx(0.0)
