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
