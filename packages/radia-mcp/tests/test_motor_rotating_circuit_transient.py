from __future__ import annotations

import copy

from radia_mcp.motor.rotating_circuit_transient_gate import (
    rotating_circuit_transient_gate,
)


def _row(index: int) -> dict:
    time_s = 0.001 * index
    angle_deg = 6.0 * 1000.0 * time_s
    phase = [2.0 + 0.1 * index, -1.0 - 0.05 * index, -1.0 - 0.05 * index]
    flux = [0.1 + 0.01 * index, -0.04, -0.06]
    branch_currents = [[10.0 + index, -10.0 - index + 1.0e-3] for _ in range(2)]
    branch_powers = [[2.0 + index, 2.0 + index + 1.0e-4] for _ in range(2)]
    components = [10.0 + index, 2.0, -0.5]
    return {
        "time_s": time_s,
        "angle_deg": angle_deg,
        "speed_rpm": 1000.0,
        "torque_nm": 1.0 + 0.1 * index,
        "phase_currents_a": phase,
        "phase_flux_linkages_wb": flux,
        "paired_branch_currents_a": branch_currents,
        "paired_branch_powers_w": branch_powers,
        "circuit_power_components_w": components,
        "reported_total_circuit_power_w": sum(components),
    }


def _summary() -> dict:
    return {
        "contract": {
            "phase_count": 3,
            "paired_branch_count": 2,
            "geometric_cycle_deg": 360.0,
            "endpoint_periodicity_policy": "require_state_match_before_fft",
            "expected_endpoint_state": "nonperiodic_transient",
            "fft_ready": False,
        },
        "rows": [_row(index) for index in range(61)],
    }


def test_rotating_circuit_gate_accepts_nonperiodic_geometric_cycle() -> None:
    result = rotating_circuit_transient_gate(_summary())
    assert result["status"] == "ok"
    assert result["endpoint_classification"] == "nonperiodic_transient"
    assert result["fft_ready"] is False


def test_rotating_circuit_gate_rejects_fft_claim_for_transient_endpoint() -> None:
    payload = _summary()
    payload["contract"]["fft_ready"] = True
    result = rotating_circuit_transient_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fft_readiness_matches_state"] is False


def test_rotating_circuit_gate_rejects_topology_and_power_corruption() -> None:
    payload = copy.deepcopy(_summary())
    payload["rows"][20]["paired_branch_currents_a"][0][1] = 5.0
    payload["rows"][20]["paired_branch_powers_w"][0][1] = -3.0
    payload["rows"][20]["reported_total_circuit_power_w"] += 4.0
    result = rotating_circuit_transient_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["paired_branch_currents_antisymmetric"] is False
    assert result["checks"]["paired_branch_powers_symmetric"] is False
    assert result["checks"]["reported_total_power_matches_component_sum"] is False
