import json
import math

from radia_mcp.motor.phase_flux_park_gate import (
    evaluate_phase_flux_park_alignment,
    phase_flux_park_alignment_gate,
)


def _balanced_flux(phase_error=0.0):
    angles = list(range(0, 181, 10))
    rows = []
    for angle in angles:
        theta = math.radians(2 * angle)
        rows.append([
            0.03 * math.cos(theta),
            0.03 * math.cos(theta - 2 * math.pi / 3 + phase_error),
            0.03 * math.cos(theta + 2 * math.pi / 3),
        ])
    return angles, rows


def test_phase_flux_park_gate_accepts_aligned_full_cycle():
    angles, rows = _balanced_flux()
    gate = evaluate_phase_flux_park_alignment(angles, rows, 2)
    assert gate["status"] == "ok"
    assert all(gate["checks"].values())
    assert json.loads(phase_flux_park_alignment_gate(json.dumps(angles), json.dumps(rows), 2))["status"] == "ok"


def test_phase_flux_park_gate_rejects_wrong_pole_pair_basis():
    angles, rows = _balanced_flux()
    gate = evaluate_phase_flux_park_alignment(angles, rows, 1)
    assert gate["status"] == "needs_attention"
    assert gate["checks"]["complete_electrical_cycle"] is False
    assert gate["checks"]["q_axis_near_zero"] is False


def test_phase_flux_park_gate_rejects_phase_error():
    angles, rows = _balanced_flux(phase_error=0.2)
    gate = evaluate_phase_flux_park_alignment(angles, rows, 2)
    assert gate["status"] == "needs_attention"
    assert gate["checks"]["q_axis_near_zero"] is False
