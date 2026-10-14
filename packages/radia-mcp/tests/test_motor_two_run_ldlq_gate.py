import copy
import json
import math

from radia_mcp.motor.server import motor_ipm_two_run_ldlq_gate
from radia_mcp.motor.two_run_ldlq_gate import evaluate_ipm_two_run_ldlq


def _summary():
    angles = list(range(0, 181, 10))
    pm_flux = []
    on_flux = []
    currents = []
    for angle in angles:
        theta = math.radians(2 * angle)
        current = [
            math.cos(theta + 3 * math.pi / 4),
            math.cos(theta + 3 * math.pi / 4 - 2 * math.pi / 3),
            math.cos(theta + 3 * math.pi / 4 + 2 * math.pi / 3),
        ]
        pm = [
            0.03 * math.cos(theta),
            0.03 * math.cos(theta - 2 * math.pi / 3),
            0.03 * math.cos(theta + 2 * math.pi / 3),
        ]
        alpha = 0.01 * (-math.sqrt(0.5))
        beta = 0.015 * math.sqrt(0.5)
        delta = [
            alpha * math.cos(theta) - beta * math.sin(theta),
            alpha * math.cos(theta - 2 * math.pi / 3) - beta * math.sin(theta - 2 * math.pi / 3),
            alpha * math.cos(theta + 2 * math.pi / 3) - beta * math.sin(theta + 2 * math.pi / 3),
        ]
        currents.append(current)
        pm_flux.append(pm)
        on_flux.append([base + inc for base, inc in zip(pm, delta)])
    return {
        "pole_pairs": 2,
        "expected_saliency": "Lq_gt_Ld",
        "pm_only": {
            "phase_order": ["A", "B", "C"],
            "mechanical_angles_deg": angles.copy(),
            "phase_flux_wb": pm_flux,
        },
        "current_on": {
            "phase_order": ["A", "B", "C"],
            "mechanical_angles_deg": angles.copy(),
            "phase_currents_a": currents,
            "phase_flux_wb": on_flux,
        },
    }


def test_two_run_ldlq_gate_accepts_pm_subtracted_ipm_pair():
    result = evaluate_ipm_two_run_ldlq(_summary())
    assert result["status"] == "ok"
    assert math.isclose(result["ld_h"], 0.01, rel_tol=1.0e-12)
    assert math.isclose(result["lq_h"], 0.015, rel_tol=1.0e-12)
    assert json.loads(motor_ipm_two_run_ldlq_gate(json.dumps(_summary())))["status"] == "ok"


def test_two_run_ldlq_gate_rejects_total_flux_without_pm_subtraction():
    bad = copy.deepcopy(_summary())
    bad["current_on"]["phase_flux_wb"] = bad["pm_only"]["phase_flux_wb"]
    result = evaluate_ipm_two_run_ldlq(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["pm_subtraction_nonzero"] is False
    assert result["checks"]["positive_finite_ld_lq"] is False


def test_two_run_ldlq_gate_rejects_angle_or_phase_mismatch():
    bad = copy.deepcopy(_summary())
    bad["current_on"]["mechanical_angles_deg"][4] += 0.5
    bad["current_on"]["phase_order"] = ["A", "C", "B"]
    result = evaluate_ipm_two_run_ldlq(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["same_angle_grid"] is False
    assert result["checks"]["canonical_phase_order"] is False
