import copy
import json

from radia_mcp.motor.server import motor_virtual_work_width_ladder_gate


def summary():
    direct = 0.136
    rows = [
        (0.5, 0.132, abs(0.132 - direct) / direct),
        (1.0, 0.1367, abs(0.1367 - direct) / direct),
        (2.0, 0.1356, abs(0.1356 - direct) / direct),
        (4.0, 0.1322, abs(0.1322 - direct) / direct),
    ]
    return {
        "excitation": {
            "phase_currents_A": {"A": 10.0, "B": -5.0, "C": -5.0},
            "expected_square_sum_A2": 150.0,
        },
        "virtual_work": [
            {
                "delta_deg": delta,
                "coenergy_derivative_torque_Nm": virtual,
                "weighted_stress_torque_Nm": direct,
                "relative_error": error,
            }
            for delta, virtual, error in rows
        ],
        "selected_virtual_work_delta_deg": 2.0,
        "mesh_element_count_relative_span": 0.008,
    }


def test_selects_interior_width_instead_of_smallest_remeshed_step():
    result = json.loads(motor_virtual_work_width_ladder_gate(json.dumps(summary())))
    assert result["status"] == "ok"
    assert result["selected_delta_deg"] == 2.0
    assert result["checks"]["smaller_width_exposes_remesh_noise"] is True
    assert result["checks"]["larger_width_exposes_truncation"] is True


def test_rejects_smallest_width_claim_and_current_imbalance():
    row = summary()
    row["selected_virtual_work_delta_deg"] = 0.5
    row["excitation"]["phase_currents_A"]["C"] = -4.0
    result = json.loads(motor_virtual_work_width_ladder_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["selected_width_is_error_minimum"] is False
    assert result["checks"]["balanced_three_phase_current_sum"] is False
    assert result["checks"]["balanced_three_phase_square_sum"] is False


def test_rejects_zero_direct_torque_even_if_differences_are_small():
    row = summary()
    for item in row["virtual_work"]:
        item["weighted_stress_torque_Nm"] = 0.0
        item["coenergy_derivative_torque_Nm"] = 0.0
        item["relative_error"] = 0.0
    result = json.loads(motor_virtual_work_width_ladder_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_direct_torque_is_nonzero_and_consistent"] is False
