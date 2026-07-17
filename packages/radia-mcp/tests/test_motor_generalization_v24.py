from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v23 import _payload_v23


def _payload_v24():
    payload = _payload_v23()
    identity = payload["artifact_identity"]
    identity[
        "loss_torque_speed_power_balance_harmonic_window_generation_identity"
    ] = {
        "power_balance_generation": "power-101",
        "torque_speed_power_balance_generation": "power-101",
        "harmonic_window_power_balance_generation": "power-101",
        "time_average_power_balance_generation": "power-101",
        "iron_loss_power_balance_generation": "power-101",
        "copper_loss_power_balance_generation": "power-101",
        "mechanical_loss_power_balance_generation": "power-101",
        "result_power_balance_generation": "power-101",
        "torque_nm": [1.0, 1.2],
        "power_balance_torque_nm": [1.0, 1.2],
        "speed_rad_s": [100.0, 100.0],
        "power_balance_speed_rad_s": [100.0, 100.0],
        "mechanical_output_w": [100.0, 120.0],
        "power_balance_mechanical_output_w": [100.0, 120.0],
        "harmonic_window_samples": [20, 120],
        "loss_harmonic_window_samples": [20, 120],
        "time_average_window_s": [0.02, 0.12],
        "loss_time_average_window_s": [0.02, 0.12],
        "iron_loss_w": [5.0, 6.0],
        "power_balance_iron_loss_w": [5.0, 6.0],
        "copper_loss_w": [3.0, 4.0],
        "power_balance_copper_loss_w": [3.0, 4.0],
        "mechanical_loss_w": [2.0, 2.0],
        "power_balance_mechanical_loss_w": [2.0, 2.0],
        "electrical_input_w": [110.0, 132.0],
        "power_balance_electrical_input_w": [110.0, 132.0],
        "power_balance_sha256": "1" * 64,
        "reported_power_balance_sha256": "1" * 64,
    }
    identity[
        "skew_slice_weight_rotor_angle_phase_periodicity_generation_identity"
    ] = {
        "skew_generation": "skew-101",
        "weight_skew_generation": "skew-101",
        "angle_skew_generation": "skew-101",
        "phase_skew_generation": "skew-101",
        "periodicity_skew_generation": "skew-101",
        "solve_skew_generation": "skew-101",
        "result_skew_generation": "skew-101",
        "slice_ids": [1, 2, 3],
        "result_slice_ids": [1, 2, 3],
        "quadrature_weights": [0.25, 0.5, 0.25],
        "result_quadrature_weights": [0.25, 0.5, 0.25],
        "rotor_angles_deg": [-5.0, 0.0, 5.0],
        "result_rotor_angles_deg": [-5.0, 0.0, 5.0],
        "current_phase_ids": ["abc@-5", "abc@0", "abc@5"],
        "result_current_phase_ids": ["abc@-5", "abc@0", "abc@5"],
        "periodic_map_ids": ["p1", "p2", "p3"],
        "result_periodic_map_ids": ["p1", "p2", "p3"],
        "slice_solve_sha256": ["2" * 64, "3" * 64, "4" * 64],
        "result_slice_solve_sha256": ["2" * 64, "3" * 64, "4" * 64],
        "slice_torque_nm": [0.8, 1.0, 1.2],
        "result_slice_torque_nm": [0.8, 1.0, 1.2],
        "weighted_torque_nm": 1.0,
        "reported_weighted_torque_nm": 1.0,
        "skew_result_sha256": "5" * 64,
        "reported_skew_result_sha256": "5" * 64,
    }
    return payload


def test_v24_public_positive_power_balance_and_skew_slice_identity():
    assert pwm_controlled_motor_loss_gate(_payload_v24())["status"] == "ok"


def test_v24_public_loss_torque_speed_power_balance_harmonic_window_generation_mismatch():
    payload = _payload_v24()
    payload["artifact_identity"][
        "loss_torque_speed_power_balance_harmonic_window_generation_identity"
    ].update(
        {
            "torque_speed_power_balance_generation": "power-100",
            "harmonic_window_power_balance_generation": "power-99",
            "time_average_power_balance_generation": "power-98",
            "iron_loss_power_balance_generation": "power-97",
            "copper_loss_power_balance_generation": "power-96",
            "mechanical_loss_power_balance_generation": "power-95",
            "power_balance_torque_nm": [1.2, 1.0],
            "power_balance_speed_rad_s": [50.0, 200.0],
            "power_balance_mechanical_output_w": [60.0, 200.0],
            "loss_harmonic_window_samples": [0, 40],
            "loss_time_average_window_s": [0.0, 0.04],
            "power_balance_iron_loss_w": [7.0, 8.0],
            "power_balance_copper_loss_w": [4.0, 5.0],
            "power_balance_mechanical_loss_w": [0.0, 0.0],
            "power_balance_electrical_input_w": [71.0, 213.0],
            "reported_power_balance_sha256": "c" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "loss_power_balance_uses_current_torque_speed_windows_and_loss_components"
    ]


def test_v24_public_skew_slice_weight_rotor_angle_phase_periodicity_generation_mismatch():
    payload = _payload_v24()
    payload["artifact_identity"][
        "skew_slice_weight_rotor_angle_phase_periodicity_generation_identity"
    ].update(
        {
            "weight_skew_generation": "skew-100",
            "angle_skew_generation": "skew-99",
            "phase_skew_generation": "skew-98",
            "periodicity_skew_generation": "skew-97",
            "solve_skew_generation": "skew-96",
            "result_quadrature_weights": [0.5, 0.5, 0.5],
            "result_rotor_angles_deg": [5.0, 0.0, -5.0],
            "result_current_phase_ids": ["acb@-5", "abc@0", "abc@5"],
            "result_periodic_map_ids": ["p3", "p2", "p1"],
            "result_slice_solve_sha256": ["d" * 64, "3" * 64, "4" * 64],
            "result_slice_torque_nm": [1.2, 1.0, 0.8],
            "reported_weighted_torque_nm": 1.5,
            "reported_skew_result_sha256": "e" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "skew_slice_result_uses_current_weights_angles_phases_and_periodicity"
    ]
