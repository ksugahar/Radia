from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v27 import _payload_v27


def _payload_v28():
    payload = _payload_v27()
    identity = payload["artifact_identity"]
    identity[
        "pwm_current_harmonic_time_electrical_angle_torque_loss_mesh_result_generation_identity"
    ] = {
        "pwm_generation": "pwm-observables-151",
        "current_pwm_generation": "pwm-observables-151",
        "time_pwm_generation": "pwm-observables-151",
        "angle_pwm_generation": "pwm-observables-151",
        "torque_pwm_generation": "pwm-observables-151",
        "loss_pwm_generation": "pwm-observables-151",
        "mesh_pwm_generation": "pwm-observables-151",
        "result_pwm_generation": "pwm-observables-151",
        "harmonic_orders": [1, 5, 7, 11],
        "result_harmonic_orders": [1, 5, 7, 11],
        "current_harmonic_a": [100.0, 8.0, 5.0, 2.0],
        "result_current_harmonic_a": [100.0, 8.0, 5.0, 2.0],
        "current_phase_deg": [0.0, -20.0, 15.0, -5.0],
        "result_current_phase_deg": [0.0, -20.0, 15.0, -5.0],
        "time_s": [0.020, 0.021, 0.022, 0.023],
        "result_time_s": [0.020, 0.021, 0.022, 0.023],
        "electrical_angle_deg": [0.0, 72.0, 144.0, 216.0],
        "result_electrical_angle_deg": [0.0, 72.0, 144.0, 216.0],
        "pole_pairs": 4,
        "result_pole_pairs": 4,
        "torque_window_s": [0.020, 0.023],
        "result_torque_window_s": [0.020, 0.023],
        "torque_average_nm": 42.0,
        "result_torque_average_nm": 42.0,
        "loss_average_w": 350.0,
        "result_loss_average_w": 350.0,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    identity[
        "skew_slice_angle_weight_frame_interpolation_torque_ripple_mesh_generation_identity"
    ] = {
        "skew_generation": "skew-average-151",
        "angle_skew_generation": "skew-average-151",
        "weight_skew_generation": "skew-average-151",
        "frame_skew_generation": "skew-average-151",
        "interpolation_skew_generation": "skew-average-151",
        "torque_skew_generation": "skew-average-151",
        "mesh_skew_generation": "skew-average-151",
        "result_skew_generation": "skew-average-151",
        "slice_angles_deg": [-3.0, 0.0, 3.0],
        "result_slice_angles_deg": [-3.0, 0.0, 3.0],
        "slice_weights": [0.25, 0.5, 0.25],
        "result_slice_weights": [0.25, 0.5, 0.25],
        "rotor_frame": "mechanical-ccw",
        "result_rotor_frame": "mechanical-ccw",
        "interpolation_rule": "periodic-cubic",
        "result_interpolation_rule": "periodic-cubic",
        "slice_torque_nm": [40.0, 44.0, 42.0],
        "result_slice_torque_nm": [40.0, 44.0, 42.0],
        "torque_average_nm": 42.5,
        "result_torque_average_nm": 42.5,
        "torque_ripple_nm": 4.0,
        "result_torque_ripple_nm": 4.0,
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v28_public_positive_pwm_observable_and_skew_average_identities():
    assert pwm_controlled_motor_loss_gate(_payload_v28())["status"] == "ok"


def test_v28_public_pwm_current_harmonic_time_alignment_electrical_angle_torque_loss_average_mismatch():
    payload = _payload_v28()
    identity = payload["artifact_identity"][
        "pwm_current_harmonic_time_electrical_angle_torque_loss_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "current_pwm_generation": "pwm-observables-150",
            "time_pwm_generation": "pwm-observables-149",
            "result_harmonic_orders": [1, 3, 5],
            "result_current_harmonic_a": [100.0, 20.0, 8.0],
            "result_current_phase_deg": [0.0, 30.0, -20.0],
            "result_time_s": [0.0, 0.001, 0.002],
            "result_electrical_angle_deg": [0.0, 60.0, 120.0],
            "result_pole_pairs": 2,
            "result_torque_window_s": [0.0, 0.002],
            "result_torque_average_nm": 38.0,
            "result_loss_average_w": 500.0,
            "result_mesh_sha256": "8" * 64,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pwm_observables_use_current_harmonics_time_angle_torque_loss_mesh_and_result"
    ]


def test_v28_public_skew_slice_angular_offset_periodic_weight_rotor_frame_torque_ripple_mismatch():
    payload = _payload_v28()
    identity = payload["artifact_identity"][
        "skew_slice_angle_weight_frame_interpolation_torque_ripple_mesh_generation_identity"
    ]
    identity.update(
        {
            "angle_skew_generation": "skew-average-150",
            "frame_skew_generation": "skew-average-149",
            "result_slice_angles_deg": [3.0, 0.0, -3.0],
            "result_slice_weights": [0.5, 0.5, 0.5],
            "result_rotor_frame": "electrical-clockwise",
            "result_interpolation_rule": "linear",
            "result_slice_torque_nm": [38.0, 45.0, 41.0],
            "result_torque_average_nm": 41.0,
            "result_torque_ripple_nm": 7.0,
            "result_mesh_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "skew_average_uses_current_angles_weights_frame_interpolation_torque_mesh_and_result"
    ]
