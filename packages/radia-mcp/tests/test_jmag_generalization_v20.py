from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_jmag_generalization_v19 import _payload_v19


def _payload_v20():
    payload = _payload_v19()
    identity = payload["artifact_identity"]
    identity["skew_slice_torque_angle_weight_periodicity_generation_identity"] = {
        "skew_generation": "skew-22",
        "torque_skew_generation": "skew-22",
        "angle_skew_generation": "skew-22",
        "weight_skew_generation": "skew-22",
        "periodicity_skew_generation": "skew-22",
        "slice_ids": [1, 2, 3],
        "torque_slice_ids": [1, 2, 3],
        "slice_angles_deg": [-5.0, 0.0, 5.0],
        "torque_slice_angles_deg": [-5.0, 0.0, 5.0],
        "quadrature_weights": [0.25, 0.5, 0.25],
        "torque_quadrature_weights": [0.25, 0.5, 0.25],
        "periodic_wrap_deg": 360.0,
        "torque_periodic_wrap_deg": 360.0,
        "skew_average_table_sha256": "1" * 64,
        "torque_skew_average_table_sha256": "1" * 64,
    }
    identity[
        "incremental_inductance_current_perturbation_phase_state_generation_identity"
    ] = {
        "operating_point_generation": "operating-point-22",
        "matrix_operating_point_generation": "operating-point-22",
        "perturbation_operating_point_generation": "operating-point-22",
        "phase_state_operating_point_generation": "operating-point-22",
        "base_solve_generation": "solve-22-base",
        "matrix_base_solve_generation": "solve-22-base",
        "perturbation_solve_generations": ["solve-22-a", "solve-22-b", "solve-22-c"],
        "matrix_perturbation_solve_generations": [
            "solve-22-a",
            "solve-22-b",
            "solve-22-c",
        ],
        "phase_names": ["a", "b", "c"],
        "matrix_phase_names": ["a", "b", "c"],
        "perturbation_currents_a": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        "matrix_perturbation_currents_a": [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        "incremental_inductance_table_sha256": "2" * 64,
        "resolved_incremental_inductance_table_sha256": "2" * 64,
    }
    return payload


def test_v20_public_positive_skew_and_incremental_inductance_identity():
    result = pwm_controlled_motor_loss_gate(_payload_v20())
    assert result["status"] == "ok"
    assert result["checks"][
        "skew_torque_uses_current_slice_angles_weights_and_periodicity"
    ]
    assert result["checks"][
        "incremental_inductance_uses_current_perturbation_phase_and_state"
    ]


def test_v20_public_skew_slice_torque_angle_weight_periodicity_generation_mismatch():
    payload = _payload_v20()
    payload["artifact_identity"][
        "skew_slice_torque_angle_weight_periodicity_generation_identity"
    ].update(
        {
            "angle_skew_generation": "skew-21",
            "weight_skew_generation": "skew-21",
            "torque_slice_ids": [3, 2, 1],
            "torque_slice_angles_deg": [5.0, 0.0, -5.0],
            "torque_quadrature_weights": [0.5, 0.25, 0.25],
            "torque_periodic_wrap_deg": 180.0,
            "torque_skew_average_table_sha256": "f" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "skew_torque_uses_current_slice_angles_weights_and_periodicity"
    ] is False


def test_v20_public_incremental_inductance_current_perturbation_phase_state_generation_mismatch():
    payload = _payload_v20()
    payload["artifact_identity"][
        "incremental_inductance_current_perturbation_phase_state_generation_identity"
    ].update(
        {
            "perturbation_operating_point_generation": "operating-point-21",
            "phase_state_operating_point_generation": "operating-point-21",
            "matrix_base_solve_generation": "solve-21-base",
            "matrix_perturbation_solve_generations": [
                "solve-22-b",
                "solve-22-c",
                "solve-22-a",
            ],
            "matrix_phase_names": ["b", "c", "a"],
            "matrix_perturbation_currents_a": [
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ],
            "resolved_incremental_inductance_table_sha256": "f" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "incremental_inductance_uses_current_perturbation_phase_and_state"
    ] is False
