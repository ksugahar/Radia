from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_jmag_generalization_v21 import _payload_v21


def _payload_v22():
    payload = _payload_v21()
    identity = payload["artifact_identity"]
    identity["motion_skew_force_harmonic_time_angle_phase_generation_identity"] = {
        "motion_study_generation": "motion-41",
        "time_motion_study_generation": "motion-41",
        "angle_motion_study_generation": "motion-41",
        "skew_motion_study_generation": "motion-41",
        "phase_motion_study_generation": "motion-41",
        "force_result_motion_study_generation": "motion-41",
        "time_s": [0.0, 0.001, 0.002],
        "force_time_s": [0.0, 0.001, 0.002],
        "mechanical_angle_deg": [0.0, 5.0, 10.0],
        "force_mechanical_angle_deg": [0.0, 5.0, 10.0],
        "skew_slice_angles_deg": [-5.0, 0.0, 5.0],
        "force_skew_slice_angles_deg": [-5.0, 0.0, 5.0],
        "slice_weights": [0.25, 0.5, 0.25],
        "force_slice_weights": [0.25, 0.5, 0.25],
        "phase_reference_deg": 30.0,
        "force_phase_reference_deg": 30.0,
        "harmonic_orders": [1, 3, 5],
        "force_harmonic_orders": [1, 3, 5],
        "force_harmonics_n": [120.0, 4.5, 1.2],
        "reported_force_harmonics_n": [120.0, 4.5, 1.2],
        "force_harmonic_table_sha256": "1" * 64,
        "resolved_force_harmonic_table_sha256": "1" * 64,
    }
    identity[
        "ipm_irreversible_demag_recoil_temperature_operating_generation_identity"
    ] = {
        "demag_study_generation": "demag-41",
        "recoil_curve_demag_study_generation": "demag-41",
        "temperature_demag_study_generation": "demag-41",
        "operating_point_demag_study_generation": "demag-41",
        "magnet_orientation_demag_study_generation": "demag-41",
        "result_demag_study_generation": "demag-41",
        "temperature_c": 120.0,
        "result_temperature_c": 120.0,
        "operating_point_id": "id=-180A;iq=240A;theta=17.5deg",
        "result_operating_point_id": "id=-180A;iq=240A;theta=17.5deg",
        "magnet_orientation_vectors": [[1.0, 0.0], [0.0, 1.0]],
        "result_magnet_orientation_vectors": [[1.0, 0.0], [0.0, 1.0]],
        "recoil_curve_sha256": "2" * 64,
        "result_recoil_curve_sha256": "2" * 64,
        "magnet_state_sha256": "3" * 64,
        "result_magnet_state_sha256": "3" * 64,
        "demag_margin_a_per_m": [175000.0, 82000.0],
        "reported_demag_margin_a_per_m": [175000.0, 82000.0],
    }
    return payload


def test_v22_public_positive_motion_force_and_demag_identity():
    assert pwm_controlled_motor_loss_gate(_payload_v22())["status"] == "ok"


def test_v22_public_motion_skew_force_harmonic_time_angle_phase_generation_mismatch():
    payload = _payload_v22()
    identity = payload["artifact_identity"][
        "motion_skew_force_harmonic_time_angle_phase_generation_identity"
    ]
    identity.update(
        {
            "time_motion_study_generation": "motion-40",
            "angle_motion_study_generation": "motion-39",
            "skew_motion_study_generation": "motion-38",
            "phase_motion_study_generation": "motion-37",
            "force_time_s": [0.0, 0.002, 0.004],
            "force_mechanical_angle_deg": [0.0, 6.0, 12.0],
            "force_slice_weights": [0.5, 0.25, 0.25],
            "force_phase_reference_deg": -30.0,
            "force_harmonic_orders": [1, 5, 7],
            "reported_force_harmonics_n": [120.0, 1.2, 0.7],
            "resolved_force_harmonic_table_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "force_harmonics_use_current_motion_skew_time_angle_and_phase"
    ]


def test_v22_public_ipm_irreversible_demag_recoil_temperature_operating_generation_mismatch():
    payload = _payload_v22()
    identity = payload["artifact_identity"][
        "ipm_irreversible_demag_recoil_temperature_operating_generation_identity"
    ]
    identity.update(
        {
            "recoil_curve_demag_study_generation": "demag-40",
            "temperature_demag_study_generation": "demag-39",
            "operating_point_demag_study_generation": "demag-38",
            "magnet_orientation_demag_study_generation": "demag-37",
            "result_temperature_c": 80.0,
            "result_operating_point_id": "id=-100A;iq=120A;theta=5deg",
            "result_magnet_orientation_vectors": [[-1.0, 0.0], [0.0, 1.0]],
            "result_recoil_curve_sha256": "b" * 64,
            "result_magnet_state_sha256": "c" * 64,
            "reported_demag_margin_a_per_m": [310000.0, 205000.0],
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "irreversible_demag_uses_current_recoil_temperature_operating_state"
    ]
