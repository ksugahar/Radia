from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v29 import _payload_v29


_PROMOTED_CASE_IDS = (
    "v30_public_ipm_dq_inductance_current_angle_frame_saturation_reciprocity_result_mismatch",
    "v30_public_srm_torque_current_position_coenergy_periodicity_phase_sequence_mismatch",
)


def _payload_v30():
    payload = _payload_v29()
    identity = payload["artifact_identity"]
    generation = "ipm-dq-171"
    identity["ipm_dq_inductance_current_angle_park_saturation_flux_derivative_reciprocity_mesh_result_identity"] = {
        "dq_generation": generation, "current_dq_generation": generation,
        "frame_dq_generation": generation, "saturation_dq_generation": generation,
        "flux_dq_generation": generation, "derivative_dq_generation": generation,
        "reciprocity_dq_generation": generation, "mesh_dq_generation": generation,
        "result_dq_generation": generation,
        "current_magnitude_a": 100.0, "result_current_magnitude_a": 100.0,
        "current_angle_electrical_deg": 30.0, "result_current_angle_electrical_deg": 30.0,
        "park_frame": "rotor_d_aligned_ccw_power_invariant",
        "result_park_frame": "rotor_d_aligned_ccw_power_invariant",
        "saturation_operating_point_a": [86.6025403784, 50.0],
        "result_saturation_operating_point_a": [86.6025403784, 50.0],
        "flux_linkage_derivative_h": [[0.003, 0.0002], [0.0002, 0.006]],
        "result_flux_linkage_derivative_h": [[0.003, 0.0002], [0.0002, 0.006]],
        "reciprocity_tolerance_h": 1.0e-9, "result_reciprocity_tolerance_h": 1.0e-9,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "srm-coenergy-171"
    torque = (2.1 - 1.9) / math.radians(2.0)
    identity["srm_torque_current_position_coenergy_periodicity_phase_sequence_mesh_result_identity"] = {
        "srm_generation": generation, "current_srm_generation": generation,
        "position_srm_generation": generation, "coenergy_srm_generation": generation,
        "periodicity_srm_generation": generation, "phase_srm_generation": generation,
        "mesh_srm_generation": generation, "result_srm_generation": generation,
        "current_a": [0.0, 25.0, 50.0], "result_current_a": [0.0, 25.0, 50.0],
        "rotor_position_mechanical_deg": [-1.0, 0.0, 1.0],
        "result_rotor_position_mechanical_deg": [-1.0, 0.0, 1.0],
        "coenergy_j_at_50a": [1.9, 2.0, 2.1], "result_coenergy_j_at_50a": [1.9, 2.0, 2.1],
        "torque_nm_at_50a": torque, "result_torque_nm_at_50a": torque,
        "sector_period_mechanical_deg": 30.0, "result_sector_period_mechanical_deg": 30.0,
        "phase_sequence": ["A", "B", "C"], "result_phase_sequence": ["A", "B", "C"],
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v30_public_positive_ipm_dq_and_srm_coenergy_identities():
    assert pwm_controlled_motor_loss_gate(_payload_v30())["status"] == "ok"


def test_v30_public_ipm_dq_inductance_current_angle_frame_saturation_reciprocity_result_mismatch():
    payload = _payload_v30()
    payload["artifact_identity"]["ipm_dq_inductance_current_angle_park_saturation_flux_derivative_reciprocity_mesh_result_identity"].update({
        "current_dq_generation": "ipm-dq-170", "mesh_dq_generation": "ipm-dq-169",
        "result_current_magnitude_a": 80.0, "result_current_angle_electrical_deg": -30.0,
        "result_park_frame": "stator_q_aligned_clockwise",
        "result_saturation_operating_point_a": [50.0, 86.6],
        "result_flux_linkage_derivative_h": [[0.004, 0.001], [-0.0005, 0.005]],
        "result_reciprocity_tolerance_h": 0.1,
        "result_mesh_sha256": "8" * 64, "accepted_result_sha256": "9" * 64,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["ipm_dq_inductance_uses_current_angle_park_frame_saturation_derivatives_reciprocity_mesh_and_result"]


def test_v30_public_srm_torque_current_position_coenergy_periodicity_phase_sequence_mismatch():
    payload = _payload_v30()
    payload["artifact_identity"]["srm_torque_current_position_coenergy_periodicity_phase_sequence_mesh_result_identity"].update({
        "current_srm_generation": "srm-coenergy-170", "periodicity_srm_generation": "srm-coenergy-169",
        "result_current_a": [0.0, 20.0, 40.0],
        "result_rotor_position_mechanical_deg": [0.0, 1.0, 2.0],
        "result_coenergy_j_at_50a": [1.9, 2.2, 2.1], "result_torque_nm_at_50a": -20.0,
        "result_sector_period_mechanical_deg": 45.0, "result_phase_sequence": ["A", "C", "B"],
        "result_mesh_sha256": "a" * 64, "accepted_result_sha256": "b" * 64,
    })
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["srm_torque_uses_current_positions_coenergy_periodicity_phase_sequence_mesh_and_result"]
