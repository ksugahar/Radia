from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import magnetic_force_method_profile_gate
from test_magnetic_force_generalization_v29 import _summary_v29

_PROMOTED_CASE_IDS = (
    "v30_public_maglev_force_stiffness_position_current_perturbation_derivative_mesh_mismatch",
    "v30_public_cogging_torque_slot_pole_period_angular_sampling_harmonic_phase_mismatch",
)

def _summary_v30():
    summary = _summary_v29(); identity = summary["artifact_identity"]
    generation = "maglev-stiffness-341"
    identity["maglev_force_stiffness_position_current_derivative_frame_mesh_result_identity"] = {
        "maglev_generation": generation, "position_maglev_generation": generation,
        "current_maglev_generation": generation, "force_maglev_generation": generation,
        "derivative_maglev_generation": generation, "frame_maglev_generation": generation,
        "mesh_maglev_generation": generation, "result_maglev_generation": generation,
        "position_m": [-0.001, 0.0, 0.001], "result_position_m": [-0.001, 0.0, 0.001],
        "current_a": [10.0, 10.0, 10.0], "result_current_a": [10.0, 10.0, 10.0],
        "force_z_n": [12.0, 10.0, 8.0], "result_force_z_n": [12.0, 10.0, 8.0],
        "stiffness_n_per_m": -2000.0, "result_stiffness_n_per_m": -2000.0,
        "coordinate_frame": "global_z_positive_up", "result_coordinate_frame": "global_z_positive_up",
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "cogging-341"
    identity["cogging_torque_slot_pole_period_origin_sampling_harmonic_phase_mesh_result_identity"] = {
        "cogging_generation": generation, "slot_cogging_generation": generation,
        "pole_cogging_generation": generation, "period_cogging_generation": generation,
        "sampling_cogging_generation": generation, "harmonic_cogging_generation": generation,
        "mesh_cogging_generation": generation, "result_cogging_generation": generation,
        "slot_count": 12, "result_slot_count": 12, "pole_count": 10, "result_pole_count": 10,
        "cogging_period_mechanical_deg": 6.0, "result_cogging_period_mechanical_deg": 6.0,
        "angular_origin_deg": 0.0, "result_angular_origin_deg": 0.0,
        "sample_angles_deg": [0.0, 1.5, 3.0, 4.5, 6.0],
        "result_sample_angles_deg": [0.0, 1.5, 3.0, 4.5, 6.0],
        "harmonic_orders": [1, 2], "result_harmonic_orders": [1, 2],
        "harmonic_phase_deg": [0.0, 90.0], "result_harmonic_phase_deg": [0.0, 90.0],
        "torque_nm": [0.0, 0.1, 0.0, -0.1, 0.0], "result_torque_nm": [0.0, 0.1, 0.0, -0.1, 0.0],
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary

def test_v30_public_positive_maglev_and_cogging_identities():
    assert magnetic_force_method_profile_gate(_summary_v30())["status"] == "ok"

def test_v30_public_maglev_force_stiffness_position_current_perturbation_derivative_mesh_mismatch():
    summary = _summary_v30()
    summary["artifact_identity"]["maglev_force_stiffness_position_current_derivative_frame_mesh_result_identity"].update({
        "position_maglev_generation": "maglev-stiffness-340", "mesh_maglev_generation": "maglev-stiffness-339",
        "result_position_m": [0.0, 0.001, 0.002], "result_current_a": [9.0, 10.0, 11.0],
        "result_force_z_n": [8.0, 10.0, 12.0], "result_stiffness_n_per_m": 2000.0,
        "result_coordinate_frame": "local_r_positive_down", "result_mesh_sha256": "7" * 64,
        "accepted_result_sha256": "8" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["maglev_stiffness_uses_fixed_current_symmetric_positions_force_derivative_frame_mesh_and_result"]

def test_v30_public_cogging_torque_slot_pole_period_angular_sampling_harmonic_phase_mismatch():
    summary = _summary_v30()
    summary["artifact_identity"]["cogging_torque_slot_pole_period_origin_sampling_harmonic_phase_mesh_result_identity"].update({
        "slot_cogging_generation": "cogging-340", "sampling_cogging_generation": "cogging-339",
        "result_slot_count": 9, "result_pole_count": 8, "result_cogging_period_mechanical_deg": 15.0,
        "result_angular_origin_deg": 2.0, "result_sample_angles_deg": [2.0, 4.0, 6.0],
        "result_harmonic_orders": [3, 7], "result_harmonic_phase_deg": [45.0, -20.0],
        "result_torque_nm": [0.2, 0.1, 0.3], "result_mesh_sha256": "9" * 64,
        "accepted_result_sha256": "a" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["cogging_torque_uses_slot_pole_period_origin_sampling_harmonics_phase_mesh_and_result"]
