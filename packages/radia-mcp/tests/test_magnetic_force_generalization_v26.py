from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import magnetic_force_method_profile_gate
from test_magnetic_force_generalization_v25 import _summary_v25


def _summary_v26():
    summary = _summary_v25()
    identity = summary["artifact_identity"]
    identity["levitation_force_displacement_gradient_stiffness_energy_derivative_frame_generation_identity"] = {
        "levitation_generation": "levitation-301", "displacement_levitation_generation": "levitation-301",
        "force_levitation_generation": "levitation-301", "energy_levitation_generation": "levitation-301",
        "gradient_levitation_generation": "levitation-301", "frame_levitation_generation": "levitation-301",
        "result_levitation_generation": "levitation-301", "displacement_m": [-0.001, 0.0, 0.001],
        "result_displacement_m": [-0.001, 0.0, 0.001], "force_n": [10.0, 0.0, -10.0],
        "result_force_n": [10.0, 0.0, -10.0], "magnetic_energy_j": [0.005, 0.0, 0.005],
        "result_magnetic_energy_j": [0.005, 0.0, 0.005],
        "negative_energy_derivative_force_n": [10.0, 0.0, -10.0],
        "result_negative_energy_derivative_force_n": [10.0, 0.0, -10.0],
        "restoring_stiffness_n_m": 10000.0, "result_restoring_stiffness_n_m": 10000.0,
        "coordinate_frame": "global_z_up", "result_coordinate_frame": "global_z_up",
        "force_sign_convention": "restoring_negative_gradient",
        "result_force_sign_convention": "restoring_negative_gradient",
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    positions = [0.0, 7.5, 15.0, 22.5, 30.0, 37.5, 45.0]
    torque = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0]
    identity["cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_identity"] = {
        "cogging_generation": "cogging-301", "position_cogging_generation": "cogging-301",
        "periodicity_cogging_generation": "cogging-301", "mesh_cogging_generation": "cogging-301",
        "interpolation_cogging_generation": "cogging-301", "reference_cogging_generation": "cogging-301",
        "result_cogging_generation": "cogging-301", "mechanical_positions_deg": positions,
        "result_mechanical_positions_deg": positions, "cogging_torque_nm": torque,
        "result_cogging_torque_nm": torque, "periodicity": 8, "result_periodicity": 8,
        "mechanical_period_deg": 45.0, "result_mechanical_period_deg": 45.0,
        "reference_angle_deg": 0.0, "result_reference_angle_deg": 0.0,
        "interpolation_method": "periodic_cubic", "result_interpolation_method": "periodic_cubic",
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v26_public_positive_levitation_and_cogging_identity():
    assert magnetic_force_method_profile_gate(_summary_v26())["status"] == "ok"


def test_v26_public_levitation_force_displacement_gradient_stiffness_energy_derivative_frame_mismatch():
    summary = _summary_v26()
    summary["artifact_identity"]["levitation_force_displacement_gradient_stiffness_energy_derivative_frame_generation_identity"].update({
        "displacement_levitation_generation": "levitation-300", "result_force_n": [-10.0, 0.0, 10.0],
        "result_magnetic_energy_j": [0.0, 0.005, 0.0],
        "result_negative_energy_derivative_force_n": [-10.0, 0.0, 10.0],
        "result_restoring_stiffness_n_m": -10000.0, "result_coordinate_frame": "local_z_down",
        "result_force_sign_convention": "positive_gradient"})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["levitation_force_stiffness_and_energy_derivative_share_displacement_frame_sign_and_generation"]


def test_v26_public_cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_mismatch():
    summary = _summary_v26()
    summary["artifact_identity"]["cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_identity"].update({
        "position_cogging_generation": "cogging-300",
        "result_mechanical_positions_deg": [5.0, 12.5, 20.0, 27.5, 35.0, 42.5],
        "result_cogging_torque_nm": [0.5, 1.0, 0.0, -1.0, 0.0, -0.5],
        "result_periodicity": 6, "result_mechanical_period_deg": 60.0,
        "result_reference_angle_deg": 5.0, "result_interpolation_method": "linear_open",
        "result_mesh_sha256": "9" * 64})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["cogging_torque_uses_current_position_periodicity_mesh_interpolation_reference_and_result"]
