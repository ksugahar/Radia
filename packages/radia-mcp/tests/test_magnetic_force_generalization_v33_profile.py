from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v32_profile import _summary_v32


_PROMOTED_CASE_IDS = (
    "v33_public_halbach_harmonic_magnetization_order_pole_pitch_phase_field_energy_force_mismatch",
    "v33_public_magnetic_bearing_force_current_displacement_stiffness_matrix_reciprocity_stability_mismatch",
)


def _summary_v33():
    summary = _summary_v32()
    identity = summary["artifact_identity"]
    generation = "halbach-harmonic-371"
    identity[
        "halbach_harmonic_magnetization_order_pitch_phase_grid_field_energy_force_geometry_owner_result_identity"
    ] = {
        "halbach_generation": generation,
        **{
            key: generation
            for key in (
                "magnetization_generation",
                "pitch_generation",
                "phase_generation",
                "grid_generation",
                "field_generation",
                "energy_generation",
                "force_generation",
                "geometry_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "magnetization_angles_deg": [0.0, 90.0, 180.0, 270.0],
        "result_magnetization_angles_deg": [0.0, 90.0, 180.0, 270.0],
        "pole_pitch_m": 0.02,
        "result_pole_pitch_m": 0.02,
        "harmonic_orders": [1, 3, 5],
        "result_harmonic_orders": [1, 3, 5],
        "harmonic_phase_deg": [0.0, 30.0, -15.0],
        "result_harmonic_phase_deg": [0.0, 30.0, -15.0],
        "sampling_grid_m": [0.0, 0.005, 0.01, 0.015, 0.02],
        "result_sampling_grid_m": [0.0, 0.005, 0.01, 0.015, 0.02],
        "field_harmonic_amplitude_t": [1.0, 0.1, 0.03],
        "result_field_harmonic_amplitude_t": [1.0, 0.1, 0.03],
        "magnetic_energy_j": 0.5,
        "result_magnetic_energy_j": 0.5,
        "force_direction": "+x",
        "result_force_direction": "+x",
        "force_n": 10.0,
        "result_force_n": 10.0,
        "geometry_sha256": "1" * 64,
        "result_geometry_sha256": "1" * 64,
        "result_owner": "halbach/case-371/harmonic-line",
        "accepted_result_owner": "halbach/case-371/harmonic-line",
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "magnetic-bearing-371"
    stiffness = [[1000.0, 50.0], [50.0, 900.0]]
    eigenvalues = [879.2893218813452, 1020.7106781186548]
    identity[
        "magnetic_bearing_force_current_displacement_stiffness_reciprocity_bias_frame_mesh_result_identity"
    ] = {
        "bearing_generation": generation,
        **{
            key: generation
            for key in (
                "current_generation",
                "displacement_generation",
                "jacobian_generation",
                "stiffness_generation",
                "reciprocity_generation",
                "bias_generation",
                "frame_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "coordinate_frame": "global_xy_right_handed",
        "result_coordinate_frame": "global_xy_right_handed",
        "bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "result_bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "bias_displacement_m": [0.0, 0.0],
        "result_bias_displacement_m": [0.0, 0.0],
        "force_current_jacobian_n_per_a": [[10.0, -10.0, 0.0, 0.0], [0.0, 0.0, 10.0, -10.0]],
        "result_force_current_jacobian_n_per_a": [[10.0, -10.0, 0.0, 0.0], [0.0, 0.0, 10.0, -10.0]],
        "force_displacement_jacobian_n_per_m": [[-1000.0, -50.0], [-50.0, -900.0]],
        "result_force_displacement_jacobian_n_per_m": [[-1000.0, -50.0], [-50.0, -900.0]],
        "stiffness_matrix_n_per_m": stiffness,
        "result_stiffness_matrix_n_per_m": [list(row) for row in stiffness],
        "stiffness_eigenvalues_n_per_m": eigenvalues,
        "result_stiffness_eigenvalues_n_per_m": list(eigenvalues),
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v33_public_positive_halbach_and_magnetic_bearing_closure():
    assert magnetic_force_method_profile_gate(_summary_v33())["status"] == "ok"


def test_v33_public_halbach_harmonic_magnetization_order_pole_pitch_phase_field_energy_force_mismatch():
    summary = _summary_v33()
    record = summary["artifact_identity"][
        "halbach_harmonic_magnetization_order_pitch_phase_grid_field_energy_force_geometry_owner_result_identity"
    ]
    record.update(
        {
            "magnetization_generation": "halbach-harmonic-370",
            "field_generation": "halbach-harmonic-369",
            "result_generation": "halbach-harmonic-368",
            "result_magnetization_angles_deg": [270.0, 180.0, 90.0, 0.0],
            "result_pole_pitch_m": 0.04,
            "result_harmonic_orders": [1, 2, 4],
            "result_harmonic_phase_deg": [180.0, 0.0, 0.0],
            "result_sampling_grid_m": [0.0, 0.01, 0.02],
            "result_field_harmonic_amplitude_t": [0.2, 0.4, 0.8],
            "result_magnetic_energy_j": -0.5,
            "result_force_direction": "-z",
            "result_force_n": -10.0,
            "result_geometry_sha256": "7" * 64,
            "accepted_result_owner": "halbach/old-case",
            "accepted_result_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "halbach_harmonics_use_current_magnetization_pitch_phase_grid_field_energy_force_geometry_owner_and_result"
    ]


def test_v33_public_magnetic_bearing_force_current_displacement_stiffness_matrix_reciprocity_stability_mismatch():
    summary = _summary_v33()
    record = summary["artifact_identity"][
        "magnetic_bearing_force_current_displacement_stiffness_reciprocity_bias_frame_mesh_result_identity"
    ]
    record.update(
        {
            "jacobian_generation": "magnetic-bearing-370",
            "stiffness_generation": "magnetic-bearing-369",
            "result_generation": "magnetic-bearing-368",
            "result_coordinate_frame": "local_yz_left_handed",
            "result_bias_currents_a": [2.0, 2.0],
            "result_bias_displacement_m": [0.001, -0.001],
            "result_force_current_jacobian_n_per_a": [[1.0, 2.0]],
            "result_force_displacement_jacobian_n_per_m": [[100.0, 200.0], [-50.0, 100.0]],
            "result_stiffness_matrix_n_per_m": [[-100.0, -200.0], [50.0, -100.0]],
            "result_stiffness_eigenvalues_n_per_m": [-100.0, -100.0],
            "result_mesh_sha256": "9" * 64,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_bearing_uses_current_force_jacobians_bias_frame_reciprocal_positive_stiffness_mesh_and_result"
    ]
