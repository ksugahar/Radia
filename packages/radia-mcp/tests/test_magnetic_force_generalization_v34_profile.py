from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v33_profile import _summary_v33


_PROMOTED_CASE_IDS = (
    "v34_public_magnetic_bearing_force_matrix_cross_coupled_stiffness_damping_stability_mismatch",
    "v34_public_moving_conductor_eddy_drag_lift_power_velocity_skin_depth_sign_mismatch",
)


def _summary_v34():
    summary = _summary_v33()
    identity = summary["artifact_identity"]

    generation = "bearing-dynamic-381"
    identity[
        "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
    ] = {
        "bearing_dynamic_generation": generation,
        **{
            key: generation
            for key in (
                "force_generation",
                "stiffness_generation",
                "damping_generation",
                "coordinate_generation",
                "reciprocity_generation",
                "stability_generation",
                "operating_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "coordinate_order": ["x", "y"],
        "result_coordinate_order": ["x", "y"],
        "displacement_perturbations_m": [[1.0e-4, 0.0], [-1.0e-4, 0.0], [0.0, 1.0e-4], [0.0, -1.0e-4]],
        "result_displacement_perturbations_m": [[1.0e-4, 0.0], [-1.0e-4, 0.0], [0.0, 1.0e-4], [0.0, -1.0e-4]],
        "force_perturbations_n": [[-0.1, -0.005], [0.1, 0.005], [-0.005, -0.09], [0.005, 0.09]],
        "result_force_perturbations_n": [[-0.1, -0.005], [0.1, 0.005], [-0.005, -0.09], [0.005, 0.09]],
        "stiffness_matrix_n_per_m": [[1000.0, 50.0], [50.0, 900.0]],
        "result_stiffness_matrix_n_per_m": [[1000.0, 50.0], [50.0, 900.0]],
        "damping_matrix_n_s_per_m": [[10.0, 2.0], [-2.0, 12.0]],
        "result_damping_matrix_n_s_per_m": [[10.0, 2.0], [-2.0, 12.0]],
        "state_eigenvalues_per_s": [[-5.0, 30.0], [-5.0, -30.0], [-6.0, 28.0], [-6.0, -28.0]],
        "result_state_eigenvalues_per_s": [[-5.0, 30.0], [-5.0, -30.0], [-6.0, 28.0], [-6.0, -28.0]],
        "operating_displacement_m": [0.0, 0.0],
        "result_operating_displacement_m": [0.0, 0.0],
        "operating_velocity_m_s": [0.0, 0.0],
        "result_operating_velocity_m_s": [0.0, 0.0],
        "bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "result_bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "bearing_mesh_sha256": "1" * 64,
        "result_bearing_mesh_sha256": "1" * 64,
        "bearing_result_owner": "bearing/case-381/linearization",
        "accepted_bearing_result_owner": "bearing/case-381/linearization",
        "bearing_result_sha256": "2" * 64,
        "accepted_bearing_result_sha256": "2" * 64,
    }

    generation = "moving-conductor-381"
    conductivity = 3.5e7
    frequency = 50.0
    skin_depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * 4.0e-7 * math.pi * conductivity))
    identity[
        "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
    ] = {
        "moving_conductor_generation": generation,
        **{
            key: generation
            for key in (
                "velocity_generation",
                "frame_generation",
                "force_generation",
                "power_generation",
                "skin_generation",
                "frequency_generation",
                "slip_generation",
                "mesh_generation",
                "owner_generation",
                "field_generation",
                "result_generation",
            )
        },
        "coordinate_frame": "global_xyz_right_handed",
        "result_coordinate_frame": "global_xyz_right_handed",
        "velocity_m_s": [10.0, 0.0, 0.0],
        "result_velocity_m_s": [10.0, 0.0, 0.0],
        "drag_force_n": [-100.0, 0.0, 0.0],
        "result_drag_force_n": [-100.0, 0.0, 0.0],
        "lift_force_n": [0.0, 20.0, 0.0],
        "result_lift_force_n": [0.0, 20.0, 0.0],
        "joule_power_w": 1000.0,
        "result_joule_power_w": 1000.0,
        "mechanical_drag_power_w": 1000.0,
        "result_mechanical_drag_power_w": 1000.0,
        "conductivity_s_m": conductivity,
        "result_conductivity_s_m": conductivity,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "excitation_frequency_hz": frequency,
        "result_excitation_frequency_hz": frequency,
        "spatial_period_m": 0.2,
        "result_spatial_period_m": 0.2,
        "slip_frequency_hz": frequency,
        "result_slip_frequency_hz": frequency,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "conductor_mesh_sha256": "3" * 64,
        "result_conductor_mesh_sha256": "3" * 64,
        "field_sha256": "4" * 64,
        "accepted_field_sha256": "4" * 64,
        "conductor_result_owner": "moving-conductor/case-381",
        "accepted_conductor_result_owner": "moving-conductor/case-381",
        "conductor_result_sha256": "5" * 64,
        "accepted_conductor_result_sha256": "5" * 64,
    }
    return summary


def test_v34_public_positive_bearing_dynamics_and_moving_conductor_closure():
    assert magnetic_force_method_profile_gate(_summary_v34())["status"] == "ok"


def test_v34_public_magnetic_bearing_force_matrix_cross_coupled_stiffness_damping_stability_mismatch():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
    ]
    record.update(
        {
            "force_generation": "bearing-dynamic-380",
            "stability_generation": "bearing-dynamic-379",
            "result_generation": "bearing-dynamic-378",
            "result_coordinate_order": ["y", "x"],
            "result_displacement_perturbations_m": [[0.0, 1.0e-3]],
            "result_force_perturbations_n": [[1.0, 1.0]],
            "result_stiffness_matrix_n_per_m": [[-1000.0, 500.0], [-50.0, -900.0]],
            "result_damping_matrix_n_s_per_m": [[-10.0, 20.0], [2.0, -12.0]],
            "result_state_eigenvalues_per_s": [[5.0, 30.0], [6.0, -28.0]],
            "result_operating_displacement_m": [0.001, -0.001],
            "result_operating_velocity_m_s": [1.0, 0.0],
            "result_bias_currents_a": [1.0, 2.0],
            "result_bearing_mesh_sha256": "a" * 64,
            "accepted_bearing_result_owner": "bearing/old",
            "accepted_bearing_result_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_bearing_dynamics_use_current_force_perturbations_stiffness_damping_coordinates_stability_operating_point_owner_and_result"
    ]


def test_v34_public_moving_conductor_eddy_drag_lift_power_velocity_skin_depth_sign_mismatch():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
    ]
    record.update(
        {
            "velocity_generation": "moving-conductor-380",
            "power_generation": "moving-conductor-379",
            "result_generation": "moving-conductor-378",
            "result_coordinate_frame": "body_left_handed",
            "result_velocity_m_s": [-10.0, 0.0, 0.0],
            "result_drag_force_n": [100.0, 0.0, 0.0],
            "result_lift_force_n": [100.0, 20.0, 0.0],
            "result_joule_power_w": -1000.0,
            "result_mechanical_drag_power_w": 500.0,
            "result_conductivity_s_m": -3.5e7,
            "result_relative_permeability": -1.0,
            "result_excitation_frequency_hz": 25.0,
            "result_spatial_period_m": 0.1,
            "result_slip_frequency_hz": -50.0,
            "result_skin_depth_m": -0.01,
            "result_conductor_mesh_sha256": "c" * 64,
            "accepted_field_sha256": "d" * 64,
            "accepted_conductor_result_owner": "moving-conductor/old",
            "accepted_conductor_result_sha256": "e" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "moving_conductor_uses_current_velocity_frame_drag_lift_power_skin_depth_slip_mesh_owner_field_and_result"
    ]


def test_v34_public_rejects_self_consistent_negative_bearing_damping():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
    ]
    damping = [[-10.0, 2.0], [-2.0, -12.0]]
    record["damping_matrix_n_s_per_m"] = damping
    record["result_damping_matrix_n_s_per_m"] = damping
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_drag_aligned_with_velocity():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
    ]
    record["drag_force_n"] = [100.0, 0.0, 0.0]
    record["result_drag_force_n"] = [100.0, 0.0, 0.0]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
