from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from radia_mcp.radia_ngsolve.server import (
    magnetic_force_method_profile_gate as mcp_magnetic_force_method_profile_gate,
)


def _summary() -> dict:
    return {
        "quantity_dimension": "3d_total",
        "force_unit": "N",
        "position_unit": "m",
        "comparison_axis": "z",
        "positions": [0.0, 0.001, 0.002, 0.003, 0.004, 0.005],
        "moving_body_element_force": [10.0, 11.0, 12.5, 14.0, 16.0, 18.0],
        "closed_surface_maxwell_stress_force": [10.2, 11.2, 12.7, 14.2, 16.2, 18.2],
        "independent_closed_surface_force": [10.1, 11.1, 12.6, 14.1, 16.1, 18.1],
        "all_body_element_force": [3.0, 3.1, 3.2, 3.4, 3.6, 3.8],
        "replay": {
            "parsed_max_abs": 0.0,
            "binary_nonlog_outputs_exact": True,
        },
    }


def _with_artifact_identity(summary: dict) -> dict:
    summary["artifact_identity"] = {
        "position_force_sample_generations": ["sweep-42"]
        * len(summary["positions"]),
        "sample_acquired_at_utc": [
            f"2026-07-16T05:00:0{index}Z"
            for index in range(len(summary["positions"]))
        ],
        "magnet_geometry": {
            "revision": "geometry-12",
            "committed_at_utc": "2026-07-16T04:59:00Z",
        },
        "demag_reference": {
            "geometry_revision": "geometry-12",
            "generated_at_utc": "2026-07-16T04:59:30Z",
        },
        "coordinate_system_binding": {
            "force_component_frame_id": "global-cartesian",
            "demag_metric_frame_id": "global-cartesian",
            "common_frame_id": "global-cartesian",
            "geometry_rotation_revision": "rotation-12",
            "force_transform_revision": "rotation-12",
            "demag_transform_revision": "rotation-12",
        },
        "force_normalization": {
            "comparison_basis": "total_3d",
            "profile_bases": {
                "moving_body_element_force": "total_3d",
                "closed_surface_maxwell_stress_force": "total_3d",
                "independent_closed_surface_force": "total_3d",
            },
            "per_length_to_total_applied": True,
        },
        "hysteresis_branch_state": {
            "observable_branch": "descending",
            "state_memory_branch": "descending",
            "tangent_branch": "descending",
            "branch_state_generation": "branch-state-12",
            "tangent_state_generation": "branch-state-12",
        },
        "remanence_frame_binding": {
            "remanence_vector_frame_id": "magnet-material-frame",
            "transform_input_frame_id": "magnet-material-frame",
            "assembly_frame_id": "assembly-global-frame",
            "transform_output_frame_id": "assembly-global-frame",
            "transformed_vector_frame_id": "assembly-global-frame",
            "transform_applied": True,
            "geometry_rotation_revision": "rotation-12",
            "remanence_transform_revision": "rotation-12",
        },
        "force_surface_body_ownership": {
            "target_body_id": "moving-magnet",
            "enclosed_body_ids": ["moving-magnet"],
            "surface_selection_generation": "force-surface-13",
            "force_integration_generation": "force-surface-13",
            "compensating_body_force_allowed": False,
        },
        "demag_branch_interpolation": {
            "operating_point_branch": "descending",
            "interpolation_source_branch": "descending",
            "branch_state_generation": "branch-state-13",
            "interpolation_state_generation": "branch-state-13",
            "bracketing_sample_ids": ["desc-17", "desc-18"],
        },
        "linear_motor_thrust_phase_identity": {
            "winding_phase_sequence": ["U", "V", "W"],
            "thrust_phase_sequence": ["U", "V", "W"],
            "winding_electrical_angle_direction": 1,
            "thrust_electrical_angle_direction": 1,
            "phase_convention_generation": "phase-convention-14",
            "thrust_observable_phase_generation": "phase-convention-14",
        },
        "demag_recoil_temperature_identity": {
            "evaluation_temperature_c": 120.0,
            "magnet_material_temperature_c": 120.0,
            "recoil_line_temperature_c": 120.0,
            "material_state_generation": "magnet-material-14",
            "recoil_line_state_generation": "magnet-material-14",
            "recoil_line_sha256": "4" * 64,
        },
        "bem_demag_surface_normal_generation_identity": {
            "active_surface_mesh_generation": "surface-mesh-15",
            "surface_element_generation": "surface-mesh-15",
            "surface_normal_generation": "surface-mesh-15",
            "demag_evaluation_surface_generation": "surface-mesh-15",
            "normal_orientation": "outward",
            "demag_kernel_normal_orientation": "outward",
            "surface_normal_sha256": "6" * 64,
            "demag_kernel_normal_sha256": "6" * 64,
        },
        "cogging_torque_periodic_sector_symmetry_identity": {
            "active_periodic_sector_count": 12,
            "torque_result_periodic_sector_count": 12,
            "symmetry_multiplier": 12.0,
            "sector_torque_scope": "one_periodic_sector",
            "reported_torque_scope": "full_machine",
            "periodic_topology_generation": "periodic-topology-15",
            "torque_result_topology_generation": "periodic-topology-15",
            "multiplier_topology_generation": "periodic-topology-15",
        },
        "bem_self_term_solid_angle_orientation_identity": {
            "active_surface_mesh_generation": "surface-mesh-16",
            "panel_generation": "surface-mesh-16",
            "panel_orientation_generation": "surface-mesh-16",
            "self_term_orientation_generation": "surface-mesh-16",
            "panel_orientation": "outward",
            "self_term_solid_angle_orientation": "outward",
            "solid_angle_sign_convention": "outward_positive",
            "self_term_sign_convention": "outward_positive",
            "panel_orientation_sha256": "a" * 64,
            "self_term_orientation_sha256": "a" * 64,
        },
        "demag_energy_force_displacement_length_unit_identity": {
            "energy_generation": "demag-energy-16",
            "displacement_generation": "displacement-grid-16",
            "force_derivative_displacement_generation": "displacement-grid-16",
            "energy_unit": "J",
            "force_unit": "N",
            "displacement_length_unit": "m",
            "displacement_scale_to_m": 1.0,
            "force_derivative_length_unit": "m",
            "force_derivative_scale_to_m": 1.0,
            "displacement_grid_sha256": "b" * 64,
            "force_derivative_grid_sha256": "b" * 64,
            "force_from_energy_convention": "negative_energy_gradient",
        },
        "bem_near_singular_quadrature_target_scale_identity": {
            "active_surface_mesh_generation": "surface-mesh-17",
            "source_panel_mesh_generation": "surface-mesh-17",
            "target_panel_mesh_generation": "surface-mesh-17",
            "target_separation_generation": "target-separation-17",
            "quadrature_target_separation_generation": "target-separation-17",
            "source_panel_id": 101,
            "target_panel_id": 205,
            "target_panel_separation_m": 1.0e-4,
            "target_panel_characteristic_length_m": 1.0e-2,
            "normalized_target_separation": 1.0e-2,
            "quadrature_normalized_target_separation": 1.0e-2,
            "quadrature_rule": "adaptive_duffy",
            "quadrature_order": 12,
            "target_pair_sha256": "1" * 64,
            "quadrature_target_pair_sha256": "1" * 64,
        },
        "force_torque_reference_origin_length_unit_identity": {
            "solve_generation": "force-solve-17",
            "force_result_generation": "force-solve-17",
            "torque_result_generation": "force-solve-17",
            "force_frame_id": "global-cartesian",
            "torque_frame_id": "global-cartesian",
            "force_unit": "N",
            "torque_unit": "N*m",
            "reference_origin_length_unit": "m",
            "force_reference_origin_length_unit": "m",
            "torque_reference_origin_length_unit": "m",
            "reference_origin_scale_to_m": 1.0,
            "force_reference_origin_scale_to_m": 1.0,
            "torque_reference_origin_scale_to_m": 1.0,
            "reference_origin_coordinates": [0.01, 0.0, -0.02],
            "force_reference_origin_coordinates": [0.01, 0.0, -0.02],
            "torque_reference_origin_coordinates": [0.01, 0.0, -0.02],
            "reference_origin_sha256": "2" * 64,
            "torque_reference_origin_sha256": "2" * 64,
        },
        "bem_solid_angle_surface_winding_identity": {
            "surface_mesh_generation": "surface-mesh-18",
            "normalized_surface_winding_generation": "surface-mesh-18",
            "solid_angle_sign_generation": "surface-mesh-18",
            "self_term_assembly_generation": "surface-mesh-18",
            "surface_winding_normalized": True,
            "solid_angle_sign_convention": "outward_positive",
            "self_term_sign_convention": "outward_positive",
            "surface_component_ids": [1, 2],
            "solid_angle_component_ids": [1, 2],
            "surface_winding_sha256": "5" * 64,
            "solid_angle_winding_sha256": "5" * 64,
        },
        "maglev_stiffness_force_displacement_generation_identity": {
            "perturbation_generation": "perturbation-18",
            "displacement_coordinate_generation": "perturbation-18",
            "force_sample_generation": "perturbation-18",
            "stiffness_derivative_generation": "perturbation-18",
            "displacement_axis": "global-z",
            "force_component_axis": "global-z",
            "displacement_unit": "m",
            "force_unit": "N",
            "stiffness_unit": "N/m",
            "displacement_sha256": "6" * 64,
            "stiffness_displacement_sha256": "6" * 64,
            "force_sample_sha256": "7" * 64,
            "stiffness_force_sample_sha256": "7" * 64,
        },
        "bem_demag_tensor_coordinate_basis_generation_identity": {
            "body_placement_generation": "placement-19",
            "surface_mesh_body_placement_generation": "placement-19",
            "demag_tensor_generation": "demag-tensor-19",
            "force_demag_tensor_generation": "demag-tensor-19",
            "demag_tensor_body_placement_generation": "placement-19",
            "body_coordinate_basis": "body-local-current",
            "tensor_coordinate_basis": "body-local-current",
            "body_basis_handedness": "right_handed",
            "tensor_basis_handedness": "right_handed",
            "body_to_global_transform_determinant": 1.0,
            "body_to_global_transform_orthogonality_error": 0.0,
            "body_to_global_transform_sha256": "1" * 64,
            "tensor_basis_transform_sha256": "1" * 64,
        },
        "magnetic_bearing_force_harmonic_phase_origin_identity": {
            "force_harmonic_generation": "bearing-harmonic-19",
            "force_sample_harmonic_generation": "bearing-harmonic-19",
            "rotor_angle_generation": "rotor-angle-19",
            "force_sample_rotor_angle_generation": "rotor-angle-19",
            "rotor_phase_origin_deg": 0.0,
            "force_harmonic_phase_origin_deg": 0.0,
            "slot_pitch_deg": 15.0,
            "force_harmonic_order": 2,
            "force_sample_harmonic_order": 2,
            "rotor_angle_basis": "mechanical_deg",
            "force_harmonic_angle_basis": "mechanical_deg",
            "fourier_sign_convention": "exp(-j*n*theta)",
            "force_harmonic_fourier_sign_convention": "exp(-j*n*theta)",
            "phase_origin_convention": "rotor_d_axis",
            "force_harmonic_phase_origin_convention": "rotor_d_axis",
            "rotor_angle_sha256": "2" * 64,
            "harmonic_input_angle_sha256": "2" * 64,
        },
        "bem_near_singular_panel_subdivision_quadrature_generation_identity": {
            "surface_generation": "surface-20",
            "near_singular_interaction_surface_generation": "surface-20",
            "panel_subdivision_surface_generation": "surface-20",
            "subdivision_generation": "subdivision-20",
            "quadrature_subdivision_generation": "subdivision-20",
            "quadrature_generation": "quadrature-20",
            "interaction_quadrature_generation": "quadrature-20",
            "interaction_ids": [101, 102, 103],
            "subdivided_interaction_ids": [101, 102, 103],
            "quadrature_orders": [12, 16, 20],
            "applied_quadrature_orders": [12, 16, 20],
            "subdivision_map_sha256": "5" * 64,
            "quadrature_input_subdivision_map_sha256": "5" * 64,
        },
        "maglev_force_coil_polarity_orientation_generation_identity": {
            "force_generation": "maglev-force-20",
            "coil_force_generation": "maglev-force-20",
            "coil_generation": "coil-20",
            "current_polarity_coil_generation": "coil-20",
            "winding_orientation_coil_generation": "coil-20",
            "force_result_coil_generation": "coil-20",
            "coil_ids": [1, 2],
            "force_coil_ids": [1, 2],
            "current_polarities": [1, -1],
            "force_current_polarities": [1, -1],
            "winding_orientations": ["clockwise", "counterclockwise"],
            "force_winding_orientations": ["clockwise", "counterclockwise"],
            "coil_orientation_map_sha256": "6" * 64,
            "force_coil_orientation_map_sha256": "6" * 64,
        },
        "bem_demag_self_term_solid_angle_orientation_generation_identity": {
            "operator_generation": "bem-operator-21",
            "result_operator_generation": "bem-operator-21",
            "boundary_generation": "boundary-21",
            "self_term_boundary_generation": "boundary-21",
            "panel_orientation_boundary_generation": "boundary-21",
            "operator_boundary_generation": "boundary-21",
            "panel_ids": [101, 102, 103],
            "self_term_panel_ids": [101, 102, 103],
            "panel_orientation_signs": [1, 1, 1],
            "self_term_orientation_signs": [1, 1, 1],
            "solid_angles_sr": [6.283185307179586] * 3,
            "applied_self_term_solid_angles_sr": [6.283185307179586] * 3,
            "orientation_convention": "outward_positive",
            "self_term_orientation_convention": "outward_positive",
            "self_term_input_sha256": "a" * 64,
            "assembled_self_term_input_sha256": "a" * 64,
        },
        "virtual_work_force_displacement_coordinate_geometry_generation_identity": {
            "force_generation": "virtual-work-force-21",
            "result_force_generation": "virtual-work-force-21",
            "geometry_generation": "geometry-21",
            "displacement_coordinate_geometry_generation": "geometry-21",
            "energy_sample_geometry_generations": ["geometry-21"] * 3,
            "force_geometry_generation": "geometry-21",
            "displacement_coordinate_generation": "displacement-21",
            "energy_sample_displacement_coordinate_generation": "displacement-21",
            "force_displacement_coordinate_generation": "displacement-21",
            "displacement_axis": "global-z",
            "force_component_axis": "global-z",
            "displacement_unit": "m",
            "energy_unit": "J",
            "force_unit": "N",
            "displacements_m": [-0.001, 0.0, 0.001],
            "force_displacements_m": [-0.001, 0.0, 0.001],
            "energy_samples_j": [1.001, 1.0, 0.999],
            "force_energy_samples_j": [1.001, 1.0, 0.999],
            "energy_sample_table_sha256": "b" * 64,
            "force_input_energy_sample_table_sha256": "b" * 64,
        },
    }
    return summary


def test_accepts_pinned_target_body_and_closed_surface_profiles() -> None:
    result = magnetic_force_method_profile_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_method_relative_difference"] < 0.05
    assert result["metrics"]["minimum_selection_scope_relative_difference"] > 0.25


def test_mcp_tool_dispatches_json() -> None:
    result = json.loads(mcp_magnetic_force_method_profile_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
    assert result["policy"] == "magnetic_force_method_profile_gate_v1"


def test_rejects_unpinned_all_body_selection() -> None:
    summary = copy.deepcopy(_summary())
    summary["all_body_element_force"] = list(summary["moving_body_element_force"])
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["selection_scope_is_materially_distinct"] is False
    assert result["checks"]["all_body_control_is_not_target_body_force"] is False


def test_rejects_method_disagreement_and_nonexact_replay() -> None:
    summary = copy.deepcopy(_summary())
    summary["closed_surface_maxwell_stress_force"][2] = 20.0
    summary["replay"]["parsed_max_abs"] = 1.0e-6
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_stress_replay_within_tolerance"] is False
    assert result["checks"]["parsed_replay_is_exact_enough"] is False


def test_rejects_profile_length_mismatch() -> None:
    summary = _summary()
    summary["all_body_element_force"].pop()
    with pytest.raises(ValueError, match="same length"):
        magnetic_force_method_profile_gate(summary)


def test_rejects_independent_surface_outlier_and_replay_drift() -> None:
    summary = copy.deepcopy(_summary())
    summary["independent_closed_surface_force"][2] *= 1.5
    summary["replay"]["parsed_max_abs"] = 1.0e-4
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["independent_stress_replay_within_tolerance"] is False
    assert result["checks"]["parsed_replay_is_exact_enough"] is False


@pytest.mark.parametrize(
    "case_id",
    ["quantity_dimension", "force_unit", "selection_control", "stress_method", "binary_replay"],
)
def test_counterfactual_curriculum90_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "quantity_dimension":
        summary["quantity_dimension"] = "2d_per_length"
    elif case_id == "force_unit":
        summary["force_unit"] = "mN"
    elif case_id == "selection_control":
        summary["all_body_element_force"] = list(summary["moving_body_element_force"])
    elif case_id == "stress_method":
        summary["closed_surface_maxwell_stress_force"][2] *= 2.0
    else:
        summary["replay"]["binary_nonlog_outputs_exact"] = False
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_generalization_v3s_rejects_unsupported_position_unit() -> None:
    summary = copy.deepcopy(_summary())
    summary["position_unit"] = "inch"
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v4_invalid_comparison_axis", "v4_position_order", "v4_element_force_sign", "v4_element_force_nonfinite", "v4_missing_target_signal"],
)
def test_counterfactual_curriculum90_v4_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v4_invalid_comparison_axis":
        summary["comparison_axis"] = "unsupported"
    elif case_id == "v4_position_order":
        summary["positions"][2] = summary["positions"][1]
    elif case_id == "v4_element_force_sign":
        summary["moving_body_element_force"][1] *= -1.0
    elif case_id == "v4_element_force_nonfinite":
        summary["moving_body_element_force"][4] = float("nan")
    else:
        summary["moving_body_element_force"] = [0.0] * len(summary["moving_body_element_force"])
    result = json.loads(mcp_magnetic_force_method_profile_gate(json.dumps(summary)))
    assert result["status"] in {"needs_attention", "invalid_input"}


def test_generalization_v5_rejects_short_target_force_profile() -> None:
    summary = copy.deepcopy(_summary())
    summary["moving_body_element_force"].pop()
    with pytest.raises(ValueError, match="same length"):
        magnetic_force_method_profile_gate(summary)


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_closed_surface_force_drift", "v6_public_force_dimension_unit_mismatch"],
)
def test_generalization_v6_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v6_public_closed_surface_force_drift":
        summary["closed_surface_maxwell_stress_force"][3] *= 1.30
    else:
        summary["force_unit"] = "N/m"
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v7_public_position_grid_unit_shadowing() -> None:
    summary = copy.deepcopy(_summary())
    summary["artifact_position_units"] = {
        "primary_positions": "m",
        "closed_surface_positions": "m",
        "independent_replay_positions": "mm",
    }
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["artifact_position_units_match_common_grid"] is False


def test_accepts_bound_sweep_and_demag_geometry_generations() -> None:
    result = magnetic_force_method_profile_gate(_with_artifact_identity(_summary()))
    assert result["status"] == "ok"
    assert result["warnings"] == []


def test_v8_public_force_position_branch_timestamp_mix() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["position_force_sample_generations"] = [
        "sweep-41" if index % 2 == 0 else "sweep-42"
        for index in range(len(summary["positions"]))
    ]
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["position_force_samples_share_one_sweep_generation"] is False


def test_v8_public_demag_reference_older_than_geometry() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["demag_reference"].update(
        {
            "geometry_revision": "geometry-11",
            "generated_at_utc": "2026-07-16T04:58:30Z",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["demag_reference_matches_current_geometry_revision"] is False


def test_v9_public_force_demag_coordinate_system_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["coordinate_system_binding"].update(
        {
            "force_component_frame_id": "rotor-local-reflected",
            "force_transform_revision": "rotation-11",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["force_and_demag_share_transformed_coordinate_system"]
        is False
    )


def test_v9_public_force_normalization_per_length_vs_total() -> None:
    summary = _with_artifact_identity(_summary())
    normalization = summary["artifact_identity"]["force_normalization"]
    normalization["profile_bases"]["moving_body_element_force"] = (
        "per_axial_length"
    )
    normalization["per_length_to_total_applied"] = False
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["force_profiles_share_total_3d_normalization"] is False


def test_v10_public_hysteresis_branch_direction_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["hysteresis_branch_state"].update(
        {
            "state_memory_branch": "ascending",
            "tangent_branch": "ascending",
            "tangent_state_generation": "branch-state-11",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hysteresis_observable_state_and_tangent_share_branch"
        ]
        is False
    )


def test_v10_public_remanence_vector_material_frame_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["remanence_frame_binding"].update(
        {
            "transform_applied": False,
            "transformed_vector_frame_id": "magnet-material-frame",
            "remanence_transform_revision": "rotation-11",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "remanence_vector_is_transformed_from_material_to_assembly_frame"
        ]
        is False
    )


def test_v11_public_force_surface_extra_magnet_compensated() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["force_surface_body_ownership"].update(
        {
            "enclosed_body_ids": ["moving-magnet", "fixed-magnet"],
            "compensating_body_force_allowed": True,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["force_surface_encloses_only_target_body"] is False


def test_v11_public_demag_operating_point_branch_interpolation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["demag_branch_interpolation"].update(
        {
            "interpolation_source_branch": "ascending",
            "interpolation_state_generation": "branch-state-12",
            "bracketing_sample_ids": ["asc-17", "asc-18"],
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "demag_operating_point_uses_active_branch_interpolation"
        ]
        is False
    )


def test_v12_public_linear_motor_thrust_phase_sequence_electrical_angle_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["linear_motor_thrust_phase_identity"].update(
        {
            "thrust_phase_sequence": ["U", "W", "V"],
            "thrust_electrical_angle_direction": -1,
            "thrust_observable_phase_generation": "phase-convention-13",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "linear_motor_thrust_uses_winding_phase_and_electrical_angle_direction"
        ]
        is False
    )


def test_v12_public_demag_recoil_line_temperature_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["demag_recoil_temperature_identity"].update(
        {
            "recoil_line_temperature_c": 20.0,
            "recoil_line_state_generation": "magnet-material-13",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "demag_recoil_line_matches_evaluation_temperature_generation"
        ]
        is False
    )


def test_v13_public_bem_demag_surface_normal_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "bem_demag_surface_normal_generation_identity"
    ].update(
        {
            "surface_normal_generation": "surface-mesh-14",
            "surface_normal_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["bem_demag_normals_match_current_surface_mesh_generation"]
        is False
    )


def test_v13_public_cogging_torque_periodic_sector_symmetry_multiplier_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "cogging_torque_periodic_sector_symmetry_identity"
    ].update(
        {
            "torque_result_periodic_sector_count": 6,
            "symmetry_multiplier": 6.0,
            "multiplier_topology_generation": "periodic-topology-14",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cogging_torque_symmetry_multiplier_matches_periodic_sector"
        ]
        is False
    )


def test_v14_public_bem_self_term_solid_angle_orientation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "bem_self_term_solid_angle_orientation_identity"
    ].update(
        {
            "self_term_solid_angle_orientation": "inward",
            "self_term_sign_convention": "inward_positive",
            "self_term_orientation_generation": "surface-mesh-15",
            "self_term_orientation_sha256": "d" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "bem_self_term_solid_angle_matches_panel_orientation_generation"
        ]
        is False
    )


def test_v14_public_demag_energy_force_displacement_length_unit_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "demag_energy_force_displacement_length_unit_identity"
    ].update(
        {
            "force_derivative_length_unit": "mm",
            "force_derivative_scale_to_m": 1.0e-3,
            "force_derivative_displacement_generation": "displacement-grid-15",
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "demag_energy_force_derivative_uses_one_displacement_length_unit"
        ]
        is False
    )


def test_v15_public_bem_near_singular_quadrature_target_scale_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "bem_near_singular_quadrature_target_scale_identity"
    ].update(
        {
            "quadrature_target_separation_generation": "target-separation-16",
            "quadrature_normalized_target_separation": 5.0e-2,
            "quadrature_target_pair_sha256": "5" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "bem_near_singular_quadrature_uses_current_target_pair_scale"
        ]
        is False
    )


def test_v15_public_force_torque_reference_origin_length_unit_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "force_torque_reference_origin_length_unit_identity"
    ].update(
        {
            "torque_reference_origin_length_unit": "mm",
            "torque_reference_origin_scale_to_m": 1.0e-3,
            "torque_reference_origin_sha256": "5" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["force_and_torque_share_reference_origin_length_unit"]
        is False
    )


def test_v16_public_bem_solid_angle_self_term_surface_winding_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"]["bem_solid_angle_surface_winding_identity"].update(
        {
            "solid_angle_sign_generation": "surface-mesh-17",
            "solid_angle_sign_convention": "inward_positive",
            "solid_angle_winding_sha256": "4" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["bem_self_term_uses_current_normalized_surface_winding"]
        is False
    )


def test_v16_public_maglev_stiffness_force_displacement_reference_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "maglev_stiffness_force_displacement_generation_identity"
    ].update(
        {
            "displacement_coordinate_generation": "perturbation-17",
            "stiffness_displacement_sha256": "4" * 64,
            "stiffness_force_sample_sha256": "4" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "maglev_stiffness_uses_one_force_displacement_perturbation_generation"
        ]
        is False
    )


def test_v17_public_bem_demag_tensor_coordinate_basis_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "bem_demag_tensor_coordinate_basis_generation_identity"
    ].update(
        {
            "demag_tensor_body_placement_generation": "placement-18",
            "tensor_coordinate_basis": "body-local-previous",
            "tensor_basis_transform_sha256": "5" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "bem_demag_tensor_uses_current_body_placement_coordinate_basis"
        ]
        is False
    )


def test_v17_public_magnetic_bearing_force_harmonic_phase_origin_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "magnetic_bearing_force_harmonic_phase_origin_identity"
    ].update(
        {
            "force_harmonic_phase_origin_deg": 15.0,
            "force_harmonic_phase_origin_convention": "next_slot_center",
            "harmonic_input_angle_sha256": "5" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "magnetic_bearing_force_harmonics_share_rotor_phase_origin"
        ]
        is False
    )


def test_v17_public_magnetic_bearing_phase_origin_is_periodic() -> None:
    summary = _with_artifact_identity(_summary())
    bearing = summary["artifact_identity"][
        "magnetic_bearing_force_harmonic_phase_origin_identity"
    ]
    bearing["force_harmonic_phase_origin_deg"] = 360.0
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "ok"


def test_v17_public_magnetic_bearing_slot_pitch_partitions_revolution() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "magnetic_bearing_force_harmonic_phase_origin_identity"
    ]["slot_pitch_deg"] = 123.456
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "magnetic_bearing_force_harmonics_share_rotor_phase_origin"
        ]
        is False
    )


def test_v17_public_demag_tensor_generation_binds_force_consumer() -> None:
    summary = _with_artifact_identity(_summary())
    summary["artifact_identity"][
        "bem_demag_tensor_coordinate_basis_generation_identity"
    ]["demag_tensor_generation"] = "stale-demag-tensor"
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "bem_demag_tensor_uses_current_body_placement_coordinate_basis"
        ]
        is False
    )


def test_v19_public_bem_demag_self_term_solid_angle_orientation_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    identity = summary["artifact_identity"][
        "bem_demag_self_term_solid_angle_orientation_generation_identity"
    ]
    identity.update(
        {
            "self_term_boundary_generation": "boundary-20",
            "panel_orientation_boundary_generation": "boundary-20",
            "self_term_panel_ids": [103, 102, 101],
            "self_term_orientation_signs": [-1, -1, -1],
            "assembled_self_term_input_sha256": "e" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "bem_demag_self_term_uses_current_boundary_orientation_generation"
        ]
        is False
    )


def test_v19_public_virtual_work_force_displacement_coordinate_generation_mismatch() -> None:
    summary = _with_artifact_identity(_summary())
    identity = summary["artifact_identity"][
        "virtual_work_force_displacement_coordinate_geometry_generation_identity"
    ]
    identity.update(
        {
            "displacement_coordinate_geometry_generation": "geometry-20",
            "energy_sample_geometry_generations": ["geometry-20"] * 3,
            "force_displacement_coordinate_generation": "displacement-20",
            "force_displacements_m": [0.001, 0.0, -0.001],
            "force_input_energy_sample_table_sha256": "e" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "virtual_work_force_uses_current_displacement_geometry_generation"
        ]
        is False
    )


def _summary_v20():
    summary = _with_artifact_identity(_summary())
    identity = summary["artifact_identity"]
    identity["demag_energy_surface_charge_normal_quadrature_generation_identity"] = {
        "energy_generation": "demag-energy-22",
        "result_energy_generation": "demag-energy-22",
        "boundary_generation": "boundary-22",
        "surface_charge_boundary_generation": "boundary-22",
        "normal_boundary_generation": "boundary-22",
        "quadrature_boundary_generation": "boundary-22",
        "energy_boundary_generation": "boundary-22",
        "panel_ids": [101, 102, 103],
        "energy_panel_ids": [101, 102, 103],
        "surface_charges_a_per_m": [1.0, -0.5, 0.25],
        "energy_surface_charges_a_per_m": [1.0, -0.5, 0.25],
        "outward_normal_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "energy_outward_normal_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "panel_quadrature_weights": [0.25, 0.5, 0.25],
        "energy_panel_quadrature_weights": [0.25, 0.5, 0.25],
        "demag_energy_input_sha256": "4" * 64,
        "assembled_demag_energy_input_sha256": "4" * 64,
    }
    identity[
        "maglev_stiffness_displacement_equilibrium_force_generation_identity"
    ] = {
        "stiffness_generation": "stiffness-22",
        "result_stiffness_generation": "stiffness-22",
        "geometry_generation": "geometry-22",
        "sample_geometry_generations": ["geometry-22"] * 3,
        "force_geometry_generations": ["geometry-22"] * 3,
        "coordinate_generation": "coordinate-22",
        "displacement_coordinate_generation": "coordinate-22",
        "force_coordinate_generation": "coordinate-22",
        "displacement_axis": "global-z",
        "force_component_axis": "global-z",
        "displacements_m": [-0.001, 0.0, 0.001],
        "force_displacements_m": [-0.001, 0.0, 0.001],
        "force_samples_n": [10.0, 0.0, -10.0],
        "stiffness_force_samples_n": [10.0, 0.0, -10.0],
        "equilibrium_sample_index": 1,
        "stiffness_equilibrium_sample_index": 1,
        "stiffness_sample_table_sha256": "5" * 64,
        "result_stiffness_sample_table_sha256": "5" * 64,
    }
    return summary


def test_v20_public_positive_demag_energy_and_maglev_stiffness_identity() -> None:
    result = magnetic_force_method_profile_gate(_summary_v20())
    assert result["status"] == "ok"
    assert result["checks"][
        "demag_energy_uses_current_surface_charge_normal_and_quadrature"
    ]
    assert result["checks"][
        "maglev_stiffness_uses_aligned_displacement_equilibrium_force_states"
    ]


def test_v20_public_demag_energy_surface_charge_normal_quadrature_generation_mismatch() -> None:
    summary = _summary_v20()
    summary["artifact_identity"][
        "demag_energy_surface_charge_normal_quadrature_generation_identity"
    ].update(
        {
            "surface_charge_boundary_generation": "boundary-21",
            "normal_boundary_generation": "boundary-21",
            "energy_panel_ids": [103, 102, 101],
            "energy_surface_charges_a_per_m": [0.25, -0.5, 1.0],
            "energy_outward_normal_sha256": ["3" * 64, "2" * 64, "1" * 64],
            "energy_panel_quadrature_weights": [0.5, 0.25, 0.25],
            "assembled_demag_energy_input_sha256": "f" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "demag_energy_uses_current_surface_charge_normal_and_quadrature"
    ] is False


def test_v20_public_maglev_stiffness_displacement_equilibrium_force_generation_mismatch() -> None:
    summary = _summary_v20()
    summary["artifact_identity"][
        "maglev_stiffness_displacement_equilibrium_force_generation_identity"
    ].update(
        {
            "sample_geometry_generations": ["geometry-21"] * 3,
            "force_coordinate_generation": "coordinate-21",
            "force_component_axis": "global-x",
            "force_displacements_m": [0.001, 0.0, -0.001],
            "stiffness_force_samples_n": [-10.0, 1.0, 10.0],
            "stiffness_equilibrium_sample_index": 0,
            "result_stiffness_sample_table_sha256": "f" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "maglev_stiffness_uses_aligned_displacement_equilibrium_force_states"
    ] is False


def _summary_v21():
    summary = _summary_v20()
    identity = summary["artifact_identity"]
    identity["bem_near_singular_distance_panel_quadrature_generation_identity"] = {
        "interaction_generation": "bem-interaction-31",
        "result_interaction_generation": "bem-interaction-31",
        "geometry_generation": "bem-geometry-31",
        "target_distance_geometry_generation": "bem-geometry-31",
        "source_panel_geometry_generation": "bem-geometry-31",
        "adaptive_quadrature_geometry_generation": "bem-geometry-31",
        "target_point_id": 501,
        "distance_target_point_id": 501,
        "source_panel_id": 101,
        "distance_source_panel_id": 101,
        "target_distance_m": 1.0e-5,
        "quadrature_target_distance_m": 1.0e-5,
        "source_panel_geometry_sha256": "1" * 64,
        "quadrature_source_panel_geometry_sha256": "1" * 64,
        "adaptive_quadrature_order": 16,
        "evaluated_quadrature_order": 16,
        "near_singular_interaction_sha256": "2" * 64,
        "assembled_near_singular_interaction_sha256": "2" * 64,
    }
    identity[
        "moving_magnet_force_position_orientation_equilibrium_generation_identity"
    ] = {
        "force_generation": "moving-force-31",
        "result_force_generation": "moving-force-31",
        "geometry_generation": "moving-geometry-31",
        "position_sample_geometry_generations": ["moving-geometry-31"] * 3,
        "orientation_geometry_generation": "moving-geometry-31",
        "equilibrium_geometry_generation": "moving-geometry-31",
        "force_geometry_generation": "moving-geometry-31",
        "position_samples_m": [
            [0.0, 0.0, -0.001],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.001],
        ],
        "force_position_samples_m": [
            [0.0, 0.0, -0.001],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.001],
        ],
        "magnet_orientation_quaternion": [1.0, 0.0, 0.0, 0.0],
        "force_orientation_quaternion": [1.0, 0.0, 0.0, 0.0],
        "force_samples_n": [
            [0.0, 0.0, 10.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -10.0],
        ],
        "differentiated_force_samples_n": [
            [0.0, 0.0, 10.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, -10.0],
        ],
        "equilibrium_sample_index": 1,
        "force_equilibrium_sample_index": 1,
        "moving_force_sample_table_sha256": "3" * 64,
        "result_moving_force_sample_table_sha256": "3" * 64,
    }
    return summary


def test_v21_public_positive_bem_and_moving_magnet_force_identity() -> None:
    result = magnetic_force_method_profile_gate(_summary_v21())
    assert result["status"] == "ok"
    assert result["checks"][
        "bem_near_singular_quadrature_uses_current_distance_and_panel_geometry"
    ]
    assert result["checks"][
        "moving_magnet_force_uses_current_position_orientation_and_equilibrium"
    ]


def test_v21_public_bem_near_singular_distance_panel_quadrature_generation_mismatch() -> None:
    summary = _summary_v21()
    summary["artifact_identity"][
        "bem_near_singular_distance_panel_quadrature_generation_identity"
    ].update(
        {
            "target_distance_geometry_generation": "bem-geometry-30",
            "source_panel_geometry_generation": "bem-geometry-29",
            "adaptive_quadrature_geometry_generation": "bem-geometry-28",
            "distance_target_point_id": 502,
            "distance_source_panel_id": 102,
            "quadrature_target_distance_m": 1.0e-3,
            "quadrature_source_panel_geometry_sha256": "a" * 64,
            "evaluated_quadrature_order": 4,
            "assembled_near_singular_interaction_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "bem_near_singular_quadrature_uses_current_distance_and_panel_geometry"
    ] is False


def test_v21_public_moving_magnet_force_position_orientation_equilibrium_generation_mismatch() -> None:
    summary = _summary_v21()
    summary["artifact_identity"][
        "moving_magnet_force_position_orientation_equilibrium_generation_identity"
    ].update(
        {
            "position_sample_geometry_generations": ["moving-geometry-30"] * 3,
            "orientation_geometry_generation": "moving-geometry-29",
            "equilibrium_geometry_generation": "moving-geometry-28",
            "force_position_samples_m": [
                [0.0, 0.0, 0.001],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, -0.001],
            ],
            "force_orientation_quaternion": [0.0, 0.0, 0.0, 1.0],
            "force_equilibrium_sample_index": 0,
            "result_moving_force_sample_table_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "moving_magnet_force_uses_current_position_orientation_and_equilibrium"
    ] is False


def _summary_v22():
    summary = _summary_v21()
    identity = summary["artifact_identity"]
    identity["motor_force_dual_lane_interface_flux_coenergy_generation_identity"] = {
        "comparison_generation": "motor-force-41",
        "lane_policy_comparison_generation": "motor-force-41",
        "interface_flux_comparison_generation": "motor-force-41",
        "coenergy_comparison_generation": "motor-force-41",
        "force_comparison_generation": "motor-force-41",
        "coupling_mesh_comparison_generation": "motor-force-41",
        "lane_ids": ["ngsolve_age", "hdiv_mmm_hcurl_eddy_bubble"],
        "result_lane_ids": ["ngsolve_age", "hdiv_mmm_hcurl_eddy_bubble"],
        "interface_normal_flux_wb": [0.031, 0.0308],
        "result_interface_normal_flux_wb": [0.031, 0.0308],
        "coenergy_j": [1.25, 1.247],
        "result_coenergy_j": [1.25, 1.247],
        "force_n": [14.2, 14.15],
        "result_force_n": [14.2, 14.15],
        "coupling_mesh_sha256": "1" * 64,
        "result_coupling_mesh_sha256": "1" * 64,
        "mixed_operator_contract_sha256": "2" * 64,
        "result_mixed_operator_contract_sha256": "2" * 64,
    }
    identity[
        "linear_motor_end_effect_translation_position_symmetry_generation_identity"
    ] = {
        "sweep_generation": "linear-motion-41",
        "position_sweep_generation": "linear-motion-41",
        "end_effect_sweep_generation": "linear-motion-41",
        "translation_frame_sweep_generation": "linear-motion-41",
        "symmetry_sweep_generation": "linear-motion-41",
        "force_result_sweep_generation": "linear-motion-41",
        "mover_positions_m": [-0.01, 0.0, 0.01],
        "result_mover_positions_m": [-0.01, 0.0, 0.01],
        "end_effect_window_m": [-0.04, 0.04],
        "result_end_effect_window_m": [-0.04, 0.04],
        "translation_frame": "stator_x",
        "result_translation_frame": "stator_x",
        "symmetry_factor": 2,
        "result_symmetry_factor": 2,
        "thrust_n": [80.0, 100.0, 78.0],
        "result_thrust_n": [80.0, 100.0, 78.0],
        "stiffness_n_per_m": [2000.0, -100.0, -2200.0],
        "result_stiffness_n_per_m": [2000.0, -100.0, -2200.0],
        "linear_motion_table_sha256": "3" * 64,
        "result_linear_motion_table_sha256": "3" * 64,
    }
    return summary


def test_v22_public_positive_canonical_dual_lane_and_linear_motor_identity() -> None:
    assert magnetic_force_method_profile_gate(_summary_v22())["status"] == "ok"


def test_v22_public_hdiv_vim_reduced_fem_interface_flux_coenergy_force_generation_mismatch() -> None:
    summary = _summary_v22()
    identity = summary["artifact_identity"][
        "motor_force_dual_lane_interface_flux_coenergy_generation_identity"
    ]
    identity.update(
        {
            "lane_policy_comparison_generation": "motor-force-40",
            "interface_flux_comparison_generation": "motor-force-39",
            "coenergy_comparison_generation": "motor-force-38",
            "force_comparison_generation": "motor-force-37",
            "result_lane_ids": ["ngsolve_age", "hdiv_vim_reduced_fem"],
            "result_interface_normal_flux_wb": [0.031, 0.026],
            "result_coenergy_j": [1.25, 1.08],
            "result_force_n": [14.2, 11.4],
            "result_coupling_mesh_sha256": "a" * 64,
            "result_mixed_operator_contract_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_force_comparison_uses_age_and_hdiv_mmm_hcurl_eddy_bubble_lanes"
    ]


def test_v22_public_linear_motor_end_effect_translation_position_symmetry_generation_mismatch() -> None:
    summary = _summary_v22()
    identity = summary["artifact_identity"][
        "linear_motor_end_effect_translation_position_symmetry_generation_identity"
    ]
    identity.update(
        {
            "position_sweep_generation": "linear-motion-40",
            "end_effect_sweep_generation": "linear-motion-39",
            "translation_frame_sweep_generation": "linear-motion-38",
            "symmetry_sweep_generation": "linear-motion-37",
            "result_mover_positions_m": [0.01, 0.0, -0.01],
            "result_end_effect_window_m": [-0.02, 0.02],
            "result_translation_frame": "mover_x",
            "result_symmetry_factor": 1,
            "result_thrust_n": [39.0, 50.0, 40.0],
            "result_stiffness_n_per_m": [-1100.0, -50.0, 1000.0],
            "result_linear_motion_table_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "linear_motor_force_uses_current_position_end_effect_frame_and_symmetry"
    ]
