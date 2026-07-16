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
