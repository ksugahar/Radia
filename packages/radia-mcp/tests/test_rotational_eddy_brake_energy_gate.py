from __future__ import annotations

import copy
import json
import math

import pytest

from radia_mcp.radia_ngsolve.rotational_eddy_brake_energy_gate import (
    rotational_eddy_brake_energy_gate as gate,
)
from radia_mcp.radia_ngsolve.server import rotational_eddy_brake_energy_gate as mcp_gate


def _summary() -> dict:
    density, radius, thickness = 7800.0, 0.1, 0.01
    inertia = 0.5 * density * math.pi * radius**4 * thickness
    times = [0.05 * index for index in range(201)]
    rate = 0.2
    omega = [100.0 * math.exp(-rate * value) for value in times]
    torque = [inertia * rate * value for value in omega]
    joule = [force * speed for force, speed in zip(torque, omega, strict=True)]
    generation = "solve-generation-42"
    frame = "laboratory-frame-generation-7"

    def replay(label: str) -> dict:
        return {
            "label": label,
            "solve_seconds": 2.0,
            "time_s": list(times),
            "angular_velocity_rad_s": list(omega),
            "braking_torque_nm": list(torque),
            "joule_loss_w": list(joule),
            "artifact_generations": {
                "time_s": generation,
                "angular_velocity_rad_s": generation,
                "braking_torque_nm": generation,
                "joule_loss_w": generation,
            },
            "artifact_coordinate_frames": {
                "time_s": frame,
                "angular_velocity_rad_s": frame,
                "braking_torque_nm": frame,
                "joule_loss_w": frame,
            },
        }

    boundary = 100
    accumulated_joule = sum(
        0.5 * (joule[index] + joule[index + 1])
        * (times[index + 1] - times[index])
        for index in range(boundary)
    )
    stored_energy = 0.5 * inertia * omega[boundary] ** 2 + 0.3

    return {
        "contract": {
            "body": "uniform_solid_conducting_disc",
            "inertia_reference": "analytic_uniform_solid_disc",
            "angular_momentum_balance": "inertia_delta_angular_velocity_plus_integrated_braking_torque_equals_zero",
            "instantaneous_power_comparison": "diagnostic_only_when_field_energy_rate_is_not_sampled_on_the_probe_grid",
            "energy_balance": "initial_kinetic_plus_magnetic_equals_final_kinetic_plus_magnetic_plus_joule",
        },
        "units": {
            "time": "s",
            "angular_velocity": "rad/s",
            "torque": "N*m",
            "power": "W",
            "inertia": "kg*m^2",
            "energy": "J",
            "density": "kg/m^3",
            "length": "m",
        },
        "disc": {
            "density_kg_m3": density,
            "radius_m": radius,
            "thickness_m": thickness,
        },
        "reported_inertia_kg_m2": inertia,
        "replays": [replay("one"), replay("two")],
        "energy_replay": {
            **replay("energy"),
            "field_energy_time_s": list(times),
            "magnetic_energy_j": [0.3 for _ in times],
            "artifact_generations": {
                "time_s": generation,
                "angular_velocity_rad_s": generation,
                "braking_torque_nm": generation,
                "joule_loss_w": generation,
                "field_energy_time_s": generation,
                "magnetic_energy_j": generation,
            },
            "artifact_coordinate_frames": {
                "time_s": frame,
                "angular_velocity_rad_s": frame,
                "braking_torque_nm": frame,
                "joule_loss_w": frame,
                "field_energy_time_s": frame,
                "magnetic_energy_j": frame,
            },
            "restart_boundaries": [
                {
                    "left_index": boundary,
                    "right_index": boundary + 1,
                    "generation_before": generation,
                    "generation_after": "solve-generation-42-restart-1",
                    "stored_energy_before_j": stored_energy,
                    "stored_energy_after_j": stored_energy,
                    "accumulated_joule_before_j": accumulated_joule,
                    "accumulated_joule_offset_after_j": accumulated_joule,
                }
            ],
        },
        "convergence_provenance": {
            "solution_generation": "solution-generation-42",
            "result_iteration_generation": "nonlinear-iteration-8",
            "convergence_table_iteration_generation": "nonlinear-iteration-8",
            "terminal_state": "converged",
            "terminal_relative_residual": 2.0e-9,
        },
        "timing_breakdown_s": {
            "attach": 0.1,
            "solve": 6.0,
            "extract": 0.1,
            "verify": 0.1,
        },
    }


def test_accepts_free_brake_with_field_energy_and_mcp_dispatch() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["kinetic_magnetic_joule_energy_closes"] is True
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_wrong_inertia_and_angular_impulse() -> None:
    summary = copy.deepcopy(_summary())
    summary["reported_inertia_kg_m2"] *= 1.25
    summary["contract"]["inertia_reference"] = "unchecked_reported_value"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["analytic_disc_inertia_matches_reported"] is False
    assert result["checks"]["angular_impulse_balance_closes"] is False


def test_rejects_missing_field_storage_and_false_power_claim() -> None:
    summary = copy.deepcopy(_summary())
    summary["contract"]["instantaneous_power_comparison"] = "always_equal_to_joule"
    summary["energy_replay"]["magnetic_energy_j"] = []
    try:
        result = gate(summary)
    except ValueError as exc:
        assert "magnetic_energy_j" in str(exc)
    else:
        assert result["status"] == "needs_attention"


def test_rejects_nonreplaying_torque_history() -> None:
    summary = copy.deepcopy(_summary())
    summary["replays"][1]["braking_torque_nm"] = list(
        summary["replays"][1]["braking_torque_nm"]
    )
    summary["replays"][1]["braking_torque_nm"][50] *= 1.2
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_replay_fields_match"] is False


def test_rejects_joule_history_inconsistent_with_motion_energy() -> None:
    summary = copy.deepcopy(_summary())
    summary["energy_replay"]["joule_loss_w"] = [
        1.4 * value for value in summary["energy_replay"]["joule_loss_w"]
    ]
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["kinetic_magnetic_joule_energy_closes"] is False


def test_rejects_isolated_magnetic_energy_spike() -> None:
    summary = copy.deepcopy(_summary())
    summary["energy_replay"]["magnetic_energy_j"][10] *= 4.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "field_energy_history_is_nonnegative_and_has_no_isolated_jump"
        ]
        is False
    )


@pytest.mark.parametrize(
    "case_id",
    ["disc_radius", "momentum_contract", "torque_replay", "torque_unit", "field_energy"],
)
def test_counterfactual_curriculum90_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "disc_radius":
        summary["disc"]["radius_m"] *= 1.2
    elif case_id == "momentum_contract":
        summary["contract"]["angular_momentum_balance"] = "unchecked"
    elif case_id == "torque_replay":
        summary["replays"][1]["braking_torque_nm"] = list(
            summary["replays"][1]["braking_torque_nm"]
        )
        summary["replays"][1]["braking_torque_nm"][25] *= 1.25
    elif case_id == "torque_unit":
        summary["units"]["torque"] = "kg"
    else:
        summary["energy_replay"]["magnetic_energy_j"][10] *= 4.0
    assert gate(summary)["status"] == "needs_attention"


def test_generalization_v3s_rejects_negative_field_energy() -> None:
    summary = copy.deepcopy(_summary())
    summary["energy_replay"]["magnetic_energy_j"][17] = -0.1
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v4_time_duplicate",
        "v4_velocity_local_rise",
        "v4_positive_braking_torque",
        "v4_negative_joule_loss",
        "v4_field_time_shift",
    ],
)
def test_counterfactual_curriculum90_v4_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v4_time_duplicate":
        summary["replays"][0]["time_s"][20] = summary["replays"][0]["time_s"][19]
    elif case_id == "v4_velocity_local_rise":
        summary["replays"][0]["angular_velocity_rad_s"][25] *= 1.5
    elif case_id == "v4_positive_braking_torque":
        summary["replays"][0]["braking_torque_nm"][20] = 1.0
    elif case_id == "v4_negative_joule_loss":
        summary["energy_replay"]["joule_loss_w"][12] = -1.0
    else:
        summary["energy_replay"]["field_energy_time_s"][30] *= 1.1
    assert gate(summary)["status"] == "needs_attention"


def test_generalization_v5_rejects_non_si_energy_unit() -> None:
    summary = copy.deepcopy(_summary())
    summary["units"]["energy"] = "mJ"
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_energy_state_drift", "v6_public_replay_joule_drift"],
)
def test_generalization_v6_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v6_public_energy_state_drift":
        summary["energy_replay"]["magnetic_energy_j"][20] *= 1.10
    else:
        summary["replays"][1]["joule_loss_w"] = list(
            summary["replays"][1]["joule_loss_w"]
        )
        summary["replays"][1]["joule_loss_w"][20] *= 1.05
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_segment_overlap_after_restart",
        "v7_public_field_loss_timebase_skew",
    ],
)
def test_generalization_v7_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v7_public_segment_overlap_after_restart":
        for replay in [*summary["replays"], summary["energy_replay"]]:
            replay["time_s"][101] = replay["time_s"][99]
        summary["energy_replay"]["field_energy_time_s"][101] = summary[
            "energy_replay"
        ]["field_energy_time_s"][99]
    else:
        for index in range(60, 140):
            summary["energy_replay"]["field_energy_time_s"][index] += 0.0125
    assert gate(summary)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_public_restart_energy_baseline_reset",
        "v8_public_loss_history_prior_solve_generation",
    ],
)
def test_generalization_v8_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v8_public_restart_energy_baseline_reset":
        summary["energy_replay"]["restart_boundaries"][0][
            "accumulated_joule_offset_after_j"
        ] = 0.0
        expected_check = "restart_energy_offsets_are_continuous"
    else:
        summary["energy_replay"]["artifact_generations"][
            "joule_loss_w"
        ] = "solve-generation-41"
        expected_check = "artifact_series_share_their_solve_generation"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected_check] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_public_torque_energy_coordinate_frame_mismatch",
        "v9_public_convergence_table_prior_nonlinear_iteration",
    ],
)
def test_generalization_v9_public(case_id: str) -> None:
    summary = copy.deepcopy(_summary())
    if case_id == "v9_public_torque_energy_coordinate_frame_mismatch":
        summary["energy_replay"]["artifact_coordinate_frames"][
            "braking_torque_nm"
        ] = "rotor-local-frame-generation-7"
        expected_check = "artifact_series_share_one_coordinate_frame"
    else:
        summary["convergence_provenance"][
            "convergence_table_iteration_generation"
        ] = "nonlinear-iteration-7"
        expected_check = "convergence_table_matches_result_iteration"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected_check] is False


def _with_v10_ownership(summary: dict) -> dict:
    summary["force_selection_identity"] = {
        "geometry_generation": "geometry-generation-12",
        "solution_geometry_generation": "geometry-generation-12",
        "integration_selection_generation": "geometry-generation-12",
        "selection_entity_digest": "d" * 64,
        "force_result_selection_digest": "d" * 64,
    }
    summary["excitation_basis_identity"] = {
        "sweep_generation": "sweep-generation-42",
        "solve_amplitude_basis": "rms",
        "extract_amplitude_basis": "rms",
        "solve_scale_to_rms": 1.0,
        "extract_scale_to_rms": 1.0,
        "torque_loss_normalization_basis": "rms_excitation",
    }
    return summary


def _with_v11_artifact_ownership(summary: dict) -> dict:
    summary = _with_v10_ownership(summary)
    summary["live_stored_force_identity"] = {
        "live_geometry_sha256": "1" * 64,
        "stored_force_geometry_sha256": "1" * 64,
        "live_solution_generation": "solution-generation-43",
        "stored_force_solution_generation": "solution-generation-43",
        "live_selection_digest": "2" * 64,
        "stored_force_selection_digest": "2" * 64,
    }
    summary["loss_partition_identity"] = {
        "aggregation_generation": "loss-aggregation-43",
        "reported_total_generation": "loss-aggregation-43",
        "partition_ownership_ids": [
            "rotor-domain",
            "stator-domain",
            "interface-layer",
        ],
        "ownership_overlap_count": 0,
        "signed_compensation_term_w": 0.0,
    }
    return summary


def _with_v12_parameter_topology_identity(summary: dict) -> dict:
    summary = _with_v11_artifact_ownership(summary)
    summary["material_property_parameter_identity"] = {
        "parameter_name": "sigma_scale",
        "parameter_value": 1.25,
        "material_evaluation_parameter_value": 1.25,
        "parameter_unit": "1",
        "material_evaluation_parameter_unit": "1",
        "parameter_generation": "parameter-row-44",
        "material_property_parameter_generation": "parameter-row-44",
    }
    summary["force_selection_topology_identity"] = {
        "geometry_rebuild_generation": "geometry-rebuild-44",
        "topology_generation": "topology-44",
        "selection_topology_generation": "topology-44",
        "force_integration_topology_generation": "topology-44",
        "selection_entity_ids": [17, 18, 19],
        "force_integration_entity_ids": [17, 18, 19],
        "selection_digest": "8" * 64,
        "force_selection_digest": "8" * 64,
    }
    return summary


def _with_v13_transform_identity(summary: dict) -> dict:
    summary = _with_v12_parameter_topology_identity(summary)
    summary["weak_form_coordinate_transform_identity"] = {
        "mapped_coordinate_generation": "mapped-coordinates-45",
        "field_coordinate_generation": "mapped-coordinates-45",
        "jacobian_coordinate_generation": "mapped-coordinates-45",
        "jacobian_orientation": "right_handed_positive",
        "integration_orientation": "right_handed_positive",
        "jacobian_orientation_sha256": "c" * 64,
        "integration_orientation_sha256": "c" * 64,
    }
    summary["time_harmonic_phasor_convention_identity"] = {
        "source_time_convention": "exp(+jwt)",
        "field_time_convention": "exp(+jwt)",
        "phase_sensitive_result_time_convention": "exp(+jwt)",
        "complex_power_formula": "0.5*V*conj(I)",
        "phase_sensitive_result_formula": "0.5*V*conj(I)",
        "phasor_generation": "harmonic-solution-45",
        "result_phasor_generation": "harmonic-solution-45",
    }
    return summary


def _with_v14_mode_ale_identity(summary: dict) -> dict:
    summary = _with_v13_transform_identity(summary)
    summary["eigenmode_mass_inner_product_identity"] = {
        "eigensolve_generation": "eigensolve-46",
        "mode_vector_generation": "eigensolve-46",
        "mode_mesh_generation": "mesh-46",
        "mass_matrix_mesh_generation": "mesh-46",
        "eigensolve_mass_matrix_generation": "mass-matrix-46",
        "normalization_mass_matrix_generation": "mass-matrix-46",
        "normalization_kind": "mass_inner_product",
        "reference_normalization_kind": "mass_inner_product",
        "mode_vector_sha256": "a" * 64,
        "normalized_mode_vector_sha256": "a" * 64,
    }
    summary["ale_material_derivative_time_level_identity"] = {
        "ale_solve_generation": "ale-solve-46",
        "field_solve_generation": "ale-solve-46",
        "mesh_velocity_solve_generation": "ale-solve-46",
        "accepted_time_level_generation": "time-level-46",
        "field_time_level_generation": "time-level-46",
        "mesh_velocity_time_level_generation": "time-level-46",
        "material_derivative_time_level_generation": "time-level-46",
        "accepted_time_index": 128,
        "field_time_index": 128,
        "mesh_velocity_time_index": 128,
        "material_derivative_time_index": 128,
        "accepted_time_grid_sha256": "b" * 64,
        "mesh_velocity_time_grid_sha256": "b" * 64,
    }
    return summary


def _with_v15_reference_jacobian_identity(summary: dict) -> dict:
    summary = _with_v14_mode_ale_identity(summary)
    summary["harmonic_reference_time_origin_identity"] = {
        "harmonic_solve_generation": "harmonic-solve-47",
        "field_phasor_generation": "harmonic-solve-47",
        "complex_power_generation": "harmonic-solve-47",
        "angular_frequency_rad_s": 3141.592653589793,
        "field_reference_time_s": 0.0,
        "complex_power_reference_time_s": 0.0,
        "field_phase_origin_sha256": "1" * 64,
        "power_phase_origin_sha256": "1" * 64,
    }
    summary["deformed_domain_integral_jacobian_identity"] = {
        "field_solve_generation": "deformed-solve-47",
        "integral_field_generation": "deformed-solve-47",
        "geometry_generation": "deformed-geometry-47",
        "integral_geometry_generation": "deformed-geometry-47",
        "volume_jacobian_geometry_generation": "deformed-geometry-47",
        "domain_selection_sha256": "2" * 64,
        "integrated_domain_selection_sha256": "2" * 64,
        "volume_jacobian_sha256": "3" * 64,
        "integral_volume_jacobian_sha256": "3" * 64,
    }
    return summary


def _with_v16_nonlinear_ale_frame_identity(summary: dict) -> dict:
    summary = _with_v15_reference_jacobian_identity(summary)
    summary["nonlinear_residual_tangent_iteration_identity"] = {
        "nonlinear_solve_generation": "nonlinear-solve-50",
        "residual_solve_generation": "nonlinear-solve-50",
        "tangent_solve_generation": "nonlinear-solve-50",
        "material_state_solve_generation": "nonlinear-solve-50",
        "residual_iteration": 12,
        "tangent_iteration": 12,
        "material_state_iteration": 12,
        "material_state_sha256": "6" * 64,
        "tangent_material_state_sha256": "6" * 64,
    }
    summary["moving_mesh_field_transfer_frame_identity"] = {
        "mesh_motion_generation": "ale-motion-50",
        "source_mesh_motion_generation": "ale-motion-50",
        "target_mesh_motion_generation": "ale-motion-50",
        "field_transfer_mesh_motion_generation": "ale-motion-50",
        "source_coordinate_frame": "material",
        "target_coordinate_frame": "material",
        "field_transfer_coordinate_frame": "material",
        "coordinate_map_sha256": "7" * 64,
        "field_transfer_coordinate_map_sha256": "7" * 64,
    }
    return summary


def _with_v17_block_scaling_port_identity(summary: dict) -> dict:
    summary = _with_v16_nonlinear_ale_frame_identity(summary)
    summary["segregated_block_residual_variable_scaling_identity"] = {
        "solver_generation": "segregated-solve-51",
        "block_residual_solver_generation": "segregated-solve-51",
        "variable_scaling_solver_generation": "segregated-solve-51",
        "block_sequence_generation": "block-sequence-51",
        "residual_block_sequence_generation": "block-sequence-51",
        "variable_scaling_block_sequence_generation": "block-sequence-51",
        "block_names": ["magnetic", "thermal"],
        "residual_block_names": ["magnetic", "thermal"],
        "variable_scaling_block_names": ["magnetic", "thermal"],
        "variable_scaling_values": [1.0, 300.0],
        "residual_variable_scaling_values": [1.0, 300.0],
        "variable_scaling_sha256": "1" * 64,
        "residual_variable_scaling_sha256": "1" * 64,
    }
    summary["modal_port_power_normalization_surface_orientation_identity"] = {
        "port_mode_generation": "port-mode-51",
        "modal_amplitude_port_mode_generation": "port-mode-51",
        "power_normalization_port_mode_generation": "port-mode-51",
        "integration_surface_mesh_generation": "port-surface-51",
        "power_normalization_surface_mesh_generation": "port-surface-51",
        "surface_orientation_mesh_generation": "port-surface-51",
        "power_normalization": "unit_forward_power",
        "modal_amplitude_normalization": "unit_forward_power",
        "surface_orientation": "outward_from_domain",
        "power_flux_normal_sign": 1,
        "surface_triangle_sha256": "2" * 64,
        "power_normalization_surface_triangle_sha256": "2" * 64,
    }
    return summary


def _with_v18_eigenmode_and_bdf_identity(summary: dict) -> dict:
    summary = _with_v17_block_scaling_port_identity(summary)
    summary["degenerate_eigenmode_subspace_tracking_basis_identity"] = {
        "eigensolve_generation": "eigensolve-52",
        "eigenvalue_cluster_generation": "eigensolve-52",
        "modal_vector_generation": "eigensolve-52",
        "tracking_basis_generation": "eigensolve-52",
        "mass_inner_product_generation": "mass-inner-product-52",
        "tracking_mass_inner_product_generation": "mass-inner-product-52",
        "cluster_mode_ids": [7, 8],
        "tracking_basis_mode_ids": [7, 8],
        "subspace_dimension": 2,
        "tracking_basis_dimension": 2,
        "modal_assurance_matrix": [[0.999, 0.001], [0.001, 0.999]],
        "tracking_modal_assurance_matrix": [[0.999, 0.001], [0.001, 0.999]],
        "eigenspace_basis_sha256": "1" * 64,
        "tracking_basis_sha256": "1" * 64,
    }
    summary["adaptive_bdf_restart_history_event_generation_identity"] = {
        "transient_generation": "adaptive-transient-52",
        "accepted_step_generation": "accepted-steps-52",
        "solution_history_step_generation": "accepted-steps-52",
        "event_restart_step_generation": "accepted-steps-52",
        "bdf_method": "bdf2",
        "solution_history_method": "bdf2",
        "history_order": 2,
        "solution_history_order": 2,
        "history_time_s": [0.0098, 0.0099, 0.01],
        "solution_history_time_s": [0.0098, 0.0099, 0.01],
        "event_id": "contact-close-52",
        "restart_event_id": "contact-close-52",
        "event_time_s": 0.01,
        "restart_event_time_s": 0.01,
        "history_state_sha256": "2" * 64,
        "restart_history_state_sha256": "2" * 64,
    }
    return summary


def test_v10_public_force_selection_generation_changed() -> None:
    summary = _with_v10_ownership(copy.deepcopy(_summary()))
    summary["force_selection_identity"].update(
        {
            "integration_selection_generation": "geometry-generation-13",
            "selection_entity_digest": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["force_integral_uses_current_geometry_selection"]
        is False
    )


def test_v10_public_excitation_peak_rms_sweep_basis_mismatch() -> None:
    summary = _with_v10_ownership(copy.deepcopy(_summary()))
    summary["excitation_basis_identity"].update(
        {
            "extract_amplitude_basis": "peak",
            "extract_scale_to_rms": 2.0**-0.5,
            "torque_loss_normalization_basis": "compensated_mixed_basis",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["sweep_excitation_uses_one_rms_basis"] is False


def test_v11_public_force_live_stored_geometry_digest_mismatch() -> None:
    summary = _with_v11_artifact_ownership(copy.deepcopy(_summary()))
    summary["live_stored_force_identity"].update(
        {
            "stored_force_geometry_sha256": "5" * 64,
            "stored_force_solution_generation": "solution-generation-42",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["stored_force_matches_live_geometry_solution_and_selection"]
        is False
    )


def test_v11_public_loss_partition_double_counts_interface() -> None:
    summary = _with_v11_artifact_ownership(copy.deepcopy(_summary()))
    summary["loss_partition_identity"].update(
        {"ownership_overlap_count": 1, "signed_compensation_term_w": -0.25}
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "loss_partitions_have_unique_ownership_without_compensation"
        ]
        is False
    )


def test_v12_public_material_property_parameter_unit_generation_mismatch() -> None:
    summary = _with_v12_parameter_topology_identity(copy.deepcopy(_summary()))
    summary["material_property_parameter_identity"].update(
        {
            "material_evaluation_parameter_unit": "mS/m",
            "material_property_parameter_generation": "parameter-row-43",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "material_property_uses_current_parameter_unit_and_generation"
        ]
        is False
    )


def test_v12_public_force_selection_geometry_topology_generation_mismatch() -> None:
    summary = _with_v12_parameter_topology_identity(copy.deepcopy(_summary()))
    summary["force_selection_topology_identity"].update(
        {
            "selection_topology_generation": "topology-43",
            "force_selection_digest": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["force_selection_matches_current_geometry_topology"]
        is False
    )


def test_v13_public_weak_form_coordinate_transform_jacobian_orientation_mismatch() -> None:
    summary = _with_v13_transform_identity(copy.deepcopy(_summary()))
    summary["weak_form_coordinate_transform_identity"].update(
        {"jacobian_coordinate_generation": "mapped-coordinates-44", "jacobian_orientation": "left_handed_negative", "jacobian_orientation_sha256": "d" * 64}
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["weak_form_uses_current_jacobian_orientation"] is False


def test_v13_public_time_harmonic_complex_field_phasor_sign_mismatch() -> None:
    summary = _with_v13_transform_identity(copy.deepcopy(_summary()))
    summary["time_harmonic_phasor_convention_identity"].update(
        {"field_time_convention": "exp(-jwt)", "phase_sensitive_result_time_convention": "exp(-jwt)"}
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["harmonic_fields_share_one_complex_time_convention"] is False


def test_v14_public_eigenmode_mass_inner_product_normalization_generation_mismatch() -> None:
    summary = _with_v14_mode_ale_identity(copy.deepcopy(_summary()))
    summary["eigenmode_mass_inner_product_identity"].update(
        {
            "mass_matrix_mesh_generation": "mesh-45",
            "normalization_mass_matrix_generation": "mass-matrix-45",
            "reference_normalization_kind": "euclidean_l2",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "eigenmodes_use_current_mass_inner_product_normalization"
        ]
        is False
    )


def test_v14_public_ale_mesh_velocity_material_time_level_mismatch() -> None:
    summary = _with_v14_mode_ale_identity(copy.deepcopy(_summary()))
    summary["ale_material_derivative_time_level_identity"].update(
        {
            "mesh_velocity_time_level_generation": "time-level-45",
            "mesh_velocity_time_index": 127,
            "mesh_velocity_time_grid_sha256": "e" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "ale_material_derivative_uses_current_mesh_velocity_time_level"
        ]
        is False
    )


def test_v15_public_positive_reference_time_and_jacobian_identity() -> None:
    result = gate(_with_v15_reference_jacobian_identity(copy.deepcopy(_summary())))
    assert result["status"] == "ok"
    assert result["checks"]["harmonic_field_and_power_share_reference_time_origin"]
    assert result["checks"]["deformed_domain_integral_uses_current_geometry_jacobian"]


def test_v15_public_harmonic_field_power_reference_time_origin_mismatch() -> None:
    summary = _with_v15_reference_jacobian_identity(copy.deepcopy(_summary()))
    summary["harmonic_reference_time_origin_identity"].update(
        {
            "complex_power_generation": "harmonic-solve-46",
            "complex_power_reference_time_s": 0.00025,
            "power_phase_origin_sha256": "6" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["harmonic_field_and_power_share_reference_time_origin"]
        is False
    )


def test_v15_public_deformed_domain_integral_geometry_jacobian_generation_mismatch() -> None:
    summary = _with_v15_reference_jacobian_identity(copy.deepcopy(_summary()))
    summary["deformed_domain_integral_jacobian_identity"].update(
        {
            "volume_jacobian_geometry_generation": "deformed-geometry-46",
            "integral_volume_jacobian_sha256": "6" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["deformed_domain_integral_uses_current_geometry_jacobian"]
        is False
    )


def test_v16_public_positive_nonlinear_and_moving_mesh_frame_identity() -> None:
    result = gate(_with_v16_nonlinear_ale_frame_identity(copy.deepcopy(_summary())))
    assert result["status"] == "ok"
    assert result["checks"]["nonlinear_residual_and_tangent_share_material_iteration"]
    assert result["checks"]["moving_mesh_field_transfer_uses_one_coordinate_frame"]


def test_v16_public_nonlinear_residual_tangent_iteration_generation_mismatch() -> None:
    summary = _with_v16_nonlinear_ale_frame_identity(copy.deepcopy(_summary()))
    summary["nonlinear_residual_tangent_iteration_identity"].update(
        {
            "tangent_solve_generation": "nonlinear-solve-49",
            "tangent_iteration": 11,
            "tangent_material_state_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "nonlinear_residual_and_tangent_share_material_iteration"
        ]
        is False
    )


def test_v16_public_moving_mesh_field_transfer_material_spatial_frame_mismatch() -> None:
    summary = _with_v16_nonlinear_ale_frame_identity(copy.deepcopy(_summary()))
    summary["moving_mesh_field_transfer_frame_identity"].update(
        {
            "target_mesh_motion_generation": "ale-motion-49",
            "target_coordinate_frame": "spatial",
            "field_transfer_coordinate_map_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["moving_mesh_field_transfer_uses_one_coordinate_frame"]
        is False
    )


def test_v17_public_positive_block_scaling_and_port_surface_identity() -> None:
    result = gate(_with_v17_block_scaling_port_identity(copy.deepcopy(_summary())))
    assert result["status"] == "ok"


def test_v17_public_segregated_solver_block_residual_variable_scaling_generation_mismatch() -> None:
    summary = _with_v17_block_scaling_port_identity(copy.deepcopy(_summary()))
    summary["segregated_block_residual_variable_scaling_identity"].update(
        {
            "variable_scaling_solver_generation": "segregated-solve-50",
            "variable_scaling_block_sequence_generation": "block-sequence-50",
            "residual_variable_scaling_values": [1.0, 293.15],
            "residual_variable_scaling_sha256": "5" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "segregated_block_residuals_use_current_variable_scaling"
        ]
        is False
    )


def test_v17_public_modal_port_power_normalization_surface_orientation_generation_mismatch() -> None:
    summary = _with_v17_block_scaling_port_identity(copy.deepcopy(_summary()))
    summary["modal_port_power_normalization_surface_orientation_identity"].update(
        {
            "power_normalization_surface_mesh_generation": "port-surface-50",
            "surface_orientation_mesh_generation": "port-surface-50",
            "surface_orientation": "inward_to_domain",
            "power_flux_normal_sign": -1,
            "power_normalization_surface_triangle_sha256": "5" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["modal_port_power_uses_current_surface_orientation"]
        is False
    )


def test_v18_public_positive_eigenmode_and_bdf_lineage() -> None:
    result = gate(_with_v18_eigenmode_and_bdf_identity(copy.deepcopy(_summary())))
    assert result["status"] == "ok"
    assert result["checks"][
        "degenerate_eigenmodes_use_current_subspace_tracking_basis"
    ]
    assert result["checks"][
        "adaptive_bdf_restart_uses_current_history_and_event_generation"
    ]


def test_v18_public_degenerate_eigenmode_subspace_tracking_basis_generation_mismatch() -> None:
    summary = _with_v18_eigenmode_and_bdf_identity(copy.deepcopy(_summary()))
    summary["degenerate_eigenmode_subspace_tracking_basis_identity"].update(
        {
            "tracking_basis_generation": "eigensolve-51",
            "tracking_mass_inner_product_generation": "mass-inner-product-51",
            "tracking_basis_mode_ids": [8, 7],
            "tracking_modal_assurance_matrix": [[0.001, 0.999], [0.999, 0.001]],
            "tracking_basis_sha256": "5" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "degenerate_eigenmodes_use_current_subspace_tracking_basis"
        ]
        is False
    )


def test_v18_public_adaptive_bdf_restart_history_event_generation_mismatch() -> None:
    summary = _with_v18_eigenmode_and_bdf_identity(copy.deepcopy(_summary()))
    summary["adaptive_bdf_restart_history_event_generation_identity"].update(
        {
            "solution_history_step_generation": "accepted-steps-51",
            "event_restart_step_generation": "accepted-steps-51",
            "solution_history_method": "bdf1",
            "solution_history_order": 1,
            "solution_history_time_s": [0.0099, 0.01],
            "restart_event_id": "contact-close-51",
            "restart_history_state_sha256": "5" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "adaptive_bdf_restart_uses_current_history_and_event_generation"
        ]
        is False
    )


def _with_v19_continuation_and_mortar_identity(summary: dict) -> dict:
    summary = _with_v18_eigenmode_and_bdf_identity(summary)
    summary["nonlinear_continuation_branch_tangent_checkpoint_identity"] = {
        "nonlinear_solve_generation": "nonlinear-solve-53",
        "continuation_state_generation": "continuation-state-53",
        "tangent_continuation_state_generation": "continuation-state-53",
        "checkpoint_continuation_state_generation": "continuation-state-53",
        "branch_id": "upper-stable-53",
        "tangent_branch_id": "upper-stable-53",
        "checkpoint_branch_id": "upper-stable-53",
        "load_parameter_name": "coil_current_A",
        "load_parameter_value": 12.5,
        "tangent_load_parameter_value": 12.5,
        "checkpoint_load_parameter_value": 12.5,
        "tangent_vector": [0.2, -0.1, 0.05],
        "checkpoint_tangent_vector": [0.2, -0.1, 0.05],
        "continuation_tangent_sha256": "5" * 64,
        "checkpoint_tangent_sha256": "5" * 64,
    }
    summary["nonconforming_mortar_projection_quadrature_mesh_identity"] = {
        "interface_generation": "mortar-interface-53",
        "source_mesh_generation": "source-mesh-53",
        "projection_source_mesh_generation": "source-mesh-53",
        "quadrature_source_mesh_generation": "source-mesh-53",
        "target_mesh_generation": "target-mesh-53",
        "projection_target_mesh_generation": "target-mesh-53",
        "quadrature_target_mesh_generation": "target-mesh-53",
        "source_trace_dof_ids": [11, 12],
        "projection_source_trace_dof_ids": [11, 12],
        "target_trace_dof_ids": [21, 22, 23],
        "projection_target_trace_dof_ids": [21, 22, 23],
        "projection_shape": [3, 2],
        "quadrature_projection_shape": [3, 2],
        "projection_operator_sha256": "6" * 64,
        "quadrature_projection_operator_sha256": "6" * 64,
    }
    return summary


def test_v19_public_positive_continuation_and_mortar_identity() -> None:
    result = gate(_with_v19_continuation_and_mortar_identity(copy.deepcopy(_summary())))
    assert result["status"] == "ok"


def test_v19_public_nonlinear_continuation_branch_tangent_load_parameter_checkpoint_mismatch() -> None:
    summary = _with_v19_continuation_and_mortar_identity(copy.deepcopy(_summary()))
    summary["nonlinear_continuation_branch_tangent_checkpoint_identity"].update(
        {
            "tangent_continuation_state_generation": "continuation-state-52",
            "checkpoint_continuation_state_generation": "continuation-state-52",
            "checkpoint_branch_id": "lower-unstable-52",
            "tangent_load_parameter_value": 12.0,
            "checkpoint_tangent_vector": [-0.2, 0.1, -0.05],
            "checkpoint_tangent_sha256": "9" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "nonlinear_continuation_uses_current_branch_tangent_checkpoint"
        ]
        is False
    )


def test_v19_public_nonconforming_mortar_interface_projection_quadrature_mesh_mismatch() -> None:
    summary = _with_v19_continuation_and_mortar_identity(copy.deepcopy(_summary()))
    summary["nonconforming_mortar_projection_quadrature_mesh_identity"].update(
        {
            "projection_source_mesh_generation": "source-mesh-52",
            "quadrature_target_mesh_generation": "target-mesh-52",
            "projection_source_trace_dof_ids": [12, 11],
            "quadrature_projection_shape": [2, 3],
            "quadrature_projection_operator_sha256": "9" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "mortar_projection_uses_current_interface_mesh_and_quadrature"
        ]
        is False
    )


def _with_v20_transfer_and_mode_tracking_identity(summary: dict) -> dict:
    summary = _with_v19_continuation_and_mortar_identity(summary)
    summary["adaptive_mesh_field_transfer_projection_conservation_identity"] = {
        "solve_generation": "adaptive-solve-62",
        "source_mesh_generation": "adaptive-mesh-61",
        "projection_source_mesh_generation": "adaptive-mesh-61",
        "conservation_source_mesh_generation": "adaptive-mesh-61",
        "target_mesh_generation": "adaptive-mesh-62",
        "projection_target_mesh_generation": "adaptive-mesh-62",
        "conservation_target_mesh_generation": "adaptive-mesh-62",
        "source_field_values": [1.0, 2.0],
        "source_integration_weights": [0.5, 0.5],
        "projected_field_values": [0.75, 1.5, 2.25],
        "target_integration_weights": [1.0 / 3.0] * 3,
        "source_conserved_integral": 1.5,
        "target_conserved_integral": 1.5,
        "projection_shape": [3, 2],
        "projection_operator_sha256": "a" * 64,
        "conservation_weight_table_sha256": "b" * 64,
        "transfer_conservation_weight_table_sha256": "b" * 64,
    }
    summary["eigenmode_phase_normalization_tracking_parameter_identity"] = {
        "parameter_table_generation": "eigen-parameter-62",
        "previous_eigensolve_generation": "eigensolve-61",
        "current_eigensolve_generation": "eigensolve-62",
        "tracker_current_eigensolve_generation": "eigensolve-62",
        "phase_anchor_current_eigensolve_generation": "eigensolve-62",
        "normalization_current_eigensolve_generation": "eigensolve-62",
        "parameter_name": "rotor_angle_deg",
        "previous_parameter_value": 10.0,
        "current_parameter_value": 12.0,
        "tracked_mode_ids": [3, 4],
        "tracker_mode_ids": [3, 4],
        "phase_anchor_dof_ids": [101, 205],
        "tracker_phase_anchor_dof_ids": [101, 205],
        "normalization_integrals": [1.0, 1.0],
        "tracker_normalization_integrals": [1.0, 1.0],
        "selected_correlation": [0.98, 0.96],
        "tracker_selected_correlation": [0.98, 0.96],
        "mode_tracking_table_sha256": "c" * 64,
        "tracker_mode_tracking_table_sha256": "c" * 64,
    }
    return summary


def test_v20_public_positive_transfer_and_mode_tracking_identity() -> None:
    result = gate(
        _with_v20_transfer_and_mode_tracking_identity(copy.deepcopy(_summary()))
    )
    assert result["status"] == "ok"


def test_v20_public_adaptive_mesh_field_transfer_projection_conservation_generation_mismatch() -> None:
    summary = _with_v20_transfer_and_mode_tracking_identity(copy.deepcopy(_summary()))
    summary["adaptive_mesh_field_transfer_projection_conservation_identity"].update(
        {
            "projection_source_mesh_generation": "adaptive-mesh-60",
            "conservation_target_mesh_generation": "adaptive-mesh-61",
            "target_integration_weights": [0.25, 0.25, 0.25],
            "target_conserved_integral": 1.125,
            "transfer_conservation_weight_table_sha256": "f" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "adaptive_field_transfer_uses_current_mesh_projection_and_conservation"
        ]
        is False
    )


def test_v20_public_eigenmode_phase_normalization_tracking_parameter_generation_mismatch() -> None:
    summary = _with_v20_transfer_and_mode_tracking_identity(copy.deepcopy(_summary()))
    summary["eigenmode_phase_normalization_tracking_parameter_identity"].update(
        {
            "tracker_current_eigensolve_generation": "eigensolve-61",
            "phase_anchor_current_eigensolve_generation": "eigensolve-61",
            "tracker_mode_ids": [4, 3],
            "tracker_phase_anchor_dof_ids": [205, 101],
            "tracker_normalization_integrals": [0.5, 2.0],
            "tracker_selected_correlation": [0.65, 0.71],
            "tracker_mode_tracking_table_sha256": "f" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "eigenmode_tracking_uses_current_phase_normalization_and_parameter_state"
        ]
        is False
    )


def _with_v21_restart_and_coupling_identity(summary: dict) -> dict:
    summary = _with_v20_transfer_and_mode_tracking_identity(summary)
    summary["parameter_sweep_branch_restart_solution_mesh_generation_identity"] = {
        "parameter_sweep_generation": "parameter-sweep-71",
        "current_branch_generation": "sweep-branch-71",
        "checkpoint_parameter_sweep_generation": "parameter-sweep-71",
        "checkpoint_branch_generation": "sweep-branch-71",
        "solution_vector_branch_generation": "sweep-branch-71",
        "continuation_state_branch_generation": "sweep-branch-71",
        "mesh_coordinate_branch_generation": "sweep-branch-71",
        "current_branch_id": 5,
        "checkpoint_branch_id": 5,
        "parameter_names": ["current_a", "rotor_angle_deg"],
        "current_parameter_tuple": [12.0, 18.0],
        "checkpoint_parameter_tuple": [12.0, 18.0],
        "solution_vector_sha256": "1" * 64,
        "checkpoint_solution_vector_sha256": "1" * 64,
        "continuation_state_sha256": "2" * 64,
        "checkpoint_continuation_state_sha256": "2" * 64,
        "mesh_coordinate_sha256": "3" * 64,
        "checkpoint_mesh_coordinate_sha256": "3" * 64,
    }
    summary[
        "multiphysics_coupling_source_frame_unit_selection_generation_identity"
    ] = {
        "coupling_generation": "multiphysics-coupling-71",
        "source_values_coupling_generation": "multiphysics-coupling-71",
        "source_frame_coupling_generation": "multiphysics-coupling-71",
        "source_unit_coupling_generation": "multiphysics-coupling-71",
        "source_selection_coupling_generation": "multiphysics-coupling-71",
        "source_coordinate_frame": "spatial",
        "assembled_coordinate_frame": "spatial",
        "source_units": ["N/m^3", "W/m^3"],
        "assembled_source_units": ["N/m^3", "W/m^3"],
        "source_selection_ids": [21, 22],
        "assembled_source_selection_ids": [21, 22],
        "source_selection_dimensions": [3, 3],
        "assembled_source_selection_dimensions": [3, 3],
        "source_values_sha256": "4" * 64,
        "assembled_source_values_sha256": "4" * 64,
        "source_selection_sha256": "5" * 64,
        "assembled_source_selection_sha256": "5" * 64,
    }
    return summary


def test_v21_public_positive_restart_and_coupling_identity() -> None:
    result = gate(_with_v21_restart_and_coupling_identity(copy.deepcopy(_summary())))
    assert result["status"] == "ok"


def test_v21_public_parameter_sweep_branch_restart_solution_mesh_generation_mismatch() -> None:
    summary = _with_v21_restart_and_coupling_identity(copy.deepcopy(_summary()))
    summary[
        "parameter_sweep_branch_restart_solution_mesh_generation_identity"
    ].update(
        {
            "checkpoint_branch_generation": "sweep-branch-70",
            "solution_vector_branch_generation": "sweep-branch-70",
            "mesh_coordinate_branch_generation": "sweep-branch-69",
            "checkpoint_branch_id": 4,
            "checkpoint_parameter_tuple": [10.0, 12.0],
            "checkpoint_solution_vector_sha256": "9" * 64,
            "checkpoint_mesh_coordinate_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "parameter_sweep_restart_uses_current_branch_solution_and_mesh"
    ]


def test_v21_public_multiphysics_coupling_source_frame_unit_selection_generation_mismatch() -> None:
    summary = _with_v21_restart_and_coupling_identity(copy.deepcopy(_summary()))
    summary[
        "multiphysics_coupling_source_frame_unit_selection_generation_identity"
    ].update(
        {
            "source_frame_coupling_generation": "multiphysics-coupling-70",
            "source_unit_coupling_generation": "multiphysics-coupling-69",
            "source_selection_coupling_generation": "multiphysics-coupling-68",
            "assembled_coordinate_frame": "material",
            "assembled_source_units": ["N/mm^3", "W/mm^3"],
            "assembled_source_selection_ids": [22, 24],
            "assembled_source_values_sha256": "b" * 64,
            "assembled_source_selection_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "multiphysics_coupling_uses_current_source_frame_units_and_selection"
    ]


def _with_v22_contact_and_acoustic_structure_identity(summary: dict) -> dict:
    summary = _with_v21_restart_and_coupling_identity(summary)
    summary["contact_active_set_friction_state_mesh_generation_identity"] = {
        "contact_generation": "contact-72",
        "active_set_contact_generation": "contact-72",
        "friction_state_contact_generation": "contact-72",
        "slave_mesh_contact_generation": "contact-72",
        "normal_orientation_contact_generation": "contact-72",
        "consistent_tangent_contact_generation": "contact-72",
        "slave_surface_ids": [31, 32],
        "active_slave_surface_ids": [31, 32],
        "normal_orientation": "slave_to_master",
        "active_set_normal_orientation": "slave_to_master",
        "friction_coefficient": 0.18,
        "active_set_friction_coefficient": 0.18,
        "slave_mesh_sha256": "1" * 64,
        "active_set_slave_mesh_sha256": "1" * 64,
        "friction_state_sha256": "2" * 64,
        "active_set_friction_state_sha256": "2" * 64,
        "consistent_tangent_sha256": "3" * 64,
        "active_set_consistent_tangent_sha256": "3" * 64,
    }
    summary["acoustic_structure_trace_impedance_order_frame_generation_identity"] = {
        "coupling_generation": "acoustic-structure-72",
        "pressure_trace_coupling_generation": "acoustic-structure-72",
        "traction_trace_coupling_generation": "acoustic-structure-72",
        "normal_frame_coupling_generation": "acoustic-structure-72",
        "impedance_coupling_generation": "acoustic-structure-72",
        "interface_mesh_coupling_generation": "acoustic-structure-72",
        "normal_frame": "acoustic_outward",
        "traction_normal_frame": "acoustic_outward",
        "pressure_to_traction_sign": -1,
        "assembled_pressure_to_traction_sign": -1,
        "impedance_order": 4,
        "assembled_impedance_order": 4,
        "pressure_trace_sha256": "4" * 64,
        "assembled_pressure_trace_sha256": "4" * 64,
        "traction_trace_sha256": "5" * 64,
        "assembled_traction_trace_sha256": "5" * 64,
        "interface_mesh_sha256": "6" * 64,
        "assembled_interface_mesh_sha256": "6" * 64,
    }
    return summary


def test_v22_public_positive_contact_and_acoustic_structure_identity() -> None:
    result = gate(
        _with_v22_contact_and_acoustic_structure_identity(copy.deepcopy(_summary()))
    )
    assert result["status"] == "ok"


def test_v22_public_contact_active_set_friction_state_mesh_generation_mismatch() -> None:
    summary = _with_v22_contact_and_acoustic_structure_identity(
        copy.deepcopy(_summary())
    )
    summary["contact_active_set_friction_state_mesh_generation_identity"].update(
        {
            "active_set_contact_generation": "contact-71",
            "friction_state_contact_generation": "contact-70",
            "slave_mesh_contact_generation": "contact-69",
            "active_slave_surface_ids": [31, 34],
            "active_set_normal_orientation": "master_to_slave",
            "active_set_friction_coefficient": 0.24,
            "active_set_slave_mesh_sha256": "b" * 64,
            "active_set_friction_state_sha256": "c" * 64,
            "active_set_consistent_tangent_sha256": "d" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "contact_active_set_uses_current_friction_state_mesh_and_tangent"
    ]


def test_v22_public_acoustic_structure_trace_impedance_order_frame_generation_mismatch() -> None:
    summary = _with_v22_contact_and_acoustic_structure_identity(
        copy.deepcopy(_summary())
    )
    summary[
        "acoustic_structure_trace_impedance_order_frame_generation_identity"
    ].update(
        {
            "traction_trace_coupling_generation": "acoustic-structure-71",
            "normal_frame_coupling_generation": "acoustic-structure-70",
            "impedance_coupling_generation": "acoustic-structure-69",
            "interface_mesh_coupling_generation": "acoustic-structure-68",
            "traction_normal_frame": "structure_outward",
            "assembled_pressure_to_traction_sign": 1,
            "assembled_impedance_order": 1,
            "assembled_pressure_trace_sha256": "e" * 64,
            "assembled_traction_trace_sha256": "f" * 64,
            "assembled_interface_mesh_sha256": "0" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "acoustic_structure_trace_uses_current_frame_impedance_and_interface"
    ]


def _with_v23_continuation_and_sweep_identity(summary: dict) -> dict:
    summary = _with_v22_contact_and_acoustic_structure_identity(summary)
    summary["nonlinear_continuation_branch_load_step_mesh_generation_identity"] = {
        "solve_generation": "continuation-81",
        "branch_solve_generation": "continuation-81",
        "load_step_solve_generation": "continuation-81",
        "tangent_state_solve_generation": "continuation-81",
        "adapted_mesh_solve_generation": "continuation-81",
        "result_solve_generation": "continuation-81",
        "branch_id": "stable-positive",
        "result_branch_id": "stable-positive",
        "load_parameters": [0.0, 0.25, 0.5, 0.75, 1.0],
        "result_load_parameters": [0.0, 0.25, 0.5, 0.75, 1.0],
        "tangent_state_sha256": ["1" * 64, "2" * 64, "3" * 64, "4" * 64, "5" * 64],
        "result_tangent_state_sha256": ["1" * 64, "2" * 64, "3" * 64, "4" * 64, "5" * 64],
        "adapted_mesh_sha256": "6" * 64,
        "result_mesh_sha256": "6" * 64,
        "continuation_table_sha256": "7" * 64,
        "result_continuation_table_sha256": "7" * 64,
    }
    summary["parametric_sequence_initial_solution_dataset_generation_identity"] = {
        "sweep_generation": "parametric-81",
        "sequence_sweep_generation": "parametric-81",
        "parameter_row_sweep_generation": "parametric-81",
        "initial_solution_sweep_generation": "parametric-81",
        "dataset_sweep_generation": "parametric-81",
        "result_sweep_generation": "parametric-81",
        "sequence_id": "continuation-from-previous-row",
        "result_sequence_id": "continuation-from-previous-row",
        "parameter_names": ["current_a", "speed_rpm"],
        "parameter_rows": [[0.0, 0.0], [5.0, 500.0], [10.0, 1000.0]],
        "result_parameter_rows": [[0.0, 0.0], [5.0, 500.0], [10.0, 1000.0]],
        "initial_solution_sha256": ["8" * 64, "9" * 64, "a" * 64],
        "result_initial_solution_sha256": ["8" * 64, "9" * 64, "a" * 64],
        "dataset_sha256": "b" * 64,
        "result_dataset_sha256": "b" * 64,
    }
    return summary


def test_v23_public_positive_continuation_and_parametric_sequence_identity() -> None:
    assert gate(_with_v23_continuation_and_sweep_identity(copy.deepcopy(_summary())))[
        "status"
    ] == "ok"


def test_v23_public_nonlinear_continuation_branch_load_step_mesh_generation_mismatch() -> None:
    summary = _with_v23_continuation_and_sweep_identity(copy.deepcopy(_summary()))
    summary["nonlinear_continuation_branch_load_step_mesh_generation_identity"].update(
        {
            "branch_solve_generation": "continuation-80",
            "load_step_solve_generation": "continuation-79",
            "tangent_state_solve_generation": "continuation-78",
            "adapted_mesh_solve_generation": "continuation-77",
            "result_branch_id": "unstable-negative",
            "result_load_parameters": [0.0, 0.5, 0.25, 0.75, 1.0],
            "result_tangent_state_sha256": ["5" * 64, "4" * 64, "3" * 64, "2" * 64, "1" * 64],
            "result_mesh_sha256": "1" * 64,
            "result_continuation_table_sha256": "2" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_continuation_result_uses_current_branch_load_tangent_and_mesh"
    ]


def test_v23_public_parametric_sweep_sequence_initial_solution_result_generation_mismatch() -> None:
    summary = _with_v23_continuation_and_sweep_identity(copy.deepcopy(_summary()))
    summary["parametric_sequence_initial_solution_dataset_generation_identity"].update(
        {
            "sequence_sweep_generation": "parametric-80",
            "parameter_row_sweep_generation": "parametric-79",
            "initial_solution_sweep_generation": "parametric-78",
            "dataset_sweep_generation": "parametric-77",
            "result_sequence_id": "independent-rows",
            "result_parameter_rows": [[10.0, 1000.0], [5.0, 500.0], [0.0, 0.0]],
            "result_initial_solution_sha256": ["a" * 64, "9" * 64, "8" * 64],
            "result_dataset_sha256": "3" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "parametric_sequence_uses_current_rows_initial_solutions_and_dataset"
    ]


def _with_v24_balance_and_degenerate_mode_identity(summary: dict) -> dict:
    summary["multiphysics_power_work_heat_balance_frame_time_generation_identity"] = {
        "balance_generation": "balance-91",
        "electromagnetic_balance_generation": "balance-91",
        "mechanical_balance_generation": "balance-91",
        "thermal_balance_generation": "balance-91",
        "stored_energy_balance_generation": "balance-91",
        "result_balance_generation": "balance-91",
        "coordinate_frame": "global-xyz",
        "result_coordinate_frame": "global-xyz",
        "time_window_s": [0.0, 0.02],
        "result_time_window_s": [0.0, 0.02],
        "electromagnetic_input_j": 120.0,
        "mechanical_work_j": 45.0,
        "heat_source_j": 60.0,
        "stored_energy_change_j": 15.0,
        "reported_balance_j": 0.0,
        "balance_table_sha256": "1" * 64,
        "result_balance_table_sha256": "1" * 64,
    }
    summary[
        "degenerate_eigenmode_subspace_normalization_phase_projection_identity"
    ] = {
        "eigensolve_generation": "eigensolve-91",
        "subspace_eigensolve_generation": "eigensolve-91",
        "normalization_eigensolve_generation": "eigensolve-91",
        "phase_eigensolve_generation": "eigensolve-91",
        "projection_eigensolve_generation": "eigensolve-91",
        "mesh_eigensolve_generation": "eigensolve-91",
        "result_eigensolve_generation": "eigensolve-91",
        "eigenvalues_hz": [1000.0, 1000.0001],
        "result_eigenvalues_hz": [1000.0, 1000.0001],
        "normalization": "unit-energy",
        "result_normalization": "unit-energy",
        "complex_phase_rad": [0.0, 0.0],
        "result_complex_phase_rad": [0.0, 0.0],
        "subspace_projection": [[1.0, 0.0], [0.0, 1.0]],
        "result_subspace_projection": [[1.0, 0.0], [0.0, 1.0]],
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "modal_table_sha256": "3" * 64,
        "result_modal_table_sha256": "3" * 64,
    }
    return summary


def test_v24_public_positive_balance_and_degenerate_mode_identity() -> None:
    result = gate(_with_v24_balance_and_degenerate_mode_identity(_summary()))
    assert result["status"] == "ok"


def test_v24_public_multiphysics_power_work_heat_balance_frame_time_generation_mismatch() -> None:
    summary = _with_v24_balance_and_degenerate_mode_identity(_summary())
    summary[
        "multiphysics_power_work_heat_balance_frame_time_generation_identity"
    ].update(
        {
            "mechanical_balance_generation": "balance-90",
            "thermal_balance_generation": "balance-89",
            "stored_energy_balance_generation": "balance-88",
            "result_coordinate_frame": "rotor-local",
            "result_time_window_s": [0.005, 0.025],
            "reported_balance_j": 12.0,
            "result_balance_table_sha256": "8" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "multiphysics_energy_balance_uses_current_frame_time_and_generation"
    ]


def test_v24_public_degenerate_eigenmode_subspace_normalization_phase_projection_generation_mismatch() -> None:
    summary = _with_v24_balance_and_degenerate_mode_identity(_summary())
    summary[
        "degenerate_eigenmode_subspace_normalization_phase_projection_identity"
    ].update(
        {
            "subspace_eigensolve_generation": "eigensolve-90",
            "normalization_eigensolve_generation": "eigensolve-89",
            "phase_eigensolve_generation": "eigensolve-88",
            "projection_eigensolve_generation": "eigensolve-87",
            "result_eigenvalues_hz": [1000.0001, 1000.0],
            "result_normalization": "unit-maximum",
            "result_complex_phase_rad": [1.5707963268, 0.0],
            "result_subspace_projection": [[0.0, 1.0], [1.0, 0.0]],
            "result_mesh_sha256": "9" * 64,
            "result_modal_table_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "degenerate_eigenmodes_use_current_subspace_normalization_phase_and_mesh"
    ]


def _with_v25_remesh_and_continuation_identity(summary: dict) -> dict:
    summary["remesh_field_projection_conservation_geometry_dataset_generation_identity"] = {
        "projection_generation": "remesh-111",
        "source_mesh_projection_generation": "remesh-111",
        "target_mesh_projection_generation": "remesh-111",
        "geometry_projection_generation": "remesh-111",
        "dataset_projection_generation": "remesh-111",
        "integral_projection_generation": "remesh-111",
        "result_projection_generation": "remesh-111",
        "source_mesh_sha256": "1" * 64,
        "projected_source_mesh_sha256": "1" * 64,
        "target_mesh_sha256": "2" * 64,
        "result_target_mesh_sha256": "2" * 64,
        "geometry_revision_sha256": "3" * 64,
        "result_geometry_revision_sha256": "3" * 64,
        "dataset_tag": "dset-remesh-2",
        "result_dataset_tag": "dset-remesh-2",
        "conserved_integral_before": 18.25,
        "conserved_integral_after": 18.25,
        "projection_map_sha256": "4" * 64,
        "result_projection_map_sha256": "4" * 64,
        "projected_field_sha256": "5" * 64,
        "result_projected_field_sha256": "5" * 64,
    }
    summary["nonlinear_continuation_load_step_branch_state_solver_generation_identity"] = {
        "continuation_generation": "continuation-111",
        "load_step_continuation_generation": "continuation-111",
        "branch_continuation_generation": "continuation-111",
        "initial_state_continuation_generation": "continuation-111",
        "solver_continuation_generation": "continuation-111",
        "result_continuation_generation": "continuation-111",
        "load_step_ids": [1, 2, 3, 4],
        "result_load_step_ids": [1, 2, 3, 4],
        "continuation_parameter": [0.25, 0.5, 0.75, 1.0],
        "result_continuation_parameter": [0.25, 0.5, 0.75, 1.0],
        "branch_id": "stable-positive-slope",
        "result_branch_id": "stable-positive-slope",
        "initial_state_sha256": "6" * 64,
        "result_initial_state_sha256": "6" * 64,
        "solver_settings_sha256": "7" * 64,
        "result_solver_settings_sha256": "7" * 64,
        "solution_table_sha256": "8" * 64,
        "result_solution_table_sha256": "8" * 64,
    }
    return summary


def test_v25_public_positive_remesh_and_continuation_identity() -> None:
    result = gate(_with_v25_remesh_and_continuation_identity(_summary()))
    assert result["status"] == "ok"


def test_v25_public_remesh_field_projection_conservation_geometry_dataset_generation_mismatch() -> None:
    summary = _with_v25_remesh_and_continuation_identity(_summary())
    summary["remesh_field_projection_conservation_geometry_dataset_generation_identity"].update(
        {
            "source_mesh_projection_generation": "remesh-110",
            "target_mesh_projection_generation": "remesh-109",
            "geometry_projection_generation": "remesh-108",
            "dataset_projection_generation": "remesh-107",
            "result_target_mesh_sha256": "c" * 64,
            "result_geometry_revision_sha256": "d" * 64,
            "result_dataset_tag": "dset-remesh-1",
            "conserved_integral_after": 16.75,
            "result_projection_map_sha256": "e" * 64,
            "result_projected_field_sha256": "f" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "remeshed_fields_use_current_projection_geometry_dataset_and_conservation"
    ]


def test_v25_public_nonlinear_continuation_load_step_branch_state_solver_generation_mismatch() -> None:
    summary = _with_v25_remesh_and_continuation_identity(_summary())
    summary["nonlinear_continuation_load_step_branch_state_solver_generation_identity"].update(
        {
            "load_step_continuation_generation": "continuation-110",
            "branch_continuation_generation": "continuation-109",
            "initial_state_continuation_generation": "continuation-108",
            "solver_continuation_generation": "continuation-107",
            "result_load_step_ids": [1, 2, 4, 5],
            "result_continuation_parameter": [0.25, 0.5, 0.9, 1.1],
            "result_branch_id": "unstable-negative-slope",
            "result_initial_state_sha256": "0" * 64,
            "result_solver_settings_sha256": "1" * 64,
            "result_solution_table_sha256": "2" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_solutions_use_current_load_steps_branch_state_and_solver"
    ]


def _with_v26_ale_and_segregated_identity(summary: dict) -> dict:
    summary["ale_moving_mesh_time_step_field_transfer_force_work_balance_generation_identity"] = {
        "ale_generation": "ale-work-131",
        "geometry_ale_generation": "ale-work-131",
        "time_step_ale_generation": "ale-work-131",
        "field_transfer_ale_generation": "ale-work-131",
        "force_ale_generation": "ale-work-131",
        "work_ale_generation": "ale-work-131",
        "result_ale_generation": "ale-work-131",
        "time_s": [0.0, 0.001, 0.002, 0.003],
        "result_time_s": [0.0, 0.001, 0.002, 0.003],
        "mesh_displacement_m": [0.01, 0.01, 0.01],
        "result_mesh_displacement_m": [0.01, 0.01, 0.01],
        "force_n": [10.0, 12.0, 14.0],
        "result_force_n": [10.0, 12.0, 14.0],
        "reported_mechanical_work_j": 0.36,
        "field_energy_change_j": 0.20,
        "dissipated_energy_j": 0.16,
        "geometry_mesh_sha256": "1" * 64,
        "result_geometry_mesh_sha256": "1" * 64,
        "field_transfer_sha256": "2" * 64,
        "result_field_transfer_sha256": "2" * 64,
        "force_work_table_sha256": "3" * 64,
        "result_force_work_table_sha256": "3" * 64,
    }
    summary["segregated_multiphysics_iteration_relaxation_residual_component_solution_generation_identity"] = {
        "solve_generation": "segregated-131",
        "iteration_solve_generation": "segregated-131",
        "relaxation_solve_generation": "segregated-131",
        "residual_solve_generation": "segregated-131",
        "component_solve_generation": "segregated-131",
        "solution_solve_generation": "segregated-131",
        "result_solve_generation": "segregated-131",
        "iteration_ids": [1, 2, 3, 4],
        "result_iteration_ids": [1, 2, 3, 4],
        "component_order": ["magnetic", "thermal", "mechanical"],
        "result_component_order": ["magnetic", "thermal", "mechanical"],
        "relaxation_factors": [0.7, 0.8, 1.0],
        "result_relaxation_factors": [0.7, 0.8, 1.0],
        "residual_norm": "l2",
        "result_residual_norm": "l2",
        "residual_history": [1.0, 0.1, 0.005, 0.00001],
        "result_residual_history": [1.0, 0.1, 0.005, 0.00001],
        "relative_tolerance": 0.0001,
        "converged": True,
        "solution_sha256": "4" * 64,
        "result_solution_sha256": "4" * 64,
    }
    return summary


def test_v26_public_positive_ale_and_segregated_identity() -> None:
    assert gate(_with_v26_ale_and_segregated_identity(_summary()))["status"] == "ok"


def test_v26_public_ale_moving_mesh_time_step_field_transfer_force_work_balance_generation_mismatch() -> None:
    summary = _with_v26_ale_and_segregated_identity(_summary())
    summary["ale_moving_mesh_time_step_field_transfer_force_work_balance_generation_identity"].update(
        {
            "geometry_ale_generation": "ale-work-130",
            "result_time_s": [0.0, 0.001, 0.0025, 0.003],
            "result_force_n": [10.0, 8.0, 14.0],
            "reported_mechanical_work_j": 0.9,
            "result_field_transfer_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ale_force_work_uses_current_geometry_time_transfer_and_energy_balance"
    ]


def test_v26_public_segregated_multiphysics_iteration_relaxation_residual_component_solution_mismatch() -> None:
    summary = _with_v26_ale_and_segregated_identity(_summary())
    summary["segregated_multiphysics_iteration_relaxation_residual_component_solution_generation_identity"].update(
        {
            "iteration_solve_generation": "segregated-130",
            "result_component_order": ["thermal", "magnetic", "mechanical"],
            "result_residual_norm": "linf",
            "result_residual_history": [1.0, 0.5, 0.8, 0.2],
            "converged": False,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "segregated_solution_uses_current_iterations_relaxation_residuals_and_components"
    ]


def _with_v27_restart_and_floquet_identity(summary: dict) -> dict:
    summary["nonlinear_state_time_integrator_tangent_load_step_restart_generation_identity"] = {
        "nonlinear_generation": "nonlinear-restart-141",
        "state_nonlinear_generation": "nonlinear-restart-141",
        "integrator_nonlinear_generation": "nonlinear-restart-141",
        "tangent_nonlinear_generation": "nonlinear-restart-141",
        "load_step_nonlinear_generation": "nonlinear-restart-141",
        "checkpoint_nonlinear_generation": "nonlinear-restart-141",
        "result_nonlinear_generation": "nonlinear-restart-141",
        "state_variable_names": ["plastic_strain", "hardening_variable"],
        "result_state_variable_names": ["plastic_strain", "hardening_variable"],
        "time_integrator": "generalized_alpha",
        "result_time_integrator": "generalized_alpha",
        "integrator_order": 2,
        "result_integrator_order": 2,
        "load_step_ids": [1, 2, 3, 4],
        "result_load_step_ids": [1, 2, 3, 4],
        "restart_time_s": 0.003,
        "result_restart_time_s": 0.003,
        "state_vector_sha256": "1" * 64,
        "result_state_vector_sha256": "1" * 64,
        "consistent_tangent_sha256": "2" * 64,
        "result_consistent_tangent_sha256": "2" * 64,
        "checkpoint_sha256": "3" * 64,
        "result_checkpoint_sha256": "3" * 64,
        "solution_sha256": "4" * 64,
        "result_solution_sha256": "4" * 64,
    }
    summary["floquet_pair_orientation_phase_wavevector_normalization_dataset_mesh_generation_identity"] = {
        "floquet_generation": "floquet-141",
        "pair_floquet_generation": "floquet-141",
        "orientation_floquet_generation": "floquet-141",
        "phase_floquet_generation": "floquet-141",
        "wavevector_floquet_generation": "floquet-141",
        "normalization_floquet_generation": "floquet-141",
        "dataset_floquet_generation": "floquet-141",
        "mesh_floquet_generation": "floquet-141",
        "result_floquet_generation": "floquet-141",
        "periodic_pair_tags": ["pair-x", "pair-y"],
        "result_periodic_pair_tags": ["pair-x", "pair-y"],
        "pair_orientation_signs": [1, -1],
        "result_pair_orientation_signs": [1, -1],
        "phase_shift_rad": [0.2, -0.1],
        "result_phase_shift_rad": [0.2, -0.1],
        "wavevector_rad_m": [20.0, -10.0, 0.0],
        "result_wavevector_rad_m": [20.0, -10.0, 0.0],
        "mode_normalization": "unit_cell_energy_1j",
        "result_mode_normalization": "unit_cell_energy_1j",
        "dataset_tag": "dset-floquet-2",
        "result_dataset_tag": "dset-floquet-2",
        "periodic_mesh_map_sha256": "5" * 64,
        "result_periodic_mesh_map_sha256": "5" * 64,
        "mode_field_sha256": "6" * 64,
        "result_mode_field_sha256": "6" * 64,
    }
    return summary


def test_v27_public_positive_restart_and_floquet_identity() -> None:
    result = gate(_with_v27_restart_and_floquet_identity(_summary()))
    assert result["status"] == "ok"


def test_v27_public_nonlinear_state_variable_time_integrator_consistent_tangent_load_step_restart_mismatch() -> None:
    summary = _with_v27_restart_and_floquet_identity(_summary())
    summary["nonlinear_state_time_integrator_tangent_load_step_restart_generation_identity"].update(
        {
            "state_nonlinear_generation": "nonlinear-restart-140",
            "result_state_variable_names": ["hardening_variable"],
            "result_time_integrator": "bdf",
            "result_load_step_ids": [1, 2, 4, 5],
            "result_checkpoint_sha256": "d" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_restart_uses_current_state_integrator_tangent_load_step_and_checkpoint"
    ]


def test_v27_public_floquet_periodic_pair_orientation_phase_wavenumber_mode_normalization_dataset_mismatch() -> None:
    summary = _with_v27_restart_and_floquet_identity(_summary())
    summary["floquet_pair_orientation_phase_wavevector_normalization_dataset_mesh_generation_identity"].update(
        {
            "orientation_floquet_generation": "floquet-139",
            "result_pair_orientation_signs": [-1, 1],
            "result_phase_shift_rad": [-0.1, 0.2],
            "result_wavevector_rad_m": [-10.0, 20.0, 0.0],
            "result_mode_normalization": "max_field_1",
            "result_dataset_tag": "dset-floquet-1",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "floquet_modes_use_current_pair_orientation_phase_wavevector_normalization_and_dataset"
    ]


def _with_v28_thermoelastic_and_field_circuit_identity(summary: dict) -> dict:
    summary = _with_v27_restart_and_floquet_identity(summary)
    generation = "thermoelastic-frequency-151"
    summary[
        "thermoelastic_frequency_reference_temperature_prestress_linearization_mesh_dataset_generation_identity"
    ] = {
        "thermoelastic_generation": generation,
        "temperature_thermoelastic_generation": generation,
        "prestress_thermoelastic_generation": generation,
        "linearization_thermoelastic_generation": generation,
        "mesh_thermoelastic_generation": generation,
        "dataset_thermoelastic_generation": generation,
        "result_thermoelastic_generation": generation,
        "reference_temperature_k": 293.15,
        "result_reference_temperature_k": 293.15,
        "frequency_grid_hz": [100.0, 200.0, 500.0, 1000.0],
        "result_frequency_grid_hz": [100.0, 200.0, 500.0, 1000.0],
        "prestress_state_sha256": "1" * 64,
        "result_prestress_state_sha256": "1" * 64,
        "linearization_state_sha256": "2" * 64,
        "result_linearization_state_sha256": "2" * 64,
        "thermal_mesh_sha256": "3" * 64,
        "result_thermal_mesh_sha256": "3" * 64,
        "structural_mesh_sha256": "4" * 64,
        "result_structural_mesh_sha256": "4" * 64,
        "dataset_tag": "dset-thermoelastic-151",
        "result_dataset_tag": "dset-thermoelastic-151",
        "frequency_response_sha256": "5" * 64,
        "accepted_frequency_response_sha256": "5" * 64,
    }
    generation = "field-circuit-151"
    summary[
        "field_circuit_coil_terminal_orientation_current_sign_gauge_power_balance_mesh_solution_generation_identity"
    ] = {
        "coupling_generation": generation,
        "terminal_coupling_generation": generation,
        "orientation_coupling_generation": generation,
        "current_sign_coupling_generation": generation,
        "gauge_coupling_generation": generation,
        "power_coupling_generation": generation,
        "mesh_coupling_generation": generation,
        "solution_coupling_generation": generation,
        "result_coupling_generation": generation,
        "coil_terminal_ids": ["coil1:p", "coil1:n"],
        "result_coil_terminal_ids": ["coil1:p", "coil1:n"],
        "terminal_orientation_signs": [1, -1],
        "result_terminal_orientation_signs": [1, -1],
        "branch_current_signs": [1, -1],
        "result_branch_current_signs": [1, -1],
        "gauge_id": "magnetic_vector_potential_coulomb",
        "result_gauge_id": "magnetic_vector_potential_coulomb",
        "circuit_branch_id": "branch-coil1",
        "result_circuit_branch_id": "branch-coil1",
        "field_complex_power_va_ri": [12.0, 3.0],
        "result_field_complex_power_va_ri": [12.0, 3.0],
        "circuit_complex_power_va_ri": [-12.0, -3.0],
        "result_circuit_complex_power_va_ri": [-12.0, -3.0],
        "coupled_mesh_sha256": "6" * 64,
        "result_coupled_mesh_sha256": "6" * 64,
        "coupled_solution_sha256": "7" * 64,
        "accepted_coupled_solution_sha256": "7" * 64,
    }
    return summary


def test_v28_public_positive_thermoelastic_and_field_circuit_identity() -> None:
    result = gate(_with_v28_thermoelastic_and_field_circuit_identity(_summary()))
    assert result["status"] == "ok"


def test_v28_public_thermoelastic_frequency_reference_temperature_prestress_linearization_mesh_dataset_mismatch() -> None:
    summary = _with_v28_thermoelastic_and_field_circuit_identity(_summary())
    summary[
        "thermoelastic_frequency_reference_temperature_prestress_linearization_mesh_dataset_generation_identity"
    ].update(
        {
            "temperature_thermoelastic_generation": "thermoelastic-frequency-150",
            "result_reference_temperature_k": 323.15,
            "result_frequency_grid_hz": [100.0, 250.0, 1000.0],
            "result_prestress_state_sha256": "c" * 64,
            "result_dataset_tag": "dset-old",
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thermoelastic_frequency_uses_current_temperature_prestress_linearization_mesh_dataset_and_result"
    ]


def test_v28_public_field_circuit_coil_terminal_orientation_current_sign_gauge_power_balance_mismatch() -> None:
    summary = _with_v28_thermoelastic_and_field_circuit_identity(_summary())
    summary[
        "field_circuit_coil_terminal_orientation_current_sign_gauge_power_balance_mesh_solution_generation_identity"
    ].update(
        {
            "terminal_coupling_generation": "field-circuit-150",
            "result_coil_terminal_ids": ["coil1:n", "coil1:p"],
            "result_terminal_orientation_signs": [-1, 1],
            "result_gauge_id": "ungauged",
            "result_circuit_complex_power_va_ri": [8.0, 1.0],
            "accepted_coupled_solution_sha256": "2" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "field_circuit_uses_current_terminals_orientation_sign_gauge_power_mesh_and_solution"
    ]
