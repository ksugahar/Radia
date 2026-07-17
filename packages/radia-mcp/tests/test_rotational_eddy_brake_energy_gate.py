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
