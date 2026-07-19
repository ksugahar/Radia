"""Solver-neutral energy gate for a freely decelerating eddy-current brake."""

from __future__ import annotations

import math
import statistics
from typing import Any


_UNITS = {
    "time": "s",
    "angular_velocity": "rad/s",
    "torque": "N*m",
    "power": "W",
    "inertia": "kg*m^2",
    "energy": "J",
    "density": "kg/m^3",
    "length": "m",
}


def _series(row: dict[str, Any], key: str) -> list[float]:
    value = row.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a list")
    result = [float(item) for item in value]
    if not result or not all(math.isfinite(item) for item in result):
        raise ValueError(f"{key} must contain finite values")
    return result


def _increasing(values: list[float]) -> bool:
    return len(values) >= 2 and all(
        right > left for left, right in zip(values, values[1:])
    )


def _integral(values: list[float], times: list[float]) -> float:
    return sum(
        0.5 * (left + right) * (t_right - t_left)
        for left, right, t_left, t_right in zip(
            values[:-1], values[1:], times[:-1], times[1:], strict=True
        )
    )


def _cumulative_integral(values: list[float], times: list[float]) -> list[float]:
    result = [0.0]
    for left, right, t_left, t_right in zip(
        values[:-1], values[1:], times[:-1], times[1:], strict=True
    ):
        result.append(result[-1] + 0.5 * (left + right) * (t_right - t_left))
    return result


def _span_error(left: list[float], right: list[float], reference: list[float]) -> float:
    span = max(reference) - min(reference)
    if len(left) != len(right) or span <= 0.0:
        return math.inf
    return max(abs(a - b) for a, b in zip(left, right, strict=True)) / span


def _parse_replay(row: Any) -> dict[str, Any]:
    if not isinstance(row, dict):
        raise ValueError("replay entries must be mappings")
    parsed = {
        key: _series(row, key)
        for key in (
            "time_s",
            "angular_velocity_rad_s",
            "braking_torque_nm",
            "joule_loss_w",
        )
    }
    parsed["solve_seconds"] = float(row.get("solve_seconds", math.nan))
    return parsed


def _artifact_generations_ok(rows: list[dict[str, Any]]) -> bool:
    evidence_present = ["artifact_generations" in row for row in rows]
    if not any(evidence_present):
        return True
    if not all(evidence_present):
        return False
    base_keys = {
        "time_s",
        "angular_velocity_rad_s",
        "braking_torque_nm",
        "joule_loss_w",
    }
    for row in rows:
        generations = row.get("artifact_generations")
        if not isinstance(generations, dict):
            return False
        required = set(base_keys)
        if "magnetic_energy_j" in row or "field_energy_time_s" in row:
            required.update({"magnetic_energy_j", "field_energy_time_s"})
        values = [generations.get(name) for name in required]
        if not all(isinstance(value, str) and value for value in values):
            return False
        if len(set(values)) != 1:
            return False
    return True


def _artifact_coordinate_frames_ok(rows: list[dict[str, Any]]) -> bool:
    evidence_present = ["artifact_coordinate_frames" in row for row in rows]
    if not any(evidence_present):
        return True
    if not all(evidence_present):
        return False
    observed: set[str] = set()
    base_keys = {
        "time_s",
        "angular_velocity_rad_s",
        "braking_torque_nm",
        "joule_loss_w",
    }
    for row in rows:
        frames = row.get("artifact_coordinate_frames")
        if not isinstance(frames, dict):
            return False
        required = set(base_keys)
        if "magnetic_energy_j" in row or "field_energy_time_s" in row:
            required.update({"magnetic_energy_j", "field_energy_time_s"})
        values = [frames.get(name) for name in required]
        if not all(isinstance(value, str) and value for value in values):
            return False
        observed.update(values)
    return len(observed) == 1


def _convergence_provenance_ok(summary: dict[str, Any]) -> bool:
    if "convergence_provenance" not in summary:
        return True
    provenance = summary.get("convergence_provenance")
    if not isinstance(provenance, dict):
        return False
    try:
        residual = float(provenance.get("terminal_relative_residual"))
    except (TypeError, ValueError):
        return False
    return (
        bool(provenance.get("solution_generation"))
        and bool(provenance.get("result_iteration_generation"))
        and provenance.get("convergence_table_iteration_generation")
        == provenance.get("result_iteration_generation")
        and provenance.get("terminal_state") == "converged"
        and math.isfinite(residual)
        and residual >= 0.0
        and residual <= 1.0e-6
    )


def _force_selection_identity_ok(summary: dict[str, Any]) -> bool:
    if "force_selection_identity" not in summary:
        return True
    identity = summary.get("force_selection_identity")
    if not isinstance(identity, dict):
        return False
    geometry = identity.get("geometry_generation")
    selection_digest = str(identity.get("selection_entity_digest") or "")
    return (
        isinstance(geometry, str)
        and bool(geometry)
        and identity.get("solution_geometry_generation") == geometry
        and identity.get("integration_selection_generation") == geometry
        and len(selection_digest) == 64
        and identity.get("force_result_selection_digest") == selection_digest
    )


def _excitation_basis_identity_ok(summary: dict[str, Any]) -> bool:
    if "excitation_basis_identity" not in summary:
        return True
    identity = summary.get("excitation_basis_identity")
    if not isinstance(identity, dict):
        return False
    try:
        solve_scale = float(identity.get("solve_scale_to_rms"))
        extract_scale = float(identity.get("extract_scale_to_rms"))
    except (TypeError, ValueError):
        return False
    return (
        bool(identity.get("sweep_generation"))
        and identity.get("solve_amplitude_basis") == "rms"
        and identity.get("extract_amplitude_basis")
        == identity.get("solve_amplitude_basis")
        and math.isclose(solve_scale, 1.0, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(extract_scale, solve_scale, rel_tol=0.0, abs_tol=1.0e-15)
        and identity.get("torque_loss_normalization_basis") == "rms_excitation"
    )


def _live_stored_force_identity_ok(summary: dict[str, Any]) -> bool:
    if "live_stored_force_identity" not in summary:
        return True
    identity = summary.get("live_stored_force_identity")
    if not isinstance(identity, dict):
        return False
    geometry_digest = str(identity.get("live_geometry_sha256") or "")
    selection_digest = str(identity.get("live_selection_digest") or "")
    solution_generation = str(identity.get("live_solution_generation") or "")
    return (
        len(geometry_digest) == 64
        and identity.get("stored_force_geometry_sha256") == geometry_digest
        and bool(solution_generation)
        and identity.get("stored_force_solution_generation") == solution_generation
        and len(selection_digest) == 64
        and identity.get("stored_force_selection_digest") == selection_digest
    )


def _loss_partition_identity_ok(summary: dict[str, Any]) -> bool:
    if "loss_partition_identity" not in summary:
        return True
    identity = summary.get("loss_partition_identity")
    if not isinstance(identity, dict):
        return False
    ownership = identity.get("partition_ownership_ids")
    if not isinstance(ownership, list):
        return False
    normalized = [str(item).strip() for item in ownership]
    try:
        overlap_count = int(identity.get("ownership_overlap_count"))
        compensation = float(identity.get("signed_compensation_term_w"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("aggregation_generation") or "")
    return (
        bool(generation)
        and identity.get("reported_total_generation") == generation
        and bool(normalized)
        and all(normalized)
        and len(set(normalized)) == len(normalized)
        and overlap_count == 0
        and math.isclose(compensation, 0.0, rel_tol=0.0, abs_tol=1.0e-15)
    )


def _material_property_parameter_identity_ok(summary: dict[str, Any]) -> bool:
    if "material_property_parameter_identity" not in summary:
        return True
    identity = summary.get("material_property_parameter_identity")
    if not isinstance(identity, dict):
        return False
    try:
        parameter_value = float(identity.get("parameter_value"))
        evaluation_value = float(identity.get("material_evaluation_parameter_value"))
    except (TypeError, ValueError):
        return False
    parameter_generation = str(identity.get("parameter_generation") or "")
    return (
        bool(identity.get("parameter_name"))
        and math.isfinite(parameter_value)
        and math.isclose(
            evaluation_value, parameter_value, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and bool(identity.get("parameter_unit"))
        and identity.get("material_evaluation_parameter_unit")
        == identity.get("parameter_unit")
        and bool(parameter_generation)
        and identity.get("material_property_parameter_generation")
        == parameter_generation
    )


def _force_selection_topology_identity_ok(summary: dict[str, Any]) -> bool:
    if "force_selection_topology_identity" not in summary:
        return True
    identity = summary.get("force_selection_topology_identity")
    if not isinstance(identity, dict):
        return False
    selection_ids = identity.get("selection_entity_ids")
    integration_ids = identity.get("force_integration_entity_ids")
    if not isinstance(selection_ids, list) or not isinstance(integration_ids, list):
        return False
    if not selection_ids or any(not isinstance(item, int) for item in selection_ids):
        return False
    selection_digest = str(identity.get("selection_digest") or "")
    topology_generation = str(identity.get("topology_generation") or "")
    return (
        bool(identity.get("geometry_rebuild_generation"))
        and bool(topology_generation)
        and identity.get("selection_topology_generation") == topology_generation
        and identity.get("force_integration_topology_generation")
        == topology_generation
        and len(set(selection_ids)) == len(selection_ids)
        and integration_ids == selection_ids
        and len(selection_digest) == 64
        and all(character in "0123456789abcdef" for character in selection_digest)
        and identity.get("force_selection_digest") == selection_digest
    )


def _weak_form_coordinate_transform_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get("weak_form_coordinate_transform_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("mapped_coordinate_generation") or "")
    digest = str(identity.get("jacobian_orientation_sha256") or "")
    return (
        bool(generation)
        and identity.get("field_coordinate_generation") == generation
        and identity.get("jacobian_coordinate_generation") == generation
        and identity.get("jacobian_orientation") == "right_handed_positive"
        and identity.get("integration_orientation") == "right_handed_positive"
        and len(digest) == 64
        and identity.get("integration_orientation_sha256") == digest
    )


def _time_harmonic_phasor_convention_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get("time_harmonic_phasor_convention_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    convention = str(identity.get("source_time_convention") or "")
    generation = str(identity.get("phasor_generation") or "")
    return (
        convention in {"exp(+jwt)", "exp(-jwt)"}
        and identity.get("field_time_convention") == convention
        and identity.get("phase_sensitive_result_time_convention") == convention
        and identity.get("complex_power_formula") == "0.5*V*conj(I)"
        and identity.get("phase_sensitive_result_formula") == "0.5*V*conj(I)"
        and bool(generation)
        and identity.get("result_phasor_generation") == generation
    )


def _eigenmode_mass_inner_product_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get("eigenmode_mass_inner_product_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    eigensolve_generation = str(identity.get("eigensolve_generation") or "")
    mesh_generation = str(identity.get("mode_mesh_generation") or "")
    mass_generation = str(
        identity.get("eigensolve_mass_matrix_generation") or ""
    )
    mode_digest = str(identity.get("mode_vector_sha256") or "")
    return (
        bool(eigensolve_generation)
        and identity.get("mode_vector_generation") == eigensolve_generation
        and bool(mesh_generation)
        and identity.get("mass_matrix_mesh_generation") == mesh_generation
        and bool(mass_generation)
        and identity.get("normalization_mass_matrix_generation")
        == mass_generation
        and identity.get("normalization_kind") == "mass_inner_product"
        and identity.get("reference_normalization_kind")
        == "mass_inner_product"
        and len(mode_digest) == 64
        and all(character in "0123456789abcdef" for character in mode_digest)
        and identity.get("normalized_mode_vector_sha256") == mode_digest
    )


def _ale_material_derivative_time_level_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("ale_material_derivative_time_level_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    solve_generation = str(identity.get("ale_solve_generation") or "")
    time_generation = str(
        identity.get("accepted_time_level_generation") or ""
    )
    time_digest = str(identity.get("accepted_time_grid_sha256") or "")
    try:
        accepted_index = int(identity.get("accepted_time_index"))
        field_index = int(identity.get("field_time_index"))
        velocity_index = int(identity.get("mesh_velocity_time_index"))
        derivative_index = int(identity.get("material_derivative_time_index"))
    except (TypeError, ValueError):
        return False
    return (
        bool(solve_generation)
        and identity.get("field_solve_generation") == solve_generation
        and identity.get("mesh_velocity_solve_generation") == solve_generation
        and bool(time_generation)
        and identity.get("field_time_level_generation") == time_generation
        and identity.get("mesh_velocity_time_level_generation") == time_generation
        and identity.get("material_derivative_time_level_generation")
        == time_generation
        and accepted_index >= 0
        and field_index == accepted_index
        and velocity_index == accepted_index
        and derivative_index == accepted_index
        and len(time_digest) == 64
        and all(character in "0123456789abcdef" for character in time_digest)
        and identity.get("mesh_velocity_time_grid_sha256") == time_digest
    )


def _harmonic_reference_time_origin_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get("harmonic_reference_time_origin_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    solve_generation = str(identity.get("harmonic_solve_generation") or "")
    phase_digest = str(identity.get("field_phase_origin_sha256") or "")
    try:
        angular_frequency = float(identity.get("angular_frequency_rad_s"))
        field_time = float(identity.get("field_reference_time_s"))
        power_time = float(identity.get("complex_power_reference_time_s"))
        phase_offset = math.remainder(
            angular_frequency * (power_time - field_time), 2.0 * math.pi
        )
    except (TypeError, ValueError):
        return False
    return (
        bool(solve_generation)
        and identity.get("field_phasor_generation") == solve_generation
        and identity.get("complex_power_generation") == solve_generation
        and math.isfinite(angular_frequency)
        and angular_frequency > 0.0
        and math.isfinite(field_time)
        and math.isfinite(power_time)
        and abs(phase_offset) <= 1.0e-12
        and len(phase_digest) == 64
        and all(character in "0123456789abcdef" for character in phase_digest)
        and identity.get("power_phase_origin_sha256") == phase_digest
    )


def _deformed_domain_integral_jacobian_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("deformed_domain_integral_jacobian_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    solve_generation = str(identity.get("field_solve_generation") or "")
    geometry_generation = str(identity.get("geometry_generation") or "")
    selection_digest = str(identity.get("domain_selection_sha256") or "")
    jacobian_digest = str(identity.get("volume_jacobian_sha256") or "")
    digests = (selection_digest, jacobian_digest)
    return (
        bool(solve_generation)
        and identity.get("integral_field_generation") == solve_generation
        and bool(geometry_generation)
        and identity.get("integral_geometry_generation") == geometry_generation
        and identity.get("volume_jacobian_geometry_generation")
        == geometry_generation
        and all(
            len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            for digest in digests
        )
        and identity.get("integrated_domain_selection_sha256")
        == selection_digest
        and identity.get("integral_volume_jacobian_sha256") == jacobian_digest
    )


def _nonlinear_residual_tangent_iteration_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("nonlinear_residual_tangent_iteration_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("nonlinear_solve_generation") or "")
    material_digest = str(identity.get("material_state_sha256") or "")
    try:
        residual_iteration = int(identity.get("residual_iteration"))
        tangent_iteration = int(identity.get("tangent_iteration"))
        material_iteration = int(identity.get("material_state_iteration"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and identity.get("residual_solve_generation") == generation
        and identity.get("tangent_solve_generation") == generation
        and identity.get("material_state_solve_generation") == generation
        and residual_iteration >= 0
        and tangent_iteration == residual_iteration
        and material_iteration == residual_iteration
        and len(material_digest) == 64
        and all(character in "0123456789abcdef" for character in material_digest)
        and identity.get("tangent_material_state_sha256") == material_digest
    )


def _moving_mesh_field_transfer_frame_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("moving_mesh_field_transfer_frame_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("mesh_motion_generation") or "")
    frame = str(identity.get("source_coordinate_frame") or "")
    map_digest = str(identity.get("coordinate_map_sha256") or "")
    return (
        bool(generation)
        and identity.get("source_mesh_motion_generation") == generation
        and identity.get("target_mesh_motion_generation") == generation
        and identity.get("field_transfer_mesh_motion_generation") == generation
        and frame in {"material", "spatial"}
        and identity.get("target_coordinate_frame") == frame
        and identity.get("field_transfer_coordinate_frame") == frame
        and len(map_digest) == 64
        and all(character in "0123456789abcdef" for character in map_digest)
        and identity.get("field_transfer_coordinate_map_sha256") == map_digest
    )


def _segregated_block_variable_scaling_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("segregated_block_residual_variable_scaling_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    block_names = identity.get("block_names")
    residual_names = identity.get("residual_block_names")
    scaling_names = identity.get("variable_scaling_block_names")
    scaling_values = identity.get("variable_scaling_values")
    residual_scaling_values = identity.get("residual_variable_scaling_values")
    if not all(
        isinstance(values, list)
        for values in (
            block_names,
            residual_names,
            scaling_names,
            scaling_values,
            residual_scaling_values,
        )
    ):
        return False
    try:
        scales = [float(value) for value in scaling_values]
        residual_scales = [float(value) for value in residual_scaling_values]
    except (TypeError, ValueError):
        return False
    solve_generation = str(identity.get("solver_generation") or "")
    sequence_generation = str(identity.get("block_sequence_generation") or "")
    digest = str(identity.get("variable_scaling_sha256") or "")
    return (
        bool(solve_generation)
        and identity.get("block_residual_solver_generation") == solve_generation
        and identity.get("variable_scaling_solver_generation") == solve_generation
        and bool(sequence_generation)
        and identity.get("residual_block_sequence_generation")
        == sequence_generation
        and identity.get("variable_scaling_block_sequence_generation")
        == sequence_generation
        and bool(block_names)
        and len(set(block_names)) == len(block_names)
        and residual_names == block_names
        and scaling_names == block_names
        and len(scales) == len(block_names)
        and all(math.isfinite(value) and value > 0.0 for value in scales)
        and residual_scales == scales
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("residual_variable_scaling_sha256") == digest
    )


def _modal_port_power_surface_orientation_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "modal_port_power_normalization_surface_orientation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    mode_generation = str(identity.get("port_mode_generation") or "")
    surface_generation = str(identity.get("integration_surface_mesh_generation") or "")
    surface_digest = str(identity.get("surface_triangle_sha256") or "")
    return (
        bool(mode_generation)
        and identity.get("modal_amplitude_port_mode_generation") == mode_generation
        and identity.get("power_normalization_port_mode_generation")
        == mode_generation
        and bool(surface_generation)
        and identity.get("power_normalization_surface_mesh_generation")
        == surface_generation
        and identity.get("surface_orientation_mesh_generation")
        == surface_generation
        and identity.get("power_normalization") == "unit_forward_power"
        and identity.get("modal_amplitude_normalization")
        == "unit_forward_power"
        and identity.get("surface_orientation") == "outward_from_domain"
        and identity.get("power_flux_normal_sign") == 1
        and len(surface_digest) == 64
        and all(character in "0123456789abcdef" for character in surface_digest)
        and identity.get("power_normalization_surface_triangle_sha256")
        == surface_digest
    )


def _degenerate_eigenmode_subspace_tracking_basis_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("degenerate_eigenmode_subspace_tracking_basis_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    mode_ids = identity.get("cluster_mode_ids")
    tracking_ids = identity.get("tracking_basis_mode_ids")
    assurance = identity.get("modal_assurance_matrix")
    tracking_assurance = identity.get("tracking_modal_assurance_matrix")
    if not all(
        isinstance(values, list)
        for values in (mode_ids, tracking_ids, assurance, tracking_assurance)
    ):
        return False
    try:
        dimension = int(identity.get("subspace_dimension"))
        tracking_dimension = int(identity.get("tracking_basis_dimension"))
        mode_values = [int(value) for value in mode_ids]
        tracking_values = [int(value) for value in tracking_ids]
        assurance_values = [[float(value) for value in row] for row in assurance]
        tracking_values_matrix = [
            [float(value) for value in row] for row in tracking_assurance
        ]
    except (TypeError, ValueError):
        return False
    solve_generation = str(identity.get("eigensolve_generation") or "")
    mass_generation = str(identity.get("mass_inner_product_generation") or "")
    digest = str(identity.get("eigenspace_basis_sha256") or "")
    matrix_shape_ok = (
        len(assurance_values) == dimension
        and len(tracking_values_matrix) == dimension
        and all(len(row) == dimension for row in assurance_values)
        and all(len(row) == dimension for row in tracking_values_matrix)
    )
    return (
        bool(solve_generation)
        and identity.get("eigenvalue_cluster_generation") == solve_generation
        and identity.get("modal_vector_generation") == solve_generation
        and identity.get("tracking_basis_generation") == solve_generation
        and bool(mass_generation)
        and identity.get("tracking_mass_inner_product_generation") == mass_generation
        and dimension > 1
        and tracking_dimension == dimension
        and len(mode_values) == dimension
        and len(set(mode_values)) == dimension
        and tracking_values == mode_values
        and matrix_shape_ok
        and all(
            math.isfinite(value) and 0.0 <= value <= 1.0
            for row in assurance_values
            for value in row
        )
        and tracking_values_matrix == assurance_values
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("tracking_basis_sha256") == digest
    )


def _adaptive_bdf_restart_history_event_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get("adaptive_bdf_restart_history_event_generation_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    history = identity.get("history_time_s")
    replay_history = identity.get("solution_history_time_s")
    if not isinstance(history, list) or not isinstance(replay_history, list):
        return False
    try:
        order = int(identity.get("history_order"))
        replay_order = int(identity.get("solution_history_order"))
        history_values = [float(value) for value in history]
        replay_values = [float(value) for value in replay_history]
        event_time = float(identity.get("event_time_s"))
        restart_time = float(identity.get("restart_event_time_s"))
    except (TypeError, ValueError):
        return False
    step_generation = str(identity.get("accepted_step_generation") or "")
    method = str(identity.get("bdf_method") or "")
    event_id = str(identity.get("event_id") or "")
    digest = str(identity.get("history_state_sha256") or "")
    return (
        bool(str(identity.get("transient_generation") or ""))
        and bool(step_generation)
        and identity.get("solution_history_step_generation") == step_generation
        and identity.get("event_restart_step_generation") == step_generation
        and method in {"bdf1", "bdf2"}
        and identity.get("solution_history_method") == method
        and order in {1, 2}
        and replay_order == order
        and ((method == "bdf1" and order == 1) or (method == "bdf2" and order == 2))
        and len(history_values) == order + 1
        and replay_values == history_values
        and all(math.isfinite(value) for value in history_values)
        and all(right > left for left, right in zip(history_values, history_values[1:]))
        and bool(event_id)
        and identity.get("restart_event_id") == event_id
        and math.isfinite(event_time)
        and math.isclose(restart_time, event_time, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(history_values[-1], event_time, rel_tol=0.0, abs_tol=1.0e-15)
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("restart_history_state_sha256") == digest
    )


def _nonlinear_continuation_branch_tangent_checkpoint_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("nonlinear_continuation_branch_tangent_checkpoint_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    tangent = identity.get("tangent_vector")
    checkpoint_tangent = identity.get("checkpoint_tangent_vector")
    if not isinstance(tangent, list) or not isinstance(checkpoint_tangent, list):
        return False
    try:
        load_value = float(identity.get("load_parameter_value"))
        tangent_load = float(identity.get("tangent_load_parameter_value"))
        checkpoint_load = float(identity.get("checkpoint_load_parameter_value"))
        tangent_values = [float(value) for value in tangent]
        checkpoint_values = [float(value) for value in checkpoint_tangent]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("continuation_state_generation") or "")
    branch = str(identity.get("branch_id") or "")
    digest = str(identity.get("continuation_tangent_sha256") or "")
    return (
        bool(str(identity.get("nonlinear_solve_generation") or ""))
        and bool(generation)
        and identity.get("tangent_continuation_state_generation") == generation
        and identity.get("checkpoint_continuation_state_generation") == generation
        and bool(branch)
        and identity.get("tangent_branch_id") == branch
        and identity.get("checkpoint_branch_id") == branch
        and bool(str(identity.get("load_parameter_name") or ""))
        and math.isfinite(load_value)
        and math.isclose(tangent_load, load_value, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(checkpoint_load, load_value, rel_tol=0.0, abs_tol=1.0e-15)
        and bool(tangent_values)
        and all(math.isfinite(value) for value in tangent_values)
        and checkpoint_values == tangent_values
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("checkpoint_tangent_sha256") == digest
    )


def _nonconforming_mortar_projection_quadrature_mesh_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get("nonconforming_mortar_projection_quadrature_mesh_identity")
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    fields = (
        identity.get("source_trace_dof_ids"),
        identity.get("projection_source_trace_dof_ids"),
        identity.get("target_trace_dof_ids"),
        identity.get("projection_target_trace_dof_ids"),
        identity.get("projection_shape"),
        identity.get("quadrature_projection_shape"),
    )
    if not all(isinstance(values, list) for values in fields):
        return False
    try:
        source_ids = [int(value) for value in fields[0]]
        projected_source_ids = [int(value) for value in fields[1]]
        target_ids = [int(value) for value in fields[2]]
        projected_target_ids = [int(value) for value in fields[3]]
        shape = [int(value) for value in fields[4]]
        quadrature_shape = [int(value) for value in fields[5]]
    except (TypeError, ValueError):
        return False
    source_generation = str(identity.get("source_mesh_generation") or "")
    target_generation = str(identity.get("target_mesh_generation") or "")
    digest = str(identity.get("projection_operator_sha256") or "")
    return (
        bool(str(identity.get("interface_generation") or ""))
        and bool(source_generation)
        and identity.get("projection_source_mesh_generation") == source_generation
        and identity.get("quadrature_source_mesh_generation") == source_generation
        and bool(target_generation)
        and identity.get("projection_target_mesh_generation") == target_generation
        and identity.get("quadrature_target_mesh_generation") == target_generation
        and bool(source_ids)
        and len(set(source_ids)) == len(source_ids)
        and projected_source_ids == source_ids
        and bool(target_ids)
        and len(set(target_ids)) == len(target_ids)
        and projected_target_ids == target_ids
        and shape == [len(target_ids), len(source_ids)]
        and quadrature_shape == shape
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("quadrature_projection_operator_sha256") == digest
    )


def _adaptive_mesh_field_transfer_conservation_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "adaptive_mesh_field_transfer_projection_conservation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    rows = (
        identity.get("source_field_values"),
        identity.get("source_integration_weights"),
        identity.get("projected_field_values"),
        identity.get("target_integration_weights"),
        identity.get("projection_shape"),
    )
    if not all(isinstance(values, list) for values in rows):
        return False
    try:
        source_values = [float(value) for value in rows[0]]
        source_weights = [float(value) for value in rows[1]]
        target_values = [float(value) for value in rows[2]]
        target_weights = [float(value) for value in rows[3]]
        shape = [int(value) for value in rows[4]]
        reported_source = float(identity.get("source_conserved_integral"))
        reported_target = float(identity.get("target_conserved_integral"))
    except (TypeError, ValueError):
        return False
    source_integral = sum(
        value * weight for value, weight in zip(source_values, source_weights)
    )
    target_integral = sum(
        value * weight for value, weight in zip(target_values, target_weights)
    )
    source_generation = str(identity.get("source_mesh_generation") or "")
    target_generation = str(identity.get("target_mesh_generation") or "")
    weight_digest = str(identity.get("conservation_weight_table_sha256") or "")
    return (
        bool(str(identity.get("solve_generation") or ""))
        and bool(source_generation)
        and identity.get("projection_source_mesh_generation") == source_generation
        and identity.get("conservation_source_mesh_generation") == source_generation
        and bool(target_generation)
        and target_generation != source_generation
        and identity.get("projection_target_mesh_generation") == target_generation
        and identity.get("conservation_target_mesh_generation") == target_generation
        and bool(source_values)
        and len(source_values) == len(source_weights)
        and bool(target_values)
        and len(target_values) == len(target_weights)
        and shape == [len(target_values), len(source_values)]
        and all(math.isfinite(value) for value in source_values + target_values)
        and all(math.isfinite(value) and value > 0.0 for value in source_weights + target_weights)
        and math.isclose(reported_source, source_integral, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(reported_target, target_integral, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(target_integral, source_integral, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and len(str(identity.get("projection_operator_sha256") or "")) == 64
        and len(weight_digest) == 64
        and all(character in "0123456789abcdef" for character in weight_digest)
        and identity.get("transfer_conservation_weight_table_sha256") == weight_digest
    )


def _eigenmode_phase_normalization_tracking_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "eigenmode_phase_normalization_tracking_parameter_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    rows = (
        identity.get("tracked_mode_ids"),
        identity.get("tracker_mode_ids"),
        identity.get("phase_anchor_dof_ids"),
        identity.get("tracker_phase_anchor_dof_ids"),
        identity.get("normalization_integrals"),
        identity.get("tracker_normalization_integrals"),
        identity.get("selected_correlation"),
        identity.get("tracker_selected_correlation"),
    )
    if not all(isinstance(values, list) for values in rows):
        return False
    try:
        modes = [int(value) for value in rows[0]]
        tracker_modes = [int(value) for value in rows[1]]
        anchors = [int(value) for value in rows[2]]
        tracker_anchors = [int(value) for value in rows[3]]
        norms = [float(value) for value in rows[4]]
        tracker_norms = [float(value) for value in rows[5]]
        correlations = [float(value) for value in rows[6]]
        tracker_correlations = [float(value) for value in rows[7]]
        previous_parameter = float(identity.get("previous_parameter_value"))
        current_parameter = float(identity.get("current_parameter_value"))
    except (TypeError, ValueError):
        return False
    current_generation = str(identity.get("current_eigensolve_generation") or "")
    digest = str(identity.get("mode_tracking_table_sha256") or "")
    return (
        bool(str(identity.get("parameter_table_generation") or ""))
        and bool(str(identity.get("previous_eigensolve_generation") or ""))
        and bool(current_generation)
        and identity.get("tracker_current_eigensolve_generation") == current_generation
        and identity.get("phase_anchor_current_eigensolve_generation") == current_generation
        and identity.get("normalization_current_eigensolve_generation") == current_generation
        and bool(str(identity.get("parameter_name") or ""))
        and math.isfinite(previous_parameter)
        and math.isfinite(current_parameter)
        and current_parameter != previous_parameter
        and bool(modes)
        and len(set(modes)) == len(modes)
        and tracker_modes == modes
        and len(anchors) == len(modes)
        and len(set(anchors)) == len(anchors)
        and tracker_anchors == anchors
        and len(norms) == len(modes)
        and all(math.isfinite(value) and value > 0.0 for value in norms)
        and tracker_norms == norms
        and len(correlations) == len(modes)
        and all(math.isfinite(value) and 0.0 <= value <= 1.0 for value in correlations)
        and tracker_correlations == correlations
        and len(digest) == 64
        and all(character in "0123456789abcdef" for character in digest)
        and identity.get("tracker_mode_tracking_table_sha256") == digest
    )


def _parameter_sweep_branch_restart_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "parameter_sweep_branch_restart_solution_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    rows = (
        identity.get("parameter_names"),
        identity.get("current_parameter_tuple"),
        identity.get("checkpoint_parameter_tuple"),
    )
    if not all(isinstance(values, list) for values in rows):
        return False
    try:
        names = [str(value) for value in rows[0]]
        current_tuple = [float(value) for value in rows[1]]
        checkpoint_tuple = [float(value) for value in rows[2]]
        current_branch_id = int(identity.get("current_branch_id"))
        checkpoint_branch_id = int(identity.get("checkpoint_branch_id"))
    except (TypeError, ValueError):
        return False
    sweep_generation = str(identity.get("parameter_sweep_generation") or "")
    branch_generation = str(identity.get("current_branch_generation") or "")
    digest_pairs = (
        ("solution_vector_sha256", "checkpoint_solution_vector_sha256"),
        ("continuation_state_sha256", "checkpoint_continuation_state_sha256"),
        ("mesh_coordinate_sha256", "checkpoint_mesh_coordinate_sha256"),
    )
    return (
        bool(sweep_generation)
        and identity.get("checkpoint_parameter_sweep_generation")
        == sweep_generation
        and bool(branch_generation)
        and identity.get("checkpoint_branch_generation") == branch_generation
        and identity.get("solution_vector_branch_generation") == branch_generation
        and identity.get("continuation_state_branch_generation")
        == branch_generation
        and identity.get("mesh_coordinate_branch_generation")
        == branch_generation
        and current_branch_id >= 0
        and checkpoint_branch_id == current_branch_id
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and len(names) == len(current_tuple) == len(checkpoint_tuple)
        and all(math.isfinite(value) for value in current_tuple + checkpoint_tuple)
        and checkpoint_tuple == current_tuple
        and all(
            len(str(identity.get(current_key) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(current_key) or "")
            )
            and identity.get(checkpoint_key) == identity.get(current_key)
            for current_key, checkpoint_key in digest_pairs
        )
    )


def _multiphysics_coupling_source_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "multiphysics_coupling_source_frame_unit_selection_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    rows = (
        identity.get("source_units"),
        identity.get("assembled_source_units"),
        identity.get("source_selection_ids"),
        identity.get("assembled_source_selection_ids"),
        identity.get("source_selection_dimensions"),
        identity.get("assembled_source_selection_dimensions"),
    )
    if not all(isinstance(values, list) for values in rows):
        return False
    try:
        source_units = [str(value) for value in rows[0]]
        assembled_units = [str(value) for value in rows[1]]
        source_ids = [int(value) for value in rows[2]]
        assembled_ids = [int(value) for value in rows[3]]
        source_dimensions = [int(value) for value in rows[4]]
        assembled_dimensions = [int(value) for value in rows[5]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("coupling_generation") or "")
    frame = str(identity.get("source_coordinate_frame") or "")
    value_digest = str(identity.get("source_values_sha256") or "")
    selection_digest = str(identity.get("source_selection_sha256") or "")
    digests = (value_digest, selection_digest)
    return (
        bool(generation)
        and identity.get("source_values_coupling_generation") == generation
        and identity.get("source_frame_coupling_generation") == generation
        and identity.get("source_unit_coupling_generation") == generation
        and identity.get("source_selection_coupling_generation") == generation
        and frame in {"material", "spatial"}
        and identity.get("assembled_coordinate_frame") == frame
        and bool(source_units)
        and all(source_units)
        and assembled_units == source_units
        and bool(source_ids)
        and all(value > 0 for value in source_ids)
        and len(set(source_ids)) == len(source_ids)
        and assembled_ids == source_ids
        and len(source_dimensions) == len(source_ids)
        and all(value in {0, 1, 2, 3} for value in source_dimensions)
        and assembled_dimensions == source_dimensions
        and all(
            len(digest) == 64
            and all(character in "0123456789abcdef" for character in digest)
            for digest in digests
        )
        and identity.get("assembled_source_values_sha256") == value_digest
        and identity.get("assembled_source_selection_sha256")
        == selection_digest
    )


def _contact_active_set_friction_state_mesh_generation_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "contact_active_set_friction_state_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    surfaces = identity.get("slave_surface_ids")
    active_surfaces = identity.get("active_slave_surface_ids")
    if not isinstance(surfaces, list) or not isinstance(active_surfaces, list):
        return False
    try:
        surface_ids = [int(value) for value in surfaces]
        active_surface_ids = [int(value) for value in active_surfaces]
        friction = float(identity.get("friction_coefficient"))
        active_friction = float(identity.get("active_set_friction_coefficient"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("contact_generation") or "")
    normal = str(identity.get("normal_orientation") or "")
    digest_pairs = (
        ("slave_mesh_sha256", "active_set_slave_mesh_sha256"),
        ("friction_state_sha256", "active_set_friction_state_sha256"),
        ("consistent_tangent_sha256", "active_set_consistent_tangent_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "active_set_contact_generation",
                "friction_state_contact_generation",
                "slave_mesh_contact_generation",
                "normal_orientation_contact_generation",
                "consistent_tangent_contact_generation",
            )
        )
        and bool(surface_ids)
        and all(value > 0 for value in surface_ids)
        and len(set(surface_ids)) == len(surface_ids)
        and active_surface_ids == surface_ids
        and normal in {"slave_to_master", "master_to_slave"}
        and identity.get("active_set_normal_orientation") == normal
        and math.isfinite(friction)
        and friction >= 0.0
        and math.isclose(active_friction, friction, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(
            len(str(identity.get(current) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(current) or "")
            )
            and identity.get(active) == identity.get(current)
            for current, active in digest_pairs
        )
    )


def _acoustic_structure_trace_impedance_order_frame_generation_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "acoustic_structure_trace_impedance_order_frame_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        sign = int(identity.get("pressure_to_traction_sign"))
        assembled_sign = int(identity.get("assembled_pressure_to_traction_sign"))
        order = int(identity.get("impedance_order"))
        assembled_order = int(identity.get("assembled_impedance_order"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("coupling_generation") or "")
    frame = str(identity.get("normal_frame") or "")
    digest_pairs = (
        ("pressure_trace_sha256", "assembled_pressure_trace_sha256"),
        ("traction_trace_sha256", "assembled_traction_trace_sha256"),
        ("interface_mesh_sha256", "assembled_interface_mesh_sha256"),
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "pressure_trace_coupling_generation",
                "traction_trace_coupling_generation",
                "normal_frame_coupling_generation",
                "impedance_coupling_generation",
                "interface_mesh_coupling_generation",
            )
        )
        and frame in {"acoustic_outward", "structure_outward"}
        and identity.get("traction_normal_frame") == frame
        and sign in {-1, 1}
        and assembled_sign == sign
        and order >= 1
        and assembled_order == order
        and all(
            len(str(identity.get(current) or "")) == 64
            and all(
                character in "0123456789abcdef"
                for character in str(identity.get(current) or "")
            )
            and identity.get(assembled) == identity.get(current)
            for current, assembled in digest_pairs
        )
    )


def _is_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        character in "0123456789abcdef" for character in value
    )


def _continuation_branch_load_mesh_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "nonlinear_continuation_branch_load_step_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    rows = (
        identity.get("load_parameters"),
        identity.get("result_load_parameters"),
        identity.get("tangent_state_sha256"),
        identity.get("result_tangent_state_sha256"),
    )
    if not all(isinstance(row, list) for row in rows):
        return False
    try:
        loads = [float(value) for value in rows[0]]
        result_loads = [float(value) for value in rows[1]]
        tangents = [str(value) for value in rows[2]]
        result_tangents = [str(value) for value in rows[3]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("solve_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "branch_solve_generation",
                "load_step_solve_generation",
                "tangent_state_solve_generation",
                "adapted_mesh_solve_generation",
                "result_solve_generation",
            )
        )
        and bool(str(identity.get("branch_id") or ""))
        and identity.get("result_branch_id") == identity.get("branch_id")
        and len(loads) >= 2
        and all(math.isfinite(value) for value in loads)
        and all(right > left for left, right in zip(loads, loads[1:]))
        and result_loads == loads
        and len(tangents) == len(loads)
        and all(_is_sha256(value) for value in tangents)
        and result_tangents == tangents
        and _is_sha256(str(identity.get("adapted_mesh_sha256") or ""))
        and identity.get("result_mesh_sha256") == identity.get("adapted_mesh_sha256")
        and _is_sha256(str(identity.get("continuation_table_sha256") or ""))
        and identity.get("result_continuation_table_sha256")
        == identity.get("continuation_table_sha256")
    )


def _parametric_sequence_initial_solution_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "parametric_sequence_initial_solution_dataset_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    rows = (
        identity.get("parameter_names"),
        identity.get("parameter_rows"),
        identity.get("result_parameter_rows"),
        identity.get("initial_solution_sha256"),
        identity.get("result_initial_solution_sha256"),
    )
    if not all(isinstance(row, list) for row in rows):
        return False
    try:
        names = [str(value) for value in rows[0]]
        parameters = [[float(value) for value in row] for row in rows[1]]
        result_parameters = [[float(value) for value in row] for row in rows[2]]
        solutions = [str(value) for value in rows[3]]
        result_solutions = [str(value) for value in rows[4]]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sweep_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "sequence_sweep_generation",
                "parameter_row_sweep_generation",
                "initial_solution_sweep_generation",
                "dataset_sweep_generation",
                "result_sweep_generation",
            )
        )
        and bool(str(identity.get("sequence_id") or ""))
        and identity.get("result_sequence_id") == identity.get("sequence_id")
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and bool(parameters)
        and all(len(row) == len(names) for row in parameters)
        and all(math.isfinite(value) for row in parameters for value in row)
        and result_parameters == parameters
        and len(solutions) == len(parameters)
        and all(_is_sha256(value) for value in solutions)
        and result_solutions == solutions
        and _is_sha256(str(identity.get("dataset_sha256") or ""))
        and identity.get("result_dataset_sha256") == identity.get("dataset_sha256")
    )


def _multiphysics_power_work_heat_balance_identity_ok(
    summary: dict[str, Any],
) -> bool:
    identity = summary.get(
        "multiphysics_power_work_heat_balance_frame_time_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        window = [float(value) for value in identity.get("time_window_s", [])]
        result_window = [
            float(value) for value in identity.get("result_time_window_s", [])
        ]
        electromagnetic = float(identity.get("electromagnetic_input_j"))
        mechanical = float(identity.get("mechanical_work_j"))
        heat = float(identity.get("heat_source_j"))
        stored = float(identity.get("stored_energy_change_j"))
        reported = float(identity.get("reported_balance_j"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("balance_generation") or "")
    frame = str(identity.get("coordinate_frame") or "")
    scale = max(abs(electromagnetic), 1.0)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "electromagnetic_balance_generation",
                "mechanical_balance_generation",
                "thermal_balance_generation",
                "stored_energy_balance_generation",
                "result_balance_generation",
            )
        )
        and bool(frame)
        and identity.get("result_coordinate_frame") == frame
        and len(window) == 2
        and all(math.isfinite(value) for value in window)
        and window[0] < window[1]
        and result_window == window
        and all(
            math.isfinite(value)
            for value in (electromagnetic, mechanical, heat, stored, reported)
        )
        and abs(electromagnetic - mechanical - heat - stored) <= 1.0e-12 * scale
        and abs(reported) <= 1.0e-12 * scale
        and _is_sha256(str(identity.get("balance_table_sha256") or ""))
        and identity.get("result_balance_table_sha256")
        == identity.get("balance_table_sha256")
    )


def _degenerate_eigenmode_subspace_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "degenerate_eigenmode_subspace_normalization_phase_projection_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    sequence_keys = (
        "eigenvalues_hz",
        "result_eigenvalues_hz",
        "complex_phase_rad",
        "result_complex_phase_rad",
        "subspace_projection",
        "result_subspace_projection",
    )
    if not all(isinstance(identity.get(key), list) for key in sequence_keys):
        return False
    try:
        eigenvalues = [float(value) for value in identity["eigenvalues_hz"]]
        result_eigenvalues = [
            float(value) for value in identity["result_eigenvalues_hz"]
        ]
        phases = [float(value) for value in identity["complex_phase_rad"]]
        result_phases = [
            float(value) for value in identity["result_complex_phase_rad"]
        ]
        projection = [
            [float(value) for value in row]
            for row in identity["subspace_projection"]
        ]
        result_projection = [
            [float(value) for value in row]
            for row in identity["result_subspace_projection"]
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("eigensolve_generation") or "")
    normalization = str(identity.get("normalization") or "")
    count = len(eigenvalues)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "subspace_eigensolve_generation",
                "normalization_eigensolve_generation",
                "phase_eigensolve_generation",
                "projection_eigensolve_generation",
                "mesh_eigensolve_generation",
                "result_eigensolve_generation",
            )
        )
        and count >= 2
        and all(math.isfinite(value) and value > 0.0 for value in eigenvalues)
        and result_eigenvalues == eigenvalues
        and bool(normalization)
        and identity.get("result_normalization") == normalization
        and len(phases) == count
        and all(math.isfinite(value) for value in phases)
        and result_phases == phases
        and len(projection) == count
        and all(len(row) == count for row in projection)
        and all(math.isfinite(value) for row in projection for value in row)
        and result_projection == projection
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("result_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("modal_table_sha256") or ""))
        and identity.get("result_modal_table_sha256")
        == identity.get("modal_table_sha256")
    )


def _remesh_field_projection_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "remesh_field_projection_conservation_geometry_dataset_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("projection_generation") or "")
    try:
        before = float(identity.get("conserved_integral_before"))
        after = float(identity.get("conserved_integral_after"))
    except (TypeError, ValueError):
        return False
    source_mesh = str(identity.get("source_mesh_sha256") or "")
    target_mesh = str(identity.get("target_mesh_sha256") or "")
    geometry = str(identity.get("geometry_revision_sha256") or "")
    projection_map = str(identity.get("projection_map_sha256") or "")
    field = str(identity.get("projected_field_sha256") or "")
    dataset = str(identity.get("dataset_tag") or "")
    conservation_scale = max(abs(before), 1.0)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "source_mesh_projection_generation",
                "target_mesh_projection_generation",
                "geometry_projection_generation",
                "dataset_projection_generation",
                "integral_projection_generation",
                "result_projection_generation",
            )
        )
        and _is_sha256(source_mesh)
        and identity.get("projected_source_mesh_sha256") == source_mesh
        and _is_sha256(target_mesh)
        and identity.get("result_target_mesh_sha256") == target_mesh
        and _is_sha256(geometry)
        and identity.get("result_geometry_revision_sha256") == geometry
        and bool(dataset)
        and identity.get("result_dataset_tag") == dataset
        and math.isfinite(before)
        and math.isfinite(after)
        and abs(after - before) <= 1.0e-12 * conservation_scale
        and _is_sha256(projection_map)
        and identity.get("result_projection_map_sha256") == projection_map
        and _is_sha256(field)
        and identity.get("result_projected_field_sha256") == field
    )


def _nonlinear_continuation_load_step_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "nonlinear_continuation_load_step_branch_state_solver_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        step_ids = [int(value) for value in identity.get("load_step_ids", [])]
        result_step_ids = [
            int(value) for value in identity.get("result_load_step_ids", [])
        ]
        parameters = [
            float(value) for value in identity.get("continuation_parameter", [])
        ]
        result_parameters = [
            float(value)
            for value in identity.get("result_continuation_parameter", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("continuation_generation") or "")
    branch = str(identity.get("branch_id") or "")
    initial_state = str(identity.get("initial_state_sha256") or "")
    solver_settings = str(identity.get("solver_settings_sha256") or "")
    solution = str(identity.get("solution_table_sha256") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "load_step_continuation_generation",
                "branch_continuation_generation",
                "initial_state_continuation_generation",
                "solver_continuation_generation",
                "result_continuation_generation",
            )
        )
        and len(step_ids) >= 2
        and all(value > 0 for value in step_ids)
        and len(set(step_ids)) == len(step_ids)
        and all(right > left for left, right in zip(step_ids, step_ids[1:]))
        and result_step_ids == step_ids
        and len(parameters) == len(step_ids)
        and all(math.isfinite(value) for value in parameters)
        and all(right > left for left, right in zip(parameters, parameters[1:]))
        and result_parameters == parameters
        and bool(branch)
        and identity.get("result_branch_id") == branch
        and _is_sha256(initial_state)
        and identity.get("result_initial_state_sha256") == initial_state
        and _is_sha256(solver_settings)
        and identity.get("result_solver_settings_sha256") == solver_settings
        and _is_sha256(solution)
        and identity.get("result_solution_table_sha256") == solution
    )


def _ale_force_work_balance_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "ale_moving_mesh_time_step_field_transfer_force_work_balance_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        time_s = [float(value) for value in identity.get("time_s", [])]
        result_time_s = [float(value) for value in identity.get("result_time_s", [])]
        displacement = [
            float(value) for value in identity.get("mesh_displacement_m", [])
        ]
        result_displacement = [
            float(value) for value in identity.get("result_mesh_displacement_m", [])
        ]
        force = [float(value) for value in identity.get("force_n", [])]
        result_force = [float(value) for value in identity.get("result_force_n", [])]
        reported_work = float(identity.get("reported_mechanical_work_j"))
        field_energy = float(identity.get("field_energy_change_j"))
        dissipation = float(identity.get("dissipated_energy_j"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("ale_generation") or "")
    mechanical_work = sum(
        force_value * displacement_value
        for force_value, displacement_value in zip(force, displacement, strict=True)
    ) if len(force) == len(displacement) else math.nan
    scale = max(abs(mechanical_work), 1.0)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "geometry_ale_generation",
                "time_step_ale_generation",
                "field_transfer_ale_generation",
                "force_ale_generation",
                "work_ale_generation",
                "result_ale_generation",
            )
        )
        and len(time_s) >= 2
        and all(math.isfinite(value) and value >= 0.0 for value in time_s)
        and all(right > left for left, right in zip(time_s, time_s[1:]))
        and result_time_s == time_s
        and len(displacement) == len(force) == len(time_s) - 1
        and all(math.isfinite(value) for value in displacement + force)
        and result_displacement == displacement
        and result_force == force
        and math.isfinite(reported_work)
        and abs(reported_work - mechanical_work) <= 1.0e-12 * scale
        and math.isfinite(field_energy)
        and math.isfinite(dissipation)
        and abs(field_energy + dissipation - mechanical_work) <= 1.0e-12 * scale
        and _is_sha256(str(identity.get("geometry_mesh_sha256") or ""))
        and identity.get("result_geometry_mesh_sha256")
        == identity.get("geometry_mesh_sha256")
        and _is_sha256(str(identity.get("field_transfer_sha256") or ""))
        and identity.get("result_field_transfer_sha256")
        == identity.get("field_transfer_sha256")
        and _is_sha256(str(identity.get("force_work_table_sha256") or ""))
        and identity.get("result_force_work_table_sha256")
        == identity.get("force_work_table_sha256")
    )


def _segregated_iteration_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "segregated_multiphysics_iteration_relaxation_residual_component_solution_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        iteration_ids = [int(value) for value in identity.get("iteration_ids", [])]
        result_iteration_ids = [
            int(value) for value in identity.get("result_iteration_ids", [])
        ]
        components = [str(value) for value in identity.get("component_order", [])]
        result_components = [
            str(value) for value in identity.get("result_component_order", [])
        ]
        relaxation = [
            float(value) for value in identity.get("relaxation_factors", [])
        ]
        result_relaxation = [
            float(value) for value in identity.get("result_relaxation_factors", [])
        ]
        residuals = [float(value) for value in identity.get("residual_history", [])]
        result_residuals = [
            float(value) for value in identity.get("result_residual_history", [])
        ]
        tolerance = float(identity.get("relative_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("solve_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "iteration_solve_generation",
                "relaxation_solve_generation",
                "residual_solve_generation",
                "component_solve_generation",
                "solution_solve_generation",
                "result_solve_generation",
            )
        )
        and len(iteration_ids) >= 2
        and iteration_ids == list(range(1, len(iteration_ids) + 1))
        and result_iteration_ids == iteration_ids
        and bool(components)
        and all(components)
        and len(set(components)) == len(components)
        and result_components == components
        and len(relaxation) == len(components)
        and all(math.isfinite(value) and 0.0 < value <= 1.0 for value in relaxation)
        and result_relaxation == relaxation
        and identity.get("residual_norm") == "l2"
        and identity.get("result_residual_norm") == identity.get("residual_norm")
        and len(residuals) == len(iteration_ids)
        and all(math.isfinite(value) and value >= 0.0 for value in residuals)
        and all(right < left for left, right in zip(residuals, residuals[1:]))
        and result_residuals == residuals
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and residuals[-1] <= tolerance
        and identity.get("converged") is True
        and _is_sha256(str(identity.get("solution_sha256") or ""))
        and identity.get("result_solution_sha256") == identity.get("solution_sha256")
    )


def _nonlinear_state_restart_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "nonlinear_state_time_integrator_tangent_load_step_restart_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        states = [str(value) for value in identity.get("state_variable_names", [])]
        result_states = [
            str(value) for value in identity.get("result_state_variable_names", [])
        ]
        load_steps = [int(value) for value in identity.get("load_step_ids", [])]
        result_load_steps = [
            int(value) for value in identity.get("result_load_step_ids", [])
        ]
        order = int(identity.get("integrator_order"))
        result_order = int(identity.get("result_integrator_order"))
        restart_time = float(identity.get("restart_time_s"))
        result_restart_time = float(identity.get("result_restart_time_s"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("nonlinear_generation") or "")
    integrator = str(identity.get("time_integrator") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "state_nonlinear_generation",
                "integrator_nonlinear_generation",
                "tangent_nonlinear_generation",
                "load_step_nonlinear_generation",
                "checkpoint_nonlinear_generation",
                "result_nonlinear_generation",
            )
        )
        and bool(states)
        and all(states)
        and len(set(states)) == len(states)
        and result_states == states
        and integrator in {"generalized_alpha", "bdf"}
        and identity.get("result_time_integrator") == integrator
        and order in {1, 2}
        and result_order == order
        and len(load_steps) >= 2
        and load_steps == list(range(load_steps[0], load_steps[0] + len(load_steps)))
        and result_load_steps == load_steps
        and math.isfinite(restart_time)
        and restart_time >= 0.0
        and result_restart_time == restart_time
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(result) == identity.get(source)
            for source, result in (
                ("state_vector_sha256", "result_state_vector_sha256"),
                ("consistent_tangent_sha256", "result_consistent_tangent_sha256"),
                ("checkpoint_sha256", "result_checkpoint_sha256"),
                ("solution_sha256", "result_solution_sha256"),
            )
        )
    )


def _floquet_pair_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "floquet_pair_orientation_phase_wavevector_normalization_dataset_mesh_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        pairs = [str(value) for value in identity.get("periodic_pair_tags", [])]
        result_pairs = [
            str(value) for value in identity.get("result_periodic_pair_tags", [])
        ]
        signs = [int(value) for value in identity.get("pair_orientation_signs", [])]
        result_signs = [
            int(value) for value in identity.get("result_pair_orientation_signs", [])
        ]
        phases = [float(value) for value in identity.get("phase_shift_rad", [])]
        result_phases = [
            float(value) for value in identity.get("result_phase_shift_rad", [])
        ]
        wavevector = [float(value) for value in identity.get("wavevector_rad_m", [])]
        result_wavevector = [
            float(value) for value in identity.get("result_wavevector_rad_m", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("floquet_generation") or "")
    normalization = str(identity.get("mode_normalization") or "")
    dataset = str(identity.get("dataset_tag") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "pair_floquet_generation",
                "orientation_floquet_generation",
                "phase_floquet_generation",
                "wavevector_floquet_generation",
                "normalization_floquet_generation",
                "dataset_floquet_generation",
                "mesh_floquet_generation",
                "result_floquet_generation",
            )
        )
        and bool(pairs)
        and all(pairs)
        and len(set(pairs)) == len(pairs)
        and result_pairs == pairs
        and len(signs) == len(pairs)
        and all(value in {-1, 1} for value in signs)
        and result_signs == signs
        and len(phases) == len(pairs)
        and all(math.isfinite(value) for value in phases)
        and result_phases == phases
        and len(wavevector) == 3
        and all(math.isfinite(value) for value in wavevector)
        and result_wavevector == wavevector
        and normalization in {"unit_cell_energy_1j", "unit_power_1w"}
        and identity.get("result_mode_normalization") == normalization
        and bool(dataset)
        and identity.get("result_dataset_tag") == dataset
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(result) == identity.get(source)
            for source, result in (
                ("periodic_mesh_map_sha256", "result_periodic_mesh_map_sha256"),
                ("mode_field_sha256", "result_mode_field_sha256"),
            )
        )
    )


def _rotating_sliding_interface_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "rotating_sliding_interface_sector_pitch_azimuth_interpolation_frame_periodicity_mesh_torque_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        pitch = float(identity.get("sector_pitch_deg"))
        result_pitch = float(identity.get("result_sector_pitch_deg"))
        count = int(identity.get("sector_count"))
        result_count = int(identity.get("result_sector_count"))
        origin = float(identity.get("azimuth_origin_deg"))
        result_origin = float(identity.get("result_azimuth_origin_deg"))
        phase = float(identity.get("periodic_phase_deg"))
        result_phase = float(identity.get("result_periodic_phase_deg"))
        samples = [
            float(value) for value in identity.get("azimuth_samples_deg", [])
        ]
        result_samples = [
            float(value) for value in identity.get("result_azimuth_samples_deg", [])
        ]
        torque = [float(value) for value in identity.get("torque_nm", [])]
        result_torque = [
            float(value) for value in identity.get("result_torque_nm", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sliding_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "sector_sliding_generation",
                "azimuth_sliding_generation",
                "interpolation_sliding_generation",
                "frame_sliding_generation",
                "periodicity_sliding_generation",
                "mesh_sliding_generation",
                "result_sliding_generation",
            )
        )
        and all(
            math.isfinite(value)
            for value in (pitch, result_pitch, origin, result_origin, phase, result_phase)
        )
        and pitch > 0.0
        and count > 1
        and math.isclose(pitch * count, 360.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_pitch == pitch
        and result_count == count
        and 0.0 <= origin < pitch
        and result_origin == origin
        and bool(str(identity.get("source_interface_tag") or ""))
        and bool(str(identity.get("target_interface_tag") or ""))
        and identity.get("source_interface_tag") != identity.get("target_interface_tag")
        and identity.get("result_source_interface_tag")
        == identity.get("source_interface_tag")
        and identity.get("result_target_interface_tag")
        == identity.get("target_interface_tag")
        and identity.get("interpolation") == "conservative_mortar_azimuth"
        and identity.get("result_interpolation") == identity.get("interpolation")
        and identity.get("rotor_frame") == "rotor_cylindrical"
        and identity.get("result_rotor_frame") == identity.get("rotor_frame")
        and result_phase == phase
        and len(samples) >= 3
        and samples[0] == origin
        and math.isclose(samples[-1] - samples[0], pitch, abs_tol=1.0e-12)
        and all(right > left for left, right in zip(samples, samples[1:]))
        and result_samples == samples
        and len(torque) == len(samples)
        and all(math.isfinite(value) for value in torque)
        and math.isclose(torque[0], torque[-1], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_torque == torque
        and _is_sha256(str(identity.get("sliding_mesh_sha256") or ""))
        and identity.get("result_sliding_mesh_sha256")
        == identity.get("sliding_mesh_sha256")
        and _is_sha256(str(identity.get("torque_result_sha256") or ""))
        and identity.get("accepted_torque_result_sha256")
        == identity.get("torque_result_sha256")
    )


def _acoustic_radiation_impedance_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "acoustic_radiation_impedance_modal_trace_reference_area_pressure_velocity_power_frequency_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        modes = [int(value) for value in identity.get("mode_indices", [])]
        result_modes = [
            int(value) for value in identity.get("result_mode_indices", [])
        ]
        area = float(identity.get("reference_area_m2"))
        result_area = float(identity.get("result_reference_area_m2"))
        frequencies = [
            float(value) for value in identity.get("frequency_grid_hz", [])
        ]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
        impedance = [
            [float(value) for value in row]
            for row in identity.get("radiation_impedance_ri", [])
        ]
        result_impedance = [
            [float(value) for value in row]
            for row in identity.get("result_radiation_impedance_ri", [])
        ]
        power = [
            float(value) for value in identity.get("outward_power_flux_w", [])
        ]
        result_power = [
            float(value)
            for value in identity.get("result_outward_power_flux_w", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("radiation_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "modal_radiation_generation",
                "trace_radiation_generation",
                "area_radiation_generation",
                "convention_radiation_generation",
                "power_radiation_generation",
                "frequency_radiation_generation",
                "result_radiation_generation",
            )
        )
        and bool(str(identity.get("modal_basis_id") or ""))
        and identity.get("result_modal_basis_id") == identity.get("modal_basis_id")
        and bool(modes)
        and all(value > 0 for value in modes)
        and len(set(modes)) == len(modes)
        and result_modes == modes
        and identity.get("trace_projection") == "l2_p1_boundary"
        and identity.get("result_trace_projection")
        == identity.get("trace_projection")
        and math.isfinite(area)
        and area > 0.0
        and result_area == area
        and identity.get("pressure_velocity_convention")
        == "outward_positive_velocity"
        and identity.get("result_pressure_velocity_convention")
        == identity.get("pressure_velocity_convention")
        and len(frequencies) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and len(impedance) == len(power) == len(frequencies)
        and all(
            len(row) == 2
            and all(math.isfinite(value) for value in row)
            and row[0] >= 0.0
            for row in impedance
        )
        and result_impedance == impedance
        and all(math.isfinite(value) and value >= 0.0 for value in power)
        and result_power == power
        and _is_sha256(str(identity.get("radiation_mesh_sha256") or ""))
        and identity.get("result_radiation_mesh_sha256")
        == identity.get("radiation_mesh_sha256")
        and _is_sha256(str(identity.get("radiation_result_sha256") or ""))
        and identity.get("accepted_radiation_result_sha256")
        == identity.get("radiation_result_sha256")
    )


def _joule_heat_energy_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "joule_heat_source_current_density_resistivity_temperature_frame_time_average_energy_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        window = [float(value) for value in identity.get("averaging_window_s", [])]
        result_window = [
            float(value) for value in identity.get("result_averaging_window_s", [])
        ]
        electric_loss = float(identity.get("electric_loss_w"))
        result_electric_loss = float(identity.get("result_electric_loss_w"))
        heat_integral = float(identity.get("heat_source_integral_w"))
        result_heat_integral = float(identity.get("result_heat_source_integral_w"))
        tolerance = float(identity.get("energy_balance_relative_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("joule_generation") or "")
    scale = max(abs(electric_loss), abs(heat_integral), 1.0e-300)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "mapping_joule_generation",
                "resistivity_joule_generation",
                "temperature_joule_generation",
                "frame_joule_generation",
                "averaging_joule_generation",
                "energy_joule_generation",
                "mesh_joule_generation",
                "result_joule_generation",
            )
        )
        and all(
            bool(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("current_density_field_id", "result_current_density_field_id"),
                ("temperature_field_id", "result_temperature_field_id"),
                ("resistivity_model_id", "result_resistivity_model_id"),
                ("source_frame", "result_source_frame"),
            )
        )
        and len(window) == 2
        and all(math.isfinite(value) for value in window)
        and window[1] > window[0] >= 0.0
        and result_window == window
        and identity.get("time_average_method")
        in {"trapezoidal_period_average", "exact_period_average"}
        and identity.get("result_time_average_method")
        == identity.get("time_average_method")
        and all(
            math.isfinite(value) and value >= 0.0
            for value in (
                electric_loss,
                result_electric_loss,
                heat_integral,
                result_heat_integral,
            )
        )
        and result_electric_loss == electric_loss
        and result_heat_integral == heat_integral
        and math.isfinite(tolerance)
        and 0.0 < tolerance < 1.0
        and abs(electric_loss - heat_integral) / scale <= tolerance
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("coupled_mesh_sha256", "result_coupled_mesh_sha256"),
                ("joule_heat_result_sha256", "accepted_joule_heat_result_sha256"),
            )
        )
    )


def _nonlinear_eigenmode_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "nonlinear_eigenmode_continuation_parameter_normalization_phase_mac_branch_eigenvalue_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        parameters = [
            float(value) for value in identity.get("continuation_parameter_values", [])
        ]
        result_parameters = [
            float(value)
            for value in identity.get("result_continuation_parameter_values", [])
        ]
        references = [
            int(value) for value in identity.get("mac_reference_branch_ids", [])
        ]
        result_references = [
            int(value)
            for value in identity.get("result_mac_reference_branch_ids", [])
        ]
        branches = [
            [int(value) for value in row]
            for row in identity.get("mode_branch_ids", [])
        ]
        result_branches = [
            [int(value) for value in row]
            for row in identity.get("result_mode_branch_ids", [])
        ]
        eigenvalues = [
            [[float(value) for value in pair] for pair in row]
            for row in identity.get("eigenvalues_ri", [])
        ]
        result_eigenvalues = [
            [[float(value) for value in pair] for pair in row]
            for row in identity.get("result_eigenvalues_ri", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("eigenmode_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "continuation_eigenmode_generation",
                "normalization_eigenmode_generation",
                "phase_eigenmode_generation",
                "mac_eigenmode_generation",
                "branch_eigenmode_generation",
                "eigenvalue_eigenmode_generation",
                "mesh_eigenmode_generation",
                "result_eigenmode_generation",
            )
        )
        and bool(str(identity.get("continuation_parameter_name") or ""))
        and identity.get("result_continuation_parameter_name")
        == identity.get("continuation_parameter_name")
        and len(parameters) >= 2
        and all(math.isfinite(value) for value in parameters)
        and all(right > left for left, right in zip(parameters, parameters[1:]))
        and result_parameters == parameters
        and identity.get("mode_normalization") in {"unit_mass", "unit_energy"}
        and identity.get("result_mode_normalization")
        == identity.get("mode_normalization")
        and bool(str(identity.get("phase_anchor_dof") or ""))
        and identity.get("result_phase_anchor_dof") == identity.get("phase_anchor_dof")
        and bool(references)
        and len(set(references)) == len(references)
        and all(value > 0 for value in references)
        and result_references == references
        and len(branches) == len(parameters)
        and all(row == references for row in branches)
        and result_branches == branches
        and len(eigenvalues) == len(parameters)
        and all(
            len(row) == len(references)
            and all(
                len(pair) == 2 and all(math.isfinite(value) for value in pair)
                for pair in row
            )
            for row in eigenvalues
        )
        and result_eigenvalues == eigenvalues
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("mac_assignment_sha256", "result_mac_assignment_sha256"),
                ("eigenmode_mesh_sha256", "result_eigenmode_mesh_sha256"),
                ("eigenmode_result_sha256", "accepted_eigenmode_result_sha256"),
            )
        )
    )


def _frequency_time_reconstruction_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "frequency_time_hermitian_spacing_window_group_delay_parseval_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        frequencies = [float(value) for value in identity.get("frequencies_hz", [])]
        result_frequencies = [
            float(value) for value in identity.get("result_frequencies_hz", [])
        ]
        spectrum = [
            [float(value) for value in pair]
            for pair in identity.get("spectrum_ri", [])
        ]
        result_spectrum = [
            [float(value) for value in pair]
            for pair in identity.get("result_spectrum_ri", [])
        ]
        spacing = float(identity.get("frequency_spacing_hz"))
        result_spacing = float(identity.get("result_frequency_spacing_hz"))
        coherent_gain = float(identity.get("window_coherent_gain"))
        result_coherent_gain = float(identity.get("result_window_coherent_gain"))
        group_delay = float(identity.get("group_delay_s"))
        result_group_delay = float(identity.get("result_group_delay_s"))
        time_origin = float(identity.get("time_origin_s"))
        result_time_origin = float(identity.get("result_time_origin_s"))
        frequency_energy = float(identity.get("frequency_energy"))
        time_energy = float(identity.get("time_energy"))
        tolerance = float(identity.get("parseval_relative_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("transform_generation") or "")
    energy_scale = max(abs(frequency_energy), abs(time_energy), 1.0e-300)
    spacing_ok = len(frequencies) >= 2 and all(
        abs((right - left) - spacing) <= 1.0e-12 * max(abs(spacing), 1.0)
        for left, right in zip(frequencies, frequencies[1:])
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "spectrum_transform_generation",
                "window_transform_generation",
                "delay_transform_generation",
                "energy_transform_generation",
                "mesh_transform_generation",
                "result_transform_generation",
            )
        )
        and all(math.isfinite(value) and value >= 0.0 for value in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and math.isfinite(spacing)
        and spacing > 0.0
        and result_spacing == spacing
        and spacing_ok
        and len(spectrum) == len(frequencies)
        and all(
            len(pair) == 2 and all(math.isfinite(value) for value in pair)
            for pair in spectrum
        )
        and spectrum[0][1] == 0.0
        and spectrum[-1][1] == 0.0
        and result_spectrum == spectrum
        and identity.get("hermitian_completion")
        == "conjugate_negative_frequencies"
        and identity.get("result_hermitian_completion")
        == identity.get("hermitian_completion")
        and identity.get("window_name") in {"hann_periodic", "none"}
        and identity.get("result_window_name") == identity.get("window_name")
        and math.isfinite(coherent_gain)
        and 0.0 < coherent_gain <= 1.0
        and result_coherent_gain == coherent_gain
        and all(
            math.isfinite(value)
            for value in (
                group_delay,
                result_group_delay,
                time_origin,
                result_time_origin,
            )
        )
        and result_group_delay == group_delay
        and result_time_origin == time_origin
        and math.isfinite(tolerance)
        and 0.0 < tolerance < 1.0
        and abs(frequency_energy - time_energy) / energy_scale <= tolerance
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("transform_mesh_sha256", "result_transform_mesh_sha256"),
                ("time_trace_sha256", "accepted_time_trace_sha256"),
            )
        )
    )


def _rotating_force_balance_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "rotating_force_virtual_work_stress_phase_frame_lever_angle_power_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        lever = [float(value) for value in identity.get("lever_arm_m", [])]
        result_lever = [
            float(value) for value in identity.get("result_lever_arm_m", [])
        ]
        angles = [
            float(value) for value in identity.get("mechanical_angles_rad", [])
        ]
        result_angles = [
            float(value)
            for value in identity.get("result_mechanical_angles_rad", [])
        ]
        virtual_work = [
            float(value) for value in identity.get("virtual_work_torque_nm", [])
        ]
        stress_tensor = [
            float(value) for value in identity.get("stress_tensor_torque_nm", [])
        ]
        torque_tolerance = float(identity.get("torque_relative_tolerance"))
        mechanical_power = float(identity.get("mechanical_power_w"))
        airgap_power = float(identity.get("airgap_power_w"))
        power_tolerance = float(identity.get("power_relative_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("force_generation") or "")
    torque_scale = max(
        max((abs(value) for value in virtual_work), default=0.0),
        max((abs(value) for value in stress_tensor), default=0.0),
        1.0e-300,
    )
    power_scale = max(abs(mechanical_power), abs(airgap_power), 1.0e-300)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "virtual_work_force_generation",
                "stress_force_generation",
                "phase_force_generation",
                "frame_force_generation",
                "angle_force_generation",
                "power_force_generation",
                "mesh_force_generation",
                "result_force_generation",
            )
        )
        and identity.get("phasor_convention") == "exp_positive_jwt_rms"
        and identity.get("result_phasor_convention")
        == identity.get("phasor_convention")
        and identity.get("coordinate_frame") == "rotor_material"
        and identity.get("result_coordinate_frame")
        == identity.get("coordinate_frame")
        and len(lever) == 3
        and all(math.isfinite(value) for value in lever)
        and result_lever == lever
        and len(angles) >= 2
        and all(math.isfinite(value) for value in angles)
        and all(right > left for left, right in zip(angles, angles[1:]))
        and result_angles == angles
        and len(virtual_work) == len(stress_tensor) == len(angles)
        and all(math.isfinite(value) for value in virtual_work + stress_tensor)
        and math.isfinite(torque_tolerance)
        and 0.0 < torque_tolerance < 1.0
        and max(
            (
                abs(left - right)
                for left, right in zip(virtual_work, stress_tensor, strict=True)
            ),
            default=math.inf,
        )
        / torque_scale
        <= torque_tolerance
        and all(math.isfinite(value) for value in (mechanical_power, airgap_power))
        and math.isfinite(power_tolerance)
        and 0.0 < power_tolerance < 1.0
        and abs(mechanical_power - airgap_power) / power_scale <= power_tolerance
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("force_mesh_sha256", "result_force_mesh_sha256"),
                ("force_result_sha256", "accepted_force_result_sha256"),
            )
        )
    )


def _nonlinear_segregated_closure_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "nonlinear_segregated_group_relaxation_residual_jacobian_continuation_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        groups = [str(value) for value in identity.get("segregated_group_order", [])]
        result_groups = [
            str(value) for value in identity.get("result_segregated_group_order", [])
        ]
        relaxation = [
            float(value) for value in identity.get("relaxation_schedule", [])
        ]
        result_relaxation = [
            float(value) for value in identity.get("result_relaxation_schedule", [])
        ]
        residuals = [float(value) for value in identity.get("residual_norms", [])]
        accepted_residuals = [
            float(value) for value in identity.get("accepted_residual_norms", [])
        ]
        residual_tolerance = float(identity.get("residual_relative_tolerance"))
        continuation = float(identity.get("continuation_parameter"))
        accepted_continuation = float(
            identity.get("accepted_continuation_parameter")
        )
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("nonlinear_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "segregated_group_generation",
                "relaxation_generation",
                "residual_generation",
                "jacobian_generation",
                "continuation_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and bool(groups)
        and len(set(groups)) == len(groups)
        and result_groups == groups
        and len(relaxation) == len(groups)
        and all(math.isfinite(value) and 0.0 < value <= 1.0 for value in relaxation)
        and result_relaxation == relaxation
        and len(residuals) >= 2
        and all(math.isfinite(value) and value >= 0.0 for value in residuals)
        and accepted_residuals == residuals
        and all(right <= left for left, right in zip(residuals, residuals[1:]))
        and math.isfinite(residual_tolerance)
        and 0.0 < residual_tolerance < 1.0
        and residuals[-1] <= residual_tolerance
        and math.isfinite(continuation)
        and accepted_continuation == continuation
        and bool(str(identity.get("continuation_unit") or ""))
        and identity.get("accepted_continuation_unit")
        == identity.get("continuation_unit")
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("jacobian_sha256", "accepted_jacobian_sha256"),
                ("nonlinear_mesh_sha256", "accepted_nonlinear_mesh_sha256"),
                ("nonlinear_solution_sha256", "accepted_nonlinear_solution_sha256"),
            )
        )
    )


def _degenerate_eigenmode_closure_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "degenerate_eigenmode_subspace_phase_normalization_participation_mass_mesh_owner_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        eigenvalues = [float(value) for value in identity.get("eigenvalues", [])]
        result_eigenvalues = [
            float(value) for value in identity.get("result_eigenvalues", [])
        ]
        anchors = [int(value) for value in identity.get("phase_anchor_dofs", [])]
        result_anchors = [
            int(value) for value in identity.get("result_phase_anchor_dofs", [])
        ]
        participation = [
            float(value) for value in identity.get("participation_factors", [])
        ]
        result_participation = [
            float(value) for value in identity.get("result_participation_factors", [])
        ]
        masses = [float(value) for value in identity.get("effective_masses_kg", [])]
        result_masses = [
            float(value) for value in identity.get("result_effective_masses_kg", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("mode_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "subspace_generation",
                "phase_generation",
                "normalization_generation",
                "participation_generation",
                "mass_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and len(eigenvalues) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in eigenvalues)
        and result_eigenvalues == eigenvalues
        and _is_sha256(str(identity.get("degenerate_subspace_sha256") or ""))
        and identity.get("result_degenerate_subspace_sha256")
        == identity.get("degenerate_subspace_sha256")
        and len(anchors) == len(eigenvalues)
        and len(set(anchors)) == len(anchors)
        and all(value >= 0 for value in anchors)
        and result_anchors == anchors
        and identity.get("normalization") == "mass_orthonormal"
        and identity.get("result_normalization") == identity.get("normalization")
        and len(participation) == len(masses) == len(eigenvalues)
        and all(math.isfinite(value) for value in participation + masses)
        and all(value >= 0.0 for value in masses)
        and result_participation == participation
        and result_masses == masses
        and all(
            math.isclose(mass, factor * factor, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for factor, mass in zip(participation, masses, strict=True)
        )
        and bool(str(identity.get("mode_owner") or ""))
        and identity.get("result_mode_owner") == identity.get("mode_owner")
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("eigenmode_mesh_sha256", "result_eigenmode_mesh_sha256"),
                ("eigenmode_result_sha256", "accepted_eigenmode_result_sha256"),
            )
        )
    )


def _contact_complementarity_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "contact_gap_pressure_active_set_friction_dissipation_normal_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        active_ids = [int(value) for value in identity.get("active_contact_ids", [])]
        result_active_ids = [
            int(value) for value in identity.get("result_active_contact_ids", [])
        ]
        gaps = [float(value) for value in identity.get("normal_gap_m", [])]
        result_gaps = [
            float(value) for value in identity.get("result_normal_gap_m", [])
        ]
        pressures = [
            float(value) for value in identity.get("normal_pressure_pa", [])
        ]
        result_pressures = [
            float(value) for value in identity.get("result_normal_pressure_pa", [])
        ]
        slips = [float(value) for value in identity.get("tangential_slip_m", [])]
        result_slips = [
            float(value) for value in identity.get("result_tangential_slip_m", [])
        ]
        tractions = [
            float(value) for value in identity.get("friction_traction_pa", [])
        ]
        result_tractions = [
            float(value) for value in identity.get("result_friction_traction_pa", [])
        ]
        areas = [float(value) for value in identity.get("contact_area_m2", [])]
        result_areas = [
            float(value) for value in identity.get("result_contact_area_m2", [])
        ]
        coefficient = float(identity.get("friction_coefficient"))
        result_coefficient = float(identity.get("result_friction_coefficient"))
        dissipation = float(identity.get("friction_dissipation_j"))
        result_dissipation = float(identity.get("result_friction_dissipation_j"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("contact_generation") or "")
    count = len(active_ids)
    recomputed_dissipation = sum(
        abs(traction * slip) * area
        for traction, slip, area in zip(tractions, slips, areas, strict=True)
    ) if count and len(tractions) == len(slips) == len(areas) == count else math.inf
    dissipation_scale = max(abs(dissipation), abs(recomputed_dissipation), 1.0e-30)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "gap_generation",
                "pressure_generation",
                "active_set_generation",
                "friction_generation",
                "dissipation_generation",
                "normal_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and bool(str(identity.get("contact_pair") or ""))
        and identity.get("result_contact_pair") == identity.get("contact_pair")
        and count > 0
        and active_ids == sorted(set(active_ids))
        and all(value > 0 for value in active_ids)
        and result_active_ids == active_ids
        and len(gaps) == len(pressures) == len(slips) == len(tractions) == len(areas) == count
        and all(math.isfinite(value) for value in gaps + pressures + slips + tractions + areas)
        and all(abs(gap) <= 1.0e-9 for gap in gaps)
        and all(pressure >= 0.0 for pressure in pressures)
        and all(abs(gap * pressure) <= 1.0e-9 * max(pressure, 1.0) for gap, pressure in zip(gaps, pressures, strict=True))
        and all(area > 0.0 for area in areas)
        and result_gaps == gaps
        and result_pressures == pressures
        and result_slips == slips
        and result_tractions == tractions
        and result_areas == areas
        and math.isfinite(coefficient)
        and coefficient >= 0.0
        and result_coefficient == coefficient
        and all(
            abs(traction) <= coefficient * pressure + 1.0e-9 * max(pressure, 1.0)
            for traction, pressure in zip(tractions, pressures, strict=True)
        )
        and math.isfinite(dissipation)
        and dissipation >= 0.0
        and abs(dissipation - recomputed_dissipation) <= 1.0e-12 * dissipation_scale
        and result_dissipation == dissipation
        and identity.get("normal_orientation") == "outward_slave_to_master"
        and identity.get("result_normal_orientation") == identity.get("normal_orientation")
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("contact_mesh_sha256", "result_contact_mesh_sha256"),
                ("contact_result_sha256", "accepted_contact_result_sha256"),
            )
        )
    )


def _field_circuit_dae_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "field_circuit_dae_charge_current_event_energy_time_dataset_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        times = [float(value) for value in identity.get("time_s", [])]
        result_times = [float(value) for value in identity.get("result_time_s", [])]
        event_time = float(identity.get("switch_event_time_s"))
        result_event_time = float(identity.get("result_switch_event_time_s"))
        charge = [float(value) for value in identity.get("charge_c", [])]
        result_charge = [float(value) for value in identity.get("result_charge_c", [])]
        integrated_current = [
            float(value) for value in identity.get("integrated_current_c", [])
        ]
        result_integrated_current = [
            float(value) for value in identity.get("result_integrated_current_c", [])
        ]
        residuals = [
            float(value) for value in identity.get("algebraic_residual_c", [])
        ]
        accepted_residuals = [
            float(value) for value in identity.get("accepted_algebraic_residual_c", [])
        ]
        tolerance = float(identity.get("algebraic_tolerance_c"))
        energy_before = float(identity.get("stored_energy_before_j"))
        result_energy_before = float(identity.get("result_stored_energy_before_j"))
        energy_after = float(identity.get("stored_energy_after_j"))
        result_energy_after = float(identity.get("result_stored_energy_after_j"))
        dissipation = float(identity.get("switch_dissipation_j"))
        result_dissipation = float(identity.get("result_switch_dissipation_j"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("dae_generation") or "")
    count = len(times)
    charge_scale = max(max((abs(value) for value in charge), default=0.0), 1.0e-30)
    energy_scale = max(abs(energy_before), 1.0e-30)
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "charge_generation",
                "current_generation",
                "event_generation",
                "energy_generation",
                "time_generation",
                "dataset_generation",
                "result_generation",
            )
        )
        and count >= 3
        and all(math.isfinite(value) and value >= 0.0 for value in times)
        and all(right > left for left, right in zip(times, times[1:]))
        and result_times == times
        and math.isfinite(event_time)
        and any(math.isclose(value, event_time, rel_tol=0.0, abs_tol=1.0e-15) for value in times)
        and result_event_time == event_time
        and identity.get("event_side") == "right_limit_after_event"
        and identity.get("result_event_side") == identity.get("event_side")
        and len(charge) == len(integrated_current) == len(residuals) == count
        and all(math.isfinite(value) for value in charge + integrated_current + residuals)
        and all(abs(left - right) <= 1.0e-12 * charge_scale for left, right in zip(charge, integrated_current, strict=True))
        and result_charge == charge
        and result_integrated_current == integrated_current
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and all(abs(value) <= tolerance for value in residuals)
        and accepted_residuals == residuals
        and identity.get("current_sign_convention") == "positive_into_field_device"
        and identity.get("result_current_sign_convention") == identity.get("current_sign_convention")
        and all(math.isfinite(value) and value >= 0.0 for value in (energy_before, energy_after, dissipation))
        and abs(energy_before - energy_after - dissipation) <= 1.0e-12 * energy_scale
        and result_energy_before == energy_before
        and result_energy_after == energy_after
        and result_dissipation == dissipation
        and bool(str(identity.get("dataset_owner") or ""))
        and identity.get("result_dataset_owner") == identity.get("dataset_owner")
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("dae_dataset_sha256", "result_dae_dataset_sha256"),
                ("dae_result_sha256", "accepted_dae_result_sha256"),
            )
        )
    )


def _thermoelastic_frequency_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "thermoelastic_frequency_reference_temperature_prestress_linearization_mesh_dataset_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        temperature = float(identity.get("reference_temperature_k"))
        result_temperature = float(identity.get("result_reference_temperature_k"))
        frequencies = [float(value) for value in identity.get("frequency_grid_hz", [])]
        result_frequencies = [
            float(value) for value in identity.get("result_frequency_grid_hz", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("thermoelastic_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "temperature_thermoelastic_generation",
                "prestress_thermoelastic_generation",
                "linearization_thermoelastic_generation",
                "mesh_thermoelastic_generation",
                "dataset_thermoelastic_generation",
                "result_thermoelastic_generation",
            )
        )
        and math.isfinite(temperature)
        and temperature > 0.0
        and result_temperature == temperature
        and len(frequencies) >= 2
        and all(math.isfinite(value) and value > 0.0 for value in frequencies)
        and all(right > left for left, right in zip(frequencies, frequencies[1:]))
        and result_frequencies == frequencies
        and bool(str(identity.get("dataset_tag") or ""))
        and identity.get("result_dataset_tag") == identity.get("dataset_tag")
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("prestress_state_sha256", "result_prestress_state_sha256"),
                ("linearization_state_sha256", "result_linearization_state_sha256"),
                ("thermal_mesh_sha256", "result_thermal_mesh_sha256"),
                ("structural_mesh_sha256", "result_structural_mesh_sha256"),
                ("frequency_response_sha256", "accepted_frequency_response_sha256"),
            )
        )
    )


def _field_circuit_power_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "field_circuit_coil_terminal_orientation_current_sign_gauge_power_balance_mesh_solution_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        terminals = [str(value) for value in identity.get("coil_terminal_ids", [])]
        result_terminals = [
            str(value) for value in identity.get("result_coil_terminal_ids", [])
        ]
        orientation = [
            int(value) for value in identity.get("terminal_orientation_signs", [])
        ]
        result_orientation = [
            int(value)
            for value in identity.get("result_terminal_orientation_signs", [])
        ]
        current_signs = [int(value) for value in identity.get("branch_current_signs", [])]
        result_current_signs = [
            int(value) for value in identity.get("result_branch_current_signs", [])
        ]
        field_power = [
            float(value) for value in identity.get("field_complex_power_va_ri", [])
        ]
        result_field_power = [
            float(value)
            for value in identity.get("result_field_complex_power_va_ri", [])
        ]
        circuit_power = [
            float(value) for value in identity.get("circuit_complex_power_va_ri", [])
        ]
        result_circuit_power = [
            float(value)
            for value in identity.get("result_circuit_complex_power_va_ri", [])
        ]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("coupling_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "terminal_coupling_generation",
                "orientation_coupling_generation",
                "current_sign_coupling_generation",
                "gauge_coupling_generation",
                "power_coupling_generation",
                "mesh_coupling_generation",
                "solution_coupling_generation",
                "result_coupling_generation",
            )
        )
        and len(terminals) == 2
        and all(terminals)
        and len(set(terminals)) == 2
        and result_terminals == terminals
        and orientation == [1, -1]
        and result_orientation == orientation
        and current_signs == [1, -1]
        and result_current_signs == current_signs
        and identity.get("gauge_id") == "magnetic_vector_potential_coulomb"
        and identity.get("result_gauge_id") == identity.get("gauge_id")
        and bool(str(identity.get("circuit_branch_id") or ""))
        and identity.get("result_circuit_branch_id") == identity.get("circuit_branch_id")
        and len(field_power) == len(circuit_power) == 2
        and all(math.isfinite(value) for value in field_power + circuit_power)
        and result_field_power == field_power
        and result_circuit_power == circuit_power
        and all(
            math.isclose(field + circuit, 0.0, rel_tol=1.0e-9, abs_tol=1.0e-12)
            for field, circuit in zip(field_power, circuit_power)
        )
        and _is_sha256(str(identity.get("coupled_mesh_sha256") or ""))
        and identity.get("result_coupled_mesh_sha256")
        == identity.get("coupled_mesh_sha256")
        and _is_sha256(str(identity.get("coupled_solution_sha256") or ""))
        and identity.get("accepted_coupled_solution_sha256")
        == identity.get("coupled_solution_sha256")
    )


def _nonlinear_arclength_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        previous = [float(value) for value in identity.get("previous_augmented_state", [])]
        result_previous = [
            float(value) for value in identity.get("result_previous_augmented_state", [])
        ]
        tangent = [float(value) for value in identity.get("predictor_tangent", [])]
        result_tangent = [
            float(value) for value in identity.get("result_predictor_tangent", [])
        ]
        step = float(identity.get("arclength_step"))
        result_step = float(identity.get("result_arclength_step"))
        predictor = [
            float(value) for value in identity.get("predictor_augmented_state", [])
        ]
        result_predictor = [
            float(value) for value in identity.get("result_predictor_augmented_state", [])
        ]
        corrected = [
            float(value) for value in identity.get("corrected_augmented_state", [])
        ]
        result_corrected = [
            float(value) for value in identity.get("result_corrected_augmented_state", [])
        ]
        residual = float(identity.get("corrected_residual_norm"))
        result_residual = float(identity.get("result_corrected_residual_norm"))
        tolerance = float(identity.get("residual_tolerance"))
        result_tolerance = float(identity.get("result_residual_tolerance"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("continuation_generation") or "")
    count = len(previous)
    tangent_norm = math.sqrt(sum(value * value for value in tangent))
    expected_predictor = [
        value + step * direction for value, direction in zip(previous, tangent, strict=True)
    ] if count and len(tangent) == count else []
    correction = [
        value - old for value, old in zip(corrected, previous, strict=True)
    ] if len(corrected) == count else []
    projected_arclength = sum(
        delta * direction for delta, direction in zip(correction, tangent, strict=True)
    ) if len(correction) == count and len(tangent) == count else math.inf
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "arclength_generation",
                "tangent_generation",
                "branch_generation",
                "turning_generation",
                "residual_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and count >= 2
        and len(tangent) == len(predictor) == len(corrected) == count
        and all(math.isfinite(value) for value in previous + tangent + predictor + corrected)
        and result_previous == previous
        and result_tangent == tangent
        and result_predictor == predictor
        and result_corrected == corrected
        and math.isclose(tangent_norm, 1.0, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isfinite(step)
        and step > 0.0
        and result_step == step
        and all(
            math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for actual, expected in zip(predictor, expected_predictor, strict=True)
        )
        and math.isclose(projected_arclength, step, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and bool(str(identity.get("branch_id") or ""))
        and identity.get("result_branch_id") == identity.get("branch_id")
        and identity.get("turning_point_side")
        == "pre_turn_positive_parameter_tangent"
        and tangent[-1] > 0.0
        and identity.get("result_turning_point_side")
        == identity.get("turning_point_side")
        and math.isfinite(residual)
        and residual >= 0.0
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and residual <= tolerance
        and result_residual == residual
        and result_tolerance == tolerance
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("continuation_mesh_sha256", "result_continuation_mesh_sha256"),
                ("continuation_result_sha256", "accepted_continuation_result_sha256"),
            )
        )
    )


def _electrochemical_conservation_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        species = [str(value) for value in identity.get("species_order", [])]
        result_species = [str(value) for value in identity.get("result_species_order", [])]
        charge = [float(value) for value in identity.get("charge_numbers", [])]
        result_charge = [float(value) for value in identity.get("result_charge_numbers", [])]
        mass = [float(value) for value in identity.get("molar_mass_basis", [])]
        result_mass = [float(value) for value in identity.get("result_molar_mass_basis", [])]
        stoichiometry = [
            float(value) for value in identity.get("reaction_stoichiometry", [])
        ]
        result_stoichiometry = [
            float(value) for value in identity.get("result_reaction_stoichiometry", [])
        ]
        extent = float(identity.get("reaction_extent_mol"))
        result_extent = float(identity.get("result_reaction_extent_mol"))
        initial = [float(value) for value in identity.get("initial_inventory_mol", [])]
        result_initial = [
            float(value) for value in identity.get("result_initial_inventory_mol", [])
        ]
        final = [float(value) for value in identity.get("final_inventory_mol", [])]
        result_final = [
            float(value) for value in identity.get("result_final_inventory_mol", [])
        ]
        boundary = [
            float(value) for value in identity.get("integrated_boundary_flux_mol", [])
        ]
        result_boundary = [
            float(value)
            for value in identity.get("result_integrated_boundary_flux_mol", [])
        ]
        current = float(identity.get("integrated_electric_current_c"))
        result_current = float(identity.get("result_integrated_electric_current_c"))
        energy_initial = float(identity.get("initial_free_energy_j"))
        result_energy_initial = float(identity.get("result_initial_free_energy_j"))
        energy_final = float(identity.get("final_free_energy_j"))
        result_energy_final = float(identity.get("result_final_free_energy_j"))
        dissipation = float(identity.get("dissipated_free_energy_j"))
        result_dissipation = float(identity.get("result_dissipated_free_energy_j"))
        times = [float(value) for value in identity.get("time_s", [])]
        result_times = [float(value) for value in identity.get("result_time_s", [])]
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("electrochemical_generation") or "")
    count = len(species)
    expected_final = [
        old + extent * coefficient + flux
        for old, coefficient, flux in zip(initial, stoichiometry, boundary, strict=True)
    ] if count and len(initial) == len(stoichiometry) == len(boundary) == count else []
    initial_mass = sum(value * weight for value, weight in zip(initial, mass, strict=True)) if len(initial) == len(mass) == count else math.inf
    final_mass = sum(value * weight for value, weight in zip(final, mass, strict=True)) if len(final) == len(mass) == count else math.inf
    boundary_mass = sum(value * weight for value, weight in zip(boundary, mass, strict=True)) if len(boundary) == len(mass) == count else math.inf
    initial_charge = sum(value * number for value, number in zip(initial, charge, strict=True)) if len(initial) == len(charge) == count else math.inf
    final_charge = sum(value * number for value, number in zip(final, charge, strict=True)) if len(final) == len(charge) == count else math.inf
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "species_generation",
                "flux_generation",
                "charge_generation",
                "mass_generation",
                "reaction_generation",
                "energy_generation",
                "time_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and count >= 2
        and all(species)
        and len(set(species)) == count
        and result_species == species
        and len(charge) == len(mass) == len(stoichiometry) == len(initial) == len(final) == len(boundary) == count
        and all(math.isfinite(value) for value in charge + mass + stoichiometry + initial + final + boundary)
        and all(value > 0.0 for value in mass)
        and all(value >= 0.0 for value in initial + final)
        and result_charge == charge
        and result_mass == mass
        and result_stoichiometry == stoichiometry
        and math.isfinite(extent)
        and extent >= 0.0
        and result_extent == extent
        and result_initial == initial
        and result_final == final
        and result_boundary == boundary
        and all(
            math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for actual, expected in zip(final, expected_final, strict=True)
        )
        and math.isclose(initial_mass + boundary_mass, final_mass, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(initial_charge, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(final_charge, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(current, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_current == current
        and all(math.isfinite(value) and value >= 0.0 for value in (energy_initial, energy_final, dissipation))
        and math.isclose(energy_initial - energy_final, dissipation, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_energy_initial == energy_initial
        and result_energy_final == energy_final
        and result_dissipation == dissipation
        and len(times) >= 2
        and all(math.isfinite(value) and value >= 0.0 for value in times)
        and all(right > left for left, right in zip(times, times[1:]))
        and result_times == times
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("electrochemical_mesh_sha256", "result_electrochemical_mesh_sha256"),
                ("electrochemical_result_sha256", "accepted_electrochemical_result_sha256"),
            )
        )
    )


def _multirate_electromechanical_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        electrical_time = [float(value) for value in identity.get("electrical_time_s", [])]
        result_electrical_time = [
            float(value) for value in identity.get("result_electrical_time_s", [])
        ]
        mechanical_time = [float(value) for value in identity.get("mechanical_time_s", [])]
        result_mechanical_time = [
            float(value) for value in identity.get("result_mechanical_time_s", [])
        ]
        event_time = float(identity.get("event_time_s"))
        result_event_time = float(identity.get("result_event_time_s"))
        electrical_energy = float(identity.get("electrical_input_energy_j"))
        result_electrical_energy = float(identity.get("result_electrical_input_energy_j"))
        mechanical_work = float(identity.get("mechanical_output_work_j"))
        result_mechanical_work = float(identity.get("result_mechanical_output_work_j"))
        dissipation = float(identity.get("dissipated_energy_j"))
        result_dissipation = float(identity.get("result_dissipated_energy_j"))
        tolerance = float(identity.get("energy_balance_tolerance_j"))
        result_tolerance = float(identity.get("result_energy_balance_tolerance_j"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("coupling_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "electrical_generation",
                "mechanical_generation",
                "event_generation",
                "timegrid_generation",
                "power_generation",
                "work_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and _increasing(electrical_time)
        and _increasing(mechanical_time)
        and all(value >= 0.0 for value in electrical_time + mechanical_time)
        and electrical_time[0] == mechanical_time[0]
        and electrical_time[-1] == mechanical_time[-1]
        and len(mechanical_time) == 2 * (len(electrical_time) - 1) + 1
        and result_electrical_time == electrical_time
        and result_mechanical_time == mechanical_time
        and event_time in electrical_time
        and event_time in mechanical_time
        and result_event_time == event_time
        and identity.get("event_interpolation_side") == "right_continuous_after_event"
        and identity.get("result_event_interpolation_side")
        == identity.get("event_interpolation_side")
        and bool(str(identity.get("substep_owner") or ""))
        and identity.get("result_substep_owner") == identity.get("substep_owner")
        and bool(str(identity.get("coordinate_frame") or ""))
        and identity.get("result_coordinate_frame") == identity.get("coordinate_frame")
        and all(
            math.isfinite(value) and value >= 0.0
            for value in (electrical_energy, mechanical_work, dissipation)
        )
        and result_electrical_energy == electrical_energy
        and result_mechanical_work == mechanical_work
        and result_dissipation == dissipation
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and abs(electrical_energy - mechanical_work - dissipation) <= tolerance
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("coupling_mesh_sha256", "result_coupling_mesh_sha256"),
                ("coupling_result_sha256", "accepted_coupling_result_sha256"),
            )
        )
    )


def _adjoint_sensitivity_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        design_scale = float(identity.get("design_scale"))
        result_design_scale = float(identity.get("result_design_scale"))
        adjoint_gradient = float(identity.get("adjoint_gradient"))
        chainrule_gradient = float(identity.get("chainrule_gradient"))
        finite_difference_gradient = float(identity.get("finite_difference_gradient"))
        tolerance = float(identity.get("gradient_tolerance"))
        result_tolerance = float(identity.get("result_gradient_tolerance"))
        perturbation = float(identity.get("fd_perturbation"))
        result_perturbation = float(identity.get("result_fd_perturbation"))
    except (TypeError, ValueError):
        return False
    generation = str(identity.get("sensitivity_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "objective_generation",
                "design_generation",
                "chainrule_generation",
                "constraint_generation",
                "fd_generation",
                "mesh_generation",
                "solution_generation",
                "result_generation",
            )
        )
        and all(
            bool(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("objective_tag", "result_objective_tag"),
                ("design_variable", "result_design_variable"),
                ("active_constraint", "result_active_constraint"),
            )
        )
        and math.isfinite(design_scale)
        and design_scale > 0.0
        and result_design_scale == design_scale
        and all(
            math.isfinite(value)
            for value in (adjoint_gradient, chainrule_gradient, finite_difference_gradient)
        )
        and math.isclose(adjoint_gradient, chainrule_gradient, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and abs(adjoint_gradient - finite_difference_gradient) <= tolerance
        and math.isfinite(perturbation)
        and 0.0 < perturbation <= 1.0e-2
        and result_perturbation == perturbation
        and all(
            _is_sha256(str(identity.get(source) or ""))
            and identity.get(target) == identity.get(source)
            for source, target in (
                ("sensitivity_mesh_sha256", "result_sensitivity_mesh_sha256"),
                ("primal_solution_sha256", "result_primal_solution_sha256"),
                ("gradient_result_sha256", "accepted_gradient_result_sha256"),
            )
        )
    )


def _magnetostatic_virtual_work_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        displacement = [float(value) for value in identity.get("displacement_m", [])]
        result_displacement = [
            float(value) for value in identity.get("result_displacement_m", [])
        ]
        coenergy = [float(value) for value in identity.get("coenergy_j", [])]
        result_coenergy = [
            float(value) for value in identity.get("result_coenergy_j", [])
        ]
        reported_central_force = float(identity.get("central_coenergy_force_n"))
        result_force = float(identity.get("result_force_n"))
        tolerance = float(identity.get("force_tolerance_n"))
        result_tolerance = float(identity.get("result_force_tolerance_n"))
    except (TypeError, ValueError):
        return False
    if len(displacement) != 3 or len(coenergy) != 3:
        return False
    central_force = (coenergy[2] - coenergy[0]) / (
        displacement[2] - displacement[0]
    )
    mesh_digests = identity.get("displaced_mesh_sha256")
    result_mesh_digests = identity.get("result_displaced_mesh_sha256")
    generation = str(identity.get("force_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "displacement_generation",
                "coenergy_generation",
                "current_generation",
                "mesh_generation",
                "frame_generation",
                "solution_generation",
                "result_generation",
            )
        )
        and _increasing(displacement)
        and displacement[1] == 0.0
        and math.isclose(
            displacement[2], -displacement[0], rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and all(math.isfinite(value) for value in displacement + coenergy)
        and result_displacement == displacement
        and result_coenergy == coenergy
        and identity.get("held_source_convention") == "constant_current"
        and identity.get("result_held_source_convention")
        == identity.get("held_source_convention")
        and identity.get("force_sign_convention") == "positive_dcoenergy_dx"
        and identity.get("result_force_sign_convention")
        == identity.get("force_sign_convention")
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and math.isclose(
            reported_central_force, central_force, rel_tol=1.0e-12, abs_tol=tolerance
        )
        and math.isclose(result_force, central_force, rel_tol=1.0e-12, abs_tol=tolerance)
        and bool(str(identity.get("coordinate_frame") or ""))
        and identity.get("result_coordinate_frame") == identity.get("coordinate_frame")
        and isinstance(mesh_digests, list)
        and len(mesh_digests) == len(displacement)
        and all(_is_sha256(str(value or "")) for value in mesh_digests)
        and result_mesh_digests == mesh_digests
        and bool(str(identity.get("force_solution_owner") or ""))
        and identity.get("result_force_solution_owner")
        == identity.get("force_solution_owner")
        and _is_sha256(str(identity.get("force_result_sha256") or ""))
        and identity.get("accepted_force_result_sha256")
        == identity.get("force_result_sha256")
    )


def _acoustic_modal_participation_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        frequencies = [float(value) for value in identity.get("mode_frequency_hz", [])]
        result_frequencies = [
            float(value) for value in identity.get("result_mode_frequency_hz", [])
        ]
        masses = [float(value) for value in identity.get("modal_mass_kg", [])]
        result_masses = [
            float(value) for value in identity.get("result_modal_mass_kg", [])
        ]
        participation = [
            float(value) for value in identity.get("participation_factor", [])
        ]
        result_participation = [
            float(value) for value in identity.get("result_participation_factor", [])
        ]
        effective_mass = [
            float(value) for value in identity.get("effective_modal_mass_kg", [])
        ]
        result_effective_mass = [
            float(value) for value in identity.get("result_effective_modal_mass_kg", [])
        ]
        damping = [float(value) for value in identity.get("damping_ratio", [])]
        result_damping = [
            float(value) for value in identity.get("result_damping_ratio", [])
        ]
        probe_factors = [
            float(value) for value in identity.get("probe_mode_factor", [])
        ]
        result_probe_factors = [
            float(value) for value in identity.get("result_probe_mode_factor", [])
        ]
        response_frequencies = [
            float(value) for value in identity.get("response_frequency_hz", [])
        ]
        result_response_frequencies = [
            float(value) for value in identity.get("result_response_frequency_hz", [])
        ]
        response = [
            [float(component) for component in value]
            for value in identity.get("probe_response_complex", [])
        ]
        result_response = [
            [float(component) for component in value]
            for value in identity.get("result_probe_response_complex", [])
        ]
        tolerance = float(identity.get("response_tolerance"))
        result_tolerance = float(identity.get("result_response_tolerance"))
    except (TypeError, ValueError):
        return False
    count = len(frequencies)
    if count < 1 or not all(
        len(values) == count
        for values in (masses, participation, effective_mass, damping, probe_factors)
    ):
        return False
    if len(response) != len(response_frequencies) or not all(
        len(value) == 2 for value in response
    ):
        return False
    reconstructed: list[complex] = []
    for frequency in response_frequencies:
        omega = 2.0 * math.pi * frequency
        value = 0.0j
        for mode_hz, zeta, factor, probe in zip(
            frequencies, damping, participation, probe_factors, strict=True
        ):
            omega_mode = 2.0 * math.pi * mode_hz
            denominator = complex(
                omega_mode**2 - omega**2,
                2.0 * zeta * omega_mode * omega,
            )
            value += probe * factor / denominator
        reconstructed.append(value)
    generation = str(identity.get("modal_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "normalization_generation",
                "mass_generation",
                "participation_generation",
                "damping_generation",
                "frequency_generation",
                "reconstruction_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and bool(str(identity.get("normalization") or ""))
        and identity.get("result_normalization") == identity.get("normalization")
        and _increasing(frequencies)
        and all(value > 0.0 and math.isfinite(value) for value in frequencies + masses)
        and all(math.isfinite(value) for value in participation + probe_factors)
        and all(0.0 <= value < 1.0 for value in damping)
        and all(
            math.isclose(actual, factor**2 * mass, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for actual, factor, mass in zip(
                effective_mass, participation, masses, strict=True
            )
        )
        and result_frequencies == frequencies
        and result_masses == masses
        and result_participation == participation
        and result_effective_mass == effective_mass
        and result_damping == damping
        and result_probe_factors == probe_factors
        and _increasing(response_frequencies)
        and all(value > 0.0 for value in response_frequencies)
        and result_response_frequencies == response_frequencies
        and all(math.isfinite(component) for value in response for component in value)
        and result_response == response
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and all(
            math.isclose(actual[0], expected.real, rel_tol=1.0e-12, abs_tol=tolerance)
            and math.isclose(actual[1], expected.imag, rel_tol=1.0e-12, abs_tol=tolerance)
            for actual, expected in zip(response, reconstructed, strict=True)
        )
        and _is_sha256(str(identity.get("modal_mesh_sha256") or ""))
        and identity.get("result_modal_mesh_sha256")
        == identity.get("modal_mesh_sha256")
        and _is_sha256(str(identity.get("modal_result_sha256") or ""))
        and identity.get("accepted_modal_result_sha256")
        == identity.get("modal_result_sha256")
    )


def _capacitance_matrix_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        terminals = [str(value) for value in identity["terminal_names"]]
        result_terminals = [str(value) for value in identity["result_terminal_names"]]
        matrix = [[float(value) for value in row] for row in identity["capacitance_matrix_f"]]
        result_matrix = [
            [float(value) for value in row]
            for row in identity["result_capacitance_matrix_f"]
        ]
        potentials = [float(value) for value in identity["terminal_potential_v"]]
        result_potentials = [
            float(value) for value in identity["result_terminal_potential_v"]
        ]
        charges = [float(value) for value in identity["terminal_charge_c"]]
        result_charges = [float(value) for value in identity["result_terminal_charge_c"]]
        energy = float(identity["stored_energy_j"])
        result_energy = float(identity["result_stored_energy_j"])
        tolerance = float(identity["reciprocity_tolerance"])
        result_tolerance = float(identity["result_reciprocity_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    size = len(matrix)
    matrix_scale = max((abs(value) for row in matrix for value in row), default=0.0)
    absolute_closure_tolerance = max(matrix_scale * tolerance, 1.0e-24)
    expected_charges = [
        sum(matrix[row][column] * potentials[column] for column in range(size))
        for row in range(size)
    ] if size and len(potentials) == size else []
    expected_energy = 0.5 * sum(
        potential * charge for potential, charge in zip(potentials, charges, strict=True)
    ) if len(potentials) == len(charges) else math.nan
    generation = str(identity.get("capacitance_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "matrix_generation",
                "charge_generation",
                "energy_generation",
                "gauge_generation",
                "reciprocity_generation",
                "terminal_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and size >= 2
        and len(terminals) == size
        and len(set(terminals)) == size
        and result_terminals == terminals
        and identity.get("reference_terminal") in terminals
        and identity.get("result_reference_terminal") == identity.get("reference_terminal")
        and all(len(row) == size for row in matrix)
        and all(math.isfinite(value) for row in matrix for value in row)
        and result_matrix == matrix
        and math.isfinite(tolerance)
        and 0.0 < tolerance <= 1.0e-6
        and result_tolerance == tolerance
        and all(
            math.isclose(
                matrix[row][column],
                matrix[column][row],
                rel_tol=0.0,
                abs_tol=absolute_closure_tolerance,
            )
            for row in range(size)
            for column in range(size)
        )
        and all(abs(sum(row)) <= absolute_closure_tolerance for row in matrix)
        and len(potentials) == len(charges) == size
        and all(math.isfinite(value) for value in potentials + charges)
        and result_potentials == potentials
        and result_charges == charges
        and all(
            math.isclose(
                actual,
                expected,
                rel_tol=1.0e-12,
                abs_tol=absolute_closure_tolerance,
            )
            for actual, expected in zip(charges, expected_charges, strict=True)
        )
        and math.isfinite(energy)
        and energy >= 0.0
        and math.isclose(
            energy,
            expected_energy,
            rel_tol=1.0e-12,
            abs_tol=absolute_closure_tolerance,
        )
        and result_energy == energy
        and bool(str(identity.get("terminal_owner") or ""))
        and identity.get("accepted_terminal_owner") == identity.get("terminal_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _thermoelastic_harmonic_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    try:
        frequency = float(identity["frequency_hz"])
        result_frequency = float(identity["result_frequency_hz"])
        heat = [float(value) for value in identity["heat_source_complex_w"]]
        result_heat = [float(value) for value in identity["result_heat_source_complex_w"]]
        phase = float(identity["heat_source_phase_rad"])
        result_phase = float(identity["result_heat_source_phase_rad"])
        temperature = [float(value) for value in identity["temperature_complex_k"]]
        result_temperature = [float(value) for value in identity["result_temperature_complex_k"]]
        displacement = [float(value) for value in identity["displacement_complex_m"]]
        result_displacement = [float(value) for value in identity["result_displacement_complex_m"]]
        work = float(identity["thermal_expansion_work_j"])
        result_work = float(identity["result_thermal_expansion_work_j"])
        loss = float(identity["mechanical_loss_j"])
        result_loss = float(identity["result_mechanical_loss_j"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(identity.get("thermoelastic_generation") or "")
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "heat_generation",
                "phase_generation",
                "temperature_generation",
                "displacement_generation",
                "work_generation",
                "loss_generation",
                "frequency_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and math.isfinite(frequency)
        and frequency > 0.0
        and result_frequency == frequency
        and len(heat) == len(temperature) == len(displacement) == 2
        and all(math.isfinite(value) for value in heat + temperature + displacement)
        and result_heat == heat
        and result_temperature == temperature
        and result_displacement == displacement
        and math.hypot(*heat) > 0.0
        and math.hypot(*temperature) > 0.0
        and math.hypot(*displacement) > 0.0
        and math.isclose(phase, math.atan2(heat[1], heat[0]), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_phase == phase
        and math.isfinite(work)
        and work > 0.0
        and result_work == work
        and math.isfinite(loss)
        and 0.0 <= loss <= work
        and result_loss == loss
        and identity.get("loss_convention") == "positive_dissipated_per_cycle"
        and identity.get("result_loss_convention") == identity.get("loss_convention")
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _restart_energy_offsets_ok(row: dict[str, Any], sample_count: int) -> bool:
    if "restart_boundaries" not in row:
        return True
    boundaries = row.get("restart_boundaries")
    if not isinstance(boundaries, list) or not boundaries:
        return False
    for boundary in boundaries:
        if not isinstance(boundary, dict):
            return False
        left = boundary.get("left_index")
        right = boundary.get("right_index")
        if (
            not isinstance(left, int)
            or not isinstance(right, int)
            or left < 0
            or right != left + 1
            or right >= sample_count
        ):
            return False
        if not all(
            isinstance(boundary.get(name), str) and boundary.get(name)
            for name in ("generation_before", "generation_after")
        ):
            return False
        try:
            stored_before = float(boundary["stored_energy_before_j"])
            stored_after = float(boundary["stored_energy_after_j"])
            accumulated_before = float(boundary["accumulated_joule_before_j"])
            accumulated_after = float(boundary["accumulated_joule_offset_after_j"])
        except (KeyError, TypeError, ValueError):
            return False
        if not all(
            math.isfinite(value)
            for value in (
                stored_before,
                stored_after,
                accumulated_before,
                accumulated_after,
            )
        ):
            return False
        if abs(accumulated_after - accumulated_before) > 1.0e-12 * max(
            abs(accumulated_before), 1.0
        ):
            return False
        if abs(stored_after - stored_before) > 0.1 * max(abs(stored_before), 1.0e-15):
            return False
    return True


def _thermoviscous_pressure_interface_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("interface_generation") or "")
    try:
        frequency = float(identity["frequency_hz"])
        result_frequency = float(identity["result_frequency_hz"])
        area = float(identity["interface_area_m2"])
        result_area = float(identity["result_interface_area_m2"])
        velocity = [float(item) for item in identity["normal_velocity_complex_m_per_s"]]
        result_velocity = [
            float(item) for item in identity["result_normal_velocity_complex_m_per_s"]
        ]
        pressure = [float(item) for item in identity["pressure_complex_pa"]]
        result_pressure = [float(item) for item in identity["result_pressure_complex_pa"]]
        interface_power = float(identity["interface_power_w"])
        result_interface_power = float(identity["result_interface_power_w"])
        viscous_loss = float(identity["viscous_loss_w"])
        result_viscous_loss = float(identity["result_viscous_loss_w"])
        thermal_loss = float(identity["thermal_loss_w"])
        result_thermal_loss = float(identity["result_thermal_loss_w"])
        outgoing_power = float(identity["outgoing_acoustic_power_w"])
        result_outgoing_power = float(identity["result_outgoing_acoustic_power_w"])
    except (KeyError, TypeError, ValueError):
        return False
    if len(velocity) != 2 or len(pressure) != 2:
        return False
    recomputed_power = 0.5 * area * (
        pressure[0] * velocity[0] + pressure[1] * velocity[1]
    )
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "velocity_generation",
                "traction_generation",
                "viscous_generation",
                "thermal_generation",
                "power_generation",
                "normal_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and all(
            math.isfinite(item)
            for item in (
                frequency,
                area,
                *velocity,
                *pressure,
                interface_power,
                viscous_loss,
                thermal_loss,
                outgoing_power,
            )
        )
        and frequency > 0.0
        and area > 0.0
        and result_frequency == frequency
        and result_area == area
        and result_velocity == velocity
        and result_pressure == pressure
        and identity.get("traction_sign") == "minus_pressure_times_outward_normal"
        and identity.get("result_traction_sign") == identity.get("traction_sign")
        and identity.get("normal_orientation")
        == "thermoviscous_to_pressure_acoustics"
        and identity.get("result_normal_orientation")
        == identity.get("normal_orientation")
        and interface_power > 0.0
        and viscous_loss >= 0.0
        and thermal_loss >= 0.0
        and outgoing_power >= 0.0
        and math.isclose(
            interface_power, recomputed_power, rel_tol=1.0e-12, abs_tol=1.0e-15
        )
        and math.isclose(
            interface_power,
            viscous_loss + thermal_loss + outgoing_power,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and result_interface_power == interface_power
        and result_viscous_loss == viscous_loss
        and result_thermal_loss == thermal_loss
        and result_outgoing_power == outgoing_power
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _piezoelectric_reciprocity_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("piezo_generation") or "")
    names = (
        "direct_coefficient_c_per_n",
        "converse_coefficient_m_per_v",
        "electric_field_v_per_m",
        "mechanical_stress_pa",
        "induced_strain",
        "induced_charge_density_c_per_m2",
        "terminal_charge_c",
        "electrical_work_j",
        "elastic_energy_j",
        "coupling_energy_j",
        "total_stored_energy_j",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        result_values = {
            name: float(identity[f"result_{name}"]) for name in names
        }
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            identity.get(key) == generation
            for key in (
                "charge_generation",
                "strain_generation",
                "reciprocity_generation",
                "electrical_energy_generation",
                "elastic_energy_generation",
                "coupling_energy_generation",
                "polarization_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and all(math.isfinite(item) for item in values.values())
        and values["direct_coefficient_c_per_n"] > 0.0
        and math.isclose(
            values["direct_coefficient_c_per_n"],
            values["converse_coefficient_m_per_v"],
            rel_tol=1.0e-12,
            abs_tol=0.0,
        )
        and math.isclose(
            values["induced_strain"],
            values["converse_coefficient_m_per_v"]
            * values["electric_field_v_per_m"],
            rel_tol=1.0e-12,
            abs_tol=1.0e-18,
        )
        and math.isclose(
            values["induced_charge_density_c_per_m2"],
            values["direct_coefficient_c_per_n"]
            * values["mechanical_stress_pa"],
            rel_tol=1.0e-12,
            abs_tol=1.0e-18,
        )
        and values["terminal_charge_c"] > 0.0
        and values["electrical_work_j"] >= 0.0
        and values["elastic_energy_j"] >= 0.0
        and values["coupling_energy_j"] >= 0.0
        and math.isclose(
            values["total_stored_energy_j"],
            values["electrical_work_j"]
            + values["elastic_energy_j"]
            - values["coupling_energy_j"],
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        and result_values == values
        and identity.get("polarization_frame") == "material_axis_3"
        and identity.get("result_polarization_frame")
        == identity.get("polarization_frame")
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _poroelastic_biot_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("poroelastic_generation") or "")
    names = (
        "biot_coefficient", "pore_pressure_pa", "volumetric_strain",
        "biot_modulus_pa", "fluid_content_increment", "permeability_m2",
        "dynamic_viscosity_pa_s", "pressure_gradient_pa_per_m",
        "darcy_flux_m_per_s", "domain_volume_m3", "time_step_s",
        "interface_traction_pa", "storage_energy_j",
        "skeleton_coupling_work_j", "fluid_dissipation_j",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    alpha = values["biot_coefficient"]
    pressure = values["pore_pressure_pa"]
    strain = values["volumetric_strain"]
    modulus = values["biot_modulus_pa"]
    permeability = values["permeability_m2"]
    viscosity = values["dynamic_viscosity_pa_s"]
    gradient = values["pressure_gradient_pa_per_m"]
    flux = values["darcy_flux_m_per_s"]
    volume = values["domain_volume_m3"]
    timestep = values["time_step_s"]
    expected_content = alpha * strain + pressure / modulus
    expected_flux = -permeability * gradient / viscosity
    expected_traction = -alpha * pressure
    expected_storage = 0.5 * pressure * pressure / modulus * volume
    expected_coupling = alpha * pressure * strain * volume
    expected_dissipation = viscosity / permeability * flux * flux * volume * timestep
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "biot_generation", "pressure_generation", "displacement_generation",
            "flux_generation", "storage_generation", "dissipation_generation",
            "interface_generation", "mesh_generation", "result_generation",
        ))
        and all(math.isfinite(item) for item in values.values())
        and 0.0 < alpha <= 1.0
        and pressure > 0.0 and strain >= 0.0 and modulus > 0.0
        and permeability > 0.0 and viscosity > 0.0
        and volume > 0.0 and timestep > 0.0
        and math.isclose(values["fluid_content_increment"], expected_content, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(flux, expected_flux, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["interface_traction_pa"], expected_traction, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["storage_energy_j"], expected_storage, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["skeleton_coupling_work_j"], expected_coupling, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["fluid_dissipation_j"], expected_dissipation, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["storage_energy_j"] >= 0.0
        and values["skeleton_coupling_work_j"] >= 0.0
        and values["fluid_dissipation_j"] >= 0.0
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and identity.get("interface_normal") == "porous_skeleton_to_free_fluid"
        and identity.get("result_interface_normal") == identity.get("interface_normal")
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _rotating_induction_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("induction_generation") or "")
    names = (
        "supply_frequency_hz", "pole_pairs", "synchronous_speed_rad_per_s",
        "rotor_speed_rad_per_s", "slip", "rotor_electrical_frequency_hz",
        "rotor_phase_current_a_rms", "rotor_phase_resistance_ohm",
        "rotor_copper_loss_w", "airgap_torque_nm", "airgap_power_w",
        "mechanical_power_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    frequency = values["supply_frequency_hz"]
    pole_pairs = values["pole_pairs"]
    synchronous = values["synchronous_speed_rad_per_s"]
    rotor_speed = values["rotor_speed_rad_per_s"]
    slip = values["slip"]
    current = values["rotor_phase_current_a_rms"]
    resistance = values["rotor_phase_resistance_ohm"]
    torque = values["airgap_torque_nm"]
    expected_synchronous = 2.0 * math.pi * frequency / pole_pairs
    expected_slip = (synchronous - rotor_speed) / synchronous
    expected_copper_loss = 3.0 * current * current * resistance
    expected_airgap_power = torque * synchronous
    expected_mechanical_power = torque * rotor_speed
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "slip_generation", "frequency_generation", "current_generation",
            "loss_generation", "torque_generation", "power_generation",
            "frame_generation", "mesh_generation", "result_generation",
        ))
        and all(math.isfinite(item) for item in values.values())
        and frequency > 0.0
        and pole_pairs > 0.0 and pole_pairs == math.floor(pole_pairs)
        and math.isclose(synchronous, expected_synchronous, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and 0.0 <= rotor_speed < synchronous
        and math.isclose(slip, expected_slip, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and 0.0 < slip < 1.0
        and math.isclose(values["rotor_electrical_frequency_hz"], slip * frequency, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and current > 0.0 and resistance > 0.0 and torque > 0.0
        and math.isclose(values["rotor_copper_loss_w"], expected_copper_loss, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["airgap_power_w"], expected_airgap_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["mechanical_power_w"], expected_mechanical_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["airgap_power_w"], values["mechanical_power_w"] + values["rotor_copper_loss_w"], rel_tol=1.0e-12, abs_tol=1.0e-10)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-12) for name in names)
        and identity.get("rotating_frame") == "rotor_mechanical_frame"
        and identity.get("result_rotating_frame") == identity.get("rotating_frame")
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _thermoacoustic_meanflow_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "thermoacoustic_meanflow_convected_wavenumber_flux_impedance_power_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("thermoacoustic_generation") or "")
    names = (
        "frequency_hz", "sound_speed_m_per_s", "mean_flow_mach",
        "mean_flow_speed_m_per_s", "downstream_wavenumber_rad_per_m",
        "upstream_wavenumber_rad_per_m", "density_kg_per_m3",
        "pressure_rms_pa", "particle_velocity_rms_m_per_s",
        "acoustic_intensity_w_per_m2", "boundary_area_m2",
        "boundary_impedance_pa_s_per_m", "boundary_flux_power_w",
        "impedance_work_w", "dissipated_power_w", "power_balance_residual_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    frequency = values["frequency_hz"]
    sound_speed = values["sound_speed_m_per_s"]
    mach = values["mean_flow_mach"]
    density = values["density_kg_per_m3"]
    pressure = values["pressure_rms_pa"]
    area = values["boundary_area_m2"]
    expected_speed = mach * sound_speed
    expected_particle_velocity = pressure / (density * sound_speed)
    expected_intensity = pressure * expected_particle_velocity
    expected_power = expected_intensity * area
    expected_residual = values["boundary_flux_power_w"] - values["dissipated_power_w"]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "meanflow_generation", "wavenumber_generation", "flux_generation",
            "impedance_generation", "power_generation", "mesh_generation",
            "result_generation",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(frequency, sound_speed, density, pressure, area) > 0.0
        and 0.0 <= mach < 1.0
        and math.isclose(values["mean_flow_speed_m_per_s"], expected_speed, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["downstream_wavenumber_rad_per_m"], 2.0 * math.pi * frequency / (sound_speed + expected_speed), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["upstream_wavenumber_rad_per_m"], 2.0 * math.pi * frequency / (sound_speed - expected_speed), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["particle_velocity_rms_m_per_s"], expected_particle_velocity, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["acoustic_intensity_w_per_m2"], expected_intensity, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["boundary_impedance_pa_s_per_m"], density * sound_speed, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["boundary_flux_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["impedance_work_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["dissipated_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["power_balance_residual_w"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and abs(expected_residual) <= 1.0e-12 * max(expected_power, 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _battery_electrothermal_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "battery_electrothermal_soc_current_heat_temperature_energy_safety_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("battery_generation") or "")
    names = (
        "capacity_c", "initial_state_of_charge", "terminal_current_a",
        "terminal_voltage_v", "time_step_s", "final_state_of_charge",
        "internal_resistance_ohm", "irreversible_heat_j", "reversible_heat_j",
        "thermal_energy_j", "electrical_energy_j", "cell_mass_kg",
        "specific_heat_j_per_kg_k", "initial_temperature_k",
        "final_temperature_k", "maximum_safe_temperature_k",
        "thermal_balance_residual_j",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    capacity = values["capacity_c"]
    initial_soc = values["initial_state_of_charge"]
    current = values["terminal_current_a"]
    voltage = values["terminal_voltage_v"]
    timestep = values["time_step_s"]
    resistance = values["internal_resistance_ohm"]
    mass = values["cell_mass_kg"]
    heat_capacity = values["specific_heat_j_per_kg_k"]
    expected_final_soc = initial_soc - current * timestep / capacity
    expected_irreversible_heat = current * current * resistance * timestep
    expected_thermal_energy = values["irreversible_heat_j"] + values["reversible_heat_j"]
    expected_temperature = values["initial_temperature_k"] + expected_thermal_energy / (mass * heat_capacity)
    expected_residual = values["thermal_energy_j"] - expected_thermal_energy
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "soc_generation", "current_generation", "heat_generation",
            "temperature_generation", "energy_generation", "safety_generation",
            "mesh_generation", "result_generation",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(capacity, current, voltage, timestep, resistance, mass, heat_capacity) > 0.0
        and 0.0 <= initial_soc <= 1.0
        and math.isclose(values["final_state_of_charge"], expected_final_soc, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and 0.0 <= values["final_state_of_charge"] <= 1.0
        and math.isclose(values["irreversible_heat_j"], expected_irreversible_heat, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["reversible_heat_j"] >= 0.0
        and math.isclose(values["thermal_energy_j"], expected_thermal_energy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["electrical_energy_j"], voltage * current * timestep, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["final_temperature_k"], expected_temperature, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["initial_temperature_k"] > 0.0
        and values["initial_temperature_k"] <= values["final_temperature_k"] <= values["maximum_safe_temperature_k"]
        and math.isclose(values["thermal_balance_residual_j"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and abs(expected_residual) <= 1.0e-12 * max(expected_thermal_energy, 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _piezoelectric_admittance_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "piezoelectric_admittance_resonance_antiresonance_coupling_energy_phase_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("piezoelectric_generation") or "")
    names = (
        "resonance_frequency_hz", "antiresonance_frequency_hz",
        "electromechanical_coupling_squared", "voltage_rms_v", "current_rms_a",
        "admittance_magnitude_s", "admittance_phase_deg",
        "real_electrical_power_w", "reactive_electrical_power_var",
        "mechanical_output_power_w", "dielectric_loss_w",
        "mechanical_stored_energy_j", "electric_stored_energy_j",
        "power_balance_residual_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    resonance = values["resonance_frequency_hz"]
    antiresonance = values["antiresonance_frequency_hz"]
    voltage = values["voltage_rms_v"]
    current = values["current_rms_a"]
    phase = math.radians(values["admittance_phase_deg"])
    apparent_power = voltage * current
    expected_coupling = 1.0 - (resonance / antiresonance) ** 2 if antiresonance > 0.0 else math.nan
    expected_residual = (
        values["real_electrical_power_w"]
        - values["mechanical_output_power_w"]
        - values["dielectric_loss_w"]
    )
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "admittance_generation", "resonance_generation", "coupling_generation",
            "phase_generation", "energy_generation", "power_generation",
            "mesh_generation", "result_generation",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(resonance, antiresonance, voltage, current) > 0.0
        and antiresonance > resonance
        and math.isclose(values["electromechanical_coupling_squared"], expected_coupling, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and 0.0 < values["electromechanical_coupling_squared"] < 1.0
        and math.isclose(values["admittance_magnitude_s"], current / voltage, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and -180.0 <= values["admittance_phase_deg"] <= 180.0
        and math.isclose(values["real_electrical_power_w"], apparent_power * math.cos(phase), rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["reactive_electrical_power_var"], apparent_power * math.sin(phase), rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["mechanical_output_power_w"] >= 0.0
        and values["dielectric_loss_w"] >= 0.0
        and values["mechanical_stored_energy_j"] >= 0.0
        and values["electric_stored_energy_j"] >= 0.0
        and math.isclose(values["power_balance_residual_w"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and abs(expected_residual) <= 1.0e-12 * max(values["real_electrical_power_w"], 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _fluidfilm_bearing_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "fluidfilm_bearing_reynolds_pressure_load_friction_temperature_power_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("fluidfilm_generation") or "")
    names = (
        "journal_radius_m", "bearing_length_m", "radial_clearance_m",
        "eccentricity_ratio", "minimum_film_thickness_m", "maximum_pressure_pa",
        "integrated_load_n", "attitude_angle_deg", "angular_speed_rad_per_s",
        "friction_torque_nm", "shaft_power_w", "viscous_dissipation_w",
        "removed_heat_w", "inlet_temperature_k", "maximum_temperature_k",
        "power_balance_residual_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    expected_thickness = values["radial_clearance_m"] * (1.0 - values["eccentricity_ratio"])
    expected_power = values["friction_torque_nm"] * values["angular_speed_rad_per_s"]
    expected_residual = expected_power - values["viscous_dissipation_w"] - values["removed_heat_w"]
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "film_generation", "pressure_generation", "load_generation",
            "friction_generation", "temperature_generation", "power_generation",
            "mesh_generation", "result_generation",
        ))
        and all(math.isfinite(item) for item in values.values())
        and min(
            values["journal_radius_m"], values["bearing_length_m"],
            values["radial_clearance_m"], values["maximum_pressure_pa"],
            values["integrated_load_n"], values["angular_speed_rad_per_s"],
            values["friction_torque_nm"], values["inlet_temperature_k"],
        ) > 0.0
        and 0.0 <= values["eccentricity_ratio"] < 1.0
        and math.isclose(values["minimum_film_thickness_m"], expected_thickness, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and -180.0 <= values["attitude_angle_deg"] <= 180.0
        and math.isclose(values["shaft_power_w"], expected_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["viscous_dissipation_w"] >= 0.0
        and values["removed_heat_w"] >= 0.0
        and values["maximum_temperature_k"] >= values["inlet_temperature_k"]
        and math.isclose(values["power_balance_residual_w"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and abs(expected_residual) <= 1.0e-12 * max(expected_power, 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _induction_heating_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "inductionheating_skin_proximity_joule_thermal_flux_temperature_energy_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("induction_generation") or "")
    names = (
        "frequency_hz", "conductivity_s_per_m", "relative_permeability",
        "skin_depth_m", "joule_loss_w", "magnetic_loss_w",
        "electromagnetic_input_power_w", "outward_thermal_flux_w",
        "ambient_temperature_k", "maximum_temperature_k", "temperature_rise_k",
        "electromagnetic_power_balance_residual_w",
        "thermal_power_balance_residual_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
        surface_current = [
            float(value) for value in identity["surface_current_density_a_per_m"]
        ]
        result_surface_current = [
            float(value)
            for value in identity["result_surface_current_density_a_per_m"]
        ]
        proximity_current = [
            float(value)
            for value in identity["proximity_current_density_a_per_m2"]
        ]
        result_proximity_current = [
            float(value)
            for value in identity["result_proximity_current_density_a_per_m2"]
        ]
    except (KeyError, TypeError, ValueError):
        return False
    omega = 2.0 * math.pi * values["frequency_hz"]
    permeability = 4.0e-7 * math.pi * values["relative_permeability"]
    expected_skin_depth = math.sqrt(
        2.0 / (omega * permeability * values["conductivity_s_per_m"])
    )
    expected_em_residual = (
        values["electromagnetic_input_power_w"]
        - values["joule_loss_w"]
        - values["magnetic_loss_w"]
    )
    expected_thermal_residual = (
        values["joule_loss_w"]
        + values["magnetic_loss_w"]
        - values["outward_thermal_flux_w"]
    )
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "skin_generation", "proximity_generation", "joule_generation",
            "thermal_generation", "temperature_generation", "energy_generation",
            "mesh_generation", "result_generation",
        ))
        and all(math.isfinite(value) for value in values.values())
        and min(
            values["frequency_hz"], values["conductivity_s_per_m"],
            values["relative_permeability"], values["skin_depth_m"],
            values["electromagnetic_input_power_w"],
            values["ambient_temperature_k"], values["maximum_temperature_k"],
        ) > 0.0
        and values["joule_loss_w"] >= 0.0
        and values["magnetic_loss_w"] >= 0.0
        and bool(surface_current) and bool(proximity_current)
        and all(math.isfinite(value) and value >= 0.0 for value in surface_current)
        and all(math.isfinite(value) and value >= 0.0 for value in proximity_current)
        and result_surface_current == surface_current
        and result_proximity_current == proximity_current
        and math.isclose(values["skin_depth_m"], expected_skin_depth, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["temperature_rise_k"], values["maximum_temperature_k"] - values["ambient_temperature_k"], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["temperature_rise_k"] >= 0.0
        and math.isclose(values["electromagnetic_power_balance_residual_w"], expected_em_residual, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["thermal_power_balance_residual_w"], expected_thermal_residual, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and abs(expected_em_residual) <= 1.0e-12 * max(values["electromagnetic_input_power_w"], 1.0e-15)
        and abs(expected_thermal_residual) <= 1.0e-12 * max(values["outward_thermal_flux_w"], 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _reacting_species_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "species_transport_reaction_diffusion_flux_massbalance_rate_temperature_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("species_generation") or "")
    names = (
        "diffusivity_m2_per_s", "temperature_k", "gas_constant_j_per_mol_k",
        "preexponential_factor_per_s", "activation_energy_j_per_mol",
        "reaction_rate_constant_per_s", "mean_concentration_mol_per_m3",
        "domain_volume_m3", "integrated_species_mol",
        "integrated_consumption_mol_per_s", "inward_boundary_flux_mol_per_s",
        "mass_balance_residual_mol_per_s",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    expected_rate = values["preexponential_factor_per_s"] * math.exp(
        -values["activation_energy_j_per_mol"]
        / (values["gas_constant_j_per_mol_k"] * values["temperature_k"])
    )
    expected_amount = (
        values["mean_concentration_mol_per_m3"] * values["domain_volume_m3"]
    )
    expected_consumption = expected_rate * expected_amount
    expected_residual = (
        values["inward_boundary_flux_mol_per_s"]
        - values["integrated_consumption_mol_per_s"]
    )
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "diffusion_generation", "reaction_generation", "flux_generation",
            "mass_generation", "temperature_generation", "mesh_generation",
            "result_generation",
        ))
        and all(math.isfinite(value) for value in values.values())
        and min(
            values["diffusivity_m2_per_s"], values["temperature_k"],
            values["gas_constant_j_per_mol_k"],
            values["preexponential_factor_per_s"],
            values["mean_concentration_mol_per_m3"], values["domain_volume_m3"],
        ) > 0.0
        and values["activation_energy_j_per_mol"] >= 0.0
        and math.isclose(values["reaction_rate_constant_per_s"], expected_rate, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["integrated_species_mol"], expected_amount, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["integrated_consumption_mol_per_s"], expected_consumption, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and values["inward_boundary_flux_mol_per_s"] >= 0.0
        and math.isclose(values["mass_balance_residual_mol_per_s"], expected_residual, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and abs(expected_residual) <= 1.0e-12 * max(expected_consumption, 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _microwave_heating_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "microwaveheating_sparameter_absorbedpower_jouleheat_temperature_energy_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("microwave_generation") or "")
    names = (
        "frequency_hz", "reference_impedance_ohm", "s11_magnitude",
        "s21_magnitude", "incident_power_w", "reflected_power_w",
        "transmitted_power_w", "absorbed_power_w", "joule_heat_w",
        "dielectric_heat_w", "electromagnetic_power_residual_w",
        "outward_thermal_flux_w", "ambient_temperature_k",
        "maximum_temperature_k", "temperature_rise_k",
        "thermal_power_residual_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    incident = values["incident_power_w"]
    expected_reflected = incident * values["s11_magnitude"] ** 2
    expected_transmitted = incident * values["s21_magnitude"] ** 2
    expected_absorbed = incident - expected_reflected - expected_transmitted
    expected_em_residual = (
        expected_absorbed - values["joule_heat_w"] - values["dielectric_heat_w"]
    )
    expected_thermal_residual = (
        values["joule_heat_w"]
        + values["dielectric_heat_w"]
        - values["outward_thermal_flux_w"]
    )
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "sparameter_generation", "power_generation", "heat_generation",
            "temperature_generation", "energy_generation", "mesh_generation",
            "result_generation",
        ))
        and all(math.isfinite(value) for value in values.values())
        and min(
            values["frequency_hz"], values["reference_impedance_ohm"],
            incident, values["ambient_temperature_k"],
            values["maximum_temperature_k"],
        ) > 0.0
        and 0.0 <= values["s11_magnitude"] <= 1.0
        and 0.0 <= values["s21_magnitude"] <= 1.0
        and values["joule_heat_w"] >= 0.0
        and values["dielectric_heat_w"] >= 0.0
        and expected_absorbed >= 0.0
        and math.isclose(values["reflected_power_w"], expected_reflected, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["transmitted_power_w"], expected_transmitted, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["absorbed_power_w"], expected_absorbed, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["electromagnetic_power_residual_w"], expected_em_residual, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(values["thermal_power_residual_w"], expected_thermal_residual, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and abs(expected_em_residual) <= 1.0e-12 * max(incident, 1.0e-15)
        and abs(expected_thermal_residual) <= 1.0e-12 * max(expected_absorbed, 1.0e-15)
        and math.isclose(values["temperature_rise_k"], values["maximum_temperature_k"] - values["ambient_temperature_k"], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and values["temperature_rise_k"] >= 0.0
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def _poroelastic_wave_identity_ok(summary: dict[str, Any]) -> bool:
    identity = summary.get(
        "poroelastic_wave_pressure_displacement_flux_dissipation_mass_energy_mesh_result_generation_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, dict):
        return False
    generation = str(identity.get("poroelastic_generation") or "")
    names = (
        "frequency_hz", "porosity", "solid_displacement_amplitude_m",
        "pore_pressure_amplitude_pa", "darcy_flux_amplitude_m_per_s",
        "pressure_displacement_phase_deg", "fluid_mass_kg",
        "fluid_mass_rate_kg_per_s", "net_inward_mass_flux_kg_per_s",
        "solid_energy_j", "fluid_energy_j", "dissipated_power_w",
        "input_power_w", "mass_balance_residual_kg_per_s",
        "energy_balance_residual_w",
    )
    try:
        values = {name: float(identity[name]) for name in names}
        results = {name: float(identity[f"result_{name}"]) for name in names}
    except (KeyError, TypeError, ValueError):
        return False
    expected_mass_residual = (
        values["net_inward_mass_flux_kg_per_s"]
        - values["fluid_mass_rate_kg_per_s"]
    )
    expected_energy_residual = (
        values["input_power_w"] - values["dissipated_power_w"]
    )
    return (
        bool(generation)
        and all(identity.get(key) == generation for key in (
            "pressure_generation", "displacement_generation", "flux_generation",
            "mass_generation", "energy_generation", "mesh_generation",
            "result_generation",
        ))
        and all(math.isfinite(value) for value in values.values())
        and values["frequency_hz"] > 0.0
        and 0.0 < values["porosity"] < 1.0
        and min(
            values["solid_displacement_amplitude_m"],
            values["pore_pressure_amplitude_pa"],
            values["darcy_flux_amplitude_m_per_s"], values["fluid_mass_kg"],
        ) >= 0.0
        and -180.0 <= values["pressure_displacement_phase_deg"] <= 180.0
        and values["solid_energy_j"] >= 0.0
        and values["fluid_energy_j"] >= 0.0
        and values["dissipated_power_w"] >= 0.0
        and values["input_power_w"] >= 0.0
        and math.isclose(values["mass_balance_residual_kg_per_s"], expected_mass_residual, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(values["energy_balance_residual_w"], expected_energy_residual, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and abs(expected_mass_residual) <= 1.0e-12 * max(abs(values["fluid_mass_rate_kg_per_s"]), 1.0e-15)
        and abs(expected_energy_residual) <= 1.0e-12 * max(values["input_power_w"], 1.0e-15)
        and all(math.isclose(results[name], values[name], rel_tol=1.0e-12, abs_tol=1.0e-15) for name in names)
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _is_sha256(str(identity.get("mesh_sha256") or ""))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _is_sha256(str(identity.get("result_sha256") or ""))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def rotational_eddy_brake_energy_gate(
    summary: dict[str, Any],
    *,
    maximum_inertia_relative_error: float = 1.0e-5,
    maximum_angular_impulse_residual: float = 0.01,
    maximum_total_energy_residual: float = 0.01,
    maximum_replay_error_over_span: float = 1.0e-6,
    maximum_field_energy_time_misalignment_s: float = 1.0e-12,
    maximum_field_energy_adjacent_jump_fraction: float = 0.1,
    maximum_field_energy_curvature_outlier_ratio: float = 50.0,
    minimum_decay_fraction: float = 0.5,
) -> dict[str, Any]:
    """Gate free rotational braking with angular momentum and full energy storage.

    The energy identity includes magnetic field storage.  Torque times speed is
    reported only as a diagnostic because it is not generally equal to Joule
    loss when the field-energy rate is unavailable.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract")
    units = summary.get("units")
    disc = summary.get("disc")
    replays = summary.get("replays")
    energy_row = summary.get("energy_replay")
    timing = summary.get("timing_breakdown_s")
    if not all(isinstance(value, dict) for value in (contract, units, disc, energy_row)):
        raise ValueError("contract, units, disc, and energy_replay must be mappings")
    if not isinstance(replays, list) or len(replays) < 2:
        raise ValueError("at least two fresh replays are required")

    density = float(disc.get("density_kg_m3", math.nan))
    radius = float(disc.get("radius_m", math.nan))
    thickness = float(disc.get("thickness_m", math.nan))
    reported_inertia = float(summary.get("reported_inertia_kg_m2", math.nan))
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (density, radius, thickness, reported_inertia)
    ):
        raise ValueError("disc dimensions, density, and inertia must be positive")
    analytic_inertia = 0.5 * density * math.pi * radius**4 * thickness
    inertia_error = abs(reported_inertia - analytic_inertia) / analytic_inertia

    parsed = [_parse_replay(row) for row in replays]
    artifact_generations_ok = _artifact_generations_ok([*replays, energy_row])
    artifact_coordinate_frames_ok = _artifact_coordinate_frames_ok(
        [*replays, energy_row]
    )
    convergence_provenance_ok = _convergence_provenance_ok(summary)
    force_selection_identity_ok = _force_selection_identity_ok(summary)
    excitation_basis_identity_ok = _excitation_basis_identity_ok(summary)
    live_stored_force_identity_ok = _live_stored_force_identity_ok(summary)
    loss_partition_identity_ok = _loss_partition_identity_ok(summary)
    material_property_parameter_identity_ok = (
        _material_property_parameter_identity_ok(summary)
    )
    force_selection_topology_identity_ok = _force_selection_topology_identity_ok(
        summary
    )
    weak_form_coordinate_transform_identity_ok = (
        _weak_form_coordinate_transform_identity_ok(summary)
    )
    time_harmonic_phasor_convention_identity_ok = (
        _time_harmonic_phasor_convention_identity_ok(summary)
    )
    eigenmode_mass_inner_product_identity_ok = (
        _eigenmode_mass_inner_product_identity_ok(summary)
    )
    ale_material_derivative_time_level_identity_ok = (
        _ale_material_derivative_time_level_identity_ok(summary)
    )
    harmonic_reference_time_origin_identity_ok = (
        _harmonic_reference_time_origin_identity_ok(summary)
    )
    deformed_domain_integral_jacobian_identity_ok = (
        _deformed_domain_integral_jacobian_identity_ok(summary)
    )
    nonlinear_residual_tangent_iteration_identity_ok = (
        _nonlinear_residual_tangent_iteration_identity_ok(summary)
    )
    moving_mesh_field_transfer_frame_identity_ok = (
        _moving_mesh_field_transfer_frame_identity_ok(summary)
    )
    segregated_block_variable_scaling_identity_ok = (
        _segregated_block_variable_scaling_identity_ok(summary)
    )
    modal_port_power_surface_orientation_identity_ok = (
        _modal_port_power_surface_orientation_identity_ok(summary)
    )
    degenerate_subspace_tracking_identity_ok = (
        _degenerate_eigenmode_subspace_tracking_basis_identity_ok(summary)
    )
    adaptive_bdf_restart_identity_ok = (
        _adaptive_bdf_restart_history_event_identity_ok(summary)
    )
    all_cardinalities = True
    all_times_increase = True
    nonnegative_dissipation = True
    monotone_decay = True
    impulse_errors: list[float] = []
    decay_fractions: list[float] = []
    power_diagnostics: list[float] = []
    for replay in parsed:
        lengths = {len(value) for key, value in replay.items() if isinstance(value, list)}
        all_cardinalities &= len(lengths) == 1 and next(iter(lengths), 0) >= 20
        times = replay["time_s"]
        omega = replay["angular_velocity_rad_s"]
        torque = replay["braking_torque_nm"]
        joule = replay["joule_loss_w"]
        all_times_increase &= _increasing(times)
        nonnegative_dissipation &= min(torque) >= 0.0 and min(joule) >= 0.0
        monotone_decay &= omega[0] > 0.0 and all(
            right <= left + 1.0e-10 * max(abs(left), 1.0)
            for left, right in zip(omega, omega[1:])
        )
        angular_scale = reported_inertia * abs(omega[0] - omega[-1])
        impulse = _cumulative_integral(torque, times)
        residual = [
            reported_inertia * (value - omega[0]) + integrated
            for value, integrated in zip(omega, impulse, strict=True)
        ]
        impulse_errors.append(
            max(abs(value) for value in residual) / angular_scale
            if angular_scale > 0.0
            else math.inf
        )
        decay_fractions.append(1.0 - omega[-1] / omega[0])
        power_scale = max(joule)
        power_diagnostics.append(
            max(abs(t * w - q) for t, w, q in zip(torque, omega, joule, strict=True))
            / power_scale
            if power_scale > 0.0
            else math.inf
        )

    reference = parsed[0]
    replay_time_errors: list[float] = []
    replay_field_errors: list[float] = []
    for replay in parsed[1:]:
        replay_time_errors.append(
            max(
                (
                    abs(a - b)
                    for a, b in zip(
                        reference["time_s"], replay["time_s"], strict=True
                    )
                ),
                default=math.inf,
            )
            if len(reference["time_s"]) == len(replay["time_s"])
            else math.inf
        )
        replay_field_errors.extend(
            _span_error(reference[key], replay[key], reference[key])
            for key in (
                "angular_velocity_rad_s",
                "braking_torque_nm",
                "joule_loss_w",
            )
        )

    energy = _parse_replay(energy_row)
    field_time = _series(energy_row, "field_energy_time_s")
    magnetic_energy = _series(energy_row, "magnetic_energy_j")
    field_energy_scale = max(abs(value) for value in magnetic_energy)
    maximum_field_energy_adjacent_jump_fraction_observed = (
        max(
            abs(right - left)
            for left, right in zip(magnetic_energy, magnetic_energy[1:])
        )
        / field_energy_scale
        if field_energy_scale > 0.0
        else math.inf
    )
    field_energy_curvatures = [
        abs(right - 2.0 * center + left)
        for left, center, right in zip(
            magnetic_energy[:-2],
            magnetic_energy[1:-1],
            magnetic_energy[2:],
            strict=True,
        )
    ]
    field_energy_curvature_baseline = statistics.median(field_energy_curvatures)
    maximum_field_energy_curvature_outlier_ratio_observed = (
        max(field_energy_curvatures)
        / max(field_energy_curvature_baseline, 1.0e-15 * field_energy_scale)
        if field_energy_curvatures and field_energy_scale > 0.0
        else 0.0
    )
    energy_times = energy["time_s"]
    energy_omega = energy["angular_velocity_rad_s"]
    energy_joule = energy["joule_loss_w"]
    energy_cardinality = (
        len({len(value) for key, value in energy.items() if isinstance(value, list)}) == 1
        and len(field_time) == len(magnetic_energy) >= 2
    )
    restart_energy_offsets_ok = _restart_energy_offsets_ok(
        energy_row, len(energy_times)
    )
    maximum_field_energy_time_misalignment_s_observed = (
        max(
            abs(field_sample - primary_sample)
            for field_sample, primary_sample in zip(
                field_time, energy_times, strict=True
            )
        )
        if len(field_time) == len(energy_times)
        else math.inf
    )
    field_time_alignment = (
        _increasing(field_time)
        and maximum_field_energy_time_misalignment_s_observed
        <= float(maximum_field_energy_time_misalignment_s)
    )
    kinetic_drop = 0.5 * reported_inertia * (
        energy_omega[0] ** 2 - energy_omega[-1] ** 2
    )
    magnetic_drop = magnetic_energy[0] - magnetic_energy[-1]
    total_stored_drop = kinetic_drop + magnetic_drop
    joule_energy = _integral(energy_joule, energy_times)
    total_energy_error = (
        abs(total_stored_drop - joule_energy) / abs(total_stored_drop)
        if total_stored_drop > 0.0
        else math.inf
    )
    energy_replay_error = max(
        _span_error(energy[key], reference[key], reference[key])
        for key in (
            "angular_velocity_rad_s",
            "braking_torque_nm",
            "joule_loss_w",
        )
    )

    timing_ok = False
    if isinstance(timing, dict) and len(timing) == 4:
        try:
            timing_ok = all(
                math.isfinite(float(value)) and float(value) >= 0.0
                for value in timing.values()
            )
        except (TypeError, ValueError):
            timing_ok = False

    checks = {
        "si_units_explicit": units == _UNITS,
        "free_brake_contract_recorded": contract.get("body")
        == "uniform_solid_conducting_disc"
        and contract.get("inertia_reference") == "analytic_uniform_solid_disc"
        and contract.get("angular_momentum_balance")
        == "inertia_delta_angular_velocity_plus_integrated_braking_torque_equals_zero",
        "instantaneous_power_is_diagnostic_only": contract.get(
            "instantaneous_power_comparison"
        )
        == "diagnostic_only_when_field_energy_rate_is_not_sampled_on_the_probe_grid",
        "full_energy_storage_contract_recorded": contract.get("energy_balance")
        == "initial_kinetic_plus_magnetic_equals_final_kinetic_plus_magnetic_plus_joule",
        "analytic_disc_inertia_matches_reported": inertia_error
        <= float(maximum_inertia_relative_error),
        "replay_series_are_complete": all_cardinalities,
        "time_axes_strictly_increase": all_times_increase,
        "braking_torque_and_joule_loss_nonnegative": nonnegative_dissipation,
        "angular_velocity_monotonically_decays": monotone_decay
        and min(decay_fractions, default=-math.inf) >= float(minimum_decay_fraction),
        "angular_impulse_balance_closes": max(impulse_errors, default=math.inf)
        <= float(maximum_angular_impulse_residual),
        "fresh_replay_time_axes_match": max(replay_time_errors, default=math.inf)
        <= 1.0e-12,
        "fresh_replay_fields_match": max(replay_field_errors, default=math.inf)
        <= float(maximum_replay_error_over_span),
        "artifact_series_share_their_solve_generation": artifact_generations_ok,
        "artifact_series_share_one_coordinate_frame": artifact_coordinate_frames_ok,
        "convergence_table_matches_result_iteration": convergence_provenance_ok,
        "force_integral_uses_current_geometry_selection": (
            force_selection_identity_ok
        ),
        "sweep_excitation_uses_one_rms_basis": excitation_basis_identity_ok,
        "stored_force_matches_live_geometry_solution_and_selection": (
            live_stored_force_identity_ok
        ),
        "loss_partitions_have_unique_ownership_without_compensation": (
            loss_partition_identity_ok
        ),
        "material_property_uses_current_parameter_unit_and_generation": (
            material_property_parameter_identity_ok
        ),
        "force_selection_matches_current_geometry_topology": (
            force_selection_topology_identity_ok
        ),
        "weak_form_uses_current_jacobian_orientation": (
            weak_form_coordinate_transform_identity_ok
        ),
        "harmonic_fields_share_one_complex_time_convention": (
            time_harmonic_phasor_convention_identity_ok
        ),
        "eigenmodes_use_current_mass_inner_product_normalization": (
            eigenmode_mass_inner_product_identity_ok
        ),
        "ale_material_derivative_uses_current_mesh_velocity_time_level": (
            ale_material_derivative_time_level_identity_ok
        ),
        "harmonic_field_and_power_share_reference_time_origin": (
            harmonic_reference_time_origin_identity_ok
        ),
        "deformed_domain_integral_uses_current_geometry_jacobian": (
            deformed_domain_integral_jacobian_identity_ok
        ),
        "nonlinear_residual_and_tangent_share_material_iteration": (
            nonlinear_residual_tangent_iteration_identity_ok
        ),
        "moving_mesh_field_transfer_uses_one_coordinate_frame": (
            moving_mesh_field_transfer_frame_identity_ok
        ),
        "segregated_block_residuals_use_current_variable_scaling": (
            segregated_block_variable_scaling_identity_ok
        ),
        "modal_port_power_uses_current_surface_orientation": (
            modal_port_power_surface_orientation_identity_ok
        ),
        "degenerate_eigenmodes_use_current_subspace_tracking_basis": (
            degenerate_subspace_tracking_identity_ok
        ),
        "adaptive_bdf_restart_uses_current_history_and_event_generation": (
            adaptive_bdf_restart_identity_ok
        ),
        "nonlinear_continuation_uses_current_branch_tangent_checkpoint": (
            _nonlinear_continuation_branch_tangent_checkpoint_identity_ok(summary)
        ),
        "mortar_projection_uses_current_interface_mesh_and_quadrature": (
            _nonconforming_mortar_projection_quadrature_mesh_identity_ok(summary)
        ),
        "adaptive_field_transfer_uses_current_mesh_projection_and_conservation": (
            _adaptive_mesh_field_transfer_conservation_identity_ok(summary)
        ),
        "eigenmode_tracking_uses_current_phase_normalization_and_parameter_state": (
            _eigenmode_phase_normalization_tracking_identity_ok(summary)
        ),
        "parameter_sweep_restart_uses_current_branch_solution_and_mesh": (
            _parameter_sweep_branch_restart_identity_ok(summary)
        ),
        "multiphysics_coupling_uses_current_source_frame_units_and_selection": (
            _multiphysics_coupling_source_identity_ok(summary)
        ),
        "contact_active_set_uses_current_friction_state_mesh_and_tangent": (
            _contact_active_set_friction_state_mesh_generation_identity_ok(summary)
        ),
        "acoustic_structure_trace_uses_current_frame_impedance_and_interface": (
            _acoustic_structure_trace_impedance_order_frame_generation_identity_ok(
                summary
            )
        ),
        "nonlinear_continuation_result_uses_current_branch_load_tangent_and_mesh": (
            _continuation_branch_load_mesh_identity_ok(summary)
        ),
        "parametric_sequence_uses_current_rows_initial_solutions_and_dataset": (
            _parametric_sequence_initial_solution_identity_ok(summary)
        ),
        "multiphysics_energy_balance_uses_current_frame_time_and_generation": (
            _multiphysics_power_work_heat_balance_identity_ok(summary)
        ),
        "degenerate_eigenmodes_use_current_subspace_normalization_phase_and_mesh": (
            _degenerate_eigenmode_subspace_identity_ok(summary)
        ),
        "remeshed_fields_use_current_projection_geometry_dataset_and_conservation": (
            _remesh_field_projection_identity_ok(summary)
        ),
        "nonlinear_solutions_use_current_load_steps_branch_state_and_solver": (
            _nonlinear_continuation_load_step_identity_ok(summary)
        ),
        "ale_force_work_uses_current_geometry_time_transfer_and_energy_balance": (
            _ale_force_work_balance_identity_ok(summary)
        ),
        "segregated_solution_uses_current_iterations_relaxation_residuals_and_components": (
            _segregated_iteration_identity_ok(summary)
        ),
        "nonlinear_restart_uses_current_state_integrator_tangent_load_step_and_checkpoint": (
            _nonlinear_state_restart_identity_ok(summary)
        ),
        "floquet_modes_use_current_pair_orientation_phase_wavevector_normalization_and_dataset": (
            _floquet_pair_identity_ok(summary)
        ),
        "thermoelastic_frequency_uses_current_temperature_prestress_linearization_mesh_dataset_and_result": (
            _thermoelastic_frequency_identity_ok(summary)
        ),
        "field_circuit_uses_current_terminals_orientation_sign_gauge_power_mesh_and_solution": (
            _field_circuit_power_identity_ok(summary)
        ),
        "rotating_sliding_interface_uses_current_sector_azimuth_interpolation_frame_periodicity_mesh_and_torque": (
            _rotating_sliding_interface_identity_ok(summary)
        ),
        "acoustic_radiation_uses_current_modes_trace_area_convention_frequency_power_and_result": (
            _acoustic_radiation_impedance_identity_ok(summary)
        ),
        "joule_heat_uses_current_mapping_resistivity_temperature_frame_average_energy_mesh_and_result": (
            _joule_heat_energy_identity_ok(summary)
        ),
        "nonlinear_eigenmodes_use_current_continuation_normalization_phase_mac_branch_eigenvalues_mesh_and_result": (
            _nonlinear_eigenmode_identity_ok(summary)
        ),
        "frequency_time_reconstruction_uses_current_hermitian_spacing_window_delay_parseval_mesh_and_result": (
            _frequency_time_reconstruction_identity_ok(summary)
        ),
        "rotating_force_uses_current_virtual_work_stress_phase_frame_angle_power_mesh_and_result": (
            _rotating_force_balance_identity_ok(summary)
        ),
        "nonlinear_segregated_solutions_use_current_groups_relaxation_residual_jacobian_continuation_mesh_and_result": (
            _nonlinear_segregated_closure_identity_ok(summary)
        ),
        "degenerate_eigenmodes_use_current_subspace_phase_normalization_participation_mass_mesh_owner_and_result": (
            _degenerate_eigenmode_closure_identity_ok(summary)
        ),
        "contact_results_satisfy_current_complementarity_active_set_friction_dissipation_normal_mesh_and_result": (
            _contact_complementarity_identity_ok(summary)
        ),
        "field_circuit_dae_results_use_current_charge_current_event_energy_time_dataset_and_result": (
            _field_circuit_dae_identity_ok(summary)
        ),
        "nonlinear_continuation_uses_current_arclength_tangent_branch_turning_residual_mesh_and_result": (
            _nonlinear_arclength_identity_ok(summary)
        ),
        "electrochemical_results_use_current_species_flux_charge_mass_reaction_energy_time_mesh_and_result": (
            _electrochemical_conservation_identity_ok(summary)
        ),
        "multirate_electromechanical_results_use_current_event_timegrids_work_power_frame_mesh_and_result": (
            _multirate_electromechanical_identity_ok(summary)
        ),
        "adjoint_sensitivities_use_current_objective_design_chainrule_constraint_fd_mesh_solution_and_result": (
            _adjoint_sensitivity_identity_ok(summary)
        ),
        "magnetostatic_force_uses_current_displacement_coenergy_source_sign_mesh_frame_owner_and_result": (
            _magnetostatic_virtual_work_identity_ok(summary)
        ),
        "acoustic_modes_use_current_normalization_mass_participation_damping_reconstruction_mesh_and_result": (
            _acoustic_modal_participation_identity_ok(summary)
        ),
        "capacitance_results_use_current_matrix_charge_energy_gauge_reciprocity_terminals_mesh_and_result": (
            _capacitance_matrix_identity_ok(summary)
        ),
        "thermoelastic_harmonics_use_current_heat_phase_temperature_displacement_work_loss_frequency_mesh_and_result": (
            _thermoelastic_harmonic_identity_ok(summary)
        ),
        "thermoviscous_interfaces_use_current_velocity_traction_dissipation_power_normal_mesh_and_result": (
            _thermoviscous_pressure_interface_identity_ok(summary)
        ),
        "piezoelectric_results_use_current_charge_strain_reciprocity_energy_polarization_mesh_and_result": (
            _piezoelectric_reciprocity_identity_ok(summary)
        ),
        "poroelastic_results_use_current_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_and_result": (
            _poroelastic_biot_identity_ok(summary)
        ),
        "rotating_induction_results_use_current_slip_frequency_current_loss_torque_power_frame_mesh_and_result": (
            _rotating_induction_identity_ok(summary)
        ),
        "thermoacoustic_results_use_current_meanflow_wavenumber_flux_impedance_power_mesh_and_result": (
            _thermoacoustic_meanflow_identity_ok(summary)
        ),
        "battery_results_use_current_soc_current_heat_temperature_energy_safety_mesh_and_result": (
            _battery_electrothermal_identity_ok(summary)
        ),
        "piezoelectric_admittance_results_use_current_resonance_coupling_phase_energy_power_mesh_and_result": (
            _piezoelectric_admittance_identity_ok(summary)
        ),
        "fluidfilm_bearing_results_use_current_film_pressure_load_friction_temperature_power_mesh_and_result": (
            _fluidfilm_bearing_identity_ok(summary)
        ),
        "induction_heating_results_use_current_skin_proximity_joule_thermal_temperature_energy_mesh_and_result": (
            _induction_heating_identity_ok(summary)
        ),
        "reacting_species_results_use_current_diffusion_rate_flux_mass_temperature_mesh_and_result": (
            _reacting_species_identity_ok(summary)
        ),
        "microwave_heating_results_use_current_sparameters_power_heat_temperature_energy_mesh_and_result": (
            _microwave_heating_identity_ok(summary)
        ),
        "poroelastic_wave_results_use_current_pressure_displacement_flux_mass_energy_mesh_and_result": (
            _poroelastic_wave_identity_ok(summary)
        ),
        "restart_energy_offsets_are_continuous": restart_energy_offsets_ok,
        "field_energy_history_is_present_and_aligned": energy_cardinality
        and field_time_alignment,
        "field_energy_history_is_nonnegative_and_has_no_isolated_jump": min(
            magnetic_energy
        )
        >= 0.0
        and maximum_field_energy_adjacent_jump_fraction_observed
        <= float(maximum_field_energy_adjacent_jump_fraction)
        and maximum_field_energy_curvature_outlier_ratio_observed
        <= float(maximum_field_energy_curvature_outlier_ratio),
        "field_energy_run_replays_primary_history": energy_replay_error
        <= float(maximum_replay_error_over_span),
        "kinetic_magnetic_joule_energy_closes": total_energy_error
        <= float(maximum_total_energy_residual),
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "rotational_eddy_brake_energy_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "analytic_inertia_kg_m2": analytic_inertia,
            "reported_inertia_kg_m2": reported_inertia,
            "inertia_relative_error": inertia_error,
            "maximum_angular_impulse_residual_relative": max(
                impulse_errors, default=math.inf
            ),
            "maximum_replay_error_over_span": max(
                replay_field_errors + [energy_replay_error], default=math.inf
            ),
            "minimum_angular_velocity_decay_fraction": min(
                decay_fractions, default=-math.inf
            ),
            "kinetic_energy_drop_j": kinetic_drop,
            "magnetic_energy_drop_j": magnetic_drop,
            "maximum_field_energy_adjacent_jump_fraction": (
                maximum_field_energy_adjacent_jump_fraction_observed
            ),
            "maximum_field_energy_time_misalignment_s": (
                maximum_field_energy_time_misalignment_s_observed
            ),
            "maximum_field_energy_curvature_outlier_ratio": (
                maximum_field_energy_curvature_outlier_ratio_observed
            ),
            "integrated_joule_loss_j": joule_energy,
            "total_energy_residual_relative": total_energy_error,
            "mechanical_joule_relative_diagnostic": max(
                power_diagnostics, default=math.inf
            ),
        },
        "tolerances": {
            "maximum_inertia_relative_error": float(maximum_inertia_relative_error),
            "maximum_angular_impulse_residual": float(
                maximum_angular_impulse_residual
            ),
            "maximum_total_energy_residual": float(maximum_total_energy_residual),
            "maximum_replay_error_over_span": float(maximum_replay_error_over_span),
            "maximum_field_energy_time_misalignment_s": float(
                maximum_field_energy_time_misalignment_s
            ),
            "maximum_field_energy_adjacent_jump_fraction": float(
                maximum_field_energy_adjacent_jump_fraction
            ),
            "maximum_field_energy_curvature_outlier_ratio": float(
                maximum_field_energy_curvature_outlier_ratio
            ),
            "minimum_decay_fraction": float(minimum_decay_fraction),
        },
        "lesson": (
            "For a freely decelerating conductor, close angular impulse and the "
            "combined kinetic-plus-magnetic energy balance. Torque times speed is "
            "only diagnostic unless the magnetic-energy rate is represented."
        ),
    }
