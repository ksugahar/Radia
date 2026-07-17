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
