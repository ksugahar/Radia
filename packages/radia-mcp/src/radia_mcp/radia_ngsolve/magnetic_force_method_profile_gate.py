"""Solver-neutral profile gate for two magnetic-force formulations."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from datetime import datetime


_DIMENSION_UNITS = {
    ("axisymmetric_total", "N"),
    ("3d_total", "N"),
    ("2d_per_length", "N/m"),
}


def _profile(value: object, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    result = [float(item) for item in value]
    if not result or not all(math.isfinite(item) for item in result):
        raise ValueError(f"{name} must contain finite values")
    return result


def _maximum_relative_difference(left: Sequence[float], right: Sequence[float]) -> float:
    return max(
        abs(a - b) / max(abs(a), abs(b), 1.0e-300)
        for a, b in zip(left, right)
    )


def _trapezoid_integral(positions: Sequence[float], values: Sequence[float]) -> float:
    return sum(
        0.5 * (left_value + right_value) * (right_x - left_x)
        for left_x, right_x, left_value, right_value in zip(
            positions, positions[1:], values, values[1:]
        )
    )


def _valid_sha256(value: object) -> bool:
    digest = str(value or "").lower()
    return len(digest) == 64 and all(
        character in "0123456789abcdef" for character in digest
    )


def _bem_panel_demag_force_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("solve_generation", "")).strip()
    panel_ids = value.get("panel_ids")
    region_ids = value.get("material_region_ids")
    frame = str(value.get("force_coordinate_frame", "")).strip()
    try:
        normals = [
            [float(component) for component in row]
            for row in value.get("outward_normals", [])
        ]
        result_normals = [
            [float(component) for component in row]
            for row in value.get("result_outward_normals", [])
        ]
        demag = [float(item) for item in value.get("demag_field_a_per_m", [])]
        result_demag = [
            float(item) for item in value.get("result_demag_field_a_per_m", [])
        ]
        forces = [
            [float(component) for component in row]
            for row in value.get("force_vectors_n", [])
        ]
        result_forces = [
            [float(component) for component in row]
            for row in value.get("result_force_vectors_n", [])
        ]
    except (TypeError, ValueError):
        return False
    digest = str(value.get("panel_force_table_sha256", "")).lower()
    panels_ok = (
        isinstance(panel_ids, list)
        and bool(panel_ids)
        and all(
            isinstance(item, int) and not isinstance(item, bool) and item > 0
            for item in panel_ids
        )
        and len(set(panel_ids)) == len(panel_ids)
    )
    normals_ok = (
        panels_ok
        and len(normals) == len(panel_ids)
        and all(
            len(row) == 3
            and all(math.isfinite(component) for component in row)
            and math.isclose(
                sum(component * component for component in row),
                1.0,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            for row in normals
        )
    )
    regions_ok = (
        isinstance(region_ids, list)
        and panels_ok
        and len(region_ids) == len(panel_ids)
        and all(
            isinstance(item, int) and not isinstance(item, bool) and item > 0
            for item in region_ids
        )
    )
    vectors_ok = (
        panels_ok
        and len(forces) == len(panel_ids)
        and all(
            len(row) == 3 and all(math.isfinite(component) for component in row)
            for row in forces
        )
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "panel_mesh_solve_generation",
                "outward_normal_solve_generation",
                "material_region_solve_generation",
                "demag_result_solve_generation",
                "force_result_solve_generation",
            )
        )
        and panels_ok
        and value.get("result_panel_ids") == panel_ids
        and normals_ok
        and result_normals == normals
        and regions_ok
        and value.get("result_material_region_ids") == region_ids
        and len(demag) == len(panel_ids)
        and all(math.isfinite(item) for item in demag)
        and result_demag == demag
        and vectors_ok
        and result_forces == forces
        and frame in {"global_xyz", "local_xyz"}
        and value.get("result_force_coordinate_frame") == frame
        and _valid_sha256(digest)
        and str(value.get("result_panel_force_table_sha256", "")).lower()
        == digest
    )


def _motor_harmonic_force_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("sweep_generation", "")).strip()
    frame = str(value.get("force_coordinate_frame", "")).strip()
    bins = value.get("harmonic_bins")
    try:
        angles = [float(item) for item in value.get("rotor_angles_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_rotor_angles_deg", [])
        ]
        phases = [float(item) for item in value.get("current_phase_deg", [])]
        result_phases = [
            float(item) for item in value.get("result_current_phase_deg", [])
        ]
        harmonics = [
            [float(component) for component in row]
            for row in value.get("force_harmonics_n", [])
        ]
        result_harmonics = [
            [float(component) for component in row]
            for row in value.get("result_force_harmonics_n", [])
        ]
    except (TypeError, ValueError):
        return False
    digest = str(value.get("harmonic_force_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "rotor_angle_sweep_generation",
                "current_phase_sweep_generation",
                "harmonic_bin_sweep_generation",
                "force_frame_sweep_generation",
                "force_result_sweep_generation",
            )
        )
        and len(angles) >= 3
        and all(math.isfinite(item) for item in angles)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and result_angles == angles
        and len(phases) == 3
        and all(math.isfinite(item) for item in phases)
        and result_phases == phases
        and isinstance(bins, list)
        and bool(bins)
        and all(
            isinstance(item, int) and not isinstance(item, bool) and item >= 0
            for item in bins
        )
        and bins == sorted(set(bins))
        and value.get("result_harmonic_bins") == bins
        and frame in {"global_xyz", "rotor_dq"}
        and value.get("result_force_coordinate_frame") == frame
        and len(harmonics) == len(bins)
        and all(
            len(row) == 2 and all(math.isfinite(component) for component in row)
            for row in harmonics
        )
        and result_harmonics == harmonics
        and _valid_sha256(digest)
        and str(value.get("result_harmonic_force_table_sha256", "")).lower()
        == digest
    )


def _maglev_force_energy_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("maglev_generation", "")).strip()
    frame = str(value.get("coordinate_frame", "")).strip()
    sign = str(value.get("force_energy_sign_convention", "")).strip()
    try:
        equilibrium_x = float(value.get("equilibrium_displacement_m"))
        result_equilibrium_x = float(value.get("result_equilibrium_displacement_m"))
        equilibrium_force = float(value.get("equilibrium_force_n"))
        result_equilibrium_force = float(value.get("result_equilibrium_force_n"))
        displacements = [
            float(item) for item in value.get("displacement_samples_m", [])
        ]
        energies = [float(item) for item in value.get("energy_samples_j", [])]
        forces = [float(item) for item in value.get("force_samples_n", [])]
        energy_forces = [
            float(item)
            for item in value.get("energy_finite_difference_force_n", [])
        ]
        stiffness = float(value.get("stiffness_n_m"))
        reported_stiffness = float(value.get("reported_stiffness_n_m"))
    except (TypeError, ValueError):
        return False
    mesh_digest = str(value.get("mesh_sha256", "")).lower()
    result_digest = str(value.get("maglev_result_sha256", "")).lower()
    finite_values = [
        equilibrium_x,
        result_equilibrium_x,
        equilibrium_force,
        result_equilibrium_force,
        stiffness,
        reported_stiffness,
        *displacements,
        *energies,
        *forces,
        *energy_forces,
    ]
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "equilibrium_maglev_generation",
                "displacement_maglev_generation",
                "energy_maglev_generation",
                "force_maglev_generation",
                "stiffness_maglev_generation",
                "coordinate_frame_maglev_generation",
                "result_maglev_generation",
            )
        )
        and len(displacements) >= 3
        and len(energies) == len(forces) == len(energy_forces) == len(displacements)
        and all(math.isfinite(item) for item in finite_values)
        and all(left < right for left, right in zip(displacements, displacements[1:]))
        and equilibrium_x in displacements
        and math.isclose(result_equilibrium_x, equilibrium_x, abs_tol=1.0e-15)
        and math.isclose(result_equilibrium_force, equilibrium_force, abs_tol=1.0e-12)
        and stiffness > 0.0
        and math.isclose(reported_stiffness, stiffness, rel_tol=1.0e-12)
        and all(
            math.isclose(energy, 0.5 * stiffness * displacement**2, rel_tol=1.0e-10, abs_tol=1.0e-14)
            for displacement, energy in zip(displacements, energies)
        )
        and all(
            math.isclose(force, -stiffness * displacement, rel_tol=1.0e-10, abs_tol=1.0e-12)
            for displacement, force in zip(displacements, forces)
        )
        and energy_forces == forces
        and math.isclose(equilibrium_force, -stiffness * equilibrium_x, rel_tol=1.0e-10, abs_tol=1.0e-12)
        and sign == "force=-dW/dx"
        and value.get("result_force_energy_sign_convention") == sign
        and frame == "global_z"
        and value.get("result_coordinate_frame") == frame
        and _valid_sha256(mesh_digest)
        and str(value.get("result_mesh_sha256", "")).lower() == mesh_digest
        and _valid_sha256(result_digest)
        and str(value.get("reported_maglev_result_sha256", "")).lower()
        == result_digest
    )


def _motor_dual_lane_alignment_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("comparison_generation", "")).strip()
    lane_ids = value.get("lane_ids")
    geometries = value.get("geometry_revision_sha256")
    excitations = value.get("excitation_table_sha256")
    frames = value.get("force_coordinate_frames")
    bins = value.get("harmonic_bins")
    angles = value.get("rotor_angles_deg")
    try:
        harmonics = [
            [[float(component) for component in row] for row in lane]
            for lane in value.get("force_harmonics_n", [])
        ]
        result_harmonics = [
            [[float(component) for component in row] for row in lane]
            for lane in value.get("result_force_harmonics_n", [])
        ]
        parsed_angles = [
            [float(item) for item in lane] for lane in angles
        ] if isinstance(angles, list) else []
    except (TypeError, ValueError):
        return False
    result_digest = str(value.get("comparison_result_sha256", "")).lower()
    bins_ok = (
        isinstance(bins, list)
        and len(bins) == 2
        and all(
            isinstance(lane, list)
            and bool(lane)
            and all(
                isinstance(item, int) and not isinstance(item, bool) and item >= 0
                for item in lane
            )
            and lane == sorted(set(lane))
            for lane in bins
        )
        and bins[0] == bins[1]
    )
    angles_ok = (
        len(parsed_angles) == 2
        and all(len(lane) >= 3 for lane in parsed_angles)
        and all(
            all(math.isfinite(item) for item in lane)
            and all(left < right for left, right in zip(lane, lane[1:]))
            for lane in parsed_angles
        )
        and parsed_angles[0] == parsed_angles[1]
    )
    harmonics_ok = (
        bins_ok
        and len(harmonics) == 2
        and all(len(lane) == len(bins[0]) for lane in harmonics)
        and all(
            len(row) == 2 and all(math.isfinite(component) for component in row)
            for lane in harmonics
            for row in lane
        )
        and harmonics[0] == harmonics[1]
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "geometry_comparison_generation",
                "excitation_comparison_generation",
                "force_frame_comparison_generation",
                "harmonic_comparison_generation",
                "rotor_angle_comparison_generation",
                "result_comparison_generation",
            )
        )
        and lane_ids == ["ngsolve-age", "hdiv-mmm-hcurl-eddy-bubble"]
        and value.get("result_lane_ids") == lane_ids
        and isinstance(geometries, list)
        and len(geometries) == 2
        and all(_valid_sha256(item) for item in geometries)
        and geometries[0] == geometries[1]
        and value.get("result_geometry_revision_sha256") == geometries
        and isinstance(excitations, list)
        and len(excitations) == 2
        and all(_valid_sha256(item) for item in excitations)
        and excitations[0] == excitations[1]
        and value.get("result_excitation_table_sha256") == excitations
        and frames == ["rotor_dq", "rotor_dq"]
        and value.get("result_force_coordinate_frames") == frames
        and bins_ok
        and value.get("result_harmonic_bins") == bins
        and angles_ok
        and value.get("result_rotor_angles_deg") == angles
        and harmonics_ok
        and result_harmonics == harmonics
        and _valid_sha256(result_digest)
        and str(value.get("reported_comparison_result_sha256", "")).lower()
        == result_digest
    )


def magnetic_force_method_profile_gate(
    summary: Mapping[str, object],
    *,
    maximum_method_relative_difference: float = 0.05,
    maximum_independent_stress_relative_difference: float = 0.02,
    minimum_selection_scope_relative_difference: float = 0.25,
    maximum_all_body_to_target_magnitude_ratio: float = 0.75,
    maximum_work_relative_difference: float = 0.05,
    maximum_parsed_replay_absolute_difference: float = 1.0e-12,
    minimum_sample_count: int = 5,
) -> dict[str, object]:
    """Compare target-body element force with closed-surface stress force.

    The gate deliberately requires an all-body element-force control. This catches
    a common false comparison where two force formulations are evaluated over
    different bodies or surfaces while being presented as method disagreement.
    """
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    if minimum_sample_count < 3:
        raise ValueError("minimum_sample_count must be at least 3")

    tolerances = {
        "maximum_method_relative_difference": float(maximum_method_relative_difference),
        "maximum_independent_stress_relative_difference": float(
            maximum_independent_stress_relative_difference
        ),
        "minimum_selection_scope_relative_difference": float(
            minimum_selection_scope_relative_difference
        ),
        "maximum_all_body_to_target_magnitude_ratio": float(
            maximum_all_body_to_target_magnitude_ratio
        ),
        "maximum_work_relative_difference": float(maximum_work_relative_difference),
        "maximum_parsed_replay_absolute_difference": float(
            maximum_parsed_replay_absolute_difference
        ),
    }
    if not all(math.isfinite(value) and value >= 0.0 for value in tolerances.values()):
        raise ValueError("all tolerances must be finite and nonnegative")

    profiles = {
        name: _profile(summary.get(name), name)
        for name in (
            "positions",
            "moving_body_element_force",
            "closed_surface_maxwell_stress_force",
            "independent_closed_surface_force",
            "all_body_element_force",
        )
    }
    lengths = {len(values) for values in profiles.values()}
    if len(lengths) != 1:
        raise ValueError("positions and all force profiles must have the same length")

    positions = profiles["positions"]
    target = profiles["moving_body_element_force"]
    stress = profiles["closed_surface_maxwell_stress_force"]
    independent_stress = profiles["independent_closed_surface_force"]
    all_body = profiles["all_body_element_force"]
    increasing_positions = all(right > left for left, right in zip(positions, positions[1:]))

    identity_value = summary.get("artifact_identity")
    identity_present = isinstance(identity_value, Mapping)
    one_sweep_generation_ok = True
    demag_reference_ok = True
    coordinate_system_binding_ok = True
    force_normalization_ok = True
    hysteresis_branch_state_ok = True
    remanence_frame_binding_ok = True
    force_surface_body_ownership_ok = True
    demag_branch_interpolation_ok = True
    linear_motor_thrust_phase_identity_ok = True
    demag_recoil_temperature_identity_ok = True
    bem_demag_surface_normal_generation_identity_ok = True
    cogging_torque_periodic_sector_symmetry_identity_ok = True
    bem_self_term_solid_angle_orientation_identity_ok = True
    demag_energy_force_displacement_length_unit_identity_ok = True
    bem_near_singular_quadrature_target_scale_identity_ok = True
    force_torque_reference_origin_length_unit_identity_ok = True
    bem_solid_angle_surface_winding_identity_ok = True
    maglev_stiffness_force_displacement_generation_identity_ok = True
    bem_demag_tensor_coordinate_basis_generation_identity_ok = True
    magnetic_bearing_force_harmonic_phase_origin_identity_ok = True
    bem_near_singular_panel_subdivision_quadrature_generation_identity_ok = True
    maglev_force_coil_polarity_orientation_generation_identity_ok = True
    bem_demag_self_term_solid_angle_orientation_generation_identity_ok = True
    virtual_work_displacement_coordinate_geometry_generation_identity_ok = True
    demag_energy_surface_charge_normal_quadrature_identity_ok = True
    maglev_stiffness_displacement_equilibrium_force_identity_ok = True
    bem_near_singular_distance_panel_quadrature_identity_ok = True
    moving_magnet_force_position_orientation_equilibrium_identity_ok = True
    motor_force_dual_lane_generation_identity_ok = True
    linear_motor_end_effect_generation_identity_ok = True
    bem_panel_demag_force_generation_identity_ok = True
    motor_harmonic_force_generation_identity_ok = True
    maglev_force_energy_generation_identity_ok = True
    motor_dual_lane_alignment_generation_identity_ok = True
    if identity_value is not None and not identity_present:
        one_sweep_generation_ok = False
        demag_reference_ok = False
        coordinate_system_binding_ok = False
        force_normalization_ok = False
        hysteresis_branch_state_ok = False
        remanence_frame_binding_ok = False
        force_surface_body_ownership_ok = False
        demag_branch_interpolation_ok = False
        linear_motor_thrust_phase_identity_ok = False
        demag_recoil_temperature_identity_ok = False
        bem_demag_surface_normal_generation_identity_ok = False
        cogging_torque_periodic_sector_symmetry_identity_ok = False
        bem_self_term_solid_angle_orientation_identity_ok = False
        demag_energy_force_displacement_length_unit_identity_ok = False
        bem_near_singular_quadrature_target_scale_identity_ok = False
        force_torque_reference_origin_length_unit_identity_ok = False
        bem_solid_angle_surface_winding_identity_ok = False
        maglev_stiffness_force_displacement_generation_identity_ok = False
        bem_demag_tensor_coordinate_basis_generation_identity_ok = False
        magnetic_bearing_force_harmonic_phase_origin_identity_ok = False
        bem_near_singular_panel_subdivision_quadrature_generation_identity_ok = False
        maglev_force_coil_polarity_orientation_generation_identity_ok = False
        bem_demag_self_term_solid_angle_orientation_generation_identity_ok = False
        virtual_work_displacement_coordinate_geometry_generation_identity_ok = False
        demag_energy_surface_charge_normal_quadrature_identity_ok = False
        maglev_stiffness_displacement_equilibrium_force_identity_ok = False
        bem_near_singular_distance_panel_quadrature_identity_ok = False
        moving_magnet_force_position_orientation_equilibrium_identity_ok = False
        motor_force_dual_lane_generation_identity_ok = False
        linear_motor_end_effect_generation_identity_ok = False
        bem_panel_demag_force_generation_identity_ok = False
        motor_harmonic_force_generation_identity_ok = False
        maglev_force_energy_generation_identity_ok = False
        motor_dual_lane_alignment_generation_identity_ok = False
    elif identity_present:
        generations = identity_value.get("position_force_sample_generations")
        timestamps = identity_value.get("sample_acquired_at_utc")

        def timestamp(value: object) -> float:
            try:
                return datetime.fromisoformat(
                    str(value).replace("Z", "+00:00")
                ).timestamp()
            except (TypeError, ValueError):
                return math.nan

        parsed_times = (
            [timestamp(value) for value in timestamps]
            if isinstance(timestamps, Sequence)
            and not isinstance(timestamps, (str, bytes))
            else []
        )
        one_sweep_generation_ok = (
            isinstance(generations, Sequence)
            and not isinstance(generations, (str, bytes))
            and len(generations) == len(positions)
            and all(bool(str(value)) for value in generations)
            and len({str(value) for value in generations if str(value)}) == 1
            and len(parsed_times) == len(positions)
            and all(math.isfinite(value) for value in parsed_times)
            and all(right > left for left, right in zip(parsed_times, parsed_times[1:]))
        )
        geometry = identity_value.get("magnet_geometry")
        reference = identity_value.get("demag_reference")
        if not isinstance(geometry, Mapping) or not isinstance(reference, Mapping):
            demag_reference_ok = False
        else:
            geometry_revision = str(geometry.get("revision", ""))
            committed_at = timestamp(geometry.get("committed_at_utc"))
            generated_at = timestamp(reference.get("generated_at_utc"))
            demag_reference_ok = (
                bool(geometry_revision)
                and str(reference.get("geometry_revision", "")) == geometry_revision
                and math.isfinite(committed_at)
                and math.isfinite(generated_at)
                and generated_at >= committed_at
            )
        coordinate_binding = identity_value.get("coordinate_system_binding")
        if coordinate_binding is not None:
            coordinate_system_binding_ok = (
                isinstance(coordinate_binding, Mapping)
                and bool(coordinate_binding.get("common_frame_id"))
                and coordinate_binding.get("force_component_frame_id")
                == coordinate_binding.get("common_frame_id")
                and coordinate_binding.get("demag_metric_frame_id")
                == coordinate_binding.get("common_frame_id")
                and bool(coordinate_binding.get("geometry_rotation_revision"))
                and coordinate_binding.get("force_transform_revision")
                == coordinate_binding.get("geometry_rotation_revision")
                and coordinate_binding.get("demag_transform_revision")
                == coordinate_binding.get("geometry_rotation_revision")
            )
        force_normalization = identity_value.get("force_normalization")
        if force_normalization is not None:
            profile_bases = (
                force_normalization.get("profile_bases")
                if isinstance(force_normalization, Mapping)
                else None
            )
            force_normalization_ok = (
                isinstance(force_normalization, Mapping)
                and force_normalization.get("comparison_basis") == "total_3d"
                and isinstance(profile_bases, Mapping)
                and set(profile_bases)
                == {
                    "moving_body_element_force",
                    "closed_surface_maxwell_stress_force",
                    "independent_closed_surface_force",
                }
                and all(value == "total_3d" for value in profile_bases.values())
                and force_normalization.get("per_length_to_total_applied") is True
            )
        branch_state = identity_value.get("hysteresis_branch_state")
        if branch_state is not None:
            observable_branch = (
                str(branch_state.get("observable_branch", ""))
                if isinstance(branch_state, Mapping)
                else ""
            )
            hysteresis_branch_state_ok = (
                isinstance(branch_state, Mapping)
                and observable_branch in {"ascending", "descending"}
                and branch_state.get("state_memory_branch") == observable_branch
                and branch_state.get("tangent_branch") == observable_branch
                and bool(branch_state.get("branch_state_generation"))
                and branch_state.get("tangent_state_generation")
                == branch_state.get("branch_state_generation")
            )
        remanence_frame = identity_value.get("remanence_frame_binding")
        if remanence_frame is not None:
            vector_frame = (
                str(remanence_frame.get("remanence_vector_frame_id", ""))
                if isinstance(remanence_frame, Mapping)
                else ""
            )
            assembly_frame = (
                str(remanence_frame.get("assembly_frame_id", ""))
                if isinstance(remanence_frame, Mapping)
                else ""
            )
            remanence_frame_binding_ok = (
                isinstance(remanence_frame, Mapping)
                and bool(vector_frame)
                and bool(assembly_frame)
                and remanence_frame.get("transform_input_frame_id") == vector_frame
                and remanence_frame.get("transform_output_frame_id")
                == assembly_frame
                and remanence_frame.get("transformed_vector_frame_id")
                == assembly_frame
                and remanence_frame.get("transform_applied")
                is (vector_frame != assembly_frame)
                and bool(remanence_frame.get("geometry_rotation_revision"))
                and remanence_frame.get("remanence_transform_revision")
                == remanence_frame.get("geometry_rotation_revision")
            )
        surface_ownership = identity_value.get("force_surface_body_ownership")
        if surface_ownership is not None:
            target_body = (
                str(surface_ownership.get("target_body_id", ""))
                if isinstance(surface_ownership, Mapping)
                else ""
            )
            enclosed = (
                surface_ownership.get("enclosed_body_ids")
                if isinstance(surface_ownership, Mapping)
                else None
            )
            enclosed_ids = (
                [str(value) for value in enclosed]
                if isinstance(enclosed, Sequence)
                and not isinstance(enclosed, (str, bytes))
                else []
            )
            selection_generation = (
                str(surface_ownership.get("surface_selection_generation", ""))
                if isinstance(surface_ownership, Mapping)
                else ""
            )
            force_surface_body_ownership_ok = (
                isinstance(surface_ownership, Mapping)
                and bool(target_body)
                and enclosed_ids == [target_body]
                and len(set(enclosed_ids)) == len(enclosed_ids)
                and bool(selection_generation)
                and surface_ownership.get("force_integration_generation")
                == selection_generation
                and surface_ownership.get("compensating_body_force_allowed")
                is False
            )
        branch_interpolation = identity_value.get("demag_branch_interpolation")
        if branch_interpolation is not None:
            operating_branch = (
                str(branch_interpolation.get("operating_point_branch", ""))
                if isinstance(branch_interpolation, Mapping)
                else ""
            )
            brackets = (
                branch_interpolation.get("bracketing_sample_ids")
                if isinstance(branch_interpolation, Mapping)
                else None
            )
            bracket_ids = (
                [str(value) for value in brackets]
                if isinstance(brackets, Sequence)
                and not isinstance(brackets, (str, bytes))
                else []
            )
            branch_generation = (
                str(branch_interpolation.get("branch_state_generation", ""))
                if isinstance(branch_interpolation, Mapping)
                else ""
            )
            demag_branch_interpolation_ok = (
                isinstance(branch_interpolation, Mapping)
                and operating_branch in {"ascending", "descending"}
                and branch_interpolation.get("interpolation_source_branch")
                == operating_branch
                and bool(branch_generation)
                and branch_interpolation.get("interpolation_state_generation")
                == branch_generation
                and len(bracket_ids) == 2
                and all(bracket_ids)
                and len(set(bracket_ids)) == 2
            )

        phase_value = identity_value.get("linear_motor_thrust_phase_identity")
        if phase_value is not None:
            phase_identity = phase_value if isinstance(phase_value, Mapping) else {}
            winding_sequence = phase_identity.get("winding_phase_sequence")
            thrust_sequence = phase_identity.get("thrust_phase_sequence")
            winding_direction = phase_identity.get(
                "winding_electrical_angle_direction"
            )
            thrust_direction = phase_identity.get(
                "thrust_electrical_angle_direction"
            )
            phase_generation = str(
                phase_identity.get("phase_convention_generation", "")
            )
            linear_motor_thrust_phase_identity_ok = (
                isinstance(winding_sequence, list)
                and len(winding_sequence) == 3
                and all(isinstance(name, str) and name for name in winding_sequence)
                and len(set(winding_sequence)) == 3
                and thrust_sequence == winding_sequence
                and type(winding_direction) is int
                and type(thrust_direction) is int
                and winding_direction in {-1, 1}
                and thrust_direction == winding_direction
                and bool(phase_generation)
                and phase_identity.get("thrust_observable_phase_generation")
                == phase_generation
            )

        recoil_value = identity_value.get("demag_recoil_temperature_identity")
        if recoil_value is not None:
            recoil_identity = (
                recoil_value if isinstance(recoil_value, Mapping) else {}
            )
            try:
                evaluation_temperature = float(
                    recoil_identity["evaluation_temperature_c"]
                )
                material_temperature = float(
                    recoil_identity["magnet_material_temperature_c"]
                )
                recoil_temperature = float(
                    recoil_identity["recoil_line_temperature_c"]
                )
            except (KeyError, TypeError, ValueError):
                evaluation_temperature = material_temperature = math.nan
                recoil_temperature = math.nan
            material_generation = str(
                recoil_identity.get("material_state_generation", "")
            )
            recoil_digest = str(recoil_identity.get("recoil_line_sha256", ""))
            demag_recoil_temperature_identity_ok = (
                all(
                    math.isfinite(value)
                    for value in (
                        evaluation_temperature,
                        material_temperature,
                        recoil_temperature,
                    )
                )
                and material_temperature == evaluation_temperature
                and recoil_temperature == evaluation_temperature
                and bool(material_generation)
                and recoil_identity.get("recoil_line_state_generation")
                == material_generation
                and len(recoil_digest) == 64
                and all(character in "0123456789abcdef" for character in recoil_digest)
            )

        normal_value = identity_value.get(
            "bem_demag_surface_normal_generation_identity"
        )
        if normal_value is not None:
            normal_identity = normal_value if isinstance(normal_value, Mapping) else {}
            mesh_generation = str(
                normal_identity.get("active_surface_mesh_generation", "")
            )
            normal_digest = str(
                normal_identity.get("surface_normal_sha256", "")
            ).lower()
            kernel_digest = str(
                normal_identity.get("demag_kernel_normal_sha256", "")
            ).lower()
            bem_demag_surface_normal_generation_identity_ok = (
                bool(mesh_generation)
                and normal_identity.get("surface_element_generation")
                == mesh_generation
                and normal_identity.get("surface_normal_generation")
                == mesh_generation
                and normal_identity.get("demag_evaluation_surface_generation")
                == mesh_generation
                and normal_identity.get("normal_orientation") == "outward"
                and normal_identity.get("demag_kernel_normal_orientation")
                == "outward"
                and len(normal_digest) == len(kernel_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in normal_digest + kernel_digest
                )
                and kernel_digest == normal_digest
            )

        symmetry_value = identity_value.get(
            "cogging_torque_periodic_sector_symmetry_identity"
        )
        if symmetry_value is not None:
            symmetry = symmetry_value if isinstance(symmetry_value, Mapping) else {}
            sector_count = symmetry.get("active_periodic_sector_count")
            result_sector_count = symmetry.get("torque_result_periodic_sector_count")
            try:
                multiplier = float(symmetry["symmetry_multiplier"])
            except (KeyError, TypeError, ValueError):
                multiplier = math.nan
            topology_generation = str(
                symmetry.get("periodic_topology_generation", "")
            )
            cogging_torque_periodic_sector_symmetry_identity_ok = (
                type(sector_count) is int
                and sector_count > 0
                and result_sector_count == sector_count
                and math.isfinite(multiplier)
                and math.isclose(
                    multiplier,
                    float(sector_count),
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                )
                and symmetry.get("sector_torque_scope") == "one_periodic_sector"
                and symmetry.get("reported_torque_scope") == "full_machine"
                and bool(topology_generation)
                and symmetry.get("torque_result_topology_generation")
                == topology_generation
                and symmetry.get("multiplier_topology_generation")
                == topology_generation
            )

        self_term_value = identity_value.get(
            "bem_self_term_solid_angle_orientation_identity"
        )
        if self_term_value is not None:
            self_term = self_term_value if isinstance(self_term_value, Mapping) else {}
            mesh_generation = str(
                self_term.get("active_surface_mesh_generation", "")
            )
            orientation = str(self_term.get("panel_orientation", ""))
            sign_convention = str(
                self_term.get("solid_angle_sign_convention", "")
            )
            orientation_digest = str(
                self_term.get("panel_orientation_sha256", "")
            ).lower()
            bem_self_term_solid_angle_orientation_identity_ok = (
                bool(mesh_generation)
                and self_term.get("panel_generation") == mesh_generation
                and self_term.get("panel_orientation_generation")
                == mesh_generation
                and self_term.get("self_term_orientation_generation")
                == mesh_generation
                and orientation in {"outward", "inward"}
                and self_term.get("self_term_solid_angle_orientation")
                == orientation
                and sign_convention in {"outward_positive", "inward_positive"}
                and self_term.get("self_term_sign_convention") == sign_convention
                and len(orientation_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in orientation_digest
                )
                and self_term.get("self_term_orientation_sha256")
                == orientation_digest
            )

        displacement_unit_value = identity_value.get(
            "demag_energy_force_displacement_length_unit_identity"
        )
        if displacement_unit_value is not None:
            displacement_unit = (
                displacement_unit_value
                if isinstance(displacement_unit_value, Mapping)
                else {}
            )
            length_units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            numeric_unit = str(
                displacement_unit.get("displacement_length_unit", "")
            )
            derivative_unit = str(
                displacement_unit.get("force_derivative_length_unit", "")
            )
            try:
                numeric_scale = float(
                    displacement_unit.get("displacement_scale_to_m")
                )
                derivative_scale = float(
                    displacement_unit.get("force_derivative_scale_to_m")
                )
            except (TypeError, ValueError):
                numeric_scale = derivative_scale = math.nan
            displacement_generation = str(
                displacement_unit.get("displacement_generation", "")
            )
            grid_digest = str(
                displacement_unit.get("displacement_grid_sha256", "")
            ).lower()
            expected_scale = length_units.get(numeric_unit)
            demag_energy_force_displacement_length_unit_identity_ok = (
                bool(str(displacement_unit.get("energy_generation", "")))
                and bool(displacement_generation)
                and displacement_unit.get(
                    "force_derivative_displacement_generation"
                )
                == displacement_generation
                and displacement_unit.get("energy_unit") == "J"
                and displacement_unit.get("force_unit") == "N"
                and expected_scale is not None
                and derivative_unit == numeric_unit
                and math.isclose(
                    numeric_scale, expected_scale, rel_tol=0.0, abs_tol=0.0
                )
                and math.isclose(
                    derivative_scale, numeric_scale, rel_tol=0.0, abs_tol=0.0
                )
                and len(grid_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in grid_digest
                )
                and displacement_unit.get("force_derivative_grid_sha256")
                == grid_digest
                and displacement_unit.get("force_from_energy_convention")
                == "negative_energy_gradient"
            )

        near_singular_value = identity_value.get(
            "bem_near_singular_quadrature_target_scale_identity"
        )
        if near_singular_value is not None:
            near_singular = (
                near_singular_value
                if isinstance(near_singular_value, Mapping)
                else {}
            )
            try:
                source_panel_id = int(near_singular.get("source_panel_id"))
                target_panel_id = int(near_singular.get("target_panel_id"))
                separation = float(
                    near_singular.get("target_panel_separation_m")
                )
                characteristic_length = float(
                    near_singular.get("target_panel_characteristic_length_m")
                )
                normalized_separation = float(
                    near_singular.get("normalized_target_separation")
                )
                quadrature_separation = float(
                    near_singular.get("quadrature_normalized_target_separation")
                )
                quadrature_order = int(near_singular.get("quadrature_order"))
            except (TypeError, ValueError):
                source_panel_id = target_panel_id = -1
                separation = characteristic_length = math.nan
                normalized_separation = quadrature_separation = math.nan
                quadrature_order = -1
            mesh_generation = str(
                near_singular.get("active_surface_mesh_generation", "")
            )
            separation_generation = str(
                near_singular.get("target_separation_generation", "")
            )
            target_digest = str(
                near_singular.get("target_pair_sha256", "")
            ).lower()
            bem_near_singular_quadrature_target_scale_identity_ok = (
                bool(mesh_generation)
                and near_singular.get("source_panel_mesh_generation")
                == mesh_generation
                and near_singular.get("target_panel_mesh_generation")
                == mesh_generation
                and bool(separation_generation)
                and near_singular.get("quadrature_target_separation_generation")
                == separation_generation
                and source_panel_id >= 0
                and target_panel_id >= 0
                and source_panel_id != target_panel_id
                and math.isfinite(separation)
                and separation > 0.0
                and math.isfinite(characteristic_length)
                and characteristic_length > 0.0
                and math.isfinite(normalized_separation)
                and math.isclose(
                    normalized_separation,
                    separation / characteristic_length,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
                and math.isfinite(quadrature_separation)
                and math.isclose(
                    quadrature_separation,
                    normalized_separation,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
                and near_singular.get("quadrature_rule") == "adaptive_duffy"
                and quadrature_order >= 2
                and len(target_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in target_digest
                )
                and str(
                    near_singular.get("quadrature_target_pair_sha256", "")
                ).lower()
                == target_digest
            )

        reference_origin_value = identity_value.get(
            "force_torque_reference_origin_length_unit_identity"
        )
        if reference_origin_value is not None:
            reference_origin = (
                reference_origin_value
                if isinstance(reference_origin_value, Mapping)
                else {}
            )
            length_units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            reference_unit = str(
                reference_origin.get("reference_origin_length_unit", "")
            )
            try:
                reference_scale = float(
                    reference_origin.get("reference_origin_scale_to_m")
                )
                force_scale = float(
                    reference_origin.get("force_reference_origin_scale_to_m")
                )
                torque_scale = float(
                    reference_origin.get("torque_reference_origin_scale_to_m")
                )
                origin = [
                    float(value)
                    for value in reference_origin.get(
                        "reference_origin_coordinates", []
                    )
                ]
            except (TypeError, ValueError):
                reference_scale = force_scale = torque_scale = math.nan
                origin = []
            solve_generation = str(reference_origin.get("solve_generation", ""))
            origin_digest = str(
                reference_origin.get("reference_origin_sha256", "")
            ).lower()
            expected_scale = length_units.get(reference_unit)
            force_torque_reference_origin_length_unit_identity_ok = (
                bool(solve_generation)
                and reference_origin.get("force_result_generation")
                == solve_generation
                and reference_origin.get("torque_result_generation")
                == solve_generation
                and reference_origin.get("force_frame_id") == "global-cartesian"
                and reference_origin.get("torque_frame_id")
                == reference_origin.get("force_frame_id")
                and reference_origin.get("force_unit") == "N"
                and reference_origin.get("torque_unit") == "N*m"
                and expected_scale is not None
                and reference_origin.get("force_reference_origin_length_unit")
                == reference_unit
                and reference_origin.get("torque_reference_origin_length_unit")
                == reference_unit
                and math.isclose(
                    reference_scale, expected_scale, rel_tol=0.0, abs_tol=0.0
                )
                and math.isclose(force_scale, reference_scale, rel_tol=0.0, abs_tol=0.0)
                and math.isclose(
                    torque_scale, reference_scale, rel_tol=0.0, abs_tol=0.0
                )
                and len(origin) == 3
                and all(math.isfinite(value) for value in origin)
                and reference_origin.get("force_reference_origin_coordinates")
                == origin
                and reference_origin.get("torque_reference_origin_coordinates")
                == origin
                and len(origin_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in origin_digest
                )
                and str(
                    reference_origin.get("torque_reference_origin_sha256", "")
                ).lower()
                == origin_digest
            )

        winding_value = identity_value.get(
            "bem_solid_angle_surface_winding_identity"
        )
        if winding_value is not None:
            winding = winding_value if isinstance(winding_value, Mapping) else {}
            mesh_generation = str(winding.get("surface_mesh_generation", ""))
            component_ids = winding.get("surface_component_ids")
            winding_digest = str(
                winding.get("surface_winding_sha256", "")
            ).lower()
            bem_solid_angle_surface_winding_identity_ok = (
                bool(mesh_generation)
                and winding.get("normalized_surface_winding_generation")
                == mesh_generation
                and winding.get("solid_angle_sign_generation") == mesh_generation
                and winding.get("self_term_assembly_generation") == mesh_generation
                and winding.get("surface_winding_normalized") is True
                and winding.get("solid_angle_sign_convention") == "outward_positive"
                and winding.get("self_term_sign_convention")
                == winding.get("solid_angle_sign_convention")
                and isinstance(component_ids, list)
                and bool(component_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in component_ids
                )
                and len(set(component_ids)) == len(component_ids)
                and winding.get("solid_angle_component_ids") == component_ids
                and len(winding_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in winding_digest
                )
                and str(winding.get("solid_angle_winding_sha256", "")).lower()
                == winding_digest
            )

        stiffness_value = identity_value.get(
            "maglev_stiffness_force_displacement_generation_identity"
        )
        if stiffness_value is not None:
            stiffness = stiffness_value if isinstance(stiffness_value, Mapping) else {}
            perturbation_generation = str(
                stiffness.get("perturbation_generation", "")
            )
            displacement_digest = str(
                stiffness.get("displacement_sha256", "")
            ).lower()
            force_digest = str(stiffness.get("force_sample_sha256", "")).lower()
            maglev_stiffness_force_displacement_generation_identity_ok = (
                bool(perturbation_generation)
                and stiffness.get("displacement_coordinate_generation")
                == perturbation_generation
                and stiffness.get("force_sample_generation")
                == perturbation_generation
                and stiffness.get("stiffness_derivative_generation")
                == perturbation_generation
                and stiffness.get("displacement_axis") == "global-z"
                and stiffness.get("force_component_axis")
                == stiffness.get("displacement_axis")
                and stiffness.get("displacement_unit") == "m"
                and stiffness.get("force_unit") == "N"
                and stiffness.get("stiffness_unit") == "N/m"
                and len(displacement_digest) == len(force_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in displacement_digest + force_digest
                )
                and str(
                    stiffness.get("stiffness_displacement_sha256", "")
                ).lower()
                == displacement_digest
                and str(
                    stiffness.get("stiffness_force_sample_sha256", "")
                ).lower()
                == force_digest
            )

        demag_tensor_value = identity_value.get(
            "bem_demag_tensor_coordinate_basis_generation_identity"
        )
        if demag_tensor_value is not None:
            demag_tensor = (
                demag_tensor_value
                if isinstance(demag_tensor_value, Mapping)
                else {}
            )
            placement_generation = str(
                demag_tensor.get("body_placement_generation", "")
            )
            tensor_generation = str(
                demag_tensor.get("demag_tensor_generation", "")
            )
            transform_digest = str(
                demag_tensor.get("body_to_global_transform_sha256", "")
            ).lower()
            try:
                transform_determinant = float(
                    demag_tensor.get("body_to_global_transform_determinant")
                )
                orthogonality_error = float(
                    demag_tensor.get("body_to_global_transform_orthogonality_error")
                )
            except (TypeError, ValueError):
                transform_determinant = orthogonality_error = math.nan
            bem_demag_tensor_coordinate_basis_generation_identity_ok = (
                bool(placement_generation)
                and demag_tensor.get("surface_mesh_body_placement_generation")
                == placement_generation
                and bool(tensor_generation)
                and demag_tensor.get("force_demag_tensor_generation")
                == tensor_generation
                and demag_tensor.get("demag_tensor_body_placement_generation")
                == placement_generation
                and demag_tensor.get("body_coordinate_basis")
                == "body-local-current"
                and demag_tensor.get("tensor_coordinate_basis")
                == demag_tensor.get("body_coordinate_basis")
                and demag_tensor.get("body_basis_handedness") == "right_handed"
                and demag_tensor.get("tensor_basis_handedness")
                == demag_tensor.get("body_basis_handedness")
                and math.isclose(
                    transform_determinant, 1.0, rel_tol=0.0, abs_tol=1.0e-12
                )
                and math.isfinite(orthogonality_error)
                and 0.0 <= orthogonality_error <= 1.0e-12
                and len(transform_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in transform_digest
                )
                and str(
                    demag_tensor.get("tensor_basis_transform_sha256", "")
                ).lower()
                == transform_digest
            )

        bearing_phase_value = identity_value.get(
            "magnetic_bearing_force_harmonic_phase_origin_identity"
        )
        if bearing_phase_value is not None:
            bearing_phase = (
                bearing_phase_value
                if isinstance(bearing_phase_value, Mapping)
                else {}
            )
            harmonic_generation = str(
                bearing_phase.get("force_harmonic_generation", "")
            )
            angle_generation = str(
                bearing_phase.get("rotor_angle_generation", "")
            )
            angle_digest = str(
                bearing_phase.get("rotor_angle_sha256", "")
            ).lower()
            try:
                rotor_origin = float(bearing_phase.get("rotor_phase_origin_deg"))
                harmonic_origin = float(
                    bearing_phase.get("force_harmonic_phase_origin_deg")
                )
                slot_pitch = float(bearing_phase.get("slot_pitch_deg"))
                harmonic_order = int(bearing_phase.get("force_harmonic_order"))
            except (TypeError, ValueError):
                rotor_origin = harmonic_origin = slot_pitch = math.nan
                harmonic_order = -1
            slot_count = (
                360.0 / slot_pitch
                if math.isfinite(slot_pitch) and slot_pitch > 0.0
                else math.nan
            )
            phase_delta = (
                math.remainder(harmonic_origin - rotor_origin, 360.0)
                if math.isfinite(rotor_origin) and math.isfinite(harmonic_origin)
                else math.nan
            )
            magnetic_bearing_force_harmonic_phase_origin_identity_ok = (
                bool(harmonic_generation)
                and bearing_phase.get("force_sample_harmonic_generation")
                == harmonic_generation
                and bool(angle_generation)
                and bearing_phase.get("force_sample_rotor_angle_generation")
                == angle_generation
                and all(
                    math.isfinite(value)
                    for value in (rotor_origin, harmonic_origin, slot_pitch)
                )
                and slot_pitch > 0.0
                and math.isfinite(slot_count)
                and math.isclose(
                    slot_count, round(slot_count), rel_tol=0.0, abs_tol=1.0e-12
                )
                and type(bearing_phase.get("force_harmonic_order")) is int
                and harmonic_order > 0
                and bearing_phase.get("force_sample_harmonic_order")
                == harmonic_order
                and math.isclose(
                    phase_delta, 0.0, rel_tol=0.0, abs_tol=1.0e-12
                )
                and bearing_phase.get("rotor_angle_basis") == "mechanical_deg"
                and bearing_phase.get("force_harmonic_angle_basis")
                == bearing_phase.get("rotor_angle_basis")
                and bearing_phase.get("fourier_sign_convention")
                == "exp(-j*n*theta)"
                and bearing_phase.get("force_harmonic_fourier_sign_convention")
                == bearing_phase.get("fourier_sign_convention")
                and bearing_phase.get("phase_origin_convention") == "rotor_d_axis"
                and bearing_phase.get("force_harmonic_phase_origin_convention")
                == bearing_phase.get("phase_origin_convention")
                and len(angle_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in angle_digest
                )
                and str(
                    bearing_phase.get("harmonic_input_angle_sha256", "")
                ).lower()
                == angle_digest
            )

        near_singular_value = identity_value.get(
            "bem_near_singular_panel_subdivision_quadrature_generation_identity"
        )
        if near_singular_value is not None:
            near_singular = (
                near_singular_value
                if isinstance(near_singular_value, Mapping)
                else {}
            )
            surface_generation = str(
                near_singular.get("surface_generation", "")
            )
            subdivision_generation = str(
                near_singular.get("subdivision_generation", "")
            )
            quadrature_generation = str(
                near_singular.get("quadrature_generation", "")
            )
            interaction_ids = near_singular.get("interaction_ids")
            quadrature_orders = near_singular.get("quadrature_orders")
            subdivision_digest = str(
                near_singular.get("subdivision_map_sha256", "")
            ).lower()
            bem_near_singular_panel_subdivision_quadrature_generation_identity_ok = (
                bool(surface_generation)
                and near_singular.get("near_singular_interaction_surface_generation")
                == surface_generation
                and near_singular.get("panel_subdivision_surface_generation")
                == surface_generation
                and bool(subdivision_generation)
                and near_singular.get("quadrature_subdivision_generation")
                == subdivision_generation
                and bool(quadrature_generation)
                and near_singular.get("interaction_quadrature_generation")
                == quadrature_generation
                and isinstance(interaction_ids, list)
                and bool(interaction_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in interaction_ids
                )
                and len(set(interaction_ids)) == len(interaction_ids)
                and near_singular.get("subdivided_interaction_ids")
                == interaction_ids
                and isinstance(quadrature_orders, list)
                and len(quadrature_orders) == len(interaction_ids)
                and all(
                    isinstance(value, int)
                    and not isinstance(value, bool)
                    and value > 0
                    for value in quadrature_orders
                )
                and near_singular.get("applied_quadrature_orders")
                == quadrature_orders
                and len(subdivision_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in subdivision_digest
                )
                and str(
                    near_singular.get(
                        "quadrature_input_subdivision_map_sha256", ""
                    )
                ).lower()
                == subdivision_digest
            )

        maglev_coil_value = identity_value.get(
            "maglev_force_coil_polarity_orientation_generation_identity"
        )
        if maglev_coil_value is not None:
            maglev_coil = (
                maglev_coil_value
                if isinstance(maglev_coil_value, Mapping)
                else {}
            )
            force_generation = str(maglev_coil.get("force_generation", ""))
            coil_generation = str(maglev_coil.get("coil_generation", ""))
            coil_ids = maglev_coil.get("coil_ids")
            current_polarities = maglev_coil.get("current_polarities")
            winding_orientations = maglev_coil.get("winding_orientations")
            orientation_digest = str(
                maglev_coil.get("coil_orientation_map_sha256", "")
            ).lower()
            maglev_force_coil_polarity_orientation_generation_identity_ok = (
                bool(force_generation)
                and maglev_coil.get("coil_force_generation") == force_generation
                and bool(coil_generation)
                and maglev_coil.get("current_polarity_coil_generation")
                == coil_generation
                and maglev_coil.get("winding_orientation_coil_generation")
                == coil_generation
                and maglev_coil.get("force_result_coil_generation")
                == coil_generation
                and isinstance(coil_ids, list)
                and bool(coil_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in coil_ids
                )
                and len(set(coil_ids)) == len(coil_ids)
                and maglev_coil.get("force_coil_ids") == coil_ids
                and isinstance(current_polarities, list)
                and len(current_polarities) == len(coil_ids)
                and all(value in (-1, 1) for value in current_polarities)
                and maglev_coil.get("force_current_polarities")
                == current_polarities
                and isinstance(winding_orientations, list)
                and len(winding_orientations) == len(coil_ids)
                and all(
                    value in {"clockwise", "counterclockwise"}
                    for value in winding_orientations
                )
                and maglev_coil.get("force_winding_orientations")
                == winding_orientations
                and len(orientation_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in orientation_digest
                )
                and str(
                    maglev_coil.get("force_coil_orientation_map_sha256", "")
                ).lower()
                == orientation_digest
            )

        self_term_generation_value = identity_value.get(
            "bem_demag_self_term_solid_angle_orientation_generation_identity"
        )
        if self_term_generation_value is not None:
            self_term_generation = (
                self_term_generation_value
                if isinstance(self_term_generation_value, Mapping)
                else {}
            )
            operator_generation = str(
                self_term_generation.get("operator_generation", "")
            )
            boundary_generation = str(
                self_term_generation.get("boundary_generation", "")
            )
            panel_ids = self_term_generation.get("panel_ids")
            orientation_signs = self_term_generation.get(
                "panel_orientation_signs"
            )
            solid_angles = self_term_generation.get("solid_angles_sr")
            input_digest = str(
                self_term_generation.get("self_term_input_sha256", "")
            ).lower()
            bem_demag_self_term_solid_angle_orientation_generation_identity_ok = (
                bool(operator_generation)
                and self_term_generation.get("result_operator_generation")
                == operator_generation
                and bool(boundary_generation)
                and self_term_generation.get("self_term_boundary_generation")
                == boundary_generation
                and self_term_generation.get(
                    "panel_orientation_boundary_generation"
                )
                == boundary_generation
                and self_term_generation.get("operator_boundary_generation")
                == boundary_generation
                and isinstance(panel_ids, list)
                and bool(panel_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in panel_ids
                )
                and len(set(panel_ids)) == len(panel_ids)
                and self_term_generation.get("self_term_panel_ids") == panel_ids
                and isinstance(orientation_signs, list)
                and len(orientation_signs) == len(panel_ids)
                and all(value in (-1, 1) for value in orientation_signs)
                and self_term_generation.get("self_term_orientation_signs")
                == orientation_signs
                and isinstance(solid_angles, list)
                and len(solid_angles) == len(panel_ids)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in solid_angles
                )
                and self_term_generation.get(
                    "applied_self_term_solid_angles_sr"
                )
                == solid_angles
                and self_term_generation.get("orientation_convention")
                == "outward_positive"
                and self_term_generation.get("self_term_orientation_convention")
                == self_term_generation.get("orientation_convention")
                and len(input_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in input_digest
                )
                and str(
                    self_term_generation.get(
                        "assembled_self_term_input_sha256", ""
                    )
                ).lower()
                == input_digest
            )

        virtual_work_value = identity_value.get(
            "virtual_work_force_displacement_coordinate_geometry_generation_identity"
        )
        if virtual_work_value is not None:
            virtual_work = (
                virtual_work_value
                if isinstance(virtual_work_value, Mapping)
                else {}
            )
            force_generation = str(virtual_work.get("force_generation", ""))
            geometry_generation = str(
                virtual_work.get("geometry_generation", "")
            )
            coordinate_generation = str(
                virtual_work.get("displacement_coordinate_generation", "")
            )
            displacements = virtual_work.get("displacements_m")
            energy_samples = virtual_work.get("energy_samples_j")
            geometry_generations = virtual_work.get(
                "energy_sample_geometry_generations"
            )
            energy_digest = str(
                virtual_work.get("energy_sample_table_sha256", "")
            ).lower()
            virtual_work_displacement_coordinate_geometry_generation_identity_ok = (
                bool(force_generation)
                and virtual_work.get("result_force_generation") == force_generation
                and bool(geometry_generation)
                and virtual_work.get(
                    "displacement_coordinate_geometry_generation"
                )
                == geometry_generation
                and virtual_work.get("force_geometry_generation")
                == geometry_generation
                and isinstance(geometry_generations, list)
                and bool(geometry_generations)
                and all(
                    generation == geometry_generation
                    for generation in geometry_generations
                )
                and bool(coordinate_generation)
                and virtual_work.get(
                    "energy_sample_displacement_coordinate_generation"
                )
                == coordinate_generation
                and virtual_work.get(
                    "force_displacement_coordinate_generation"
                )
                == coordinate_generation
                and virtual_work.get("displacement_axis") == "global-z"
                and virtual_work.get("force_component_axis")
                == virtual_work.get("displacement_axis")
                and virtual_work.get("displacement_unit") == "m"
                and virtual_work.get("energy_unit") == "J"
                and virtual_work.get("force_unit") == "N"
                and isinstance(displacements, list)
                and len(displacements) >= 3
                and len(displacements) == len(geometry_generations)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in displacements
                )
                and all(
                    right > left
                    for left, right in zip(displacements, displacements[1:])
                )
                and virtual_work.get("force_displacements_m") == displacements
                and isinstance(energy_samples, list)
                and len(energy_samples) == len(displacements)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in energy_samples
                )
                and virtual_work.get("force_energy_samples_j") == energy_samples
                and len(energy_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in energy_digest
                )
                and str(
                    virtual_work.get(
                        "force_input_energy_sample_table_sha256", ""
                    )
                ).lower()
                == energy_digest
            )

        demag_energy_value = identity_value.get(
            "demag_energy_surface_charge_normal_quadrature_generation_identity"
        )
        if demag_energy_value is not None:
            demag_energy = (
                demag_energy_value
                if isinstance(demag_energy_value, Mapping)
                else {}
            )
            energy_generation = str(
                demag_energy.get("energy_generation", "")
            ).strip()
            boundary_generation = str(
                demag_energy.get("boundary_generation", "")
            ).strip()
            panel_ids = demag_energy.get("panel_ids")
            charges = demag_energy.get("surface_charges_a_per_m")
            normal_digests = demag_energy.get("outward_normal_sha256")
            weights = demag_energy.get("panel_quadrature_weights")
            digest = str(
                demag_energy.get("demag_energy_input_sha256", "")
            ).lower()
            demag_energy_surface_charge_normal_quadrature_identity_ok = (
                bool(energy_generation)
                and demag_energy.get("result_energy_generation")
                == energy_generation
                and bool(boundary_generation)
                and all(
                    demag_energy.get(key) == boundary_generation
                    for key in (
                        "surface_charge_boundary_generation",
                        "normal_boundary_generation",
                        "quadrature_boundary_generation",
                        "energy_boundary_generation",
                    )
                )
                and isinstance(panel_ids, list)
                and bool(panel_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in panel_ids
                )
                and len(set(panel_ids)) == len(panel_ids)
                and demag_energy.get("energy_panel_ids") == panel_ids
                and isinstance(charges, list)
                and len(charges) == len(panel_ids)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in charges
                )
                and demag_energy.get("energy_surface_charges_a_per_m") == charges
                and isinstance(normal_digests, list)
                and len(normal_digests) == len(panel_ids)
                and all(
                    isinstance(value, str)
                    and len(value) == 64
                    and all(character in "0123456789abcdef" for character in value)
                    for value in normal_digests
                )
                and demag_energy.get("energy_outward_normal_sha256")
                == normal_digests
                and isinstance(weights, list)
                and len(weights) == len(panel_ids)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    and float(value) >= 0.0
                    for value in weights
                )
                and demag_energy.get("energy_panel_quadrature_weights") == weights
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    demag_energy.get("assembled_demag_energy_input_sha256", "")
                ).lower()
                == digest
            )

        stiffness_value = identity_value.get(
            "maglev_stiffness_displacement_equilibrium_force_generation_identity"
        )
        if stiffness_value is not None:
            stiffness = stiffness_value if isinstance(stiffness_value, Mapping) else {}
            stiffness_generation = str(
                stiffness.get("stiffness_generation", "")
            ).strip()
            geometry_generation = str(
                stiffness.get("geometry_generation", "")
            ).strip()
            coordinate_generation = str(
                stiffness.get("coordinate_generation", "")
            ).strip()
            sample_geometry = stiffness.get("sample_geometry_generations")
            force_geometry = stiffness.get("force_geometry_generations")
            displacements = stiffness.get("displacements_m")
            forces = stiffness.get("force_samples_n")
            digest = str(
                stiffness.get("stiffness_sample_table_sha256", "")
            ).lower()
            try:
                equilibrium_index = int(stiffness.get("equilibrium_sample_index"))
                result_equilibrium_index = int(
                    stiffness.get("stiffness_equilibrium_sample_index")
                )
            except (TypeError, ValueError):
                equilibrium_index = result_equilibrium_index = -1
            maglev_stiffness_displacement_equilibrium_force_identity_ok = (
                bool(stiffness_generation)
                and stiffness.get("result_stiffness_generation")
                == stiffness_generation
                and bool(geometry_generation)
                and isinstance(sample_geometry, list)
                and len(sample_geometry) >= 3
                and all(value == geometry_generation for value in sample_geometry)
                and force_geometry == sample_geometry
                and bool(coordinate_generation)
                and stiffness.get("displacement_coordinate_generation")
                == coordinate_generation
                and stiffness.get("force_coordinate_generation")
                == coordinate_generation
                and stiffness.get("displacement_axis") in {
                    "global-x",
                    "global-y",
                    "global-z",
                }
                and stiffness.get("force_component_axis")
                == stiffness.get("displacement_axis")
                and isinstance(displacements, list)
                and len(displacements) == len(sample_geometry)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in displacements
                )
                and all(
                    right > left
                    for left, right in zip(displacements, displacements[1:])
                )
                and stiffness.get("force_displacements_m") == displacements
                and isinstance(forces, list)
                and len(forces) == len(displacements)
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in forces
                )
                and stiffness.get("stiffness_force_samples_n") == forces
                and 0 <= equilibrium_index < len(displacements)
                and result_equilibrium_index == equilibrium_index
                and abs(float(displacements[equilibrium_index])) <= 1.0e-15
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    stiffness.get("result_stiffness_sample_table_sha256", "")
                ).lower()
                == digest
            )

        near_singular_value = identity_value.get(
            "bem_near_singular_distance_panel_quadrature_generation_identity"
        )
        if near_singular_value is not None:
            near_singular = (
                near_singular_value
                if isinstance(near_singular_value, Mapping)
                else {}
            )
            interaction_generation = str(
                near_singular.get("interaction_generation", "")
            ).strip()
            geometry_generation = str(
                near_singular.get("geometry_generation", "")
            ).strip()
            target_id = near_singular.get("target_point_id")
            panel_id = near_singular.get("source_panel_id")
            panel_digest = str(
                near_singular.get("source_panel_geometry_sha256", "")
            ).lower()
            interaction_digest = str(
                near_singular.get("near_singular_interaction_sha256", "")
            ).lower()
            try:
                distance = float(near_singular.get("target_distance_m"))
                quadrature_distance = float(
                    near_singular.get("quadrature_target_distance_m")
                )
            except (TypeError, ValueError):
                distance = quadrature_distance = math.nan
            quadrature_order = near_singular.get("adaptive_quadrature_order")
            bem_near_singular_distance_panel_quadrature_identity_ok = (
                bool(interaction_generation)
                and near_singular.get("result_interaction_generation")
                == interaction_generation
                and bool(geometry_generation)
                and all(
                    near_singular.get(key) == geometry_generation
                    for key in (
                        "target_distance_geometry_generation",
                        "source_panel_geometry_generation",
                        "adaptive_quadrature_geometry_generation",
                    )
                )
                and isinstance(target_id, int)
                and not isinstance(target_id, bool)
                and target_id > 0
                and near_singular.get("distance_target_point_id") == target_id
                and isinstance(panel_id, int)
                and not isinstance(panel_id, bool)
                and panel_id > 0
                and near_singular.get("distance_source_panel_id") == panel_id
                and math.isfinite(distance)
                and distance > 0.0
                and quadrature_distance == distance
                and len(panel_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in panel_digest
                )
                and str(
                    near_singular.get(
                        "quadrature_source_panel_geometry_sha256", ""
                    )
                ).lower()
                == panel_digest
                and isinstance(quadrature_order, int)
                and not isinstance(quadrature_order, bool)
                and quadrature_order > 0
                and near_singular.get("evaluated_quadrature_order")
                == quadrature_order
                and len(interaction_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in interaction_digest
                )
                and str(
                    near_singular.get(
                        "assembled_near_singular_interaction_sha256", ""
                    )
                ).lower()
                == interaction_digest
            )

        moving_force_value = identity_value.get(
            "moving_magnet_force_position_orientation_equilibrium_generation_identity"
        )
        if moving_force_value is not None:
            moving_force = (
                moving_force_value
                if isinstance(moving_force_value, Mapping)
                else {}
            )
            force_generation = str(
                moving_force.get("force_generation", "")
            ).strip()
            geometry_generation = str(
                moving_force.get("geometry_generation", "")
            ).strip()
            sample_generations = moving_force.get(
                "position_sample_geometry_generations"
            )
            moving_positions = moving_force.get("position_samples_m")
            orientation = moving_force.get("magnet_orientation_quaternion")
            moving_force_samples = moving_force.get("force_samples_n")
            equilibrium_index = moving_force.get("equilibrium_sample_index")
            digest = str(
                moving_force.get("moving_force_sample_table_sha256", "")
            ).lower()

            def finite_vectors(value: object, count: int) -> bool:
                return (
                    isinstance(value, list)
                    and len(value) == count
                    and all(
                        isinstance(row, list)
                        and len(row) == 3
                        and all(
                            isinstance(component, (int, float))
                            and not isinstance(component, bool)
                            and math.isfinite(float(component))
                            for component in row
                        )
                        for row in value
                    )
                )

            positions_ok = finite_vectors(
                moving_positions, len(sample_generations)
            ) if isinstance(sample_generations, list) else False
            force_samples_ok = finite_vectors(
                moving_force_samples, len(sample_generations)
            ) if isinstance(sample_generations, list) else False
            orientation_ok = (
                isinstance(orientation, list)
                and len(orientation) == 4
                and all(
                    isinstance(value, (int, float))
                    and not isinstance(value, bool)
                    and math.isfinite(float(value))
                    for value in orientation
                )
                and math.isclose(
                    sum(float(value) ** 2 for value in orientation),
                    1.0,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                )
            )
            moving_magnet_force_position_orientation_equilibrium_identity_ok = (
                bool(force_generation)
                and moving_force.get("result_force_generation") == force_generation
                and bool(geometry_generation)
                and isinstance(sample_generations, list)
                and len(sample_generations) >= 3
                and all(
                    value == geometry_generation for value in sample_generations
                )
                and all(
                    moving_force.get(key) == geometry_generation
                    for key in (
                        "orientation_geometry_generation",
                        "equilibrium_geometry_generation",
                        "force_geometry_generation",
                    )
                )
                and positions_ok
                and moving_force.get("force_position_samples_m")
                == moving_positions
                and orientation_ok
                and moving_force.get("force_orientation_quaternion") == orientation
                and force_samples_ok
                and moving_force.get("differentiated_force_samples_n")
                == moving_force_samples
                and isinstance(equilibrium_index, int)
                and not isinstance(equilibrium_index, bool)
                and 0 <= equilibrium_index < len(sample_generations)
                and moving_force.get("force_equilibrium_sample_index")
                == equilibrium_index
                and sum(
                    float(value) ** 2
                    for value in moving_positions[equilibrium_index]
                )
                <= 1.0e-30
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    moving_force.get("result_moving_force_sample_table_sha256", "")
                ).lower()
                == digest
            )

        dual_lane_value = identity_value.get(
            "motor_force_dual_lane_interface_flux_coenergy_generation_identity"
        )
        if dual_lane_value is not None:
            dual_lane = (
                dual_lane_value if isinstance(dual_lane_value, Mapping) else {}
            )
            generation = str(
                dual_lane.get("comparison_generation", "")
            ).strip()
            lane_ids = dual_lane.get("lane_ids")
            try:
                interface_flux = [
                    float(value)
                    for value in dual_lane.get("interface_normal_flux_wb", [])
                ]
                result_interface_flux = [
                    float(value)
                    for value in dual_lane.get(
                        "result_interface_normal_flux_wb", []
                    )
                ]
                coenergy = [
                    float(value) for value in dual_lane.get("coenergy_j", [])
                ]
                result_coenergy = [
                    float(value)
                    for value in dual_lane.get("result_coenergy_j", [])
                ]
                force = [float(value) for value in dual_lane.get("force_n", [])]
                result_force = [
                    float(value) for value in dual_lane.get("result_force_n", [])
                ]
            except (TypeError, ValueError):
                interface_flux = result_interface_flux = []
                coenergy = result_coenergy = []
                force = result_force = []
            mesh_digest = str(
                dual_lane.get("coupling_mesh_sha256", "")
            ).lower()
            operator_digest = str(
                dual_lane.get("mixed_operator_contract_sha256", "")
            ).lower()
            canonical_lanes = ["ngsolve_age", "hdiv_mmm_hcurl_eddy_bubble"]
            motor_force_dual_lane_generation_identity_ok = (
                bool(generation)
                and all(
                    dual_lane.get(key) == generation
                    for key in (
                        "lane_policy_comparison_generation",
                        "interface_flux_comparison_generation",
                        "coenergy_comparison_generation",
                        "force_comparison_generation",
                        "coupling_mesh_comparison_generation",
                    )
                )
                and lane_ids == canonical_lanes
                and dual_lane.get("result_lane_ids") == canonical_lanes
                and len(interface_flux) == len(canonical_lanes)
                and all(math.isfinite(value) for value in interface_flux)
                and result_interface_flux == interface_flux
                and len(coenergy) == len(canonical_lanes)
                and all(math.isfinite(value) for value in coenergy)
                and result_coenergy == coenergy
                and len(force) == len(canonical_lanes)
                and all(math.isfinite(value) for value in force)
                and result_force == force
                and len(mesh_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in mesh_digest
                )
                and str(dual_lane.get("result_coupling_mesh_sha256", "")).lower()
                == mesh_digest
                and len(operator_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in operator_digest
                )
                and str(
                    dual_lane.get("result_mixed_operator_contract_sha256", "")
                ).lower()
                == operator_digest
            )

        linear_motion_value = identity_value.get(
            "linear_motor_end_effect_translation_position_symmetry_generation_identity"
        )
        if linear_motion_value is not None:
            linear_motion = (
                linear_motion_value
                if isinstance(linear_motion_value, Mapping)
                else {}
            )
            generation = str(
                linear_motion.get("sweep_generation", "")
            ).strip()
            frame = str(linear_motion.get("translation_frame", "")).strip()
            symmetry_factor = linear_motion.get("symmetry_factor")
            try:
                positions_m = [
                    float(value)
                    for value in linear_motion.get("mover_positions_m", [])
                ]
                result_positions_m = [
                    float(value)
                    for value in linear_motion.get("result_mover_positions_m", [])
                ]
                window_m = [
                    float(value)
                    for value in linear_motion.get("end_effect_window_m", [])
                ]
                result_window_m = [
                    float(value)
                    for value in linear_motion.get("result_end_effect_window_m", [])
                ]
                thrust_n = [
                    float(value) for value in linear_motion.get("thrust_n", [])
                ]
                result_thrust_n = [
                    float(value)
                    for value in linear_motion.get("result_thrust_n", [])
                ]
                stiffness = [
                    float(value)
                    for value in linear_motion.get("stiffness_n_per_m", [])
                ]
                result_stiffness = [
                    float(value)
                    for value in linear_motion.get(
                        "result_stiffness_n_per_m", []
                    )
                ]
            except (TypeError, ValueError):
                positions_m = result_positions_m = []
                window_m = result_window_m = []
                thrust_n = result_thrust_n = []
                stiffness = result_stiffness = []
            digest = str(
                linear_motion.get("linear_motion_table_sha256", "")
            ).lower()
            linear_motor_end_effect_generation_identity_ok = (
                bool(generation)
                and all(
                    linear_motion.get(key) == generation
                    for key in (
                        "position_sweep_generation",
                        "end_effect_sweep_generation",
                        "translation_frame_sweep_generation",
                        "symmetry_sweep_generation",
                        "force_result_sweep_generation",
                    )
                )
                and len(positions_m) >= 3
                and all(math.isfinite(value) for value in positions_m)
                and all(
                    right > left for left, right in zip(positions_m, positions_m[1:])
                )
                and result_positions_m == positions_m
                and len(window_m) == 2
                and all(math.isfinite(value) for value in window_m)
                and window_m[0] < min(positions_m)
                and window_m[1] > max(positions_m)
                and result_window_m == window_m
                and bool(frame)
                and linear_motion.get("result_translation_frame") == frame
                and isinstance(symmetry_factor, int)
                and not isinstance(symmetry_factor, bool)
                and symmetry_factor > 0
                and linear_motion.get("result_symmetry_factor") == symmetry_factor
                and len(thrust_n) == len(positions_m)
                and all(math.isfinite(value) for value in thrust_n)
                and result_thrust_n == thrust_n
                and len(stiffness) == len(positions_m)
                and all(math.isfinite(value) for value in stiffness)
                and result_stiffness == stiffness
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    linear_motion.get("result_linear_motion_table_sha256", "")
                ).lower()
                == digest
            )

        bem_panel_demag_force_generation_identity_ok = (
            _bem_panel_demag_force_identity_ok(
                identity_value.get(
                    "bem_panel_normal_material_region_demag_force_generation_identity"
                )
            )
        )
        motor_harmonic_force_generation_identity_ok = (
            _motor_harmonic_force_identity_ok(
                identity_value.get(
                    "motor_harmonic_rotor_angle_current_phase_force_frame_generation_identity"
                )
            )
        )
        maglev_force_energy_generation_identity_ok = (
            _maglev_force_energy_identity_ok(
                identity_value.get(
                    "maglev_force_stiffness_equilibrium_energy_finite_difference_generation_identity"
                )
            )
        )
        motor_dual_lane_alignment_generation_identity_ok = (
            _motor_dual_lane_alignment_identity_ok(
                identity_value.get(
                    "motor_dual_lane_geometry_excitation_force_frame_harmonic_alignment_generation_identity"
                )
            )
        )

    method_difference = _maximum_relative_difference(target, stress)
    independent_stress_difference = _maximum_relative_difference(stress, independent_stress)
    selection_differences = [
        abs(all_value - target_value)
        / max(abs(all_value), abs(target_value), 1.0e-300)
        for all_value, target_value in zip(all_body, target)
    ]
    all_body_ratios = [
        abs(all_value) / max(abs(target_value), 1.0e-300)
        for all_value, target_value in zip(all_body, target)
    ]
    target_integral = _trapezoid_integral(positions, target)
    stress_integral = _trapezoid_integral(positions, stress)
    work_difference = abs(target_integral - stress_integral) / max(
        abs(target_integral), abs(stress_integral), 1.0e-300
    )
    same_sign = all(a * b > 0.0 for a, b in zip(target, stress))
    same_trend = all(
        (right_a - left_a) * (right_b - left_b) >= 0.0
        for left_a, right_a, left_b, right_b in zip(
            target, target[1:], stress, stress[1:]
        )
    )

    replay = summary.get("replay")
    replay = replay if isinstance(replay, Mapping) else {}
    position_unit = str(summary.get("position_unit") or "")
    artifact_position_units = summary.get("artifact_position_units")
    artifact_position_units = (
        artifact_position_units if isinstance(artifact_position_units, Mapping) else None
    )
    artifact_units_match = artifact_position_units is None or (
        bool(artifact_position_units)
        and all(str(unit) == position_unit for unit in artifact_position_units.values())
    )
    parsed_replay = float(replay.get("parsed_max_abs", math.inf))
    checks = {
        "sample_count_sufficient": len(positions) >= minimum_sample_count,
        "positions_strictly_increase": increasing_positions,
        "dimension_and_force_unit_consistent": (
            str(summary.get("quantity_dimension") or ""),
            str(summary.get("force_unit") or ""),
        )
        in _DIMENSION_UNITS,
        "position_unit_recorded": position_unit in {"m", "mm"},
        "artifact_position_units_match_common_grid": artifact_units_match,
        "comparison_axis_recorded": summary.get("comparison_axis") in {"x", "y", "z"},
        "target_methods_share_sign": same_sign,
        "target_methods_share_stepwise_trend": same_trend,
        "target_method_closure_within_tolerance": method_difference
        <= tolerances["maximum_method_relative_difference"],
        "independent_stress_replay_within_tolerance": independent_stress_difference
        <= tolerances["maximum_independent_stress_relative_difference"],
        "selection_scope_is_materially_distinct": min(selection_differences)
        >= tolerances["minimum_selection_scope_relative_difference"],
        "all_body_control_is_not_target_body_force": max(all_body_ratios)
        <= tolerances["maximum_all_body_to_target_magnitude_ratio"],
        "force_position_integrals_close": work_difference
        <= tolerances["maximum_work_relative_difference"],
        "parsed_replay_is_exact_enough": math.isfinite(parsed_replay)
        and 0.0 <= parsed_replay
        <= tolerances["maximum_parsed_replay_absolute_difference"],
        "binary_nonlog_outputs_replay_exact": replay.get(
            "binary_nonlog_outputs_exact"
        )
        is True,
        "position_force_samples_share_one_sweep_generation": one_sweep_generation_ok,
        "demag_reference_matches_current_geometry_revision": demag_reference_ok,
        "force_and_demag_share_transformed_coordinate_system": (
            coordinate_system_binding_ok
        ),
        "force_profiles_share_total_3d_normalization": force_normalization_ok,
        "hysteresis_observable_state_and_tangent_share_branch": (
            hysteresis_branch_state_ok
        ),
        "remanence_vector_is_transformed_from_material_to_assembly_frame": (
            remanence_frame_binding_ok
        ),
        "force_surface_encloses_only_target_body": force_surface_body_ownership_ok,
        "demag_operating_point_uses_active_branch_interpolation": (
            demag_branch_interpolation_ok
        ),
        "linear_motor_thrust_uses_winding_phase_and_electrical_angle_direction": (
            linear_motor_thrust_phase_identity_ok
        ),
        "demag_recoil_line_matches_evaluation_temperature_generation": (
            demag_recoil_temperature_identity_ok
        ),
        "bem_demag_normals_match_current_surface_mesh_generation": (
            bem_demag_surface_normal_generation_identity_ok
        ),
        "cogging_torque_symmetry_multiplier_matches_periodic_sector": (
            cogging_torque_periodic_sector_symmetry_identity_ok
        ),
        "bem_self_term_solid_angle_matches_panel_orientation_generation": (
            bem_self_term_solid_angle_orientation_identity_ok
        ),
        "demag_energy_force_derivative_uses_one_displacement_length_unit": (
            demag_energy_force_displacement_length_unit_identity_ok
        ),
        "bem_near_singular_quadrature_uses_current_target_pair_scale": (
            bem_near_singular_quadrature_target_scale_identity_ok
        ),
        "force_and_torque_share_reference_origin_length_unit": (
            force_torque_reference_origin_length_unit_identity_ok
        ),
        "bem_self_term_uses_current_normalized_surface_winding": (
            bem_solid_angle_surface_winding_identity_ok
        ),
        "maglev_stiffness_uses_one_force_displacement_perturbation_generation": (
            maglev_stiffness_force_displacement_generation_identity_ok
        ),
        "bem_demag_tensor_uses_current_body_placement_coordinate_basis": (
            bem_demag_tensor_coordinate_basis_generation_identity_ok
        ),
        "magnetic_bearing_force_harmonics_share_rotor_phase_origin": (
            magnetic_bearing_force_harmonic_phase_origin_identity_ok
        ),
        "bem_near_singular_quadrature_uses_current_panel_subdivision": (
            bem_near_singular_panel_subdivision_quadrature_generation_identity_ok
        ),
        "maglev_force_uses_current_coil_polarity_and_winding_orientation": (
            maglev_force_coil_polarity_orientation_generation_identity_ok
        ),
        "bem_demag_self_term_uses_current_boundary_orientation_generation": (
            bem_demag_self_term_solid_angle_orientation_generation_identity_ok
        ),
        "virtual_work_force_uses_current_displacement_geometry_generation": (
            virtual_work_displacement_coordinate_geometry_generation_identity_ok
        ),
        "demag_energy_uses_current_surface_charge_normal_and_quadrature": (
            demag_energy_surface_charge_normal_quadrature_identity_ok
        ),
        "maglev_stiffness_uses_aligned_displacement_equilibrium_force_states": (
            maglev_stiffness_displacement_equilibrium_force_identity_ok
        ),
        "bem_near_singular_quadrature_uses_current_distance_and_panel_geometry": (
            bem_near_singular_distance_panel_quadrature_identity_ok
        ),
        "moving_magnet_force_uses_current_position_orientation_and_equilibrium": (
            moving_magnet_force_position_orientation_equilibrium_identity_ok
        ),
        "motor_force_comparison_uses_age_and_hdiv_mmm_hcurl_eddy_bubble_lanes": (
            motor_force_dual_lane_generation_identity_ok
        ),
        "linear_motor_force_uses_current_position_end_effect_frame_and_symmetry": (
            linear_motor_end_effect_generation_identity_ok
        ),
        "bem_demag_force_uses_current_panels_normals_materials_and_results": (
            bem_panel_demag_force_generation_identity_ok
        ),
        "motor_force_harmonics_use_current_angles_phases_bins_and_frame": (
            motor_harmonic_force_generation_identity_ok
        ),
        "maglev_force_stiffness_energy_and_equilibrium_share_current_state": (
            maglev_force_energy_generation_identity_ok
        ),
        "motor_dual_lanes_share_geometry_excitation_frame_harmonics_and_angles": (
            motor_dual_lane_alignment_generation_identity_ok
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "magnetic_force_method_profile_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "warnings": [] if identity_present else ["artifact_identity_not_recorded"],
        "metrics": {
            "sample_count": len(positions),
            "maximum_method_relative_difference": method_difference,
            "maximum_independent_stress_relative_difference": independent_stress_difference,
            "minimum_selection_scope_relative_difference": min(selection_differences),
            "maximum_all_body_to_target_magnitude_ratio": max(all_body_ratios),
            "target_force_position_integral": target_integral,
            "stress_force_position_integral": stress_integral,
            "force_position_integral_relative_difference": work_difference,
            "parsed_replay_maximum_absolute_difference": parsed_replay,
        },
        "tolerances": tolerances,
        "lesson": (
            "Force-method closure is meaningful only when the target body, closed "
            "stress surface, comparison axis, units, and dimensional convention are "
            "pinned. Bind every artifact's position grid to the same declared unit; equal "
            "numeric coordinates do not prove equal physical positions. Keep an all-body "
            "force as a negative control so a selection-scope "
            "error cannot masquerade as disagreement between force formulations. Bind "
            "hysteresis observables to their active branch and transform material-frame "
            "remanence vectors into the assembly frame before reuse."
        ),
    }
