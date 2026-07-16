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
            bem_demag_tensor_coordinate_basis_generation_identity_ok = (
                bool(placement_generation)
                and demag_tensor.get("surface_mesh_body_placement_generation")
                == placement_generation
                and bool(tensor_generation)
                and demag_tensor.get("demag_tensor_body_placement_generation")
                == placement_generation
                and demag_tensor.get("body_coordinate_basis")
                == "body-local-current"
                and demag_tensor.get("tensor_coordinate_basis")
                == demag_tensor.get("body_coordinate_basis")
                and demag_tensor.get("body_basis_handedness") == "right_handed"
                and demag_tensor.get("tensor_basis_handedness")
                == demag_tensor.get("body_basis_handedness")
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
            except (TypeError, ValueError):
                rotor_origin = harmonic_origin = slot_pitch = math.nan
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
                and math.isclose(
                    harmonic_origin, rotor_origin, rel_tol=0.0, abs_tol=1.0e-12
                )
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
