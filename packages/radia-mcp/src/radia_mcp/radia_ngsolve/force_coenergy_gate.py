"""Virtual-work consistency gate for displacement-force sweeps."""
from __future__ import annotations

import math


def force_coenergy_displacement_gate(
    positions_m,
    coenergy_j,
    forces_along_displacement_n,
    *,
    energy_kind: str = "constant_current_coenergy",
    max_central_relative_error: float = 0.02,
    min_sample_count: int = 5,
    artifact_identity: dict | None = None,
):
    """Compare direct force with the central derivative of magnetic coenergy.

    The caller must project the direct force onto the increasing displacement
    coordinate before calling this gate.  Endpoints are reported using one-sided
    differences but are not part of the acceptance metric.
    """
    x = [float(value) for value in positions_m]
    w = [float(value) for value in coenergy_j]
    force = [float(value) for value in forces_along_displacement_n]
    if not (len(x) == len(w) == len(force)):
        raise ValueError("positions, coenergy, and force must have the same length")
    if min_sample_count < 5:
        raise ValueError("min_sample_count must be >= 5")
    if max_central_relative_error < 0.0:
        raise ValueError("max_central_relative_error must be >= 0")

    identity_present = isinstance(artifact_identity, dict)
    force_snapshot_ok = True
    mesh_family_ok = True
    displacement_unit_ok = True
    force_frame_ok = True
    force_normalization_ok = True
    force_body_selection_ok = True
    virtual_work_constraint_basis_ok = True
    eddy_loss_harmonic_basis_ok = True
    axisymmetric_force_measure_ok = True
    eddy_loss_material_frequency_ok = True
    weighted_stress_mask_mesh_identity_ok = True
    complex_current_phasor_basis_identity_ok = True
    axisymmetric_force_radius_jacobian_coordinate_identity_ok = True
    nonlinear_bh_interpolation_extrapolation_identity_ok = True
    weighted_stress_mask_material_interface_identity_ok = True
    axisymmetric_coil_voltage_measure_identity_ok = True
    nonlinear_energy_coenergy_bh_iteration_identity_ok = True
    virtual_displacement_force_geometry_field_identity_ok = True
    if artifact_identity is not None and not identity_present:
        force_snapshot_ok = False
        mesh_family_ok = False
        displacement_unit_ok = False
        force_frame_ok = False
        force_normalization_ok = False
        force_body_selection_ok = False
        virtual_work_constraint_basis_ok = False
        eddy_loss_harmonic_basis_ok = False
        axisymmetric_force_measure_ok = False
        eddy_loss_material_frequency_ok = False
        weighted_stress_mask_mesh_identity_ok = False
        complex_current_phasor_basis_identity_ok = False
        axisymmetric_force_radius_jacobian_coordinate_identity_ok = False
        nonlinear_bh_interpolation_extrapolation_identity_ok = False
        weighted_stress_mask_material_interface_identity_ok = False
        axisymmetric_coil_voltage_measure_identity_ok = False
        nonlinear_energy_coenergy_bh_iteration_identity_ok = False
        virtual_displacement_force_geometry_field_identity_ok = False
    elif identity_present:
        direct = artifact_identity.get("direct_force_snapshot")
        derivative = artifact_identity.get("coenergy_derivative_snapshot")
        if not isinstance(direct, dict) or not isinstance(derivative, dict):
            force_snapshot_ok = False
        else:
            direct_step = str(direct.get("load_step_id", ""))
            derivative_step = str(derivative.get("load_step_id", ""))
            try:
                direct_time = float(direct["time_s"])
                derivative_time = float(derivative["time_s"])
            except (KeyError, TypeError, ValueError):
                direct_time = math.nan
                derivative_time = math.nan
            force_snapshot_ok = (
                bool(direct_step)
                and direct_step == derivative_step
                and math.isfinite(direct_time)
                and math.isfinite(derivative_time)
                and direct_time == derivative_time
            )
        generations = artifact_identity.get("coenergy_mesh_family_generations")
        mesh_family_ok = (
            isinstance(generations, list)
            and len(generations) == len(x)
            and all(isinstance(value, str) and bool(value) for value in generations)
            and len(set(generations)) == 1
        )
        displacement_axis = artifact_identity.get("displacement_axis")
        if displacement_axis is not None:
            displacement_unit_ok = (
                isinstance(displacement_axis, dict)
                and displacement_axis.get("numeric_unit") == "m"
                and displacement_axis.get("derivative_unit") == "m"
                and displacement_axis.get("scale_to_si") == 1.0
            )
        force_frame = artifact_identity.get("force_frame")
        if force_frame is not None:
            direct_axis = (
                force_frame.get("direct_axis")
                if isinstance(force_frame, dict)
                else None
            )
            derivative_axis = (
                force_frame.get("derivative_axis") if isinstance(force_frame, dict) else None
            )
            axes_are_finite = (
                isinstance(direct_axis, list)
                and isinstance(derivative_axis, list)
                and len(direct_axis) == len(derivative_axis) == 3
                and all(
                    isinstance(value, (int, float)) and math.isfinite(float(value))
                    for value in direct_axis + derivative_axis
                )
            )
            force_frame_ok = (
                isinstance(force_frame, dict)
                and bool(force_frame.get("direct_frame_id"))
                and force_frame.get("direct_frame_id")
                == force_frame.get("derivative_frame_id")
                and axes_are_finite
                and [float(value) for value in direct_axis]
                == [float(value) for value in derivative_axis]
                and force_frame.get("reflection_applied") is True
            )
        normalization = artifact_identity.get("force_normalization")
        if normalization is not None:
            force_normalization_ok = (
                isinstance(normalization, dict)
                and normalization.get("formulation") == "axisymmetric"
                and normalization.get("solver_result_scope") == "total_3d_force"
                and normalization.get("reported_result_scope")
                == normalization.get("solver_result_scope")
                and normalization.get("revolution_factor_application_count") == 0
            )
        selection = artifact_identity.get("force_body_selection")
        if selection is not None:
            target_groups = (
                selection.get("target_group_ids")
                if isinstance(selection, dict)
                else None
            )
            selected_groups = (
                selection.get("weighted_stress_selected_group_ids")
                if isinstance(selection, dict)
                else None
            )
            roles = selection.get("material_roles") if isinstance(selection, dict) else None
            excluded_air = (
                selection.get("excluded_air_group_ids")
                if isinstance(selection, dict)
                else None
            )
            force_body_selection_ok = (
                isinstance(selection, dict)
                and target_groups == [1]
                and selected_groups == target_groups
                and isinstance(roles, dict)
                and roles.get("1") == "magnetic_body"
                and excluded_air == [0]
                and roles.get("0") == "air"
                and bool(selection.get("selection_generation"))
            )
        constraint_basis = artifact_identity.get("virtual_work_constraint_basis")
        if constraint_basis is not None:
            virtual_work_constraint_basis_ok = (
                isinstance(constraint_basis, dict)
                and constraint_basis.get("direct_force_constraint") == "fixed_current"
                and constraint_basis.get("coenergy_derivative_constraint")
                == "fixed_current"
                and bool(constraint_basis.get("current_control_generation"))
                and constraint_basis.get("derivative_control_generation")
                == constraint_basis.get("current_control_generation")
                and constraint_basis.get("flux_constraint_active") is False
            )
        harmonic_basis = artifact_identity.get("eddy_loss_harmonic_basis")
        if harmonic_basis is not None:
            amplitude_frequencies = (
                harmonic_basis.get("harmonic_frequency_hz")
                if isinstance(harmonic_basis, dict)
                else None
            )
            material_frequencies = (
                harmonic_basis.get("skin_depth_state_frequency_hz")
                if isinstance(harmonic_basis, dict)
                else None
            )
            frequency_rows_valid = (
                isinstance(amplitude_frequencies, list)
                and isinstance(material_frequencies, list)
                and bool(amplitude_frequencies)
                and len(amplitude_frequencies) == len(material_frequencies)
                and all(
                    isinstance(value, (int, float))
                    and math.isfinite(float(value))
                    and float(value) > 0.0
                    for value in amplitude_frequencies + material_frequencies
                )
            )
            eddy_loss_harmonic_basis_ok = (
                isinstance(harmonic_basis, dict)
                and frequency_rows_valid
                and [float(value) for value in amplitude_frequencies]
                == [float(value) for value in material_frequencies]
                and all(
                    right > left
                    for left, right in zip(
                        amplitude_frequencies, amplitude_frequencies[1:]
                    )
                )
                and bool(harmonic_basis.get("amplitude_basis_id"))
                and harmonic_basis.get("material_state_basis_id")
                == harmonic_basis.get("amplitude_basis_id")
                and bool(harmonic_basis.get("solve_generation"))
                and harmonic_basis.get("material_state_solve_generation")
                == harmonic_basis.get("solve_generation")
            )
        force_measure = artifact_identity.get("axisymmetric_force_measure_identity")
        if force_measure is not None:
            axisymmetric_force_measure_ok = (
                isinstance(force_measure, dict)
                and force_measure.get("formulation") == "axisymmetric"
                and force_measure.get("integration_measure") == "2*pi*r*dr*dz"
                and force_measure.get("reference_integration_measure")
                == force_measure.get("integration_measure")
                and force_measure.get("radius_weighting_basis_id")
                == "axisymmetric-radius-weighted-v1"
                and force_measure.get("force_result_basis_id")
                == force_measure.get("radius_weighting_basis_id")
                and force_measure.get("radius_coordinate_frame") == "cylindrical-rz"
                and force_measure.get("force_component_frame") == "global-z"
                and bool(force_measure.get("solve_generation"))
                and force_measure.get("integration_solve_generation")
                == force_measure.get("solve_generation")
            )
        material_frequency = artifact_identity.get(
            "eddy_loss_material_frequency_identity"
        )
        if material_frequency is not None:
            try:
                field_frequency = float(
                    material_frequency.get("field_solution_frequency_hz")
                )
                loss_frequency = float(
                    material_frequency.get("loss_evaluation_frequency_hz")
                )
                material_frequency_hz = float(
                    material_frequency.get("material_state_frequency_hz")
                )
            except (AttributeError, TypeError, ValueError):
                field_frequency = math.nan
                loss_frequency = math.nan
                material_frequency_hz = math.nan
            assignment_generation = (
                material_frequency.get("material_assignment_generation")
                if isinstance(material_frequency, dict)
                else None
            )
            conductivity_digest = str(
                material_frequency.get("conductivity_sha256", "")
                if isinstance(material_frequency, dict)
                else ""
            ).lower()
            lamination_digest = str(
                material_frequency.get("lamination_state_sha256", "")
                if isinstance(material_frequency, dict)
                else ""
            ).lower()
            eddy_loss_material_frequency_ok = (
                isinstance(material_frequency, dict)
                and math.isfinite(field_frequency)
                and field_frequency > 0.0
                and loss_frequency == field_frequency
                and material_frequency_hz == field_frequency
                and bool(assignment_generation)
                and material_frequency.get("conductivity_material_generation")
                == assignment_generation
                and material_frequency.get("lamination_material_generation")
                == assignment_generation
                and bool(material_frequency.get("solve_generation"))
                and material_frequency.get("material_state_solve_generation")
                == material_frequency.get("solve_generation")
                and len(conductivity_digest) == 64
                and len(lamination_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in conductivity_digest + lamination_digest
                )
            )
        weighted_mask = artifact_identity.get("weighted_stress_mask_mesh_identity")
        if weighted_mask is not None:
            mask_digest = str(
                weighted_mask.get("weighted_mask_sha256", "")
                if isinstance(weighted_mask, dict)
                else ""
            ).lower()
            force_mask_digest = str(
                weighted_mask.get("force_mask_sha256", "")
                if isinstance(weighted_mask, dict)
                else ""
            ).lower()
            mesh_generations = (
                [
                    weighted_mask.get("active_air_mesh_generation"),
                    weighted_mask.get("field_solution_mesh_generation"),
                    weighted_mask.get("weighted_mask_mesh_generation"),
                    weighted_mask.get("force_integration_mesh_generation"),
                ]
                if isinstance(weighted_mask, dict)
                else []
            )
            weighted_stress_mask_mesh_identity_ok = (
                isinstance(weighted_mask, dict)
                and all(isinstance(value, str) and bool(value) for value in mesh_generations)
                and len(set(mesh_generations)) == 1
                and weighted_mask.get("mask_basis") == "nodal_weighting_function"
                and weighted_mask.get("force_method") == "weighted_stress_tensor"
                and len(mask_digest) == len(force_mask_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in mask_digest + force_mask_digest
                )
                and force_mask_digest == mask_digest
            )
        phasor_basis = artifact_identity.get("complex_current_phasor_basis_identity")
        if phasor_basis is not None:
            try:
                source_scale = float(phasor_basis.get("source_scale_to_rms"))
                field_scale = float(phasor_basis.get("field_scale_to_rms"))
                result_scale = float(phasor_basis.get("force_loss_scale_to_rms"))
            except (AttributeError, TypeError, ValueError):
                source_scale = math.nan
                field_scale = math.nan
                result_scale = math.nan
            source_basis = (
                phasor_basis.get("source_current_basis")
                if isinstance(phasor_basis, dict)
                else None
            )
            expected_scale = {
                "rms_phasor": 1.0,
                "peak_phasor": 1.0 / math.sqrt(2.0),
            }.get(source_basis)
            complex_current_phasor_basis_identity_ok = (
                isinstance(phasor_basis, dict)
                and expected_scale is not None
                and phasor_basis.get("field_current_basis") == source_basis
                and phasor_basis.get("force_loss_current_basis") == source_basis
                and all(
                    math.isfinite(value)
                    and math.isclose(value, expected_scale, rel_tol=1.0e-12, abs_tol=1.0e-15)
                    for value in (source_scale, field_scale, result_scale)
                )
                and phasor_basis.get("complex_time_convention")
                in {"exp(+jwt)", "exp(-jwt)"}
                and phasor_basis.get("result_time_convention")
                == phasor_basis.get("complex_time_convention")
                and bool(phasor_basis.get("solve_generation"))
                and phasor_basis.get("result_generation")
                == phasor_basis.get("solve_generation")
            )

        radius_jacobian = artifact_identity.get(
            "axisymmetric_force_radius_jacobian_coordinate_identity"
        )
        if radius_jacobian is not None:
            length_units = {"m": 1.0, "cm": 1.0e-2, "mm": 1.0e-3}
            radius_unit = (
                str(radius_jacobian.get("radius_length_unit", "")).strip()
                if isinstance(radius_jacobian, dict)
                else ""
            )
            stress_unit = (
                str(radius_jacobian.get("stress_coordinate_length_unit", "")).strip()
                if isinstance(radius_jacobian, dict)
                else ""
            )
            radius_digest = str(
                radius_jacobian.get("radius_coordinate_sha256", "")
                if isinstance(radius_jacobian, dict)
                else ""
            ).lower()
            try:
                radius_scale = float(radius_jacobian.get("radius_scale_to_m"))
                stress_scale = float(
                    radius_jacobian.get("stress_coordinate_scale_to_m")
                )
            except (AttributeError, TypeError, ValueError):
                radius_scale = math.nan
                stress_scale = math.nan
            coordinate_generation = (
                radius_jacobian.get("coordinate_generation")
                if isinstance(radius_jacobian, dict)
                else None
            )
            solution_generation = (
                radius_jacobian.get("field_solution_generation")
                if isinstance(radius_jacobian, dict)
                else None
            )
            expected_scale = length_units.get(radius_unit)
            axisymmetric_force_radius_jacobian_coordinate_identity_ok = (
                isinstance(radius_jacobian, dict)
                and bool(solution_generation)
                and radius_jacobian.get("stress_field_solution_generation")
                == solution_generation
                and bool(coordinate_generation)
                and radius_jacobian.get("stress_coordinate_generation")
                == coordinate_generation
                and radius_jacobian.get("radius_jacobian_coordinate_generation")
                == coordinate_generation
                and radius_jacobian.get("force_integration_coordinate_generation")
                == coordinate_generation
                and radius_jacobian.get("radius_coordinate_frame") == "cylindrical-rz"
                and radius_jacobian.get("stress_coordinate_frame")
                == radius_jacobian.get("radius_coordinate_frame")
                and expected_scale is not None
                and stress_unit == radius_unit
                and math.isclose(radius_scale, expected_scale, rel_tol=0.0, abs_tol=0.0)
                and math.isclose(stress_scale, radius_scale, rel_tol=0.0, abs_tol=0.0)
                and len(radius_digest) == 64
                and all(character in "0123456789abcdef" for character in radius_digest)
                and radius_jacobian.get("force_radius_coordinate_sha256")
                == radius_digest
                and radius_jacobian.get("integration_measure") == "2*pi*r*dr*dz"
            )

        bh_interpolation = artifact_identity.get(
            "nonlinear_bh_interpolation_extrapolation_identity"
        )
        if bh_interpolation is not None:
            table_digest = str(
                bh_interpolation.get("bh_table_sha256", "")
                if isinstance(bh_interpolation, dict)
                else ""
            ).lower()
            material_generation = (
                bh_interpolation.get("material_generation")
                if isinstance(bh_interpolation, dict)
                else None
            )
            interpolation_method = (
                bh_interpolation.get("interpolation_method")
                if isinstance(bh_interpolation, dict)
                else None
            )
            endpoint_branch = (
                bh_interpolation.get("endpoint_extrapolation_branch")
                if isinstance(bh_interpolation, dict)
                else None
            )
            solve_generation = (
                bh_interpolation.get("solve_generation")
                if isinstance(bh_interpolation, dict)
                else None
            )
            nonlinear_bh_interpolation_extrapolation_identity_ok = (
                isinstance(bh_interpolation, dict)
                and bool(material_generation)
                and bh_interpolation.get("bh_table_material_generation")
                == material_generation
                and bh_interpolation.get("field_solution_material_generation")
                == material_generation
                and len(table_digest) == 64
                and all(character in "0123456789abcdef" for character in table_digest)
                and bh_interpolation.get("field_bh_table_sha256") == table_digest
                and interpolation_method in {"monotone_piecewise_linear", "cubic_hermite"}
                and bh_interpolation.get("field_interpolation_method")
                == interpolation_method
                and endpoint_branch in {"last_segment_slope", "constant_mu0"}
                and bh_interpolation.get("field_endpoint_extrapolation_branch")
                == endpoint_branch
                and bh_interpolation.get("evaluation_region")
                in {"inside_table", "lower_endpoint_extrapolation", "upper_endpoint_extrapolation"}
                and bool(solve_generation)
                and bh_interpolation.get("field_state_solve_generation")
                == solve_generation
            )

        material_interface = artifact_identity.get(
            "weighted_stress_mask_material_interface_identity"
        )
        if material_interface is not None:
            solve_generation = (
                material_interface.get("field_solution_generation")
                if isinstance(material_interface, dict)
                else None
            )
            material_generation = (
                material_interface.get("material_generation")
                if isinstance(material_interface, dict)
                else None
            )
            topology_generation = (
                material_interface.get("mesh_topology_generation")
                if isinstance(material_interface, dict)
                else None
            )
            integration_material_ids = (
                material_interface.get("integration_material_ids")
                if isinstance(material_interface, dict)
                else None
            )
            mask_material_ids = (
                material_interface.get("mask_material_ids")
                if isinstance(material_interface, dict)
                else None
            )
            interface_face_ids = (
                material_interface.get("material_interface_face_ids")
                if isinstance(material_interface, dict)
                else None
            )
            excluded_face_ids = (
                material_interface.get("mask_excluded_interface_face_ids")
                if isinstance(material_interface, dict)
                else None
            )
            mask_digest = str(
                material_interface.get("mask_sha256", "")
                if isinstance(material_interface, dict)
                else ""
            ).lower()
            weighted_stress_mask_material_interface_identity_ok = (
                isinstance(material_interface, dict)
                and bool(solve_generation)
                and material_interface.get("stress_mask_solution_generation")
                == solve_generation
                and bool(material_generation)
                and material_interface.get("mask_material_generation")
                == material_generation
                and bool(topology_generation)
                and material_interface.get("mask_topology_generation")
                == topology_generation
                and isinstance(integration_material_ids, list)
                and bool(integration_material_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in integration_material_ids
                )
                and len(set(integration_material_ids))
                == len(integration_material_ids)
                and mask_material_ids == integration_material_ids
                and isinstance(interface_face_ids, list)
                and bool(interface_face_ids)
                and all(
                    isinstance(value, int) and not isinstance(value, bool)
                    for value in interface_face_ids
                )
                and len(set(interface_face_ids)) == len(interface_face_ids)
                and excluded_face_ids == interface_face_ids
                and len(mask_digest) == 64
                and all(character in "0123456789abcdef" for character in mask_digest)
                and str(material_interface.get("stress_mask_sha256", "")).lower()
                == mask_digest
            )

        coil_voltage = artifact_identity.get(
            "axisymmetric_coil_voltage_measure_identity"
        )
        if coil_voltage is not None:
            solve_generation = (
                coil_voltage.get("solve_generation")
                if isinstance(coil_voltage, dict)
                else None
            )
            coordinate_generation = (
                coil_voltage.get("radius_coordinate_generation")
                if isinstance(coil_voltage, dict)
                else None
            )
            radius_digest = str(
                coil_voltage.get("radius_coordinate_sha256", "")
                if isinstance(coil_voltage, dict)
                else ""
            ).lower()
            factor_count = (
                coil_voltage.get("two_pi_radius_factor_count")
                if isinstance(coil_voltage, dict)
                else None
            )
            axisymmetric_coil_voltage_measure_identity_ok = (
                isinstance(coil_voltage, dict)
                and bool(solve_generation)
                and coil_voltage.get("winding_voltage_generation")
                == solve_generation
                and bool(coordinate_generation)
                and coil_voltage.get("voltage_radius_coordinate_generation")
                == coordinate_generation
                and coil_voltage.get("potential_voltage_basis") == "per_radian"
                and coil_voltage.get("reported_voltage_basis") == "total_3d"
                and coil_voltage.get("integration_measure") == "2*pi*r*dr*dz"
                and isinstance(factor_count, int)
                and not isinstance(factor_count, bool)
                and factor_count == 1
                and len(radius_digest) == 64
                and all(character in "0123456789abcdef" for character in radius_digest)
                and str(
                    coil_voltage.get("voltage_radius_coordinate_sha256", "")
                ).lower()
                == radius_digest
            )

        nonlinear_energy = artifact_identity.get(
            "nonlinear_energy_coenergy_bh_iteration_identity"
        )
        if nonlinear_energy is not None:
            field_generation = (
                nonlinear_energy.get("field_solve_generation")
                if isinstance(nonlinear_energy, dict)
                else None
            )
            iteration = (
                nonlinear_energy.get("nonlinear_iteration")
                if isinstance(nonlinear_energy, dict)
                else None
            )
            branch_generation = (
                nonlinear_energy.get("bh_branch_generation")
                if isinstance(nonlinear_energy, dict)
                else None
            )
            state_digest = str(
                nonlinear_energy.get("bh_state_sha256", "")
                if isinstance(nonlinear_energy, dict)
                else ""
            ).lower()
            nonlinear_energy_coenergy_bh_iteration_identity_ok = (
                isinstance(nonlinear_energy, dict)
                and bool(field_generation)
                and nonlinear_energy.get("energy_field_solve_generation")
                == field_generation
                and nonlinear_energy.get("coenergy_field_solve_generation")
                == field_generation
                and isinstance(iteration, int)
                and not isinstance(iteration, bool)
                and iteration >= 0
                and nonlinear_energy.get("energy_nonlinear_iteration") == iteration
                and nonlinear_energy.get("coenergy_nonlinear_iteration") == iteration
                and bool(branch_generation)
                and nonlinear_energy.get("energy_bh_branch_generation")
                == branch_generation
                and nonlinear_energy.get("coenergy_bh_branch_generation")
                == branch_generation
                and len(state_digest) == 64
                and all(character in "0123456789abcdef" for character in state_digest)
                and str(nonlinear_energy.get("energy_bh_state_sha256", "")).lower()
                == state_digest
                and str(nonlinear_energy.get("coenergy_bh_state_sha256", "")).lower()
                == state_digest
            )

        virtual_displacement = artifact_identity.get(
            "virtual_displacement_force_geometry_field_identity"
        )
        if virtual_displacement is not None:
            base_geometry = (
                virtual_displacement.get("base_geometry_generation")
                if isinstance(virtual_displacement, dict)
                else None
            )
            perturbed_geometry = (
                virtual_displacement.get("perturbed_geometry_generation")
                if isinstance(virtual_displacement, dict)
                else None
            )
            base_solve = (
                virtual_displacement.get("base_field_solve_generation")
                if isinstance(virtual_displacement, dict)
                else None
            )
            perturbed_solve = (
                virtual_displacement.get("perturbed_field_solve_generation")
                if isinstance(virtual_displacement, dict)
                else None
            )
            base_digest = str(
                virtual_displacement.get("base_field_sha256", "")
                if isinstance(virtual_displacement, dict)
                else ""
            ).lower()
            perturbed_digest = str(
                virtual_displacement.get("perturbed_field_sha256", "")
                if isinstance(virtual_displacement, dict)
                else ""
            ).lower()
            try:
                displacement_step = float(virtual_displacement["displacement_step_m"])
                field_displacement_step = float(
                    virtual_displacement["field_displacement_step_m"]
                )
            except (KeyError, TypeError, ValueError):
                displacement_step = math.nan
                field_displacement_step = math.nan
            virtual_displacement_force_geometry_field_identity_ok = (
                isinstance(virtual_displacement, dict)
                and bool(virtual_displacement.get("force_evaluation_generation"))
                and bool(base_geometry)
                and bool(perturbed_geometry)
                and base_geometry != perturbed_geometry
                and virtual_displacement.get("base_field_geometry_generation")
                == base_geometry
                and virtual_displacement.get("perturbed_field_geometry_generation")
                == perturbed_geometry
                and bool(base_solve)
                and virtual_displacement.get("force_base_field_solve_generation")
                == base_solve
                and bool(perturbed_solve)
                and virtual_displacement.get("force_perturbed_field_solve_generation")
                == perturbed_solve
                and math.isfinite(displacement_step)
                and displacement_step > 0.0
                and field_displacement_step == displacement_step
                and len(base_digest) == len(perturbed_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in base_digest + perturbed_digest
                )
                and str(
                    virtual_displacement.get("force_base_field_sha256", "")
                ).lower()
                == base_digest
                and str(
                    virtual_displacement.get("force_perturbed_field_sha256", "")
                ).lower()
                == perturbed_digest
            )

    finite = all(math.isfinite(value) for value in x + w + force)
    increasing = finite and all(right > left for left, right in zip(x, x[1:]))
    rows = []
    central_errors = []
    if finite and increasing and len(x) >= 2:
        for index in range(len(x)):
            if index == 0:
                derivative = (w[1] - w[0]) / (x[1] - x[0])
                stencil = "forward"
            elif index == len(x) - 1:
                derivative = (w[-1] - w[-2]) / (x[-1] - x[-2])
                stencil = "backward"
            else:
                derivative = (w[index + 1] - w[index - 1]) / (x[index + 1] - x[index - 1])
                stencil = "central"
            scale = max(abs(force[index]), abs(derivative), 1.0e-30)
            relative_error = abs(derivative - force[index]) / scale
            if stencil == "central":
                central_errors.append(relative_error)
            rows.append(
                {
                    "index": index,
                    "position_m": x[index],
                    "direct_force_N": force[index],
                    "coenergy_derivative_force_N": derivative,
                    "stencil": stencil,
                    "relative_error": relative_error,
                }
            )

    max_error = max(central_errors) if central_errors else math.inf
    checks = {
        "sample_count_sufficient": len(x) >= min_sample_count,
        "all_finite": finite,
        "positions_strictly_increase": increasing,
        "constant_current_coenergy_recorded": energy_kind == "constant_current_coenergy",
        "coenergy_nontrivial": finite and bool(w) and max(w) > min(w),
        "central_rows_available": len(central_errors) >= 3,
        "central_virtual_work_matches_direct_force": max_error <= max_central_relative_error,
        "force_and_coenergy_share_load_step_snapshot": force_snapshot_ok,
        "coenergy_stencil_uses_one_mesh_family_generation": mesh_family_ok,
        "displacement_axis_uses_one_si_unit": displacement_unit_ok,
        "force_vectors_share_transformed_frame": force_frame_ok,
        "axisymmetric_force_is_already_total_3d": force_normalization_ok,
        "weighted_stress_selects_only_target_magnetic_body": force_body_selection_ok,
        "force_and_coenergy_use_same_fixed_current_constraint": (
            virtual_work_constraint_basis_ok
        ),
        "eddy_loss_harmonics_share_frequency_and_material_basis": (
            eddy_loss_harmonic_basis_ok
        ),
        "axisymmetric_force_uses_radius_weighted_measure": (
            axisymmetric_force_measure_ok
        ),
        "eddy_loss_uses_current_frequency_and_material_generation": (
            eddy_loss_material_frequency_ok
        ),
        "weighted_stress_mask_matches_current_air_mesh_generation": (
            weighted_stress_mask_mesh_identity_ok
        ),
        "complex_current_force_and_loss_share_phasor_basis": (
            complex_current_phasor_basis_identity_ok
        ),
        "axisymmetric_force_radius_jacobian_uses_current_coordinates": (
            axisymmetric_force_radius_jacobian_coordinate_identity_ok
        ),
        "nonlinear_bh_state_uses_one_interpolation_and_extrapolation_branch": (
            nonlinear_bh_interpolation_extrapolation_identity_ok
        ),
        "weighted_stress_mask_excludes_current_material_interfaces": (
            weighted_stress_mask_material_interface_identity_ok
        ),
        "axisymmetric_coil_voltage_applies_two_pi_radius_once": (
            axisymmetric_coil_voltage_measure_identity_ok
        ),
        "nonlinear_energy_and_coenergy_share_current_bh_iteration": (
            nonlinear_energy_coenergy_bh_iteration_identity_ok
        ),
        "virtual_displacement_force_uses_paired_geometry_field_generations": (
            virtual_displacement_force_geometry_field_identity_ok
        ),
    }
    return {
        "policy": "force_coenergy_displacement_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(x),
        "central_sample_count": len(central_errors),
        "max_central_relative_error": max_error,
        "mean_central_relative_error": (
            sum(central_errors) / len(central_errors) if central_errors else None
        ),
        "endpoint_errors_are_diagnostic_only": True,
        "checks": checks,
        "warnings": [] if identity_present else ["artifact_identity_not_recorded"],
        "rows": rows,
        "lesson": (
            "At fixed current, direct force projected onto the displacement axis "
            "must match dW'/dx. Use central differences for the acceptance gate; "
            "one-sided endpoint errors are diagnostics only."
        ),
    }
