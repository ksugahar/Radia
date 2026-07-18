"""Virtual-work consistency gate for displacement-force sweeps."""
from __future__ import annotations

import math

from .magnetic_force_method_profile_gate import (
    _airgap_stress_harmonic_torque_identity_ok,
    _laminated_core_loss_identity_ok,
)


def _valid_sha256(value):
    digest = str(value or "").lower()
    return len(digest) == 64 and all(
        character in "0123456789abcdef" for character in digest
    )


def _nonlinear_force_operating_point_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("solve_generation", "")).strip()
    branch_id = str(value.get("branch_id", "")).strip()
    try:
        current = float(value.get("operating_point_current_a"))
        force_current = float(value.get("force_operating_point_current_a"))
        flux_density = [
            float(item) for item in value.get("operating_point_flux_density_t", [])
        ]
        force_flux_density = [
            float(item)
            for item in value.get("force_operating_point_flux_density_t", [])
        ]
        force = [float(item) for item in value.get("force_n", [])]
        reported_force = [
            float(item) for item in value.get("reported_force_n", [])
        ]
    except (TypeError, ValueError):
        return False
    permeability_digest = str(value.get("permeability_state_sha256", "")).lower()
    mesh_digest = str(value.get("force_mesh_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "branch_solve_generation",
                "operating_point_solve_generation",
                "permeability_solve_generation",
                "force_mesh_solve_generation",
                "force_result_solve_generation",
            )
        )
        and bool(branch_id)
        and value.get("force_branch_id") == branch_id
        and math.isfinite(current)
        and math.isclose(force_current, current, rel_tol=0.0, abs_tol=1.0e-15)
        and bool(flux_density)
        and all(math.isfinite(item) for item in flux_density)
        and force_flux_density == flux_density
        and _valid_sha256(permeability_digest)
        and value.get("force_permeability_state_sha256") == permeability_digest
        and _valid_sha256(mesh_digest)
        and value.get("integrated_force_mesh_sha256") == mesh_digest
        and bool(force)
        and all(math.isfinite(item) for item in force)
        and reported_force == force
    )


def _sliding_band_harmonic_torque_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("sweep_generation", "")).strip()
    try:
        angles = [float(item) for item in value.get("rotor_angles_deg", [])]
        torque_angles = [
            float(item) for item in value.get("torque_rotor_angles_deg", [])
        ]
        samples = [float(item) for item in value.get("torque_samples_nm", [])]
        harmonic_samples = [
            float(item) for item in value.get("harmonic_torque_samples_nm", [])
        ]
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        reported_orders = [
            int(item) for item in value.get("reported_harmonic_orders", [])
        ]
        amplitudes = [
            float(item) for item in value.get("harmonic_amplitudes_nm", [])
        ]
        reported_amplitudes = [
            float(item)
            for item in value.get("reported_harmonic_amplitudes_nm", [])
        ]
    except (TypeError, ValueError):
        return False
    mesh_digest = str(value.get("airgap_mesh_sha256", "")).lower()
    current_digest = str(value.get("phase_current_table_sha256", "")).lower()
    sample_digest = str(value.get("torque_sample_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "angle_sweep_generation",
                "airgap_mesh_sweep_generation",
                "phase_current_sweep_generation",
                "torque_sample_sweep_generation",
                "harmonic_sweep_generation",
            )
        )
        and len(angles) >= 3
        and all(math.isfinite(item) for item in angles)
        and all(left < right for left, right in zip(angles, angles[1:]))
        and torque_angles == angles
        and _valid_sha256(mesh_digest)
        and value.get("torque_airgap_mesh_sha256") == mesh_digest
        and _valid_sha256(current_digest)
        and value.get("torque_phase_current_table_sha256") == current_digest
        and len(samples) == len(angles)
        and all(math.isfinite(item) for item in samples)
        and harmonic_samples == samples
        and bool(orders)
        and all(item >= 0 for item in orders)
        and len(set(orders)) == len(orders)
        and reported_orders == orders
        and len(amplitudes) == len(orders)
        and all(math.isfinite(item) and item >= 0.0 for item in amplitudes)
        and reported_amplitudes == amplitudes
        and _valid_sha256(sample_digest)
        and value.get("harmonic_sample_table_sha256") == sample_digest
    )


def _weighted_stress_energy_derivative_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("force_generation", "")).strip()
    mesh_digest = str(value.get("weighted_stress_mesh_sha256", "")).lower()
    result_digest = str(value.get("force_result_sha256", "")).lower()
    try:
        weighted_force = [
            float(item) for item in value.get("weighted_stress_force_n", [])
        ]
        derivative_force = [
            float(item) for item in value.get("energy_derivative_force_n", [])
        ]
    except (TypeError, ValueError):
        return False
    frame_id = str(value.get("displacement_frame_id", "")).strip()
    displacement_unit = str(value.get("displacement_unit", "")).strip()
    force_unit = str(value.get("force_unit", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "weighted_stress_force_generation",
                "energy_derivative_force_generation",
                "mesh_force_generation",
                "displacement_frame_force_generation",
                "unit_force_generation",
                "result_force_generation",
            )
        )
        and _valid_sha256(mesh_digest)
        and value.get("energy_derivative_mesh_sha256") == mesh_digest
        and bool(frame_id)
        and value.get("energy_derivative_displacement_frame_id") == frame_id
        and displacement_unit == "m"
        and value.get("energy_derivative_displacement_unit") == displacement_unit
        and force_unit == "N"
        and value.get("energy_derivative_force_unit") == force_unit
        and bool(weighted_force)
        and all(math.isfinite(item) for item in weighted_force)
        and derivative_force == weighted_force
        and _valid_sha256(result_digest)
        and value.get("energy_derivative_result_sha256") == result_digest
    )


def _axisymmetric_revolved_energy_force_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("axisymmetric_generation", "")).strip()
    measure = str(value.get("jacobian_measure", "")).strip().lower()
    derham_id = str(value.get("derham_sequence_id", "")).strip()
    try:
        angle = float(value.get("revolution_angle_deg"))
        axisymmetric_energy = float(value.get("axisymmetric_energy_j"))
        revolved_energy = float(value.get("revolved_energy_j"))
        axisymmetric_force = [
            float(item) for item in value.get("axisymmetric_force_n", [])
        ]
        revolved_force = [
            float(item) for item in value.get("revolved_force_n", [])
        ]
    except (TypeError, ValueError):
        return False
    digest_pairs = (
        ("field_state_sha256", "revolved_field_state_sha256"),
        ("material_map_sha256", "revolved_material_map_sha256"),
        ("axisymmetric_mesh_sha256", "revolved_source_mesh_sha256"),
        ("result_sha256", "revolved_result_sha256"),
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "jacobian_axisymmetric_generation",
                "field_axisymmetric_generation",
                "material_axisymmetric_generation",
                "mesh_axisymmetric_generation",
                "revolved_result_axisymmetric_generation",
            )
        )
        and measure == "2*pi*r"
        and str(value.get("revolved_jacobian_measure", "")).strip().lower()
        == measure
        and all(
            _valid_sha256(value.get(source))
            and value.get(result) == value.get(source)
            for source, result in digest_pairs
        )
        and math.isfinite(angle)
        and math.isclose(angle, 360.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isfinite(axisymmetric_energy)
        and math.isfinite(revolved_energy)
        and math.isclose(
            revolved_energy,
            axisymmetric_energy,
            rel_tol=1.0e-12,
            abs_tol=1.0e-18,
        )
        and bool(axisymmetric_force)
        and all(math.isfinite(item) for item in axisymmetric_force)
        and revolved_force == axisymmetric_force
        and bool(derham_id)
        and value.get("revolved_derham_sequence_id") == derham_id
    )


def _nonlinear_bh_incremental_force_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("nonlinear_generation", "")).strip()
    material_ids = [
        str(item).strip() for item in value.get("nonlinear_material_ids", [])
    ]
    result_material_ids = [
        str(item).strip()
        for item in value.get("result_nonlinear_material_ids", [])
    ]
    branch_id = str(value.get("branch_id", "")).strip()
    try:
        current = float(value.get("load_current_a"))
        result_current = float(value.get("result_load_current_a"))
        energy = float(value.get("magnetic_energy_j"))
        result_energy = float(value.get("result_magnetic_energy_j"))
        coenergy = float(value.get("magnetic_coenergy_j"))
        result_coenergy = float(value.get("result_magnetic_coenergy_j"))
        force = [float(item) for item in value.get("incremental_force_n", [])]
        result_force = [
            float(item) for item in value.get("result_incremental_force_n", [])
        ]
    except (TypeError, ValueError):
        return False
    digest_pairs = (
        ("bh_curve_sha256", "result_bh_curve_sha256"),
        ("material_map_sha256", "result_material_map_sha256"),
        ("incremental_state_sha256", "result_incremental_state_sha256"),
        ("mesh_sha256", "result_mesh_sha256"),
        ("result_sha256", "accepted_result_sha256"),
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "bh_curve_nonlinear_generation",
                "material_nonlinear_generation",
                "branch_nonlinear_generation",
                "incremental_state_nonlinear_generation",
                "mesh_nonlinear_generation",
                "energy_nonlinear_generation",
                "coenergy_nonlinear_generation",
                "force_nonlinear_generation",
                "result_nonlinear_generation",
            )
        )
        and bool(material_ids)
        and all(material_ids)
        and len(set(material_ids)) == len(material_ids)
        and result_material_ids == material_ids
        and bool(branch_id)
        and value.get("result_branch_id") == branch_id
        and math.isfinite(current)
        and math.isclose(result_current, current, rel_tol=0.0, abs_tol=1.0e-15)
        and all(
            _valid_sha256(value.get(source))
            and value.get(result) == value.get(source)
            for source, result in digest_pairs
        )
        and math.isfinite(energy)
        and energy >= 0.0
        and math.isclose(result_energy, energy, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isfinite(coenergy)
        and coenergy >= 0.0
        and math.isclose(
            result_coenergy, coenergy, rel_tol=1.0e-12, abs_tol=1.0e-18
        )
        and bool(force)
        and all(math.isfinite(item) for item in force)
        and result_force == force
    )


def _open_boundary_decay_multipole_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("boundary_generation", "")).strip()
    boundary_type = str(value.get("boundary_type", "")).strip()
    try:
        source_radius = float(value.get("source_radius_m"))
        result_source_radius = float(value.get("result_source_radius_m"))
        outer_radius = float(value.get("outer_radius_m"))
        result_outer_radius = float(value.get("result_outer_radius_m"))
        multipole_order = int(value.get("multipole_order"))
        result_multipole_order = int(value.get("result_multipole_order"))
        radii = [float(item) for item in value.get("decay_sample_radii_m", [])]
        result_radii = [
            float(item) for item in value.get("result_decay_sample_radii_m", [])
        ]
        flux_density = [
            float(item) for item in value.get("decay_flux_density_t", [])
        ]
        result_flux_density = [
            float(item) for item in value.get("result_decay_flux_density_t", [])
        ]
    except (TypeError, ValueError):
        return False
    digest_pairs = (
        ("material_map_sha256", "result_material_map_sha256"),
        ("mesh_sha256", "result_mesh_sha256"),
        ("multipole_moment_sha256", "result_multipole_moment_sha256"),
        ("result_sha256", "accepted_result_sha256"),
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "domain_boundary_generation",
                "mesh_boundary_generation",
                "material_boundary_generation",
                "multipole_boundary_generation",
                "decay_boundary_generation",
                "result_boundary_generation",
            )
        )
        and boundary_type == "asymptotic_multipole"
        and value.get("result_boundary_type") == boundary_type
        and math.isfinite(source_radius)
        and source_radius > 0.0
        and math.isclose(
            result_source_radius, source_radius, rel_tol=0.0, abs_tol=1.0e-18
        )
        and math.isfinite(outer_radius)
        and outer_radius > source_radius
        and math.isclose(
            result_outer_radius, outer_radius, rel_tol=0.0, abs_tol=1.0e-18
        )
        and multipole_order >= 1
        and result_multipole_order == multipole_order
        and all(
            _valid_sha256(value.get(source))
            and value.get(result) == value.get(source)
            for source, result in digest_pairs
        )
        and len(radii) >= 3
        and all(math.isfinite(item) for item in radii)
        and all(left < right for left, right in zip(radii, radii[1:]))
        and radii[0] > source_radius
        and radii[-1] <= outer_radius
        and result_radii == radii
        and len(flux_density) == len(radii)
        and all(math.isfinite(item) and item >= 0.0 for item in flux_density)
        and all(
            left >= right for left, right in zip(flux_density, flux_density[1:])
        )
        and result_flux_density == flux_density
    )


def _weighted_stress_tensor_closure_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("force_generation", "")).strip()
    body_groups = value.get("body_group_ids", [])
    mask_nodes = value.get("weighted_mask_node_ids", [])
    region = str(value.get("integration_region", "")).strip()
    try:
        body_groups = [int(item) for item in body_groups]
        result_body_groups = [int(item) for item in value.get("result_body_group_ids", [])]
        mask_nodes = [int(item) for item in mask_nodes]
        result_mask_nodes = [int(item) for item in value.get("result_weighted_mask_node_ids", [])]
        force = [float(item) for item in value.get("weighted_force_n", [])]
        result_force = [float(item) for item in value.get("result_weighted_force_n", [])]
        torque = [float(item) for item in value.get("weighted_torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_weighted_torque_nm", [])]
        energy = float(value.get("magnetic_energy_j"))
        result_energy = float(value.get("result_magnetic_energy_j"))
        coenergy = float(value.get("magnetic_coenergy_j"))
        result_coenergy = float(value.get("result_magnetic_coenergy_j"))
    except (TypeError, ValueError):
        return False
    digest_pairs = (
        ("mesh_sha256", "result_mesh_sha256"),
        ("mask_sha256", "result_mask_sha256"),
        ("integration_region_sha256", "result_integration_region_sha256"),
        ("result_sha256", "accepted_result_sha256"),
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "mask_force_generation", "mesh_force_generation", "region_force_generation",
            "energy_force_generation", "torque_force_generation", "result_force_generation"))
        and bool(body_groups) and all(item > 0 for item in body_groups)
        and len(set(body_groups)) == len(body_groups) and result_body_groups == body_groups
        and bool(mask_nodes) and all(item > 0 for item in mask_nodes)
        and len(set(mask_nodes)) == len(mask_nodes) and result_mask_nodes == mask_nodes
        and bool(region) and value.get("result_integration_region") == region
        and all(_valid_sha256(value.get(source)) and value.get(result) == value.get(source)
                for source, result in digest_pairs)
        and len(force) == 3 and all(math.isfinite(item) for item in force) and result_force == force
        and len(torque) == 3 and all(math.isfinite(item) for item in torque) and result_torque == torque
        and math.isfinite(energy) and energy >= 0.0
        and math.isclose(result_energy, energy, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isfinite(coenergy) and coenergy >= 0.0
        and math.isclose(result_coenergy, coenergy, rel_tol=1.0e-12, abs_tol=1.0e-18)
    )


def _axisymmetric_planar_normalization_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("normalization_generation", "")).strip()
    try:
        depth = float(value.get("planar_depth_m"))
        result_depth = float(value.get("result_planar_depth_m"))
        planar_per_m = float(value.get("planar_force_n_per_m"))
        result_planar_per_m = float(value.get("result_planar_force_n_per_m"))
        planar_total = float(value.get("planar_total_force_n"))
        result_planar_total = float(value.get("result_planar_total_force_n"))
        radius = float(value.get("axisymmetric_radius_m"))
        result_radius = float(value.get("result_axisymmetric_radius_m"))
        meridian_force = float(value.get("axisymmetric_meridian_force_n_per_rad"))
        result_meridian_force = float(value.get("result_axisymmetric_meridian_force_n_per_rad"))
        axisymmetric_total = float(value.get("axisymmetric_total_force_n"))
        result_axisymmetric_total = float(value.get("result_axisymmetric_total_force_n"))
    except (TypeError, ValueError):
        return False
    finite_values = (
        depth, result_depth, planar_per_m, result_planar_per_m, planar_total,
        result_planar_total, radius, result_radius, meridian_force,
        result_meridian_force, axisymmetric_total, result_axisymmetric_total,
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "planar_depth_normalization_generation", "radius_normalization_generation",
            "coordinate_normalization_generation", "unit_normalization_generation",
            "mesh_normalization_generation", "result_normalization_generation"))
        and all(math.isfinite(item) for item in finite_values)
        and depth > 0.0 and result_depth == depth
        and result_planar_per_m == planar_per_m
        and math.isclose(planar_total, planar_per_m * depth, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(result_planar_total, planar_total, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and radius > 0.0 and result_radius == radius
        and result_meridian_force == meridian_force
        and math.isclose(axisymmetric_total, 2.0 * math.pi * meridian_force,
                         rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(result_axisymmetric_total, axisymmetric_total,
                         rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(axisymmetric_total, planar_total, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and value.get("radius_measure_convention") == "2*pi*r"
        and value.get("result_radius_measure_convention") == "2*pi*r"
        and value.get("coordinate_convention") == "r_z_right_handed"
        and value.get("result_coordinate_convention") == "r_z_right_handed"
        and value.get("force_unit") == "N_total_3d"
        and value.get("result_force_unit") == "N_total_3d"
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
    )


def _nonlinear_minor_loop_force_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("nonlinear_generation", "")).strip()
    try:
        state = [float(item) for item in value.get("state_point_am", [])]
        result_state = [float(item) for item in value.get("result_state_point_am", [])]
        coenergy = float(value.get("magnetic_coenergy_j"))
        result_coenergy = float(value.get("result_magnetic_coenergy_j"))
        force = [float(item) for item in value.get("force_n", [])]
        result_force = [float(item) for item in value.get("result_force_n", [])]
    except (TypeError, ValueError):
        return False
    branch = str(value.get("bh_branch", ""))
    interpolation = str(value.get("interpolation_rule", ""))
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "branch_nonlinear_generation", "interpolation_nonlinear_generation",
            "state_nonlinear_generation", "coenergy_nonlinear_generation",
            "mesh_nonlinear_generation", "force_nonlinear_generation",
            "result_nonlinear_generation"))
        and branch in {"ascending-minor-loop", "descending-minor-loop"}
        and value.get("result_bh_branch") == branch
        and interpolation == "monotone-cubic-h"
        and value.get("result_interpolation_rule") == interpolation
        and len(state) == 2 and all(math.isfinite(item) for item in state)
        and result_state == state
        and math.isfinite(coenergy) and coenergy >= 0.0
        and math.isclose(result_coenergy, coenergy, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and len(force) == 2 and all(math.isfinite(item) for item in force)
        and result_force == force
        and all(_valid_sha256(value.get(source)) and value.get(result) == value.get(source)
                for source, result in (
                    ("bh_table_sha256", "result_bh_table_sha256"),
                    ("mesh_sha256", "result_mesh_sha256"),
                    ("solution_sha256", "accepted_solution_sha256")))
    )


def _harmonic_eddy_loss_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("eddy_generation", "")).strip()
    try:
        conductivity = float(value.get("conductivity_s_m"))
        result_conductivity = float(value.get("result_conductivity_s_m"))
        permeability = float(value.get("relative_permeability"))
        result_permeability = float(value.get("result_relative_permeability"))
        frequency = float(value.get("frequency_hz"))
        result_frequency = float(value.get("result_frequency_hz"))
        skin_depth = float(value.get("skin_depth_m"))
        result_skin_depth = float(value.get("result_skin_depth_m"))
        elements = float(value.get("minimum_elements_per_skin_depth"))
        result_elements = float(value.get("result_minimum_elements_per_skin_depth"))
        loss = float(value.get("joule_loss_w"))
        result_loss = float(value.get("result_joule_loss_w"))
    except (TypeError, ValueError):
        return False
    expected_skin_depth = math.sqrt(
        2.0 / (2.0 * math.pi * frequency * 4.0e-7 * math.pi * permeability * conductivity)
    ) if conductivity > 0.0 and permeability > 0.0 and frequency > 0.0 else math.nan
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "phasor_eddy_generation", "conductivity_eddy_generation", "skin_eddy_generation",
            "frequency_eddy_generation", "loss_eddy_generation", "mesh_eddy_generation",
            "result_eddy_generation"))
        and value.get("phasor_convention") == "exp(+jwt)"
        and value.get("result_phasor_convention") == "exp(+jwt)"
        and math.isfinite(conductivity) and conductivity > 0.0 and result_conductivity == conductivity
        and math.isfinite(permeability) and permeability > 0.0 and result_permeability == permeability
        and math.isfinite(frequency) and frequency > 0.0 and result_frequency == frequency
        and math.isfinite(skin_depth) and skin_depth > 0.0 and result_skin_depth == skin_depth
        and math.isclose(skin_depth, expected_skin_depth, rel_tol=1.0e-12, abs_tol=1.0e-18)
        and math.isfinite(elements) and elements >= 3.0 and result_elements == elements
        and math.isfinite(loss) and loss >= 0.0 and result_loss == loss
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("loss_result_sha256"))
        and value.get("accepted_loss_result_sha256") == value.get("loss_result_sha256")
    )


def _axisymmetric_aphi_force_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("axisymmetric_generation", "")).strip()
    try:
        regions = [int(item) for item in value.get("region_ids", [])]
        result_regions = [int(item) for item in value.get("result_region_ids", [])]
        displacement = [float(item) for item in value.get("displacement_m", [])]
        coenergy = [float(item) for item in value.get("magnetic_coenergy_j", [])]
        result_coenergy = [float(item) for item in value.get("result_magnetic_coenergy_j", [])]
        force = float(value.get("force_from_energy_n"))
        result_force = float(value.get("result_force_from_energy_n"))
    except (TypeError, ValueError):
        return False
    derivative = (
        (coenergy[2] - coenergy[0]) / (displacement[2] - displacement[0])
        if len(displacement) == len(coenergy) == 3 and displacement[2] != displacement[0]
        else math.nan
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "aphi_axisymmetric_generation", "weight_axisymmetric_generation",
            "region_axisymmetric_generation", "energy_axisymmetric_generation",
            "force_axisymmetric_generation", "mesh_axisymmetric_generation",
            "solution_axisymmetric_generation", "result_axisymmetric_generation"))
        and value.get("formulation") == "Aphi" and value.get("result_formulation") == "Aphi"
        and value.get("radial_weighting") == "2*pi*r"
        and value.get("result_radial_weighting") == "2*pi*r"
        and bool(regions) and all(item > 0 for item in regions)
        and len(set(regions)) == len(regions) and result_regions == regions
        and len(displacement) == 3 and all(math.isfinite(item) for item in displacement)
        and displacement[0] < displacement[1] < displacement[2]
        and len(coenergy) == 3 and all(math.isfinite(item) and item >= 0.0 for item in coenergy)
        and result_coenergy == coenergy
        and value.get("force_axis") == "z" and value.get("result_force_axis") == "z"
        and math.isfinite(force) and result_force == force
        and math.isclose(force, derivative, rel_tol=1.0e-10, abs_tol=1.0e-12)
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("solution_sha256"))
        and value.get("accepted_solution_sha256") == value.get("solution_sha256")
    )


def _permanent_magnet_operating_point_identity_ok(value):
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("magnet_generation", "")).strip()
    try:
        recoil_mu = float(value.get("recoil_relative_permeability"))
        result_recoil_mu = float(value.get("result_recoil_relative_permeability"))
        remanence = float(value.get("remanence_t"))
        result_remanence = float(value.get("result_remanence_t"))
        temperature = float(value.get("magnet_temperature_c"))
        result_temperature = float(value.get("result_magnet_temperature_c"))
        operating_point = [float(item) for item in value.get("operating_point_bh", [])]
        result_operating_point = [float(item) for item in value.get("result_operating_point_bh", [])]
        margin = float(value.get("demag_margin_a_m"))
        result_margin = float(value.get("result_demag_margin_a_m"))
        force = [float(item) for item in value.get("force_n", [])]
        result_force = [float(item) for item in value.get("result_force_n", [])]
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "recoil_magnet_generation", "temperature_magnet_generation",
            "operating_point_magnet_generation", "frame_magnet_generation",
            "demag_magnet_generation", "force_magnet_generation",
            "mesh_magnet_generation", "result_magnet_generation"))
        and math.isfinite(recoil_mu) and recoil_mu > 0.0 and result_recoil_mu == recoil_mu
        and math.isfinite(remanence) and remanence > 0.0 and result_remanence == remanence
        and math.isfinite(temperature) and result_temperature == temperature
        and len(operating_point) == 2 and all(math.isfinite(item) for item in operating_point)
        and operating_point[0] > 0.0 and result_operating_point == operating_point
        and math.isfinite(margin) and margin > 0.0 and result_margin == margin
        and _valid_sha256(value.get("magnetization_frame_sha256"))
        and value.get("result_magnetization_frame_sha256") == value.get("magnetization_frame_sha256")
        and len(force) == 2 and all(math.isfinite(item) for item in force) and result_force == force
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


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
    axisymmetric_weighted_stress_force_mask_mesh_identity_ok = True
    lorentz_current_density_orientation_identity_ok = True
    circuit_current_phasor_convention_identity_ok = True
    incremental_permeability_operating_point_identity_ok = True
    weighted_stress_air_mask_nodal_weight_identity_ok = True
    sliding_band_periodic_angle_rotor_position_identity_ok = True
    coenergy_torque_angle_difference_remesh_state_identity_ok = True
    axisymmetric_henrotte_hodge_radius_weight_coordinate_identity_ok = True
    weighted_stress_force_mask_material_mesh_generation_identity_ok = True
    harmonic_loss_phase_frequency_lamination_generation_identity_ok = True
    axisymmetric_force_energy_normalization_generation_identity_ok = True
    nonlinear_incremental_force_branch_generation_identity_ok = True
    nonlinear_force_operating_point_identity_ok = True
    sliding_band_harmonic_torque_identity_ok = True
    weighted_stress_energy_derivative_identity_ok = True
    axisymmetric_revolved_energy_force_identity_ok = True
    nonlinear_bh_incremental_force_identity_ok = True
    open_boundary_decay_multipole_identity_ok = True
    weighted_stress_tensor_closure_identity_ok = True
    axisymmetric_planar_normalization_identity_ok = True
    nonlinear_minor_loop_force_identity_ok = True
    harmonic_eddy_loss_identity_ok = True
    axisymmetric_aphi_force_identity_ok = True
    permanent_magnet_operating_point_identity_ok = True
    airgap_stress_harmonic_torque_identity_ok = True
    laminated_core_loss_identity_ok = True
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
        axisymmetric_weighted_stress_force_mask_mesh_identity_ok = False
        lorentz_current_density_orientation_identity_ok = False
        circuit_current_phasor_convention_identity_ok = False
        incremental_permeability_operating_point_identity_ok = False
        weighted_stress_air_mask_nodal_weight_identity_ok = False
        sliding_band_periodic_angle_rotor_position_identity_ok = False
        coenergy_torque_angle_difference_remesh_state_identity_ok = False
        axisymmetric_henrotte_hodge_radius_weight_coordinate_identity_ok = False
        weighted_stress_force_mask_material_mesh_generation_identity_ok = False
        harmonic_loss_phase_frequency_lamination_generation_identity_ok = False
        axisymmetric_force_energy_normalization_generation_identity_ok = False
        nonlinear_incremental_force_branch_generation_identity_ok = False
        nonlinear_force_operating_point_identity_ok = False
        sliding_band_harmonic_torque_identity_ok = False
        weighted_stress_energy_derivative_identity_ok = False
        axisymmetric_revolved_energy_force_identity_ok = False
        nonlinear_bh_incremental_force_identity_ok = False
        open_boundary_decay_multipole_identity_ok = False
        weighted_stress_tensor_closure_identity_ok = False
        axisymmetric_planar_normalization_identity_ok = False
        nonlinear_minor_loop_force_identity_ok = False
        harmonic_eddy_loss_identity_ok = False
        axisymmetric_aphi_force_identity_ok = False
        permanent_magnet_operating_point_identity_ok = False
        airgap_stress_harmonic_torque_identity_ok = False
        laminated_core_loss_identity_ok = False
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

        axisymmetric_mask = artifact_identity.get(
            "axisymmetric_weighted_stress_force_mask_mesh_identity"
        )
        if axisymmetric_mask is not None:
            current_mask = artifact_identity.get("weighted_stress_mask_mesh_identity")
            force_measure = artifact_identity.get("axisymmetric_force_measure_identity")
            radius_identity = artifact_identity.get(
                "axisymmetric_force_radius_jacobian_coordinate_identity"
            )
            mesh_generation = (
                axisymmetric_mask.get("force_mesh_generation")
                if isinstance(axisymmetric_mask, dict)
                else None
            )
            mask_solution_generation = (
                axisymmetric_mask.get("mask_solution_generation")
                if isinstance(axisymmetric_mask, dict)
                else None
            )
            mask_digest = str(
                axisymmetric_mask.get("mask_field_sha256", "")
                if isinstance(axisymmetric_mask, dict)
                else ""
            ).lower()
            axisymmetric_weighted_stress_force_mask_mesh_identity_ok = (
                isinstance(axisymmetric_mask, dict)
                and isinstance(current_mask, dict)
                and isinstance(force_measure, dict)
                and isinstance(radius_identity, dict)
                and bool(mesh_generation)
                and mesh_generation
                == current_mask.get("active_air_mesh_generation")
                and axisymmetric_mask.get("mask_solve_mesh_generation")
                == mesh_generation
                and axisymmetric_mask.get(
                    "weighted_stress_integral_mesh_generation"
                )
                == mesh_generation
                and bool(mask_solution_generation)
                and axisymmetric_mask.get("force_mask_solution_generation")
                == mask_solution_generation
                and axisymmetric_mask.get("axisymmetric_measure")
                == "2*pi*r*dr*dz"
                and axisymmetric_mask.get("force_integral_measure")
                == axisymmetric_mask.get("axisymmetric_measure")
                and axisymmetric_mask.get("axisymmetric_measure")
                == force_measure.get("integration_measure")
                and axisymmetric_mask.get("axisymmetric_measure")
                == radius_identity.get("integration_measure")
                and len(mask_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in mask_digest
                )
                and str(
                    axisymmetric_mask.get("force_mask_field_sha256", "")
                ).lower()
                == mask_digest
                and mask_digest
                == str(current_mask.get("weighted_mask_sha256", "")).lower()
            )

        lorentz_orientation = artifact_identity.get(
            "lorentz_force_current_density_out_of_plane_orientation_identity"
        )
        if lorentz_orientation is not None:
            solve_generation = (
                lorentz_orientation.get("field_solve_generation")
                if isinstance(lorentz_orientation, dict)
                else None
            )
            current_digest = str(
                lorentz_orientation.get("current_density_sha256", "")
                if isinstance(lorentz_orientation, dict)
                else ""
            ).lower()
            orientation_sign = (
                lorentz_orientation.get("orientation_sign")
                if isinstance(lorentz_orientation, dict)
                else None
            )
            magnetic_flux_digest = str(
                lorentz_orientation.get("magnetic_flux_density_sha256", "")
                if isinstance(lorentz_orientation, dict)
                else ""
            ).lower()
            force_digest = str(
                lorentz_orientation.get("lorentz_force_sha256", "")
                if isinstance(lorentz_orientation, dict)
                else ""
            ).lower()
            lorentz_current_density_orientation_identity_ok = (
                isinstance(lorentz_orientation, dict)
                and bool(solve_generation)
                and lorentz_orientation.get(
                    "current_density_field_solve_generation"
                )
                == solve_generation
                and lorentz_orientation.get("lorentz_force_field_solve_generation")
                == solve_generation
                and lorentz_orientation.get("coordinate_plane") == "xy"
                and lorentz_orientation.get("out_of_plane_axis") == "+z"
                and lorentz_orientation.get("current_density_component") == "Jz"
                and lorentz_orientation.get("current_density_positive_axis") == "+z"
                and lorentz_orientation.get("lorentz_cross_product") == "J_cross_B"
                and lorentz_orientation.get("coordinate_handedness")
                == "right_handed"
                and lorentz_orientation.get("magnetic_flux_density_frame")
                == "global_cartesian"
                and lorentz_orientation.get("force_component_frame")
                == "global_cartesian"
                and lorentz_orientation.get("force_component_formula")
                == ["Fx=-Jz*By", "Fy=Jz*Bx", "Fz=0"]
                and type(orientation_sign) is int
                and orientation_sign == 1
                and lorentz_orientation.get("force_orientation_sign")
                == orientation_sign
                and len(current_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in current_digest
                )
                and str(
                    lorentz_orientation.get("force_current_density_sha256", "")
                ).lower()
                == current_digest
                and len(magnetic_flux_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in magnetic_flux_digest
                )
                and str(
                    lorentz_orientation.get(
                        "force_magnetic_flux_density_sha256", ""
                    )
                ).lower()
                == magnetic_flux_digest
                and len(force_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in force_digest
                )
                and str(
                    lorentz_orientation.get("result_force_sha256", "")
                ).lower()
                == force_digest
            )

        circuit_phasor = artifact_identity.get(
            "circuit_current_peak_rms_phasor_phase_generation_identity"
        )
        if circuit_phasor is not None:
            circuit_phasor = circuit_phasor if isinstance(circuit_phasor, dict) else {}
            circuit_generation = str(
                circuit_phasor.get("circuit_generation", "")
            ).strip()
            phasor_generation = str(
                circuit_phasor.get("phasor_generation", "")
            ).strip()
            phasor_digest = str(
                circuit_phasor.get("current_phasor_sha256", "")
            ).lower()
            try:
                current_amplitude = float(circuit_phasor.get("current_amplitude_a"))
                response_amplitude = float(
                    circuit_phasor.get("field_response_current_amplitude_a")
                )
                phase_value = float(circuit_phasor.get("phase_value"))
                response_phase = float(
                    circuit_phasor.get("field_response_phase_value")
                )
            except (TypeError, ValueError):
                current_amplitude = response_amplitude = math.nan
                phase_value = response_phase = math.nan
            circuit_current_phasor_convention_identity_ok = (
                bool(circuit_generation)
                and circuit_phasor.get("field_response_circuit_generation")
                == circuit_generation
                and bool(phasor_generation)
                and circuit_phasor.get("circuit_current_phasor_generation")
                == phasor_generation
                and circuit_phasor.get("field_response_phasor_generation")
                == phasor_generation
                and circuit_phasor.get("current_amplitude_convention")
                in {"peak", "rms"}
                and circuit_phasor.get("field_response_amplitude_convention")
                == circuit_phasor.get("current_amplitude_convention")
                and math.isfinite(current_amplitude)
                and current_amplitude > 0.0
                and response_amplitude == current_amplitude
                and circuit_phasor.get("phase_unit") in {"deg", "rad"}
                and circuit_phasor.get("field_response_phase_unit")
                == circuit_phasor.get("phase_unit")
                and math.isfinite(phase_value)
                and response_phase == phase_value
                and len(phasor_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in phasor_digest
                )
                and str(
                    circuit_phasor.get("field_response_current_phasor_sha256", "")
                ).lower()
                == phasor_digest
            )

        incremental_mu = artifact_identity.get(
            "incremental_permeability_bh_operating_point_iteration_identity"
        )
        if incremental_mu is not None:
            incremental_mu = incremental_mu if isinstance(incremental_mu, dict) else {}
            solve_generation = str(
                incremental_mu.get("nonlinear_solve_generation", "")
            ).strip()
            iteration_generation = str(
                incremental_mu.get("operating_point_iteration_generation", "")
            ).strip()
            operating_digest = str(
                incremental_mu.get("operating_point_sha256", "")
            ).lower()
            try:
                operating_iteration = int(
                    incremental_mu.get("operating_point_iteration", -1)
                )
                incremental_iteration = int(
                    incremental_mu.get("incremental_permeability_iteration", -1)
                )
                force_iteration = int(
                    incremental_mu.get("force_sensitivity_iteration", -1)
                )
                operating_b = [
                    float(value) for value in incremental_mu.get("operating_point_b_t", [])
                ]
                incremental_b = [
                    float(value)
                    for value in incremental_mu.get("incremental_permeability_b_t", [])
                ]
            except (TypeError, ValueError):
                operating_iteration = incremental_iteration = force_iteration = -1
                operating_b = []
                incremental_b = []
            incremental_permeability_operating_point_identity_ok = (
                bool(solve_generation)
                and incremental_mu.get("operating_point_solve_generation")
                == solve_generation
                and incremental_mu.get("incremental_permeability_solve_generation")
                == solve_generation
                and incremental_mu.get("force_sensitivity_solve_generation")
                == solve_generation
                and bool(iteration_generation)
                and incremental_mu.get(
                    "incremental_permeability_iteration_generation"
                )
                == iteration_generation
                and incremental_mu.get("force_sensitivity_iteration_generation")
                == iteration_generation
                and operating_iteration >= 0
                and incremental_iteration == operating_iteration
                and force_iteration == operating_iteration
                and bool(operating_b)
                and all(math.isfinite(value) for value in operating_b)
                and incremental_b == operating_b
                and len(operating_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in operating_digest
                )
                and str(
                    incremental_mu.get(
                        "incremental_permeability_operating_point_sha256", ""
                    )
                ).lower()
                == operating_digest
                and str(
                    incremental_mu.get("force_sensitivity_operating_point_sha256", "")
                ).lower()
                == operating_digest
            )

        air_mask_weight = artifact_identity.get(
            "weighted_stress_air_mask_nodal_weight_mesh_generation_identity"
        )
        if air_mask_weight is not None:
            air_mask_weight = air_mask_weight if isinstance(air_mask_weight, dict) else {}
            mesh_generation = str(
                air_mask_weight.get("field_mesh_generation", "")
            ).strip()
            mask_digest = str(air_mask_weight.get("air_mask_sha256", "")).lower()
            weight_digest = str(
                air_mask_weight.get("nodal_weight_sha256", "")
            ).lower()
            try:
                air_regions = [
                    int(value) for value in air_mask_weight.get("air_region_ids", [])
                ]
                mask_regions = [
                    int(value)
                    for value in air_mask_weight.get("mask_air_region_ids", [])
                ]
                weight_nodes = [
                    int(value)
                    for value in air_mask_weight.get("nodal_weight_node_ids", [])
                ]
                force_nodes = [
                    int(value)
                    for value in air_mask_weight.get("force_weight_node_ids", [])
                ]
            except (TypeError, ValueError):
                air_regions = mask_regions = weight_nodes = force_nodes = []
            weighted_stress_air_mask_nodal_weight_identity_ok = (
                bool(mesh_generation)
                and air_mask_weight.get("air_mask_mesh_generation")
                == mesh_generation
                and air_mask_weight.get("nodal_weight_mesh_generation")
                == mesh_generation
                and air_mask_weight.get("force_integral_mesh_generation")
                == mesh_generation
                and bool(air_regions)
                and len(set(air_regions)) == len(air_regions)
                and mask_regions == air_regions
                and bool(weight_nodes)
                and len(set(weight_nodes)) == len(weight_nodes)
                and force_nodes == weight_nodes
                and len(mask_digest) == 64
                and all(character in "0123456789abcdef" for character in mask_digest)
                and str(air_mask_weight.get("force_air_mask_sha256", "")).lower()
                == mask_digest
                and len(weight_digest) == 64
                and all(character in "0123456789abcdef" for character in weight_digest)
                and str(
                    air_mask_weight.get("force_nodal_weight_sha256", "")
                ).lower()
                == weight_digest
            )

        sliding_band = artifact_identity.get(
            "sliding_band_periodic_angle_rotor_position_generation_identity"
        )
        if sliding_band is not None:
            sliding_band = sliding_band if isinstance(sliding_band, dict) else {}
            rotor_generation = str(
                sliding_band.get("rotor_position_generation", "")
            ).strip()
            angle_generation = str(
                sliding_band.get("periodic_angle_generation", "")
            ).strip()
            map_digest = str(
                sliding_band.get("sliding_band_map_sha256", "")
            ).lower()
            try:
                rotor_angle = float(sliding_band.get("rotor_angle_deg"))
                mapped_rotor_angle = float(
                    sliding_band.get("sliding_band_rotor_angle_deg")
                )
                angle_pairs = [
                    [float(value) for value in row]
                    for row in sliding_band.get("periodic_angle_pairs_deg", [])
                ]
                torque_pairs = [
                    [float(value) for value in row]
                    for row in sliding_band.get(
                        "torque_periodic_angle_pairs_deg", []
                    )
                ]
            except (TypeError, ValueError):
                rotor_angle = mapped_rotor_angle = math.nan
                angle_pairs = torque_pairs = []
            sliding_band_periodic_angle_rotor_position_identity_ok = (
                bool(rotor_generation)
                and sliding_band.get("sliding_band_rotor_position_generation")
                == rotor_generation
                and sliding_band.get("torque_rotor_position_generation")
                == rotor_generation
                and bool(angle_generation)
                and sliding_band.get("sliding_band_periodic_angle_generation")
                == angle_generation
                and sliding_band.get("torque_periodic_angle_generation")
                == angle_generation
                and math.isfinite(rotor_angle)
                and mapped_rotor_angle == rotor_angle
                and bool(angle_pairs)
                and all(
                    len(row) == 2 and all(math.isfinite(value) for value in row)
                    for row in angle_pairs
                )
                and torque_pairs == angle_pairs
                and len(map_digest) == 64
                and all(character in "0123456789abcdef" for character in map_digest)
                and str(
                    sliding_band.get("torque_sliding_band_map_sha256", "")
                ).lower()
                == map_digest
            )

        torque_states = artifact_identity.get(
            "coenergy_torque_angle_difference_remesh_state_generation_identity"
        )
        if torque_states is not None:
            torque_states = torque_states if isinstance(torque_states, dict) else {}
            torque_generation = str(
                torque_states.get("torque_generation", "")
            ).strip()
            table_generation = str(
                torque_states.get("state_table_generation", "")
            ).strip()
            angle_generation = str(
                torque_states.get("angle_spacing_generation", "")
            ).strip()
            solve_generations = [
                str(value)
                for value in torque_states.get("coenergy_solve_generations", [])
            ]
            mesh_generations = [
                str(value)
                for value in torque_states.get("mesh_remap_solve_generations", [])
            ]
            excitation_generations = [
                str(value)
                for value in torque_states.get("excitation_solve_generations", [])
            ]
            angle_state_generations = [
                str(value)
                for value in torque_states.get("angle_state_solve_generations", [])
            ]
            try:
                angles = [
                    float(value) for value in torque_states.get("angles_deg", [])
                ]
                derivative_angles = [
                    float(value)
                    for value in torque_states.get("derivative_angles_deg", [])
                ]
            except (TypeError, ValueError):
                angles = derivative_angles = []
            digest = str(
                torque_states.get("coenergy_state_table_sha256", "")
            ).lower()
            coenergy_torque_angle_difference_remesh_state_identity_ok = (
                bool(torque_generation)
                and torque_states.get("derivative_torque_generation")
                == torque_generation
                and bool(table_generation)
                and torque_states.get("derivative_state_table_generation")
                == table_generation
                and len(solve_generations) >= 3
                and all(solve_generations)
                and len(set(solve_generations)) == len(solve_generations)
                and mesh_generations == solve_generations
                and excitation_generations == solve_generations
                and angle_state_generations == solve_generations
                and len(angles) == len(solve_generations)
                and all(math.isfinite(value) for value in angles)
                and all(left < right for left, right in zip(angles, angles[1:]))
                and derivative_angles == angles
                and bool(angle_generation)
                and torque_states.get("derivative_angle_spacing_generation")
                == angle_generation
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    torque_states.get("derivative_state_table_sha256", "")
                ).lower()
                == digest
            )

        hodge_radius = artifact_identity.get(
            "axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_identity"
        )
        if hodge_radius is not None:
            hodge_radius = hodge_radius if isinstance(hodge_radius, dict) else {}
            mesh_generation = str(
                hodge_radius.get("mesh_geometry_generation", "")
            ).strip()
            try:
                node_ids = [int(value) for value in hodge_radius.get("node_ids", [])]
                field_node_ids = [
                    int(value) for value in hodge_radius.get("field_node_ids", [])
                ]
                radius = [float(value) for value in hodge_radius.get("radius_m", [])]
                hodge_weights = [
                    float(value)
                    for value in hodge_radius.get("hodge_radius_weight_m", [])
                ]
                cylindrical_radius = [
                    float(value)
                    for value in hodge_radius.get("cylindrical_r_coordinate_m", [])
                ]
            except (TypeError, ValueError):
                node_ids = field_node_ids = []
                radius = hodge_weights = cylindrical_radius = []
            weight_digest = str(
                hodge_radius.get("radius_weight_table_sha256", "")
            ).lower()
            coordinate_digest = str(
                hodge_radius.get("coordinate_table_sha256", "")
            ).lower()
            axisymmetric_henrotte_hodge_radius_weight_coordinate_identity_ok = (
                bool(mesh_generation)
                and hodge_radius.get("field_mesh_geometry_generation")
                == mesh_generation
                and hodge_radius.get("radius_weight_mesh_geometry_generation")
                == mesh_generation
                and hodge_radius.get(
                    "cylindrical_coordinate_mesh_geometry_generation"
                )
                == mesh_generation
                and bool(node_ids)
                and len(set(node_ids)) == len(node_ids)
                and field_node_ids == node_ids
                and len(radius) == len(node_ids)
                and all(math.isfinite(value) and value >= 0.0 for value in radius)
                and hodge_weights == radius
                and cylindrical_radius == radius
                and len(weight_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in weight_digest
                )
                and str(
                    hodge_radius.get("force_radius_weight_table_sha256", "")
                ).lower()
                == weight_digest
                and len(coordinate_digest) == 64
                and all(
                    character in "0123456789abcdef"
                    for character in coordinate_digest
                )
                and str(
                    hodge_radius.get("field_coordinate_table_sha256", "")
                ).lower()
                == coordinate_digest
            )

        weighted_force = artifact_identity.get(
            "weighted_stress_force_mask_material_mesh_generation_identity"
        )
        if weighted_force is not None:
            weighted_force = weighted_force if isinstance(weighted_force, dict) else {}
            solve_generation = str(
                weighted_force.get("solve_generation", "")
            ).strip()
            try:
                body_groups = [
                    int(value) for value in weighted_force.get("body_group_ids", [])
                ]
                mask_groups = [
                    int(value)
                    for value in weighted_force.get("mask_body_group_ids", [])
                ]
                force_values = [
                    float(value)
                    for value in weighted_force.get("weighted_force_n", [])
                ]
                reported_force = [
                    float(value)
                    for value in weighted_force.get("reported_force_n", [])
                ]
            except (TypeError, ValueError):
                body_groups = mask_groups = []
                force_values = reported_force = []
            material_labels = [
                str(value).strip()
                for value in weighted_force.get("material_labels", [])
            ]
            resolved_labels = [
                str(value).strip()
                for value in weighted_force.get("resolved_material_labels", [])
            ]
            mask_digest = str(
                weighted_force.get("mask_table_sha256", "")
            ).lower()
            field_digest = str(
                weighted_force.get("mesh_field_sha256", "")
            ).lower()
            weighted_stress_force_mask_material_mesh_generation_identity_ok = (
                bool(solve_generation)
                and all(
                    weighted_force.get(key) == solve_generation
                    for key in (
                        "weighted_stress_mask_solve_generation",
                        "material_label_solve_generation",
                        "mesh_field_solve_generation",
                        "force_integral_solve_generation",
                    )
                )
                and bool(body_groups)
                and len(set(body_groups)) == len(body_groups)
                and mask_groups == body_groups
                and len(material_labels) == len(body_groups)
                and all(material_labels)
                and resolved_labels == material_labels
                and bool(force_values)
                and all(math.isfinite(value) for value in force_values)
                and reported_force == force_values
                and len(mask_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in mask_digest
                )
                and str(weighted_force.get("force_mask_table_sha256", "")).lower()
                == mask_digest
                and len(field_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in field_digest
                )
                and str(weighted_force.get("force_mesh_field_sha256", "")).lower()
                == field_digest
            )

        harmonic_loss = artifact_identity.get(
            "harmonic_loss_phase_frequency_lamination_generation_identity"
        )
        if harmonic_loss is not None:
            harmonic_loss = harmonic_loss if isinstance(harmonic_loss, dict) else {}
            generation = str(
                harmonic_loss.get("analysis_generation", "")
            ).strip()
            phase = str(harmonic_loss.get("phase_convention", "")).strip()
            try:
                frequency = float(harmonic_loss.get("frequency_hz"))
                loss_frequency = float(harmonic_loss.get("loss_frequency_hz"))
                coefficients = [
                    [float(value) for value in row]
                    for row in harmonic_loss.get("material_loss_coefficients", [])
                ]
                loss_coefficients = [
                    [float(value) for value in row]
                    for row in harmonic_loss.get("loss_material_coefficients", [])
                ]
            except (TypeError, ValueError):
                frequency = loss_frequency = math.nan
                coefficients = loss_coefficients = []
            orientations = [
                str(value).strip()
                for value in harmonic_loss.get("lamination_orientations", [])
            ]
            loss_orientations = [
                str(value).strip()
                for value in harmonic_loss.get("loss_lamination_orientations", [])
            ]
            table_digest = str(
                harmonic_loss.get("material_loss_table_sha256", "")
            ).lower()
            harmonic_loss_phase_frequency_lamination_generation_identity_ok = (
                bool(generation)
                and all(
                    harmonic_loss.get(key) == generation
                    for key in (
                        "phase_convention_analysis_generation",
                        "frequency_analysis_generation",
                        "lamination_analysis_generation",
                        "material_loss_analysis_generation",
                        "loss_result_analysis_generation",
                    )
                )
                and phase in {"exp(+jwt)", "exp(-jwt)"}
                and harmonic_loss.get("loss_phase_convention") == phase
                and math.isfinite(frequency)
                and frequency > 0.0
                and math.isclose(
                    loss_frequency, frequency, rel_tol=1.0e-12, abs_tol=1.0e-15
                )
                and bool(orientations)
                and all(orientations)
                and loss_orientations == orientations
                and len(coefficients) == len(orientations)
                and all(
                    row and all(math.isfinite(value) for value in row)
                    for row in coefficients
                )
                and loss_coefficients == coefficients
                and len(table_digest) == 64
                and all(
                    character in "0123456789abcdef" for character in table_digest
                )
                and str(harmonic_loss.get("loss_material_table_sha256", "")).lower()
                == table_digest
            )

        axisym_normalization = artifact_identity.get(
            "axisymmetric_force_energy_measure_depth_coordinate_generation_identity"
        )
        if axisym_normalization is not None:
            axisym_normalization = (
                axisym_normalization if isinstance(axisym_normalization, dict) else {}
            )
            generation = str(
                axisym_normalization.get("solve_generation", "")
            ).strip()
            problem_type = str(
                axisym_normalization.get("problem_type", "")
            ).strip()
            measure = str(
                axisym_normalization.get("measure_convention", "")
            ).strip()
            coordinate = str(
                axisym_normalization.get("coordinate_convention", "")
            ).strip()
            try:
                depth = float(axisym_normalization.get("planar_depth_m"))
                result_depth = float(
                    axisym_normalization.get("result_planar_depth_m")
                )
                values = [
                    float(value)
                    for value in axisym_normalization.get("force_energy_values", [])
                ]
                reported_values = [
                    float(value)
                    for value in axisym_normalization.get(
                        "reported_force_energy_values", []
                    )
                ]
            except (TypeError, ValueError):
                depth = result_depth = math.nan
                values = reported_values = []
            digest = str(
                axisym_normalization.get("normalization_table_sha256", "")
            ).lower()
            axisymmetric_force_energy_normalization_generation_identity_ok = (
                bool(generation)
                and all(
                    axisym_normalization.get(key) == generation
                    for key in (
                        "measure_solve_generation",
                        "depth_solve_generation",
                        "coordinate_solve_generation",
                        "result_solve_generation",
                    )
                )
                and problem_type in {"planar", "axisymmetric"}
                and axisym_normalization.get("result_problem_type") == problem_type
                and measure in {"planar_depth", "2*pi*r"}
                and axisym_normalization.get("result_measure_convention") == measure
                and ((problem_type == "axisymmetric" and measure == "2*pi*r") or (problem_type == "planar" and measure == "planar_depth"))
                and math.isfinite(depth)
                and depth > 0.0
                and math.isclose(result_depth, depth, rel_tol=0.0, abs_tol=1.0e-15)
                and coordinate in {"x_y", "r_z"}
                and axisym_normalization.get("result_coordinate_convention")
                == coordinate
                and ((problem_type == "axisymmetric" and coordinate == "r_z") or (problem_type == "planar" and coordinate == "x_y"))
                and bool(values)
                and all(math.isfinite(value) for value in values)
                and reported_values == values
                and len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest)
                and str(
                    axisym_normalization.get(
                        "result_normalization_table_sha256", ""
                    )
                ).lower()
                == digest
            )

        incremental_force = artifact_identity.get(
            "nonlinear_incremental_mu_force_branch_perturbation_generation_identity"
        )
        if incremental_force is not None:
            incremental_force = (
                incremental_force if isinstance(incremental_force, dict) else {}
            )
            generation = str(
                incremental_force.get("operating_point_generation", "")
            ).strip()
            branch_id = str(incremental_force.get("branch_id", "")).strip()
            try:
                perturbation = float(
                    incremental_force.get("perturbation_current_a")
                )
                force_perturbation = float(
                    incremental_force.get("force_perturbation_current_a")
                )
                force_values = [
                    float(value)
                    for value in incremental_force.get("incremental_force_n", [])
                ]
                reported_force = [
                    float(value)
                    for value in incremental_force.get(
                        "reported_incremental_force_n", []
                    )
                ]
            except (TypeError, ValueError):
                perturbation = force_perturbation = math.nan
                force_values = reported_force = []
            mu_digest = str(
                incremental_force.get("differential_mu_sha256", "")
            ).lower()
            state_digest = str(
                incremental_force.get("incremental_state_sha256", "")
            ).lower()
            nonlinear_incremental_force_branch_generation_identity_ok = (
                bool(generation)
                and all(
                    incremental_force.get(key) == generation
                    for key in (
                        "branch_operating_point_generation",
                        "differential_mu_operating_point_generation",
                        "perturbation_operating_point_generation",
                        "force_operating_point_generation",
                    )
                )
                and bool(branch_id)
                and incremental_force.get("force_branch_id") == branch_id
                and math.isfinite(perturbation)
                and perturbation != 0.0
                and math.isclose(
                    force_perturbation,
                    perturbation,
                    rel_tol=0.0,
                    abs_tol=1.0e-15,
                )
                and len(mu_digest) == 64
                and all(character in "0123456789abcdef" for character in mu_digest)
                and str(
                    incremental_force.get("force_differential_mu_sha256", "")
                ).lower()
                == mu_digest
                and bool(force_values)
                and all(math.isfinite(value) for value in force_values)
                and reported_force == force_values
                and len(state_digest) == 64
                and all(character in "0123456789abcdef" for character in state_digest)
                and str(
                    incremental_force.get("force_incremental_state_sha256", "")
                ).lower()
                == state_digest
            )

        nonlinear_force_operating_point_identity_ok = (
            _nonlinear_force_operating_point_identity_ok(
                artifact_identity.get(
                    "nonlinear_bh_branch_operating_point_force_mesh_generation_identity"
                )
            )
        )
        sliding_band_harmonic_torque_identity_ok = (
            _sliding_band_harmonic_torque_identity_ok(
                artifact_identity.get(
                    "sliding_band_angle_mesh_harmonic_torque_generation_identity"
                )
            )
        )
        weighted_stress_energy_derivative_identity_ok = (
            _weighted_stress_energy_derivative_identity_ok(
                artifact_identity.get(
                    "weighted_stress_energy_derivative_force_mesh_frame_unit_generation_identity"
                )
            )
        )
        axisymmetric_revolved_energy_force_identity_ok = (
            _axisymmetric_revolved_energy_force_identity_ok(
                artifact_identity.get(
                    "axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_identity"
                )
            )
        )
        nonlinear_bh_incremental_force_identity_ok = (
            _nonlinear_bh_incremental_force_identity_ok(
                artifact_identity.get(
                    "nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_generation_identity"
                )
            )
        )
        open_boundary_decay_multipole_identity_ok = (
            _open_boundary_decay_multipole_identity_ok(
                artifact_identity.get(
                    "open_boundary_domain_decay_multipole_moment_material_generation_identity"
                )
            )
        )
        weighted_stress_tensor_closure_identity_ok = (
            _weighted_stress_tensor_closure_identity_ok(
                artifact_identity.get(
                    "weighted_stress_tensor_mask_mesh_region_force_torque_energy_generation_identity"
                )
            )
        )
        axisymmetric_planar_normalization_identity_ok = (
            _axisymmetric_planar_normalization_identity_ok(
                artifact_identity.get(
                    "axisymmetric_planar_depth_two_pi_r_force_normalization_coordinate_unit_generation_identity"
                )
            )
        )
        nonlinear_minor_loop_force_identity_ok = _nonlinear_minor_loop_force_identity_ok(
            artifact_identity.get(
                "nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_generation_identity"
            )
        )
        harmonic_eddy_loss_identity_ok = _harmonic_eddy_loss_identity_ok(
            artifact_identity.get(
                "harmonic_eddy_phasor_conductivity_skin_depth_frequency_loss_mesh_generation_identity"
            )
        )
        axisymmetric_aphi_force_identity_ok = _axisymmetric_aphi_force_identity_ok(
            artifact_identity.get(
                "axisymmetric_aphi_radial_weight_region_energy_force_mesh_solution_generation_identity"
            )
        )
        permanent_magnet_operating_point_identity_ok = (
            _permanent_magnet_operating_point_identity_ok(
                artifact_identity.get(
                    "permanent_magnet_recoil_temperature_operating_point_demag_force_generation_identity"
                )
            )
        )
        airgap_stress_harmonic_torque_identity_ok = (
            _airgap_stress_harmonic_torque_identity_ok(
                artifact_identity.get(
                    "airgap_stress_harmonic_sector_periodicity_origin_sampling_alias_radius_torque_generation_identity"
                )
            )
        )
        laminated_core_loss_identity_ok = _laminated_core_loss_identity_ok(
            artifact_identity.get(
                "laminated_core_hysteresis_eddy_excess_frequency_flux_lamination_volume_result_generation_identity"
            )
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
        "axisymmetric_weighted_stress_force_uses_current_mask_mesh": (
            axisymmetric_weighted_stress_force_mask_mesh_identity_ok
        ),
        "planar_lorentz_force_uses_current_jz_out_of_plane_orientation": (
            lorentz_current_density_orientation_identity_ok
        ),
        "circuit_current_and_field_response_share_phasor_convention_generation": (
            circuit_current_phasor_convention_identity_ok
        ),
        "incremental_permeability_and_force_use_current_bh_operating_point": (
            incremental_permeability_operating_point_identity_ok
        ),
        "weighted_stress_air_mask_and_nodal_weights_use_current_mesh": (
            weighted_stress_air_mask_nodal_weight_identity_ok
        ),
        "sliding_band_angles_and_torque_use_current_rotor_position": (
            sliding_band_periodic_angle_rotor_position_identity_ok
        ),
        "coenergy_torque_states_share_remesh_excitation_and_angle_generations": (
            coenergy_torque_angle_difference_remesh_state_identity_ok
        ),
        "axisymmetric_hodge_radius_weights_use_current_mesh_coordinates": (
            axisymmetric_henrotte_hodge_radius_weight_coordinate_identity_ok
        ),
        "weighted_stress_force_uses_current_mask_materials_and_mesh_field": (
            weighted_stress_force_mask_material_mesh_generation_identity_ok
        ),
        "harmonic_loss_uses_current_phase_frequency_lamination_and_material_data": (
            harmonic_loss_phase_frequency_lamination_generation_identity_ok
        ),
        "axisymmetric_force_energy_uses_current_measure_depth_and_coordinates": (
            axisymmetric_force_energy_normalization_generation_identity_ok
        ),
        "nonlinear_incremental_force_uses_current_branch_mu_and_perturbation": (
            nonlinear_incremental_force_branch_generation_identity_ok
        ),
        "nonlinear_force_uses_current_bh_branch_operating_point_mu_and_mesh": (
            nonlinear_force_operating_point_identity_ok
        ),
        "sliding_band_harmonics_use_current_angles_mesh_currents_and_samples": (
            sliding_band_harmonic_torque_identity_ok
        ),
        "weighted_stress_and_energy_derivative_share_mesh_frame_units_and_generation": (
            weighted_stress_energy_derivative_identity_ok
        ),
        "axisymmetric_revolved_energy_force_share_2pir_hodge_field_material_and_mesh": (
            axisymmetric_revolved_energy_force_identity_ok
        ),
        "nonlinear_bh_incremental_force_uses_current_branch_material_mesh_energy_and_coenergy": (
            nonlinear_bh_incremental_force_identity_ok
        ),
        "open_boundary_uses_current_domain_decay_multipole_material_and_mesh": (
            open_boundary_decay_multipole_identity_ok
        ),
        "weighted_stress_tensor_uses_current_mask_mesh_region_force_torque_and_energy": (
            weighted_stress_tensor_closure_identity_ok
        ),
        "axisymmetric_and_planar_force_share_depth_two_pi_r_coordinates_units_and_mesh": (
            axisymmetric_planar_normalization_identity_ok
        ),
        "nonlinear_force_uses_current_minor_loop_branch_interpolation_state_coenergy_and_mesh": (
            nonlinear_minor_loop_force_identity_ok
        ),
        "harmonic_eddy_loss_uses_current_phasor_conductivity_skin_depth_frequency_and_mesh": (
            harmonic_eddy_loss_identity_ok
        ),
        "axisymmetric_force_uses_current_aphi_weight_regions_energy_axis_mesh_and_solution": (
            axisymmetric_aphi_force_identity_ok
        ),
        "permanent_magnet_force_uses_current_recoil_temperature_operating_point_frame_demag_and_mesh": (
            permanent_magnet_operating_point_identity_ok
        ),
        "airgap_torque_uses_current_sector_sampling_alias_harmonics_geometry_mesh_and_result": (
            airgap_stress_harmonic_torque_identity_ok
        ),
        "laminated_core_loss_uses_current_frequency_flux_lamination_volume_components_and_result": (
            laminated_core_loss_identity_ok
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
