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


def _bem_demag_surface_material_frame_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("solve_generation", "")).strip()
    surface_ids = value.get("surface_ids")
    orientation = str(value.get("surface_orientation", "")).strip()
    frame = str(value.get("coordinate_frame", "")).strip()
    material_id = value.get("material_region_id")
    try:
        normals = [
            [float(component) for component in row]
            for row in value.get("outward_normals", [])
        ]
        result_normals = [
            [float(component) for component in row]
            for row in value.get("result_outward_normals", [])
        ]
        magnetization = [
            float(component)
            for component in value.get("magnetization_vector_a_per_m", [])
        ]
        result_magnetization = [
            float(component)
            for component in value.get("result_magnetization_vector_a_per_m", [])
        ]
        volume = float(value.get("body_volume_m3"))
        result_volume = float(value.get("result_body_volume_m3"))
    except (TypeError, ValueError):
        return False
    ids_ok = (
        isinstance(surface_ids, list)
        and bool(surface_ids)
        and all(
            isinstance(item, int) and not isinstance(item, bool) and item > 0
            for item in surface_ids
        )
        and len(set(surface_ids)) == len(surface_ids)
    )
    normals_ok = (
        ids_ok
        and len(normals) == len(surface_ids)
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
    surface_digest = str(value.get("surface_mesh_sha256", "")).lower()
    result_digest = str(value.get("demag_result_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "surface_mesh_solve_generation",
                "surface_orientation_solve_generation",
                "magnetization_solve_generation",
                "body_volume_solve_generation",
                "material_region_solve_generation",
                "coordinate_frame_solve_generation",
                "result_solve_generation",
            )
        )
        and ids_ok
        and value.get("result_surface_ids") == surface_ids
        and orientation == "outward_from_magnet"
        and value.get("result_surface_orientation") == orientation
        and normals_ok
        and result_normals == normals
        and len(magnetization) == 3
        and all(math.isfinite(component) for component in magnetization)
        and sum(component * component for component in magnetization) > 0.0
        and result_magnetization == magnetization
        and math.isfinite(volume)
        and volume > 0.0
        and math.isclose(result_volume, volume, rel_tol=0.0, abs_tol=0.0)
        and isinstance(material_id, int)
        and not isinstance(material_id, bool)
        and material_id > 0
        and value.get("result_material_region_id") == material_id
        and frame == "global_xyz"
        and value.get("result_coordinate_frame") == frame
        and _valid_sha256(surface_digest)
        and str(value.get("result_surface_mesh_sha256", "")).lower()
        == surface_digest
        and _valid_sha256(result_digest)
        and str(value.get("reported_demag_result_sha256", "")).lower()
        == result_digest
    )


def _linear_motor_thrust_ripple_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("sweep_generation", "")).strip()
    frame = str(value.get("force_coordinate_frame", "")).strip()
    order = value.get("sample_order")
    try:
        period = float(value.get("mechanical_period_m"))
        result_period = float(value.get("result_mechanical_period_m"))
        positions = [float(item) for item in value.get("mover_positions_m", [])]
        result_positions = [
            float(item) for item in value.get("result_mover_positions_m", [])
        ]
        phases = [float(item) for item in value.get("excitation_phase_deg", [])]
        result_phases = [
            float(item) for item in value.get("result_excitation_phase_deg", [])
        ]
        thrust = [float(item) for item in value.get("thrust_samples_n", [])]
        result_thrust = [
            float(item) for item in value.get("result_thrust_samples_n", [])
        ]
        ripple = float(value.get("thrust_ripple_peak_to_peak_n"))
        reported_ripple = float(
            value.get("reported_thrust_ripple_peak_to_peak_n")
        )
    except (TypeError, ValueError):
        return False
    digest = str(value.get("thrust_table_sha256", "")).lower()
    sample_count = len(positions)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "period_sweep_generation",
                "position_sweep_generation",
                "phase_sweep_generation",
                "force_frame_sweep_generation",
                "sample_order_sweep_generation",
                "result_sweep_generation",
            )
        )
        and math.isfinite(period)
        and period > 0.0
        and math.isclose(result_period, period, rel_tol=0.0, abs_tol=0.0)
        and sample_count >= 4
        and all(math.isfinite(item) for item in positions)
        and all(left < right for left, right in zip(positions, positions[1:]))
        and math.isclose(positions[-1] - positions[0], period, rel_tol=1.0e-12)
        and result_positions == positions
        and len(phases) == sample_count
        and all(math.isfinite(item) for item in phases)
        and all(left < right for left, right in zip(phases, phases[1:]))
        and math.isclose(phases[-1] - phases[0], 360.0, abs_tol=1.0e-12)
        and result_phases == phases
        and order == list(range(sample_count))
        and value.get("result_sample_order") == order
        and frame == "global_x"
        and value.get("result_force_coordinate_frame") == frame
        and len(thrust) == sample_count
        and all(math.isfinite(item) for item in thrust)
        and math.isclose(thrust[0], thrust[-1], rel_tol=1.0e-12)
        and result_thrust == thrust
        and math.isclose(ripple, max(thrust) - min(thrust), rel_tol=1.0e-12)
        and math.isclose(reported_ripple, ripple, rel_tol=1.0e-12)
        and _valid_sha256(digest)
        and str(value.get("result_thrust_table_sha256", "")).lower() == digest
    )


def _levitation_gradient_energy_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("levitation_generation", "")).strip()
    try:
        displacement = [float(item) for item in value.get("displacement_m", [])]
        result_displacement = [float(item) for item in value.get("result_displacement_m", [])]
        force = [float(item) for item in value.get("force_n", [])]
        result_force = [float(item) for item in value.get("result_force_n", [])]
        energy = [float(item) for item in value.get("magnetic_energy_j", [])]
        result_energy = [float(item) for item in value.get("result_magnetic_energy_j", [])]
        derivative_force = [float(item) for item in value.get("negative_energy_derivative_force_n", [])]
        result_derivative_force = [float(item) for item in value.get("result_negative_energy_derivative_force_n", [])]
        stiffness = float(value.get("restoring_stiffness_n_m"))
        result_stiffness = float(value.get("result_restoring_stiffness_n_m"))
    except (TypeError, ValueError):
        return False
    count = len(displacement)
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "displacement_levitation_generation", "force_levitation_generation",
            "energy_levitation_generation", "gradient_levitation_generation",
            "frame_levitation_generation", "result_levitation_generation"))
        and count >= 3 and len(force) == len(energy) == len(derivative_force) == count
        and all(math.isfinite(item) for item in displacement + force + energy + derivative_force)
        and all(left < right for left, right in zip(displacement, displacement[1:]))
        and result_displacement == displacement and result_force == force
        and result_energy == energy and result_derivative_force == derivative_force
        and math.isfinite(stiffness) and stiffness > 0.0 and result_stiffness == stiffness
        and all(math.isclose(item, -stiffness * x, rel_tol=1.0e-12, abs_tol=1.0e-12)
                for x, item in zip(displacement, force))
        and all(math.isclose(item, 0.5 * stiffness * x * x, rel_tol=1.0e-12, abs_tol=1.0e-15)
                for x, item in zip(displacement, energy))
        and derivative_force == force
        and value.get("coordinate_frame") == "global_z_up"
        and value.get("result_coordinate_frame") == "global_z_up"
        and value.get("force_sign_convention") == "restoring_negative_gradient"
        and value.get("result_force_sign_convention") == "restoring_negative_gradient"
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _cogging_periodic_interpolation_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("cogging_generation", "")).strip()
    try:
        positions = [float(item) for item in value.get("mechanical_positions_deg", [])]
        result_positions = [float(item) for item in value.get("result_mechanical_positions_deg", [])]
        torque = [float(item) for item in value.get("cogging_torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_cogging_torque_nm", [])]
        periodicity = int(value.get("periodicity"))
        result_periodicity = int(value.get("result_periodicity"))
        period = float(value.get("mechanical_period_deg"))
        result_period = float(value.get("result_mechanical_period_deg"))
        reference = float(value.get("reference_angle_deg"))
        result_reference = float(value.get("result_reference_angle_deg"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "position_cogging_generation", "periodicity_cogging_generation",
            "mesh_cogging_generation", "interpolation_cogging_generation",
            "reference_cogging_generation", "result_cogging_generation"))
        and len(positions) >= 4 and len(torque) == len(positions)
        and all(math.isfinite(item) for item in positions + torque)
        and all(left < right for left, right in zip(positions, positions[1:]))
        and result_positions == positions and result_torque == torque
        and periodicity > 0 and result_periodicity == periodicity
        and math.isclose(period, 360.0 / periodicity, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_period == period and result_reference == reference
        and math.isclose(positions[0], reference, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(positions[-1] - positions[0], period, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(torque[0], torque[-1], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and value.get("interpolation_method") == "periodic_cubic"
        and value.get("result_interpolation_method") == "periodic_cubic"
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _bem_panel_self_term_energy_force_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("panel_generation", "")).strip()
    try:
        panel_ids = [int(item) for item in value.get("panel_ids", [])]
        result_panel_ids = [int(item) for item in value.get("result_panel_ids", [])]
        normals = [
            [float(component) for component in row]
            for row in value.get("outward_unit_normals", [])
        ]
        result_normals = [
            [float(component) for component in row]
            for row in value.get("result_outward_unit_normals", [])
        ]
        areas = [float(item) for item in value.get("panel_area_m2", [])]
        result_areas = [float(item) for item in value.get("result_panel_area_m2", [])]
        magnetization = [
            [float(component) for component in row]
            for row in value.get("magnetization_a_m", [])
        ]
        result_magnetization = [
            [float(component) for component in row]
            for row in value.get("result_magnetization_a_m", [])
        ]
        displacement = [float(item) for item in value.get("displacement_m", [])]
        result_displacement = [
            float(item) for item in value.get("result_displacement_m", [])
        ]
        energy = [float(item) for item in value.get("magnetic_energy_j", [])]
        result_energy = [
            float(item) for item in value.get("result_magnetic_energy_j", [])
        ]
        force = [
            float(item)
            for item in value.get("negative_energy_derivative_force_n", [])
        ]
        result_force = [float(item) for item in value.get("result_force_n", [])]
    except (TypeError, ValueError):
        return False
    count = len(panel_ids)
    closure = [
        sum(area * normal[axis] for area, normal in zip(areas, normals))
        for axis in range(3)
    ] if count and len(areas) == len(normals) == count else [math.inf] * 3
    symmetric_displacement = (
        len(displacement) == 3
        and displacement[0] < 0.0
        and displacement[1] == 0.0
        and math.isclose(displacement[2], -displacement[0], rel_tol=0.0, abs_tol=1.0e-15)
    )
    stiffness = (
        2.0 * energy[2] / (displacement[2] * displacement[2])
        if symmetric_displacement and len(energy) == 3 and displacement[2] != 0.0
        else math.nan
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "orientation_panel_generation",
                "magnetization_panel_generation",
                "self_term_panel_generation",
                "energy_panel_generation",
                "force_panel_generation",
                "mesh_panel_generation",
                "result_panel_generation",
            )
        )
        and count >= 4
        and all(item > 0 for item in panel_ids)
        and len(set(panel_ids)) == count
        and result_panel_ids == panel_ids
        and len(normals) == len(areas) == len(magnetization) == count
        and all(
            len(row) == 3
            and all(math.isfinite(item) for item in row)
            and math.isclose(sum(item * item for item in row), 1.0, rel_tol=1.0e-12)
            for row in normals
        )
        and result_normals == normals
        and all(math.isfinite(item) and item > 0.0 for item in areas)
        and result_areas == areas
        and all(abs(item) <= 1.0e-12 * sum(areas) for item in closure)
        and all(
            len(row) == 3 and all(math.isfinite(item) for item in row)
            for row in magnetization
        )
        and result_magnetization == magnetization
        and value.get("magnetization_frame") == "global-cartesian"
        and value.get("result_magnetization_frame") == "global-cartesian"
        and value.get("singular_self_term") == "analytic-solid-angle"
        and value.get("result_singular_self_term") == "analytic-solid-angle"
        and symmetric_displacement
        and len(energy) == len(force) == 3
        and all(math.isfinite(item) for item in displacement + energy + force)
        and result_displacement == displacement
        and result_energy == energy
        and math.isfinite(stiffness)
        and stiffness > 0.0
        and all(
            math.isclose(item, 0.5 * stiffness * x * x, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for x, item in zip(displacement, energy)
        )
        and all(
            math.isclose(item, -stiffness * x, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for x, item in zip(displacement, force)
        )
        and result_force == force
        and _valid_sha256(value.get("panel_mesh_sha256"))
        and value.get("result_panel_mesh_sha256") == value.get("panel_mesh_sha256")
        and _valid_sha256(value.get("force_result_sha256"))
        and value.get("accepted_force_result_sha256") == value.get("force_result_sha256")
    )


def _motor_reduced_basis_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("reduced_generation", "")).strip()
    snapshot_ids = value.get("snapshot_ids")
    try:
        basis_dimension = int(value.get("basis_dimension"))
        result_basis_dimension = int(value.get("result_basis_dimension"))
        points = [
            [float(component) for component in row]
            for row in value.get("snapshot_operating_points", [])
        ]
        result_points = [
            [float(component) for component in row]
            for row in value.get("result_snapshot_operating_points", [])
        ]
        query = [float(item) for item in value.get("query_operating_point", [])]
        result_query = [
            float(item) for item in value.get("result_query_operating_point", [])
        ]
        weights = [float(item) for item in value.get("interpolation_weights", [])]
        result_weights = [
            float(item) for item in value.get("result_interpolation_weights", [])
        ]
        torque = [float(item) for item in value.get("snapshot_torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_snapshot_torque_nm", [])]
        reduced_torque = float(value.get("reduced_torque_nm"))
        result_reduced_torque = float(value.get("result_reduced_torque_nm"))
        residual = float(value.get("relative_residual"))
        result_residual = float(value.get("result_relative_residual"))
        residual_limit = float(value.get("accepted_relative_residual"))
    except (TypeError, ValueError):
        return False
    count = len(snapshot_ids) if isinstance(snapshot_ids, list) else 0
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "basis_reduced_generation",
                "snapshot_reduced_generation",
                "operating_point_reduced_generation",
                "weight_reduced_generation",
                "torque_reduced_generation",
                "residual_reduced_generation",
                "result_reduced_generation",
            )
        )
        and count >= 2
        and all(isinstance(item, str) and item.strip() for item in snapshot_ids)
        and len(set(snapshot_ids)) == count
        and value.get("result_snapshot_ids") == snapshot_ids
        and 1 <= basis_dimension <= count
        and result_basis_dimension == basis_dimension
        and len(points) == len(weights) == len(torque) == count
        and all(len(row) == 3 and all(math.isfinite(item) for item in row) for row in points)
        and result_points == points
        and len(query) == 3
        and all(math.isfinite(item) for item in query)
        and result_query == query
        and all(math.isfinite(item) and item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_weights == weights
        and all(math.isfinite(item) for item in torque)
        and result_torque == torque
        and math.isclose(
            reduced_torque,
            sum(weight * item for weight, item in zip(weights, torque)),
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and result_reduced_torque == reduced_torque
        and math.isfinite(residual)
        and 0.0 <= residual <= residual_limit
        and result_residual == residual
        and _valid_sha256(value.get("basis_sha256"))
        and value.get("loaded_basis_sha256") == value.get("basis_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _maglev_force_stiffness_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("stiffness_generation", "")).strip()
    try:
        displacement = [float(item) for item in value.get("displacement_m", [])]
        result_displacement = [
            float(item) for item in value.get("result_displacement_m", [])
        ]
        step = float(value.get("displacement_step_m"))
        result_step = float(value.get("result_displacement_step_m"))
        force = [float(item) for item in value.get("force_n", [])]
        result_force = [float(item) for item in value.get("result_force_n", [])]
        stiffness = float(value.get("stiffness_n_m"))
        result_stiffness = float(value.get("result_stiffness_n_m"))
    except (TypeError, ValueError):
        return False
    slopes = [
        -(right_force - left_force) / (right_x - left_x)
        for left_x, right_x, left_force, right_force in zip(
            displacement, displacement[1:], force, force[1:]
        )
    ] if len(displacement) == len(force) and len(displacement) >= 3 else []
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "displacement_stiffness_generation",
                "coordinate_stiffness_generation",
                "geometry_stiffness_generation",
                "mesh_stiffness_generation",
                "force_stiffness_generation",
                "derivative_stiffness_generation",
                "solution_stiffness_generation",
                "result_stiffness_generation",
            )
        )
        and len(displacement) >= 3
        and len(force) == len(displacement)
        and all(math.isfinite(item) for item in displacement + force)
        and all(right > left for left, right in zip(displacement, displacement[1:]))
        and math.isfinite(step)
        and step > 0.0
        and all(
            math.isclose(right - left, step, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for left, right in zip(displacement, displacement[1:])
        )
        and result_displacement == displacement
        and math.isclose(result_step, step, rel_tol=0.0, abs_tol=1.0e-15)
        and value.get("coordinate_direction") == "global-z-positive"
        and value.get("result_coordinate_direction") == "global-z-positive"
        and result_force == force
        and value.get("derivative_convention")
        == "stiffness-equals-negative-force-derivative"
        and value.get("result_derivative_convention")
        == "stiffness-equals-negative-force-derivative"
        and slopes
        and all(math.isfinite(item) and item > 0.0 for item in slopes)
        and all(
            math.isclose(item, stiffness, rel_tol=1.0e-12, abs_tol=1.0e-9)
            for item in slopes
        )
        and math.isclose(result_stiffness, stiffness, rel_tol=0.0, abs_tol=1.0e-9)
        and _valid_sha256(value.get("geometry_sha256"))
        and value.get("result_geometry_sha256") == value.get("geometry_sha256")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("solution_sha256"))
        and value.get("accepted_solution_sha256") == value.get("solution_sha256")
    )


def _motor_coenergy_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("torque_generation", "")).strip()
    phase_order = value.get("phase_order")
    try:
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        currents = [
            [float(component) for component in row]
            for row in value.get("phase_current_harmonic_a", [])
        ]
        result_currents = [
            [float(component) for component in row]
            for row in value.get("result_phase_current_harmonic_a", [])
        ]
        phases = [float(item) for item in value.get("current_phase_deg", [])]
        result_phases = [
            float(item) for item in value.get("result_current_phase_deg", [])
        ]
        angles = [
            float(item) for item in value.get("rotor_mechanical_angle_deg", [])
        ]
        result_angles = [
            float(item)
            for item in value.get("result_rotor_mechanical_angle_deg", [])
        ]
        coenergy = [float(item) for item in value.get("coenergy_j", [])]
        result_coenergy = [float(item) for item in value.get("result_coenergy_j", [])]
        torque = [float(item) for item in value.get("torque_nm", [])]
        result_torque = [float(item) for item in value.get("result_torque_nm", [])]
    except (TypeError, ValueError):
        return False
    derived_torque = [
        (right_energy - left_energy) / math.radians(right_angle - left_angle)
        for left_angle, right_angle, left_energy, right_energy in zip(
            angles, angles[1:], coenergy, coenergy[1:]
        )
    ] if len(angles) == len(coenergy) and len(angles) >= 3 else []
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "winding_torque_generation",
                "harmonic_torque_generation",
                "current_torque_generation",
                "phase_torque_generation",
                "angle_torque_generation",
                "coenergy_torque_generation",
                "mesh_torque_generation",
                "result_torque_generation",
            )
        )
        and phase_order == ["U", "V", "W"]
        and value.get("result_phase_order") == phase_order
        and bool(orders)
        and all(item > 0 for item in orders)
        and len(set(orders)) == len(orders)
        and result_orders == orders
        and len(currents) == len(orders)
        and all(
            len(row) == 3
            and all(math.isfinite(item) for item in row)
            and math.isclose(sum(row), 0.0, rel_tol=0.0, abs_tol=1.0e-12)
            for row in currents
        )
        and result_currents == currents
        and len(phases) == 3
        and all(math.isfinite(item) for item in phases)
        and result_phases == phases
        and len(angles) >= 3
        and len(coenergy) == len(angles)
        and all(math.isfinite(item) for item in angles + coenergy)
        and all(right > left for left, right in zip(angles, angles[1:]))
        and result_angles == angles
        and result_coenergy == coenergy
        and value.get("torque_convention") == "positive-coenergy-angle-derivative"
        and value.get("result_torque_convention")
        == "positive-coenergy-angle-derivative"
        and len(torque) == len(angles) - 1
        and all(math.isfinite(item) for item in torque)
        and all(
            math.isclose(item, derived, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for item, derived in zip(torque, derived_torque)
        )
        and result_torque == torque
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _airgap_stress_harmonic_torque_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("stress_generation", "")).strip()
    try:
        pitch = float(value.get("sector_pitch_deg"))
        result_pitch = float(value.get("result_sector_pitch_deg"))
        sectors = int(value.get("sector_count"))
        result_sectors = int(value.get("result_sector_count"))
        origin = float(value.get("angular_origin_deg"))
        result_origin = float(value.get("result_angular_origin_deg"))
        samples = int(value.get("angular_sample_count"))
        result_samples = int(value.get("result_angular_sample_count"))
        sector_samples = int(value.get("sector_sample_count"))
        result_sector_samples = int(value.get("result_sector_sample_count"))
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        harmonics = [float(item) for item in value.get("torque_harmonics_nm", [])]
        result_harmonics = [float(item) for item in value.get("result_torque_harmonics_nm", [])]
        cutoff = int(value.get("alias_cutoff_order"))
        result_cutoff = int(value.get("result_alias_cutoff_order"))
        radius = float(value.get("airgap_radius_m"))
        result_radius = float(value.get("result_airgap_radius_m"))
        length = float(value.get("axial_length_m"))
        result_length = float(value.get("result_axial_length_m"))
        torque = float(value.get("torque_nm"))
        result_torque = float(value.get("result_torque_nm"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "harmonic_stress_generation", "sector_stress_generation",
            "sampling_stress_generation", "alias_stress_generation",
            "geometry_stress_generation", "mesh_stress_generation",
            "result_stress_generation"))
        and all(math.isfinite(item) for item in (pitch, origin, radius, length, torque))
        and pitch > 0.0 and sectors > 1
        and math.isclose(pitch * sectors, 360.0, rel_tol=0.0, abs_tol=1.0e-12)
        and result_pitch == pitch and result_sectors == sectors and result_origin == origin
        and samples > 0 and samples % sectors == 0
        and sector_samples == samples // sectors
        and result_samples == samples and result_sector_samples == sector_samples
        and bool(orders) and orders == sorted(set(orders)) and orders[0] == 0
        and all(order >= 0 and order % sectors in {0, sectors // 2} for order in orders)
        and result_orders == orders
        and len(harmonics) == len(orders) and all(math.isfinite(item) for item in harmonics)
        and result_harmonics == harmonics
        and value.get("alias_filter") == "truncate_below_nyquist"
        and value.get("result_alias_filter") == value.get("alias_filter")
        and 0 <= cutoff < samples // 2 and cutoff >= max(orders)
        and result_cutoff == cutoff
        and radius > 0.0 and length > 0.0
        and result_radius == radius and result_length == length
        and math.isclose(torque, sum(harmonics), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_torque == torque
        and _valid_sha256(value.get("airgap_mesh_sha256"))
        and value.get("result_airgap_mesh_sha256") == value.get("airgap_mesh_sha256")
        and _valid_sha256(value.get("torque_result_sha256"))
        and value.get("accepted_torque_result_sha256") == value.get("torque_result_sha256")
    )


def _laminated_core_loss_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("loss_generation", "")).strip()
    names = (
        "frequency_hz", "peak_flux_density_t", "lamination_thickness_m",
        "magnetic_volume_m3", "hysteresis_coefficient", "hysteresis_exponent",
        "eddy_coefficient", "excess_coefficient", "hysteresis_loss_w",
        "eddy_loss_w", "excess_loss_w", "total_core_loss_w",
    )
    try:
        data = {name: float(value.get(name)) for name in names}
        result = {name: float(value.get(f"result_{name}")) for name in names}
    except (TypeError, ValueError):
        return False
    frequency = data["frequency_hz"]
    flux = data["peak_flux_density_t"]
    thickness = data["lamination_thickness_m"]
    volume = data["magnetic_volume_m3"]
    expected_hysteresis = data["hysteresis_coefficient"] * frequency * flux ** data["hysteresis_exponent"] * volume
    expected_eddy = data["eddy_coefficient"] * frequency**2 * flux**2 * thickness**2 * volume
    expected_excess = data["excess_coefficient"] * frequency**1.5 * flux**1.5 * volume
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "hysteresis_loss_generation", "eddy_loss_generation", "excess_loss_generation",
            "frequency_loss_generation", "flux_loss_generation", "lamination_loss_generation",
            "volume_loss_generation", "result_loss_generation"))
        and all(math.isfinite(item) and item > 0.0 for item in data.values())
        and result == data
        and math.isclose(data["hysteresis_loss_w"], expected_hysteresis, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(data["eddy_loss_w"], expected_eddy, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(data["excess_loss_w"], expected_excess, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(data["total_core_loss_w"], expected_hysteresis + expected_eddy + expected_excess, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and _valid_sha256(value.get("material_sha256"))
        and value.get("result_material_sha256") == value.get("material_sha256")
        and _valid_sha256(value.get("loss_result_sha256"))
        and value.get("accepted_loss_result_sha256") == value.get("loss_result_sha256")
    )


def _magnet_demag_volume_fraction_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    try:
        recoil = float(value.get("recoil_relative_permeability"))
        result_recoil = float(value.get("result_recoil_relative_permeability"))
        knee = float(value.get("knee_field_a_m"))
        result_knee = float(value.get("result_knee_field_a_m"))
        temperature = float(value.get("temperature_c"))
        result_temperature = float(value.get("result_temperature_c"))
        element_ids = [int(item) for item in value.get("element_ids", [])]
        result_element_ids = [int(item) for item in value.get("result_element_ids", [])]
        local_field = [float(item) for item in value.get("local_recoil_axis_field_a_m", [])]
        result_local_field = [float(item) for item in value.get("result_local_recoil_axis_field_a_m", [])]
        volumes = [float(item) for item in value.get("element_volumes_m3", [])]
        result_volumes = [float(item) for item in value.get("result_element_volumes_m3", [])]
        magnet_volume = float(value.get("magnet_volume_m3"))
        result_magnet_volume = float(value.get("result_magnet_volume_m3"))
        fraction = float(value.get("irreversible_volume_fraction"))
        result_fraction = float(value.get("result_irreversible_volume_fraction"))
    except (TypeError, ValueError):
        return False
    mask = value.get("irreversible_mask")
    result_mask = value.get("result_irreversible_mask")
    derived_mask = [item <= knee for item in local_field]
    count = len(element_ids)
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "recoil_demag_generation", "knee_demag_generation",
            "temperature_demag_generation", "field_demag_generation",
            "mask_demag_generation", "volume_demag_generation",
            "mesh_demag_generation", "result_demag_generation"))
        and math.isfinite(recoil) and recoil > 0.0
        and math.isclose(result_recoil, recoil, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isfinite(knee) and knee < 0.0
        and math.isclose(result_knee, knee, rel_tol=0.0, abs_tol=1.0e-9)
        and math.isfinite(temperature)
        and math.isclose(result_temperature, temperature, rel_tol=0.0, abs_tol=1.0e-12)
        and count > 0 and len(set(element_ids)) == count
        and all(item > 0 for item in element_ids)
        and len(local_field) == len(volumes) == count
        and all(math.isfinite(item) for item in local_field)
        and all(math.isfinite(item) and item > 0.0 for item in volumes)
        and result_element_ids == element_ids
        and result_local_field == local_field
        and result_volumes == volumes
        and isinstance(mask, list) and len(mask) == count
        and all(isinstance(item, bool) for item in mask)
        and mask == derived_mask and result_mask == mask
        and math.isclose(magnet_volume, sum(volumes), rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(result_magnet_volume, magnet_volume, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(
            fraction,
            sum(volume for volume, failed in zip(volumes, mask) if failed) / magnet_volume,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and math.isclose(result_fraction, fraction, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("material_state_sha256"))
        and value.get("result_material_state_sha256") == value.get("material_state_sha256")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _linear_motor_wave_end_effect_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("linear_motor_generation", "")).strip()
    phase_sequence = value.get("phase_sequence")
    try:
        pitch = float(value.get("pole_pitch_m"))
        result_pitch = float(value.get("result_pole_pitch_m"))
        positions = [float(item) for item in value.get("position_m", [])]
        result_positions = [float(item) for item in value.get("result_position_m", [])]
        end_effect = [float(item) for item in value.get("end_effect_factor", [])]
        result_end_effect = [float(item) for item in value.get("result_end_effect_factor", [])]
        force = [float(item) for item in value.get("force_n", [])]
        result_force = [float(item) for item in value.get("result_force_n", [])]
        mean_force = float(value.get("mean_force_n"))
        result_mean_force = float(value.get("result_mean_force_n"))
        ripple = float(value.get("force_ripple_peak_to_peak_n"))
        result_ripple = float(value.get("result_force_ripple_peak_to_peak_n"))
    except (TypeError, ValueError):
        return False
    count = len(positions)
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "end_effect_linear_motor_generation", "phase_linear_motor_generation",
            "wave_linear_motor_generation", "pitch_linear_motor_generation",
            "position_linear_motor_generation", "force_linear_motor_generation",
            "ripple_linear_motor_generation", "result_linear_motor_generation"))
        and phase_sequence == ["U", "V", "W"]
        and value.get("result_phase_sequence") == phase_sequence
        and value.get("traveling_wave_direction") == "global-x-positive"
        and value.get("result_traveling_wave_direction") == value.get("traveling_wave_direction")
        and math.isfinite(pitch) and pitch > 0.0
        and math.isclose(result_pitch, pitch, rel_tol=0.0, abs_tol=1.0e-15)
        and count >= 5 and len(end_effect) == len(force) == count
        and all(math.isfinite(item) for item in positions + end_effect + force)
        and all(right > left for left, right in zip(positions, positions[1:]))
        and math.isclose(positions[-1] - positions[0], pitch, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(0.0 < item <= 1.0 for item in end_effect)
        and result_positions == positions and result_end_effect == end_effect
        and result_force == force
        and math.isclose(force[0], force[-1], rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(mean_force, sum(force) / count, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(result_mean_force, mean_force, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(ripple, max(force) - min(force), rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(result_ripple, ripple, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _maglev_stiffness_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("maglev_generation", "")).strip()
    try:
        position = [float(item) for item in value.get("position_m", [])]
        currents = [float(item) for item in value.get("current_a", [])]
        force = [float(item) for item in value.get("force_z_n", [])]
        stiffness = float(value.get("stiffness_n_per_m"))
    except (TypeError, ValueError):
        return False
    derivative = (
        (force[2] - force[0]) / (position[2] - position[0])
        if len(force) == 3 and len(position) == 3 and position[2] != position[0]
        else math.nan
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "position_maglev_generation", "current_maglev_generation", "force_maglev_generation",
            "derivative_maglev_generation", "frame_maglev_generation", "mesh_maglev_generation",
            "result_maglev_generation"))
        and len(position) == 3 and position[0] < position[1] < position[2]
        and math.isclose(position[1], 0.0, abs_tol=1.0e-15)
        and math.isclose(position[2], -position[0], rel_tol=1.0e-12)
        and value.get("result_position_m") == position
        and len(currents) == 3 and all(math.isfinite(item) for item in currents)
        and len(set(currents)) == 1 and value.get("result_current_a") == currents
        and len(force) == 3 and all(math.isfinite(item) for item in force)
        and value.get("result_force_z_n") == force
        and math.isfinite(stiffness) and stiffness < 0.0
        and math.isclose(stiffness, derivative, rel_tol=1.0e-10)
        and value.get("result_stiffness_n_per_m") == stiffness
        and value.get("coordinate_frame") == "global_z_positive_up"
        and value.get("result_coordinate_frame") == value.get("coordinate_frame")
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _cogging_torque_sampling_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("cogging_generation", "")).strip()
    try:
        slots = int(value.get("slot_count")); poles = int(value.get("pole_count"))
        period = float(value.get("cogging_period_mechanical_deg"))
        origin = float(value.get("angular_origin_deg"))
        angles = [float(item) for item in value.get("sample_angles_deg", [])]
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        phases = [float(item) for item in value.get("harmonic_phase_deg", [])]
        torque = [float(item) for item in value.get("torque_nm", [])]
    except (TypeError, ValueError):
        return False
    expected_period = 360.0 / math.lcm(slots, poles) if slots > 0 and poles > 0 else math.nan
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "slot_cogging_generation", "pole_cogging_generation", "period_cogging_generation",
            "sampling_cogging_generation", "harmonic_cogging_generation", "mesh_cogging_generation",
            "result_cogging_generation"))
        and slots > 0 and value.get("result_slot_count") == slots
        and poles > 0 and poles % 2 == 0 and value.get("result_pole_count") == poles
        and math.isclose(period, expected_period, rel_tol=1.0e-12)
        and value.get("result_cogging_period_mechanical_deg") == period
        and math.isfinite(origin) and value.get("result_angular_origin_deg") == origin
        and len(angles) >= 5 and angles[0] == origin
        and all(left < right for left, right in zip(angles, angles[1:]))
        and math.isclose(angles[-1] - angles[0], period, rel_tol=1.0e-12)
        and value.get("result_sample_angles_deg") == angles
        and bool(orders) and all(item > 0 for item in orders) and len(set(orders)) == len(orders)
        and value.get("result_harmonic_orders") == orders
        and len(phases) == len(orders) and all(math.isfinite(item) for item in phases)
        and value.get("result_harmonic_phase_deg") == phases
        and len(torque) == len(angles) and all(math.isfinite(item) for item in torque)
        and value.get("result_torque_nm") == torque
        and math.isclose(torque[0], torque[-1], rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _bem_near_singular_force_identity_ok(value: object) -> bool:
    if value is None: return True
    if not isinstance(value, dict): return False
    generation = str(value.get("bem_generation", "")).strip()
    try:
        gap = float(value.get("gap_m")); panel = float(value.get("panel_size_m")); order = int(value.get("quadrature_order"))
        source_normal = [float(item) for item in value.get("source_normal", [])]; target_normal = [float(item) for item in value.get("target_normal", [])]
        force_source = [float(item) for item in value.get("force_on_source_n", [])]; force_target = [float(item) for item in value.get("force_on_target_n", [])]
        residual = float(value.get("action_reaction_residual_n"))
    except (TypeError, ValueError): return False
    return (
        bool(generation) and all(value.get(key) == generation for key in ("gap_bem_generation", "quadrature_bem_generation", "normal_bem_generation", "order_bem_generation", "force_bem_generation", "reciprocity_bem_generation", "geometry_bem_generation", "result_bem_generation"))
        and 0.0 < gap < panel and value.get("result_gap_m") == gap and value.get("result_panel_size_m") == panel
        and value.get("quadrature_policy") == "gap_adaptive_duffy" and value.get("result_quadrature_policy") == value.get("quadrature_policy")
        and order >= 8 and value.get("result_quadrature_order") == order
        and len(source_normal) == len(target_normal) == 3 and value.get("result_source_normal") == source_normal and value.get("result_target_normal") == target_normal
        and all(math.isclose(a, -b, abs_tol=1.0e-12) for a, b in zip(source_normal, target_normal))
        and value.get("source_target_order") == ["body_a", "body_b"] and value.get("result_source_target_order") == value.get("source_target_order")
        and len(force_source) == len(force_target) == 3 and value.get("result_force_on_source_n") == force_source and value.get("result_force_on_target_n") == force_target
        and all(math.isclose(a, -b, abs_tol=1.0e-12) for a, b in zip(force_source, force_target))
        and math.isclose(residual, 0.0, abs_tol=1.0e-12) and value.get("result_action_reaction_residual_n") == residual
        and _valid_sha256(value.get("geometry_sha256")) and value.get("result_geometry_sha256") == value.get("geometry_sha256")
        and _valid_sha256(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _hysteresis_minor_loop_identity_ok(value: object) -> bool:
    if value is None: return True
    if not isinstance(value, dict): return False
    generation = str(value.get("loop_generation", "")).strip()
    try:
        time = [float(item) for item in value.get("time_s", [])]; drive = [float(item) for item in value.get("drive_h_a_per_m", [])]
        reversals = [int(item) for item in value.get("reversal_indices", [])]; remanence = float(value.get("remanence_t")); energy = float(value.get("loop_energy_j_per_m3"))
    except (TypeError, ValueError): return False
    return (
        bool(generation) and all(value.get(key) == generation for key in ("state_loop_generation", "reversal_loop_generation", "memory_loop_generation", "remanence_loop_generation", "energy_loop_generation", "time_loop_generation", "material_loop_generation", "result_loop_generation"))
        and _valid_sha256(value.get("initial_state_sha256")) and value.get("result_initial_state_sha256") == value.get("initial_state_sha256")
        and len(time) == len(drive) >= 5 and all(right > left for left, right in zip(time, time[1:])) and value.get("result_time_s") == time and value.get("result_drive_h_a_per_m") == drive
        and reversals == [1, 2, 3] and value.get("result_reversal_indices") == reversals
        and value.get("return_point_memory_closed") is True and value.get("result_return_point_memory_closed") is True
        and math.isfinite(remanence) and remanence >= 0.0 and value.get("result_remanence_t") == remanence
        and math.isfinite(energy) and energy > 0.0 and value.get("result_loop_energy_j_per_m3") == energy
        and bool(str(value.get("material_owner", "")).strip()) and value.get("result_material_owner") == value.get("material_owner")
        and _valid_sha256(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _maglev_equilibrium_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("equilibrium_generation", "")).strip()
    try:
        displacements = [float(item) for item in value.get("displacement_samples_m", [])]
        forces = [float(item) for item in value.get("magnetic_force_samples_n", [])]
        derivative = float(value.get("force_derivative_n_per_m"))
        stiffness = float(value.get("vertical_stiffness_n_per_m"))
        gravity = float(value.get("gravity_force_n"))
    except (TypeError, ValueError):
        return False
    if len(displacements) != 3 or len(forces) != 3:
        return False
    spacing_left = displacements[1] - displacements[0]
    spacing_right = displacements[2] - displacements[1]
    if not (
        all(math.isfinite(item) for item in (*displacements, *forces, derivative, stiffness, gravity))
        and spacing_left > 0.0
        and math.isclose(spacing_left, spacing_right, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(displacements[1], 0.0, rel_tol=0.0, abs_tol=1.0e-15)
    ):
        return False
    central_derivative = (forces[2] - forces[0]) / (
        displacements[2] - displacements[0]
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "force_equilibrium_generation",
                "frame_equilibrium_generation",
                "displacement_equilibrium_generation",
                "derivative_equilibrium_generation",
                "stiffness_equilibrium_generation",
                "gravity_equilibrium_generation",
                "mesh_equilibrium_generation",
                "result_equilibrium_generation",
            )
        )
        and value.get("force_sign_convention") == "positive_up"
        and value.get("result_force_sign_convention") == "positive_up"
        and value.get("displacement_frame") == "global_z_up"
        and value.get("result_displacement_frame") == "global_z_up"
        and value.get("result_displacement_samples_m") == displacements
        and value.get("result_magnetic_force_samples_n") == forces
        and value.get("derivative_stencil") == "symmetric_central_difference"
        and value.get("result_derivative_stencil") == "symmetric_central_difference"
        and math.isclose(derivative, central_derivative, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and value.get("result_force_derivative_n_per_m") == derivative
        and stiffness > 0.0
        and math.isclose(stiffness, -derivative, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and value.get("result_vertical_stiffness_n_per_m") == stiffness
        and math.isclose(forces[1] + gravity, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and value.get("result_gravity_force_n") == gravity
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _bem_surface_charge_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("bem_generation", "")).strip()
    try:
        net_charge = float(value.get("net_surface_charge"))
        charge_tolerance = float(value.get("charge_balance_tolerance"))
        source_normal = [float(item) for item in value.get("source_normal", [])]
        target_normal = [float(item) for item in value.get("target_normal", [])]
        energy = float(value.get("field_energy_j"))
        reciprocity_residual = float(value.get("reciprocity_residual"))
        reciprocity_tolerance = float(value.get("reciprocity_tolerance"))
    except (TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "charge_bem_generation",
                "gauge_bem_generation",
                "normal_bem_generation",
                "energy_bem_generation",
                "reciprocity_bem_generation",
                "geometry_bem_generation",
                "owner_bem_generation",
                "result_bem_generation",
            )
        )
        and math.isfinite(net_charge)
        and math.isfinite(charge_tolerance)
        and charge_tolerance >= 0.0
        and abs(net_charge) <= charge_tolerance
        and value.get("result_net_surface_charge") == net_charge
        and value.get("gauge_reference") == "mean_zero_scalar_potential"
        and value.get("result_gauge_reference") == value.get("gauge_reference")
        and len(source_normal) == len(target_normal) == 3
        and math.isclose(sum(item * item for item in source_normal), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(sum(item * item for item in target_normal), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and all(math.isclose(a, -b, rel_tol=0.0, abs_tol=1.0e-12) for a, b in zip(source_normal, target_normal))
        and value.get("result_source_normal") == source_normal
        and value.get("result_target_normal") == target_normal
        and math.isfinite(energy)
        and energy >= 0.0
        and value.get("result_field_energy_j") == energy
        and math.isfinite(reciprocity_residual)
        and math.isfinite(reciprocity_tolerance)
        and reciprocity_tolerance >= 0.0
        and abs(reciprocity_residual) <= reciprocity_tolerance
        and value.get("result_reciprocity_residual") == reciprocity_residual
        and _valid_sha256(value.get("geometry_sha256"))
        and value.get("result_geometry_sha256") == value.get("geometry_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _halbach_harmonic_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("halbach_generation", "")).strip()
    try:
        angles = [float(item) for item in value.get("magnetization_angles_deg", [])]
        result_angles = [
            float(item) for item in value.get("result_magnetization_angles_deg", [])
        ]
        pitch = float(value.get("pole_pitch_m"))
        result_pitch = float(value.get("result_pole_pitch_m"))
        orders = [int(item) for item in value.get("harmonic_orders", [])]
        result_orders = [int(item) for item in value.get("result_harmonic_orders", [])]
        phases = [float(item) for item in value.get("harmonic_phase_deg", [])]
        result_phases = [
            float(item) for item in value.get("result_harmonic_phase_deg", [])
        ]
        grid = [float(item) for item in value.get("sampling_grid_m", [])]
        result_grid = [
            float(item) for item in value.get("result_sampling_grid_m", [])
        ]
        amplitudes = [
            float(item) for item in value.get("field_harmonic_amplitude_t", [])
        ]
        result_amplitudes = [
            float(item)
            for item in value.get("result_field_harmonic_amplitude_t", [])
        ]
        energy = float(value.get("magnetic_energy_j"))
        result_energy = float(value.get("result_magnetic_energy_j"))
        force = float(value.get("force_n"))
        result_force = float(value.get("result_force_n"))
    except (TypeError, ValueError):
        return False
    wrapped_steps = [
        (angles[(index + 1) % len(angles)] - angles[index]) % 360.0
        for index in range(len(angles))
    ] if angles else []
    expected_step = 360.0 / len(angles) if angles else math.nan
    direction = value.get("force_direction")
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "magnetization_generation",
                "pitch_generation",
                "phase_generation",
                "grid_generation",
                "field_generation",
                "energy_generation",
                "force_generation",
                "geometry_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and len(angles) >= 4
        and all(math.isfinite(item) for item in angles)
        and all(
            math.isclose(step, expected_step, rel_tol=0.0, abs_tol=1.0e-12)
            for step in wrapped_steps
        )
        and result_angles == angles
        and math.isfinite(pitch)
        and pitch > 0.0
        and math.isclose(result_pitch, pitch, rel_tol=0.0, abs_tol=1.0e-15)
        and bool(orders)
        and orders == sorted(set(orders))
        and all(item > 0 and item % 2 == 1 for item in orders)
        and result_orders == orders
        and len(phases) == len(orders)
        and all(math.isfinite(item) for item in phases)
        and result_phases == phases
        and len(grid) >= 5
        and all(math.isfinite(item) for item in grid)
        and all(left < right for left, right in zip(grid, grid[1:]))
        and math.isclose(grid[0], 0.0, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(grid[-1], pitch, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and result_grid == grid
        and len(amplitudes) == len(orders)
        and all(math.isfinite(item) and item >= 0.0 for item in amplitudes)
        and all(
            left >= right for left, right in zip(amplitudes, amplitudes[1:])
        )
        and result_amplitudes == amplitudes
        and math.isfinite(energy)
        and energy >= 0.0
        and math.isclose(result_energy, energy, rel_tol=0.0, abs_tol=1.0e-12)
        and direction in {"+x", "-x"}
        and value.get("result_force_direction") == direction
        and math.isfinite(force)
        and ((direction == "+x" and force >= 0.0) or (direction == "-x" and force <= 0.0))
        and math.isclose(result_force, force, rel_tol=0.0, abs_tol=1.0e-12)
        and _valid_sha256(value.get("geometry_sha256"))
        and value.get("result_geometry_sha256") == value.get("geometry_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _magnetic_bearing_linearization_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("bearing_generation", "")).strip()
    try:
        currents = [float(item) for item in value.get("bias_currents_a", [])]
        result_currents = [
            float(item) for item in value.get("result_bias_currents_a", [])
        ]
        displacement = [
            float(item) for item in value.get("bias_displacement_m", [])
        ]
        result_displacement = [
            float(item) for item in value.get("result_bias_displacement_m", [])
        ]
        current_jacobian = [
            [float(item) for item in row]
            for row in value.get("force_current_jacobian_n_per_a", [])
        ]
        result_current_jacobian = [
            [float(item) for item in row]
            for row in value.get("result_force_current_jacobian_n_per_a", [])
        ]
        displacement_jacobian = [
            [float(item) for item in row]
            for row in value.get("force_displacement_jacobian_n_per_m", [])
        ]
        result_displacement_jacobian = [
            [float(item) for item in row]
            for row in value.get("result_force_displacement_jacobian_n_per_m", [])
        ]
        stiffness = [
            [float(item) for item in row]
            for row in value.get("stiffness_matrix_n_per_m", [])
        ]
        result_stiffness = [
            [float(item) for item in row]
            for row in value.get("result_stiffness_matrix_n_per_m", [])
        ]
        eigenvalues = [
            float(item) for item in value.get("stiffness_eigenvalues_n_per_m", [])
        ]
        result_eigenvalues = [
            float(item)
            for item in value.get("result_stiffness_eigenvalues_n_per_m", [])
        ]
    except (TypeError, ValueError):
        return False
    matrices_are_2d = (
        len(current_jacobian) == 2
        and all(len(row) == len(currents) for row in current_jacobian)
        and len(displacement_jacobian) == 2
        and all(len(row) == 2 for row in displacement_jacobian)
        and len(stiffness) == 2
        and all(len(row) == 2 for row in stiffness)
    )
    if not matrices_are_2d:
        return False
    trace = stiffness[0][0] + stiffness[1][1]
    discriminant = (stiffness[0][0] - stiffness[1][1]) ** 2 + 4.0 * stiffness[0][1] * stiffness[1][0]
    derived_eigenvalues = sorted(
        ((trace - math.sqrt(discriminant)) / 2.0, (trace + math.sqrt(discriminant)) / 2.0)
    ) if discriminant >= 0.0 else []
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "current_generation",
                "displacement_generation",
                "jacobian_generation",
                "stiffness_generation",
                "reciprocity_generation",
                "bias_generation",
                "frame_generation",
                "mesh_generation",
                "result_generation",
            )
        )
        and value.get("coordinate_frame") == "global_xy_right_handed"
        and value.get("result_coordinate_frame") == value.get("coordinate_frame")
        and len(currents) >= 2
        and all(math.isfinite(item) for item in currents)
        and result_currents == currents
        and len(displacement) == 2
        and all(math.isfinite(item) for item in displacement)
        and result_displacement == displacement
        and all(math.isfinite(item) for row in current_jacobian for item in row)
        and result_current_jacobian == current_jacobian
        and all(
            math.isfinite(item) for row in displacement_jacobian for item in row
        )
        and result_displacement_jacobian == displacement_jacobian
        and all(math.isfinite(item) for row in stiffness for item in row)
        and all(
            math.isclose(stiffness[row][column], -displacement_jacobian[row][column], rel_tol=1.0e-12, abs_tol=1.0e-12)
            for row in range(2)
            for column in range(2)
        )
        and math.isclose(stiffness[0][1], stiffness[1][0], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_stiffness == stiffness
        and len(eigenvalues) == 2
        and all(math.isfinite(item) and item > 0.0 for item in eigenvalues)
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-9)
            for observed, expected in zip(sorted(eigenvalues), derived_eigenvalues)
        )
        and result_eigenvalues == eigenvalues
        and _valid_sha256(value.get("mesh_sha256"))
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _valid_sha256(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _magnetic_bearing_dynamic_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("bearing_dynamic_generation", "")).strip()
    try:
        displacement = [
            [float(item) for item in row]
            for row in value.get("displacement_perturbations_m", [])
        ]
        force = [
            [float(item) for item in row]
            for row in value.get("force_perturbations_n", [])
        ]
        stiffness = [
            [float(item) for item in row]
            for row in value.get("stiffness_matrix_n_per_m", [])
        ]
        damping = [
            [float(item) for item in row]
            for row in value.get("damping_matrix_n_s_per_m", [])
        ]
        eigenvalues = [
            [float(item) for item in row]
            for row in value.get("state_eigenvalues_per_s", [])
        ]
        operating_displacement = [
            float(item) for item in value.get("operating_displacement_m", [])
        ]
        operating_velocity = [
            float(item) for item in value.get("operating_velocity_m_s", [])
        ]
        currents = [float(item) for item in value.get("bias_currents_a", [])]
    except (TypeError, ValueError):
        return False
    if (
        len(stiffness) != 2
        or any(len(row) != 2 for row in stiffness)
        or len(damping) != 2
        or any(len(row) != 2 for row in damping)
    ):
        return False
    stiffness_positive = (
        stiffness[0][0] > 0.0
        and stiffness[1][1] > 0.0
        and stiffness[0][0] * stiffness[1][1]
        - stiffness[0][1] * stiffness[1][0]
        > 0.0
    )
    symmetric_damping_cross = 0.5 * (damping[0][1] + damping[1][0])
    damping_positive = (
        damping[0][0] > 0.0
        and damping[1][1] > 0.0
        and damping[0][0] * damping[1][1]
        - symmetric_damping_cross * symmetric_damping_cross
        > 0.0
    )
    derived_force = [
        [
            -sum(stiffness[row][column] * delta[column] for column in range(2))
            for row in range(2)
        ]
        for delta in displacement
        if len(delta) == 2
    ]
    mirrored_fields = (
        "coordinate_order",
        "displacement_perturbations_m",
        "force_perturbations_n",
        "stiffness_matrix_n_per_m",
        "damping_matrix_n_s_per_m",
        "state_eigenvalues_per_s",
        "operating_displacement_m",
        "operating_velocity_m_s",
        "bias_currents_a",
        "bearing_mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "force_generation",
                "stiffness_generation",
                "damping_generation",
                "coordinate_generation",
                "reciprocity_generation",
                "stability_generation",
                "operating_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and value.get("coordinate_order") == ["x", "y"]
        and len(displacement) == len(force) >= 4
        and all(len(row) == 2 for row in displacement + force)
        and all(math.isfinite(item) for row in displacement + force for item in row)
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for observed_row, expected_row in zip(force, derived_force)
            for observed, expected in zip(observed_row, expected_row)
        )
        and all(math.isfinite(item) for row in stiffness + damping for item in row)
        and math.isclose(stiffness[0][1], stiffness[1][0], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and stiffness_positive
        and damping_positive
        and len(eigenvalues) == 4
        and all(len(row) == 2 for row in eigenvalues)
        and all(math.isfinite(item) for row in eigenvalues for item in row)
        and all(row[0] < 0.0 for row in eigenvalues)
        and len(operating_displacement) == len(operating_velocity) == 2
        and all(
            math.isfinite(item)
            for item in operating_displacement + operating_velocity + currents
        )
        and len(currents) >= 2
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored_fields)
        and _valid_sha256(value.get("bearing_mesh_sha256"))
        and bool(str(value.get("bearing_result_owner", "")).strip())
        and value.get("accepted_bearing_result_owner")
        == value.get("bearing_result_owner")
        and _valid_sha256(value.get("bearing_result_sha256"))
        and value.get("accepted_bearing_result_sha256")
        == value.get("bearing_result_sha256")
    )


def _moving_conductor_drag_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("moving_conductor_generation", "")).strip()
    try:
        velocity = [float(item) for item in value.get("velocity_m_s", [])]
        drag = [float(item) for item in value.get("drag_force_n", [])]
        lift = [float(item) for item in value.get("lift_force_n", [])]
        joule_power = float(value.get("joule_power_w"))
        mechanical_power = float(value.get("mechanical_drag_power_w"))
        conductivity = float(value.get("conductivity_s_m"))
        relative_permeability = float(value.get("relative_permeability"))
        frequency = float(value.get("excitation_frequency_hz"))
        spatial_period = float(value.get("spatial_period_m"))
        slip_frequency = float(value.get("slip_frequency_hz"))
        skin_depth = float(value.get("skin_depth_m"))
    except (TypeError, ValueError):
        return False
    if len(velocity) != 3 or len(drag) != 3 or len(lift) != 3:
        return False
    speed = math.sqrt(sum(item * item for item in velocity))
    drag_work = -sum(force * motion for force, motion in zip(drag, velocity))
    lift_work = sum(force * motion for force, motion in zip(lift, velocity))
    expected_skin_depth = (
        math.sqrt(
            2.0
            / (
                2.0
                * math.pi
                * frequency
                * 4.0e-7
                * math.pi
                * relative_permeability
                * conductivity
            )
        )
        if frequency > 0.0 and conductivity > 0.0 and relative_permeability > 0.0
        else math.nan
    )
    mirrored_fields = (
        "coordinate_frame",
        "velocity_m_s",
        "drag_force_n",
        "lift_force_n",
        "joule_power_w",
        "mechanical_drag_power_w",
        "conductivity_s_m",
        "relative_permeability",
        "excitation_frequency_hz",
        "spatial_period_m",
        "slip_frequency_hz",
        "skin_depth_m",
        "conductor_mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "velocity_generation",
                "frame_generation",
                "force_generation",
                "power_generation",
                "skin_generation",
                "frequency_generation",
                "slip_generation",
                "mesh_generation",
                "owner_generation",
                "field_generation",
                "result_generation",
            )
        )
        and value.get("coordinate_frame") == "global_xyz_right_handed"
        and all(math.isfinite(item) for item in velocity + drag + lift)
        and speed > 0.0
        and sum(force * motion for force, motion in zip(drag, velocity)) < 0.0
        and math.isclose(lift_work, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and all(math.isfinite(item) and item > 0.0 for item in (joule_power, mechanical_power, conductivity, relative_permeability, frequency, spatial_period, slip_frequency, skin_depth))
        and math.isclose(mechanical_power, drag_work, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(joule_power, mechanical_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(slip_frequency, speed / spatial_period, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(frequency, slip_frequency, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(skin_depth, expected_skin_depth, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored_fields)
        and _valid_sha256(value.get("conductor_mesh_sha256"))
        and _valid_sha256(value.get("field_sha256"))
        and value.get("accepted_field_sha256") == value.get("field_sha256")
        and bool(str(value.get("conductor_result_owner", "")).strip())
        and value.get("accepted_conductor_result_owner")
        == value.get("conductor_result_owner")
        and _valid_sha256(value.get("conductor_result_sha256"))
        and value.get("accepted_conductor_result_sha256")
        == value.get("conductor_result_sha256")
    )


def _magnetic_bearing_bias_sweep_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("bearing_generation", "")).strip()
    try:
        bias_current = float(value.get("bias_current_a"))
        displacement = [
            float(item) for item in value.get("displacement_samples_m", [])
        ]
        force = [float(item) for item in value.get("force_x_samples_n", [])]
        stiffness = [
            [float(item) for item in row]
            for row in value.get("stiffness_matrix_n_m", [])
        ]
    except (TypeError, ValueError):
        return False
    if (
        len(displacement) != 3
        or len(force) != 3
        or len(stiffness) != 2
        or any(len(row) != 2 for row in stiffness)
    ):
        return False
    span = displacement[2] - displacement[0]
    derived_stiffness = -(force[2] - force[0]) / span if span > 0.0 else math.nan
    determinant = stiffness[0][0] * stiffness[1][1] - stiffness[0][1] * stiffness[1][0]
    mirrored_fields = (
        "bias_current_a",
        "displacement_samples_m",
        "force_x_samples_n",
        "stiffness_matrix_n_m",
        "coordinate_frame",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "bias_generation",
                "displacement_generation",
                "force_generation",
                "stiffness_generation",
                "crosscoupling_generation",
                "frame_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and math.isfinite(bias_current)
        and bias_current > 0.0
        and all(math.isfinite(item) for item in displacement + force)
        and displacement[0] < displacement[1] < displacement[2]
        and math.isclose(displacement[1], 0.0, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(displacement[0], -displacement[2], rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(force[1], 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(force[0], -force[2], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(math.isfinite(item) for row in stiffness for item in row)
        and stiffness[0][0] > 0.0
        and stiffness[1][1] > 0.0
        and determinant > 0.0
        and math.isclose(stiffness[0][1], stiffness[1][0], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(stiffness[0][0], derived_stiffness, rel_tol=1.0e-12, abs_tol=1.0e-9)
        and value.get("coordinate_frame") == "global_xyz_right_handed"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored_fields)
        and bool(str(value.get("bearing_owner", "")).strip())
        and value.get("accepted_bearing_owner") == value.get("bearing_owner")
        and _valid_sha256(value.get("bearing_result_sha256"))
        and value.get("accepted_bearing_result_sha256")
        == value.get("bearing_result_sha256")
    )


def _pm_demag_recoil_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    try:
        reference_temperature = float(value.get("reference_temperature_c"))
        operating_temperature = float(value.get("operating_temperature_c"))
        remanence_reference = float(value.get("remanence_reference_t"))
        coefficient = float(value.get("remanence_temperature_coefficient_per_c"))
        remanence_temperature = float(value.get("temperature_adjusted_remanence_t"))
        recoil_mu = float(value.get("recoil_relative_permeability"))
        knee_field = float(value.get("knee_field_a_m"))
        h_points = [float(item) for item in value.get("loadline_h_a_m", [])]
        b_points = [float(item) for item in value.get("loadline_b_t", [])]
        irreversible_loss = float(value.get("irreversible_flux_loss_fraction"))
    except (TypeError, ValueError):
        return False
    expected_remanence = remanence_reference * (
        1.0 + coefficient * (operating_temperature - reference_temperature)
    )
    expected_b = [
        remanence_temperature + 4.0e-7 * math.pi * recoil_mu * field
        for field in h_points
    ]
    knee_crossed = bool(h_points) and min(h_points) <= knee_field
    mirrored_fields = (
        "reference_temperature_c",
        "operating_temperature_c",
        "remanence_reference_t",
        "remanence_temperature_coefficient_per_c",
        "temperature_adjusted_remanence_t",
        "recoil_relative_permeability",
        "knee_field_a_m",
        "loadline_h_a_m",
        "loadline_b_t",
        "knee_crossed",
        "irreversible_flux_loss_fraction",
        "field_orientation",
        "mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "recoil_generation",
                "knee_generation",
                "loadline_generation",
                "temperature_generation",
                "irreversible_generation",
                "orientation_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and all(
            math.isfinite(item)
            for item in (
                reference_temperature,
                operating_temperature,
                remanence_reference,
                coefficient,
                remanence_temperature,
                recoil_mu,
                knee_field,
                irreversible_loss,
            )
        )
        and remanence_reference > 0.0
        and recoil_mu > 0.0
        and coefficient <= 0.0
        and math.isclose(remanence_temperature, expected_remanence, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and len(h_points) >= 3
        and len(b_points) == len(h_points)
        and all(math.isfinite(item) for item in h_points + b_points)
        and all(right > left for left, right in zip(h_points, h_points[1:]))
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for observed, expected in zip(b_points, expected_b)
        )
        and value.get("knee_crossed") is knee_crossed
        and 0.0 <= irreversible_loss <= 1.0
        and ((knee_crossed and irreversible_loss > 0.0) or (not knee_crossed and irreversible_loss == 0.0))
        and value.get("field_orientation") == "magnetization_antiparallel_h"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored_fields)
        and _valid_sha256(value.get("mesh_sha256"))
        and bool(str(value.get("demag_owner", "")).strip())
        and value.get("accepted_demag_owner") == value.get("demag_owner")
        and _valid_sha256(value.get("demag_result_sha256"))
        and value.get("accepted_demag_result_sha256") == value.get("demag_result_sha256")
    )


def _magnetic_gear_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("magnetic_gear_generation", "")).strip()
    try:
        high_poles = int(value.get("high_speed_pole_pairs"))
        low_poles = int(value.get("low_speed_pole_pairs"))
        modulator_poles = int(value.get("modulator_pole_count"))
        harmonic = int(value.get("transmitted_harmonic_order"))
        high_torque = float(value.get("high_speed_torque_nm"))
        low_torque = float(value.get("low_speed_torque_nm"))
        high_speed = float(value.get("high_speed_angular_velocity_rad_s"))
        low_speed = float(value.get("low_speed_angular_velocity_rad_s"))
        high_phase = float(value.get("high_speed_harmonic_phase_rad"))
        low_phase = float(value.get("low_speed_harmonic_phase_rad"))
        modulator_phase = float(value.get("modulator_phase_rad"))
        transmitted_phase = float(value.get("transmitted_phase_rad"))
    except (TypeError, ValueError):
        return False
    numeric = (
        high_torque,
        low_torque,
        high_speed,
        low_speed,
        high_phase,
        low_phase,
        modulator_phase,
        transmitted_phase,
    )
    ratio = low_poles / high_poles if high_poles > 0 else math.nan
    phase_error = math.remainder(
        transmitted_phase - (high_phase - low_phase + modulator_phase),
        2.0 * math.pi,
    )
    mirrored_fields = (
        "high_speed_pole_pairs",
        "low_speed_pole_pairs",
        "modulator_pole_count",
        "transmitted_harmonic_order",
        "high_speed_torque_nm",
        "low_speed_torque_nm",
        "high_speed_angular_velocity_rad_s",
        "low_speed_angular_velocity_rad_s",
        "high_speed_harmonic_phase_rad",
        "low_speed_harmonic_phase_rad",
        "modulator_phase_rad",
        "transmitted_phase_rad",
        "coordinate_frame",
        "gear_mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "pole_generation",
                "harmonic_generation",
                "torque_generation",
                "phase_generation",
                "power_generation",
                "frame_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        )
        and all(
            isinstance(value.get(key), int) and not isinstance(value.get(key), bool)
            for key in (
                "high_speed_pole_pairs",
                "low_speed_pole_pairs",
                "modulator_pole_count",
                "transmitted_harmonic_order",
            )
        )
        and high_poles > 0
        and low_poles > high_poles
        and modulator_poles == high_poles + low_poles
        and harmonic == low_poles
        and all(math.isfinite(item) for item in numeric)
        and high_torque * low_torque < 0.0
        and high_speed > 0.0
        and low_speed > 0.0
        and math.isclose(abs(low_torque / high_torque), ratio, rel_tol=1.0e-12)
        and math.isclose(high_speed / low_speed, ratio, rel_tol=1.0e-12)
        and math.isclose(
            high_torque * high_speed + low_torque * low_speed,
            0.0,
            rel_tol=0.0,
            abs_tol=1.0e-9,
        )
        and math.isclose(phase_error, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and value.get("coordinate_frame") == "global_xyz_right_handed"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored_fields)
        and _valid_sha256(value.get("gear_mesh_sha256"))
        and bool(str(value.get("gear_result_owner", "")).strip())
        and value.get("accepted_gear_result_owner") == value.get("gear_result_owner")
        and _valid_sha256(value.get("gear_result_sha256"))
        and value.get("accepted_gear_result_sha256") == value.get("gear_result_sha256")
    )


def _demag_bem_charge_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("demag_bem_generation", "")).strip()
    try:
        areas = [float(item) for item in value.get("panel_areas_m2", [])]
        charges = [
            float(item) for item in value.get("surface_charge_density_a_m", [])
        ]
        charge_integral = float(value.get("surface_charge_integral_a_m"))
        normals = [
            [float(item) for item in row] for row in value.get("outward_normals", [])
        ]
        jump = [float(item) for item in value.get("normal_field_jump_a_m", [])]
        radii = [float(item) for item in value.get("farfield_radius_m", [])]
        potential = [float(item) for item in value.get("farfield_potential_a", [])]
        field = [float(item) for item in value.get("farfield_field_a_m", [])]
        energy = float(value.get("magnetic_energy_j"))
    except (TypeError, ValueError):
        return False
    count = len(areas)
    neutral_integral = sum(area * charge for area, charge in zip(areas, charges))
    potential_scale = [item * radius**2 for item, radius in zip(potential, radii)]
    field_scale = [item * radius**3 for item, radius in zip(field, radii)]
    mirrored_fields = (
        "panel_areas_m2",
        "surface_charge_density_a_m",
        "surface_charge_integral_a_m",
        "outward_normals",
        "outward_orientation_verified",
        "normal_field_jump_a_m",
        "farfield_radius_m",
        "farfield_potential_a",
        "farfield_field_a_m",
        "magnetic_energy_j",
        "boundary_mesh_sha256",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "charge_generation",
                "normal_generation",
                "jump_generation",
                "farfield_generation",
                "energy_generation",
                "mesh_generation",
                "owner_generation",
                "solution_generation",
            )
        )
        and count >= 4
        and len(charges) == len(normals) == len(jump) == count
        and all(area > 0.0 and math.isfinite(area) for area in areas)
        and all(math.isfinite(charge) for charge in charges)
        and all(len(row) == 3 for row in normals)
        and all(math.isfinite(item) for row in normals for item in row)
        and all(
            math.isclose(
                sum(item * item for item in row), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12
            )
            for row in normals
        )
        and value.get("outward_orientation_verified") is True
        and all(
            math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for observed, expected in zip(jump, charges)
        )
        and math.isclose(charge_integral, neutral_integral, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(neutral_integral, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and len(radii) == len(potential) == len(field) >= 3
        and all(math.isfinite(item) and item > 0.0 for item in radii + potential + field)
        and all(right > left for left, right in zip(radii, radii[1:]))
        and all(right < left for left, right in zip(potential, potential[1:]))
        and all(right < left for left, right in zip(field, field[1:]))
        and all(
            math.isclose(item, potential_scale[0], rel_tol=1.0e-12, abs_tol=1.0e-12)
            for item in potential_scale
        )
        and all(
            math.isclose(item, field_scale[0], rel_tol=1.0e-12, abs_tol=1.0e-12)
            for item in field_scale
        )
        and math.isfinite(energy)
        and energy > 0.0
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored_fields)
        and _valid_sha256(value.get("boundary_mesh_sha256"))
        and bool(str(value.get("demag_solution_owner", "")).strip())
        and value.get("accepted_demag_solution_owner") == value.get("demag_solution_owner")
        and _valid_sha256(value.get("demag_solution_sha256"))
        and value.get("accepted_demag_solution_sha256") == value.get("demag_solution_sha256")
    )


def _maglev_dynamic_stiffness_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("maglev_generation", "")).strip()
    try:
        bias = float(value.get("bias_current_a"))
        gap = float(value.get("equilibrium_gap_m"))
        equilibrium_force = float(value.get("equilibrium_force_n"))
        load = float(value.get("supported_load_n"))
        frequency = float(value.get("excitation_frequency_hz"))
        stiffness = [float(item) for item in value.get("complex_stiffness_n_m", [])]
        damping = float(value.get("viscous_damping_n_s_m"))
        displacement = [float(item) for item in value.get("displacement_phasor_m", [])]
        force = [float(item) for item in value.get("force_phasor_n", [])]
        phase = float(value.get("force_displacement_phase_rad"))
    except (TypeError, ValueError):
        return False
    fields = (
        "bias_current_a", "equilibrium_gap_m", "equilibrium_force_n",
        "supported_load_n", "excitation_frequency_hz", "complex_stiffness_n_m",
        "viscous_damping_n_s_m", "displacement_phasor_m", "force_phasor_n",
        "force_displacement_phase_rad", "coordinate_frame",
    )
    values = [bias, gap, equilibrium_force, load, frequency, damping, phase]
    if (
        len(stiffness) != 2 or len(displacement) != 2 or len(force) != 2
        or not all(math.isfinite(item) for item in values + stiffness + displacement + force)
    ):
        return False
    omega = 2.0 * math.pi * frequency
    expected_force = [
        stiffness[0] * displacement[0] - stiffness[1] * displacement[1],
        stiffness[0] * displacement[1] + stiffness[1] * displacement[0],
    ]
    expected_phase = math.atan2(force[1], force[0]) - math.atan2(
        displacement[1], displacement[0]
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "bias_generation", "equilibrium_generation", "frequency_generation",
            "stiffness_generation", "damping_generation", "force_generation",
            "displacement_generation", "phase_generation", "frame_generation",
            "owner_generation", "result_generation"))
        and bias > 0.0 and gap > 0.0 and frequency > 0.0
        and equilibrium_force > 0.0
        and math.isclose(equilibrium_force, load, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and stiffness[0] > 0.0 and damping >= 0.0
        and math.isclose(stiffness[1], omega * damping, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.hypot(*displacement) > 0.0
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
                for item, expected in zip(force, expected_force))
        and math.isclose(phase, expected_phase, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and value.get("coordinate_frame") == "global_z_up_force_positive"
        and all(value.get(f"result_{field}") == value.get(field) for field in fields)
        and bool(str(value.get("maglev_owner", "")).strip())
        and value.get("accepted_maglev_owner") == value.get("maglev_owner")
        and _valid_sha256(value.get("maglev_result_sha256"))
        and value.get("accepted_maglev_result_sha256")
        == value.get("maglev_result_sha256")
    )


def _bem_demag_reciprocity_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("demag_generation", "")).strip()
    try:
        energy_12 = float(value.get("interaction_energy_12_j"))
        energy_21 = float(value.get("interaction_energy_21_j"))
        field_12 = [float(item) for item in value.get("field_1_due_2_a_m", [])]
        field_21 = [float(item) for item in value.get("field_2_due_1_a_m", [])]
        magnetization_1 = [float(item) for item in value.get("magnetization_1_a_m", [])]
        magnetization_2 = [float(item) for item in value.get("magnetization_2_a_m", [])]
        volumes = [float(item) for item in value.get("region_volumes_m3", [])]
    except (TypeError, ValueError):
        return False
    vectors = (field_12, field_21, magnetization_1, magnetization_2)
    mirrored = (
        "interaction_energy_12_j", "interaction_energy_21_j",
        "field_1_due_2_a_m", "field_2_due_1_a_m", "magnetization_1_a_m",
        "magnetization_2_a_m", "surface_orientation", "region_volumes_m3",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "reciprocity_generation", "energy_generation", "field_generation",
            "magnetization_generation", "surface_generation", "volume_generation",
            "mesh_generation", "solution_generation", "result_generation"))
        and all(len(vector) == 3 for vector in vectors)
        and len(volumes) == 2
        and all(math.isfinite(item) for vector in vectors for item in vector)
        and all(math.isfinite(item) and item > 0.0 for item in volumes)
        and math.isfinite(energy_12) and math.isfinite(energy_21)
        and math.isclose(energy_12, energy_21, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and sum(a * b for a, b in zip(magnetization_1, field_12)) < 0.0
        and sum(a * b for a, b in zip(magnetization_2, field_21)) < 0.0
        and value.get("surface_orientation") == "outward_right_handed"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and bool(str(value.get("solution_owner", "")).strip())
        and value.get("accepted_solution_owner") == value.get("solution_owner")
        and _valid_sha256(value.get("demag_result_sha256"))
        and value.get("accepted_demag_result_sha256")
        == value.get("demag_result_sha256")
    )


def _eddy_maglev_power_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("maglev_generation", "")).strip()
    try:
        velocity = float(value.get("plate_velocity_m_s"))
        pole_pitch = float(value.get("pole_pitch_m"))
        frequency = float(value.get("excitation_frequency_hz"))
        conductivity = float(value.get("plate_conductivity_s_m"))
        relative_permeability = float(value.get("relative_permeability"))
        skin_depth = float(value.get("skin_depth_m"))
        lift = float(value.get("lift_force_n"))
        drag = float(value.get("drag_force_n"))
        joule_loss = float(value.get("joule_loss_w"))
        drag_power = float(value.get("mechanical_drag_power_w"))
        residual = float(value.get("power_balance_residual_w"))
        tolerance = float(value.get("power_tolerance_w"))
    except (TypeError, ValueError):
        return False
    mirrored = (
        "plate_velocity_m_s", "pole_pitch_m", "excitation_frequency_hz",
        "plate_conductivity_s_m", "relative_permeability", "skin_depth_m",
        "lift_force_n", "drag_force_n", "joule_loss_w",
        "mechanical_drag_power_w", "power_balance_residual_w",
        "power_tolerance_w",
    )
    values = (
        velocity, pole_pitch, frequency, conductivity, relative_permeability,
        skin_depth, lift, drag, joule_loss, drag_power, residual, tolerance,
    )
    if not all(math.isfinite(item) for item in values):
        return False
    expected_frequency = velocity / pole_pitch if pole_pitch > 0.0 else math.nan
    expected_skin_depth = (
        math.sqrt(
            2.0
            / (
                2.0
                * math.pi
                * frequency
                * (4.0e-7 * math.pi)
                * relative_permeability
                * conductivity
            )
        )
        if frequency > 0.0 and conductivity > 0.0 and relative_permeability > 0.0
        else math.nan
    )
    expected_drag_power = drag * velocity
    expected_residual = drag_power - joule_loss
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "velocity_generation", "frequency_generation",
            "conductivity_generation", "skin_generation", "force_generation",
            "loss_generation", "power_generation", "mesh_generation",
            "owner_generation", "result_generation"))
        and velocity > 0.0 and pole_pitch > 0.0 and frequency > 0.0
        and conductivity > 0.0 and relative_permeability > 0.0
        and skin_depth > 0.0 and lift > 0.0 and drag > 0.0
        and joule_loss >= 0.0 and drag_power >= 0.0 and tolerance >= 0.0
        and math.isclose(frequency, expected_frequency, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(skin_depth, expected_skin_depth, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(drag_power, expected_drag_power, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(joule_loss, drag_power, rel_tol=1.0e-12, abs_tol=tolerance)
        and math.isclose(residual, expected_residual, rel_tol=0.0, abs_tol=1.0e-12)
        and abs(residual) <= tolerance
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _valid_sha256(value.get("maglev_result_sha256"))
        and value.get("accepted_maglev_result_sha256")
        == value.get("maglev_result_sha256")
    )


def _pm_coupling_energy_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    pole_pairs_value = value.get("pole_pairs")
    try:
        angle = float(value.get("relative_angle_rad"))
        pole_pairs = int(pole_pairs_value)
        period = float(value.get("pole_period_rad"))
        delta = float(value.get("angle_perturbation_rad"))
        energy_minus = float(value.get("energy_minus_j"))
        energy_center = float(value.get("energy_center_j"))
        energy_plus = float(value.get("energy_plus_j"))
        periodic_energy = float(value.get("periodic_energy_j"))
        derivative_torque = float(value.get("energy_derivative_torque_nm"))
        driver_torque = float(value.get("driver_torque_nm"))
        driven_torque = float(value.get("driven_torque_nm"))
    except (TypeError, ValueError):
        return False
    mirrored = (
        "relative_angle_rad", "pole_pairs", "pole_period_rad",
        "angle_perturbation_rad", "energy_minus_j", "energy_center_j",
        "energy_plus_j", "periodic_energy_j", "energy_derivative_torque_nm",
        "driver_torque_nm", "driven_torque_nm", "torque_frame",
    )
    values = (
        angle, period, delta, energy_minus, energy_center, energy_plus,
        periodic_energy, derivative_torque, driver_torque, driven_torque,
    )
    if not all(math.isfinite(item) for item in values):
        return False
    expected_period = 2.0 * math.pi / pole_pairs if pole_pairs > 0 else math.nan
    expected_torque = -(energy_plus - energy_minus) / (2.0 * delta) if delta > 0.0 else math.nan
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "angle_generation", "periodicity_generation", "energy_generation",
            "derivative_generation", "torque_generation", "reaction_generation",
            "frame_generation", "mesh_generation", "owner_generation",
            "result_generation"))
        and isinstance(pole_pairs_value, int)
        and not isinstance(pole_pairs_value, bool)
        and pole_pairs > 0 and delta > 0.0
        and math.isclose(period, expected_period, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(periodic_energy, energy_center, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(derivative_torque, expected_torque, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(driver_torque, derivative_torque, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(driven_torque, -driver_torque, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(driver_torque + driven_torque, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and value.get("torque_frame") == "relative_angle_driver_positive"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _valid_sha256(value.get("coupling_result_sha256"))
        and value.get("accepted_coupling_result_sha256")
        == value.get("coupling_result_sha256")
    )


def _thin_conductor_surface_impedance_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("thin_generation", "")).strip()
    try:
        frequency = float(value.get("frequency_hz"))
        conductivity = float(value.get("conductivity_s_m"))
        mu_r = float(value.get("relative_permeability"))
        skin_depth = float(value.get("skin_depth_m"))
        impedance = [float(item) for item in value.get("surface_impedance_ohm", [])]
        current = [float(item) for item in value.get("sheet_current_peak_a_m", [])]
        field_jump = [float(item) for item in value.get("tangential_field_jump_peak_a_m", [])]
        area = float(value.get("surface_area_m2"))
        loss = float(value.get("joule_loss_w"))
        reactive = float(value.get("reactive_power_var"))
    except (TypeError, ValueError):
        return False
    mu = 4.0e-7 * math.pi * mu_r
    expected_skin = math.sqrt(2.0 / (2.0 * math.pi * frequency * mu * conductivity)) if frequency > 0.0 and mu > 0.0 and conductivity > 0.0 else math.nan
    expected_resistance = math.sqrt(math.pi * frequency * mu / conductivity) if frequency > 0.0 and mu > 0.0 and conductivity > 0.0 else math.nan
    current_norm_sq = sum(item * item for item in current)
    expected_loss = 0.5 * expected_resistance * current_norm_sq * area
    mirrored = (
        "frequency_hz", "conductivity_s_m", "relative_permeability",
        "skin_depth_m", "surface_impedance_ohm", "sheet_current_peak_a_m",
        "tangential_field_jump_peak_a_m", "surface_area_m2", "joule_loss_w",
        "reactive_power_var",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "impedance_generation", "skin_generation", "current_generation",
            "field_generation", "power_generation", "surface_generation",
            "owner_generation", "result_generation"))
        and frequency > 0.0 and conductivity > 0.0 and mu_r > 0.0
        and math.isfinite(skin_depth) and skin_depth > 0.0
        and math.isclose(skin_depth, expected_skin, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and len(impedance) == 2
        and math.isclose(impedance[0], expected_resistance, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(impedance[1], expected_resistance, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and len(current) == len(field_jump) == 2
        and all(math.isfinite(item) for item in current + field_jump)
        and current == field_jump and current_norm_sq > 0.0
        and math.isfinite(area) and area > 0.0
        and math.isclose(loss, expected_loss, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(reactive, expected_loss, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("surface_owner", "")).startswith("surface:")
        and value.get("accepted_surface_owner") == value.get("surface_owner")
        and _valid_sha256(value.get("thin_result_sha256"))
        and value.get("accepted_thin_result_sha256") == value.get("thin_result_sha256")
    )


def _magnetic_gear_action_reaction_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("gear_generation", "")).strip()
    try:
        high_poles = int(value.get("high_speed_pole_pairs"))
        low_poles = int(value.get("low_speed_pole_pairs"))
        modulators = int(value.get("modulator_segment_count"))
        harmonic = int(value.get("working_harmonic_order"))
        phase = float(value.get("mechanical_phase_rad"))
        high_speed = float(value.get("high_speed_rad_s"))
        low_speed = float(value.get("low_speed_rad_s"))
        ratio = float(value.get("gear_ratio"))
        high_torque = float(value.get("high_speed_torque_nm"))
        low_torque = float(value.get("low_speed_torque_nm"))
        reaction = float(value.get("modulator_reaction_torque_nm"))
        residual = float(value.get("power_balance_residual_w"))
    except (TypeError, ValueError):
        return False
    mirrored = (
        "high_speed_pole_pairs", "low_speed_pole_pairs",
        "modulator_segment_count", "working_harmonic_order",
        "mechanical_phase_rad", "high_speed_rad_s", "low_speed_rad_s",
        "gear_ratio", "high_speed_torque_nm", "low_speed_torque_nm",
        "modulator_reaction_torque_nm", "power_balance_residual_w",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "harmonic_generation", "pole_generation", "phase_generation",
            "ratio_generation", "torque_generation", "reaction_generation",
            "power_generation", "owner_generation", "result_generation"))
        and high_poles > 0 and low_poles > 0
        and modulators == high_poles + low_poles and harmonic == modulators
        and math.isfinite(phase) and 0.0 <= phase < 2.0 * math.pi / modulators
        and all(math.isfinite(item) for item in (high_speed, low_speed, ratio, high_torque, low_torque, reaction, residual))
        and high_speed != 0.0 and low_speed != 0.0
        and math.isclose(high_poles * high_speed + low_poles * low_speed, 0.0, rel_tol=0.0, abs_tol=1.0e-12)
        and math.isclose(ratio, low_speed / high_speed, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(ratio, -high_poles / low_poles, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and high_torque > 0.0 and low_torque > 0.0
        and math.isclose(reaction, -(high_torque + low_torque), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(residual, high_torque * high_speed + low_torque * low_speed, rel_tol=0.0, abs_tol=1.0e-12)
        and abs(residual) <= 1.0e-12
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("model_owner", "")).startswith("gear:")
        and value.get("accepted_model_owner") == value.get("model_owner")
        and _valid_sha256(value.get("gear_result_sha256"))
        and value.get("accepted_gear_result_sha256") == value.get("gear_result_sha256")
    )


def _multilayer_shield_closure_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("shield_generation", "")).strip()
    try:
        permeability = [float(item) for item in value.get("relative_permeability", [])]
        thickness = [float(item) for item in value.get("layer_thickness_m", [])]
        radii = [float(item) for item in value.get("layer_mean_radius_m", [])]
        factors = [float(item) for item in value.get("layer_shielding_factor", [])]
        flux = [float(item) for item in value.get("interface_normal_flux_t", [])]
        external_field = float(value.get("external_field_t"))
        attenuation = float(value.get("attenuation_factor"))
        leakage = float(value.get("leakage_field_t"))
        volume = float(value.get("cavity_volume_m3"))
        energy = float(value.get("stored_energy_j"))
    except (TypeError, ValueError):
        return False
    expected_factors = [
        1.0 + (mu_r - 1.0) * layer_thickness / (2.0 * radius)
        for mu_r, layer_thickness, radius in zip(permeability, thickness, radii)
    ]
    expected_attenuation = math.prod(expected_factors)
    expected_leakage = external_field / expected_attenuation
    expected_energy = expected_leakage**2 * volume / (2.0 * 4.0e-7 * math.pi)
    mirrored = (
        "relative_permeability", "layer_thickness_m", "layer_mean_radius_m",
        "layer_shielding_factor", "interface_normal_flux_t", "external_field_t",
        "attenuation_factor", "leakage_field_t", "cavity_volume_m3",
        "stored_energy_j",
    )
    count = len(permeability)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "material_generation", "thickness_generation", "geometry_generation",
                "flux_generation", "attenuation_generation", "field_generation",
                "energy_generation", "owner_generation", "result_generation",
            )
        )
        and count >= 2 and len(thickness) == len(radii) == len(factors) == count
        and len(flux) == count + 1
        and all(math.isfinite(item) for item in permeability + thickness + radii + factors + flux)
        and all(item > 1.0 for item in permeability)
        and all(item > 0.0 for item in thickness + radii)
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-12) for item, expected in zip(factors, expected_factors))
        and external_field > 0.0 and attenuation > 1.0 and leakage > 0.0
        and math.isclose(attenuation, expected_attenuation, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(leakage, expected_leakage, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(math.isclose(item, leakage, rel_tol=1.0e-12, abs_tol=1.0e-15) for item in flux)
        and volume > 0.0 and energy > 0.0
        and math.isclose(energy, expected_energy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("geometry_owner", "")).startswith("geometry:")
        and value.get("accepted_geometry_owner") == value.get("geometry_owner")
        and _valid_sha256(value.get("shield_result_sha256"))
        and value.get("accepted_shield_result_sha256") == value.get("shield_result_sha256")
    )


def _transformer_energy_force_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("transformer_generation", "")).strip()
    try:
        matrix = [[float(item) for item in row] for row in value.get("inductance_matrix_h", [])]
        leakage = float(value.get("primary_leakage_inductance_h"))
        currents = [float(item) for item in value.get("winding_currents_a", [])]
        linkages = [float(item) for item in value.get("flux_linkages_wb_turn", [])]
        reciprocity = float(value.get("reciprocity_residual_h"))
        coenergy = float(value.get("coenergy_j"))
        mutual_gradient = float(value.get("mutual_inductance_gradient_h_per_m"))
        force = float(value.get("force_n"))
    except (TypeError, ValueError):
        return False
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix) or len(currents) != 2 or len(linkages) != 2:
        return False
    l_primary, mutual_12 = matrix[0]
    mutual_21, l_secondary = matrix[1]
    if l_primary <= 0.0 or l_secondary <= 0.0:
        return False
    expected_leakage = l_primary - mutual_12**2 / l_secondary
    expected_linkages = [
        l_primary * currents[0] + mutual_12 * currents[1],
        mutual_21 * currents[0] + l_secondary * currents[1],
    ]
    expected_coenergy = 0.5 * sum(current * linkage for current, linkage in zip(currents, expected_linkages))
    expected_force = mutual_gradient * currents[0] * currents[1]
    mirrored = (
        "inductance_matrix_h", "primary_leakage_inductance_h",
        "winding_currents_a", "flux_linkages_wb_turn", "reciprocity_residual_h",
        "coenergy_j", "mutual_inductance_gradient_h_per_m", "force_n",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "inductance_generation", "leakage_generation", "flux_generation",
                "reciprocity_generation", "energy_generation", "force_generation",
                "winding_generation", "result_generation",
            )
        )
        and all(math.isfinite(item) for row in matrix for item in row)
        and all(math.isfinite(item) for item in currents + linkages)
        and math.isclose(mutual_12, mutual_21, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and l_primary * l_secondary - mutual_12**2 >= -1.0e-15
        and leakage >= 0.0
        and math.isclose(leakage, expected_leakage, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(math.isclose(item, expected, rel_tol=1.0e-12, abs_tol=1.0e-15) for item, expected in zip(linkages, expected_linkages))
        and math.isfinite(reciprocity) and abs(reciprocity) <= 1.0e-12
        and math.isfinite(coenergy) and coenergy >= 0.0
        and math.isclose(coenergy, expected_coenergy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isfinite(mutual_gradient) and math.isfinite(force)
        and math.isclose(force, expected_force, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("winding_owner", "")).startswith("winding:")
        and value.get("accepted_winding_owner") == value.get("winding_owner")
        and _valid_sha256(value.get("transformer_result_sha256"))
        and value.get("accepted_transformer_result_sha256") == value.get("transformer_result_sha256")
    )


def _maglev_equilibrium_energy_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("maglev_generation", "")).strip()
    try:
        air_gap = float(value.get("air_gap_m"))
        positions = [float(item) for item in value.get("sample_position_m", [])]
        force = [float(item) for item in value.get("force_n", [])]
        equilibrium_position = float(value.get("equilibrium_position_m"))
        equilibrium_force = float(value.get("equilibrium_force_n"))
        force_gradient = float(value.get("force_gradient_n_per_m"))
        stiffness = float(value.get("stiffness_n_per_m"))
        energy = [float(item) for item in value.get("potential_energy_j", [])]
        energy_curvature = float(value.get("energy_curvature_j_per_m2"))
    except (TypeError, ValueError):
        return False
    if len(positions) != 3 or len(force) != 3 or len(energy) != 3:
        return False
    equilibrium_index = min(
        range(len(positions)), key=lambda index: abs(positions[index] - equilibrium_position)
    )
    expected_force = [
        -stiffness * (position - equilibrium_position) for position in positions
    ]
    energy_offset = energy[equilibrium_index]
    expected_energy = [
        energy_offset + 0.5 * stiffness * (position - equilibrium_position) ** 2
        for position in positions
    ]
    expected_gradient = (force[2] - force[0]) / (positions[2] - positions[0])
    left_step = positions[1] - positions[0]
    right_step = positions[2] - positions[1]
    expected_curvature = (
        2.0
        * (
            (energy[2] - energy[1]) / right_step
            - (energy[1] - energy[0]) / left_step
        )
        / (left_step + right_step)
    )
    mirrored = (
        "air_gap_m", "sample_position_m", "force_n", "equilibrium_position_m",
        "equilibrium_force_n", "force_gradient_n_per_m", "stiffness_n_per_m",
        "potential_energy_j", "energy_curvature_j_per_m2", "stability",
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "gap_generation", "position_generation", "force_generation",
                "gradient_generation", "stiffness_generation", "energy_generation",
                "stability_generation", "geometry_generation", "result_generation",
            )
        )
        and all(math.isfinite(item) for item in positions + force + energy)
        and all(
            math.isfinite(item)
            for item in (
                air_gap, equilibrium_position, equilibrium_force, force_gradient,
                stiffness, energy_curvature,
            )
        )
        and air_gap > 0.0
        and max(abs(position - equilibrium_position) for position in positions) < air_gap
        and all(left < right for left, right in zip(positions, positions[1:]))
        and math.isclose(positions[equilibrium_index], equilibrium_position, rel_tol=0.0, abs_tol=1.0e-15)
        and math.isclose(force[equilibrium_index], equilibrium_force, rel_tol=0.0, abs_tol=1.0e-12)
        and abs(equilibrium_force) <= 1.0e-12
        and all(
            math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
            for actual, expected in zip(force, expected_force)
        )
        and math.isclose(force_gradient, expected_gradient, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(stiffness, -force_gradient, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and stiffness > 0.0
        and all(
            math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for actual, expected in zip(energy, expected_energy)
        )
        and energy[equilibrium_index] == min(energy)
        and math.isclose(energy_curvature, expected_curvature, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(energy_curvature, stiffness, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and value.get("stability") == "stable"
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and str(value.get("geometry_owner", "")).startswith("geometry:")
        and value.get("accepted_geometry_owner") == value.get("geometry_owner")
        and _valid_sha256(value.get("maglev_result_sha256"))
        and value.get("accepted_maglev_result_sha256") == value.get("maglev_result_sha256")
    )


def _eddy_shield_frequency_identity_ok(value: object) -> bool:
    if value is None:
        return True
    if not isinstance(value, dict):
        return False
    generation = str(value.get("eddy_shield_generation", "")).strip()
    fields = (
        "frequency_hz", "angular_frequency_rad_s", "relative_permeability",
        "conductivity_s_per_m", "skin_depth_m", "shield_thickness_m",
        "attenuation_factor", "phase_lag_rad", "incident_field_t",
        "transmitted_field_t", "surface_resistance_ohm", "shield_area_m2",
        "eddy_loss_w", "shield_volume_m3", "stored_energy_j",
    )
    try:
        numbers = {field: float(value.get(field)) for field in fields}
    except (TypeError, ValueError):
        return False
    mu0 = 4.0e-7 * math.pi
    permeability = mu0 * numbers["relative_permeability"]
    expected_omega = 2.0 * math.pi * numbers["frequency_hz"]
    expected_skin_depth = math.sqrt(
        2.0
        / (
            numbers["angular_frequency_rad_s"]
            * permeability
            * numbers["conductivity_s_per_m"]
        )
    )
    expected_attenuation = math.exp(
        -numbers["shield_thickness_m"] / numbers["skin_depth_m"]
    )
    expected_phase = -numbers["shield_thickness_m"] / numbers["skin_depth_m"]
    expected_transmitted = numbers["incident_field_t"] * numbers["attenuation_factor"]
    expected_surface_resistance = 1.0 / (
        numbers["conductivity_s_per_m"] * numbers["skin_depth_m"]
    )
    magnetic_field = numbers["incident_field_t"] / permeability
    expected_loss = (
        0.5
        * numbers["surface_resistance_ohm"]
        * magnetic_field**2
        * numbers["shield_area_m2"]
    )
    expected_volume = numbers["shield_area_m2"] * numbers["shield_thickness_m"]
    expected_energy = (
        numbers["transmitted_field_t"] ** 2
        * numbers["shield_volume_m3"]
        / (2.0 * permeability)
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "frequency_generation", "material_generation", "skin_depth_generation",
                "thickness_generation", "phase_generation", "field_generation",
                "loss_generation", "energy_generation", "geometry_generation",
                "result_generation",
            )
        )
        and all(math.isfinite(item) for item in numbers.values())
        and numbers["frequency_hz"] > 0.0
        and math.isclose(numbers["angular_frequency_rad_s"], expected_omega, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and numbers["relative_permeability"] > 0.0
        and numbers["conductivity_s_per_m"] > 0.0
        and math.isclose(numbers["skin_depth_m"], expected_skin_depth, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and numbers["shield_thickness_m"] > 0.0
        and math.isclose(numbers["attenuation_factor"], expected_attenuation, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and 0.0 < numbers["attenuation_factor"] < 1.0
        and math.isclose(numbers["phase_lag_rad"], expected_phase, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and numbers["incident_field_t"] > 0.0
        and math.isclose(numbers["transmitted_field_t"], expected_transmitted, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(numbers["surface_resistance_ohm"], expected_surface_resistance, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and numbers["shield_area_m2"] > 0.0
        and math.isclose(numbers["eddy_loss_w"], expected_loss, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(numbers["shield_volume_m3"], expected_volume, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and math.isclose(numbers["stored_energy_j"], expected_energy, rel_tol=1.0e-12, abs_tol=1.0e-15)
        and all(value.get(f"result_{field}") == value.get(field) for field in fields)
        and str(value.get("geometry_owner", "")).startswith("geometry:")
        and value.get("accepted_geometry_owner") == value.get("geometry_owner")
        and _valid_sha256(value.get("eddy_shield_result_sha256"))
        and value.get("accepted_eddy_shield_result_sha256") == value.get("eddy_shield_result_sha256")
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
    bem_demag_surface_material_frame_generation_identity_ok = True
    linear_motor_thrust_ripple_generation_identity_ok = True
    levitation_gradient_energy_identity_ok = True
    cogging_periodic_interpolation_identity_ok = True
    bem_panel_self_term_energy_force_identity_ok = True
    motor_reduced_basis_torque_identity_ok = True
    maglev_force_stiffness_identity_ok = True
    motor_coenergy_torque_identity_ok = True
    airgap_stress_harmonic_torque_identity_ok = True
    laminated_core_loss_identity_ok = True
    magnet_demag_volume_fraction_identity_ok = True
    linear_motor_wave_end_effect_identity_ok = True
    maglev_stiffness_identity_ok = True
    cogging_torque_sampling_identity_ok = True
    bem_near_singular_force_identity_ok = True
    hysteresis_minor_loop_identity_ok = True
    maglev_equilibrium_closure_identity_ok = True
    bem_surface_charge_closure_identity_ok = True
    halbach_harmonic_closure_identity_ok = True
    magnetic_bearing_linearization_identity_ok = True
    magnetic_bearing_dynamic_identity_ok = True
    moving_conductor_drag_identity_ok = True
    magnetic_gear_identity_ok = True
    demag_bem_charge_identity_ok = True
    magnetic_bearing_bias_sweep_identity_ok = True
    pm_demag_recoil_identity_ok = True
    maglev_dynamic_stiffness_identity_ok = True
    bem_demag_reciprocity_identity_ok = True
    eddy_maglev_power_identity_ok = True
    pm_coupling_energy_identity_ok = True
    thin_conductor_surface_impedance_identity_ok = True
    magnetic_gear_action_reaction_identity_ok = True
    multilayer_shield_closure_identity_ok = True
    transformer_energy_force_identity_ok = True
    maglev_equilibrium_energy_identity_ok = True
    eddy_shield_frequency_identity_ok = True
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
        bem_demag_surface_material_frame_generation_identity_ok = False
        linear_motor_thrust_ripple_generation_identity_ok = False
        levitation_gradient_energy_identity_ok = False
        cogging_periodic_interpolation_identity_ok = False
        bem_panel_self_term_energy_force_identity_ok = False
        motor_reduced_basis_torque_identity_ok = False
        maglev_force_stiffness_identity_ok = False
        motor_coenergy_torque_identity_ok = False
        airgap_stress_harmonic_torque_identity_ok = False
        laminated_core_loss_identity_ok = False
        magnet_demag_volume_fraction_identity_ok = False
        linear_motor_wave_end_effect_identity_ok = False
        maglev_stiffness_identity_ok = False
        cogging_torque_sampling_identity_ok = False
        bem_near_singular_force_identity_ok = False
        hysteresis_minor_loop_identity_ok = False
        maglev_equilibrium_closure_identity_ok = False
        bem_surface_charge_closure_identity_ok = False
        halbach_harmonic_closure_identity_ok = False
        magnetic_bearing_linearization_identity_ok = False
        magnetic_bearing_dynamic_identity_ok = False
        moving_conductor_drag_identity_ok = False
        magnetic_gear_identity_ok = False
        demag_bem_charge_identity_ok = False
        magnetic_bearing_bias_sweep_identity_ok = False
        pm_demag_recoil_identity_ok = False
        maglev_dynamic_stiffness_identity_ok = False
        bem_demag_reciprocity_identity_ok = False
        eddy_maglev_power_identity_ok = False
        pm_coupling_energy_identity_ok = False
        thin_conductor_surface_impedance_identity_ok = False
        magnetic_gear_action_reaction_identity_ok = False
        multilayer_shield_closure_identity_ok = False
        transformer_energy_force_identity_ok = False
        maglev_equilibrium_energy_identity_ok = False
        eddy_shield_frequency_identity_ok = False
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
        bem_demag_surface_material_frame_generation_identity_ok = (
            _bem_demag_surface_material_frame_identity_ok(
                identity_value.get(
                    "bem_demag_surface_orientation_magnetization_volume_material_frame_generation_identity"
                )
            )
        )
        linear_motor_thrust_ripple_generation_identity_ok = (
            _linear_motor_thrust_ripple_identity_ok(
                identity_value.get(
                    "linear_motor_thrust_ripple_period_position_phase_frame_generation_identity"
                )
            )
        )
        levitation_gradient_energy_identity_ok = _levitation_gradient_energy_identity_ok(
            identity_value.get(
                "levitation_force_displacement_gradient_stiffness_energy_derivative_frame_generation_identity"
            )
        )
        cogging_periodic_interpolation_identity_ok = _cogging_periodic_interpolation_identity_ok(
            identity_value.get(
                "cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_identity"
            )
        )
        bem_panel_self_term_energy_force_identity_ok = (
            _bem_panel_self_term_energy_force_identity_ok(
                identity_value.get(
                    "bem_panel_orientation_magnetization_frame_self_term_energy_force_generation_identity"
                )
            )
        )
        motor_reduced_basis_torque_identity_ok = _motor_reduced_basis_torque_identity_ok(
            identity_value.get(
                "motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_generation_identity"
            )
        )
        maglev_force_stiffness_identity_ok = _maglev_force_stiffness_identity_ok(
            identity_value.get(
                "maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_generation_identity"
            )
        )
        motor_coenergy_torque_identity_ok = _motor_coenergy_torque_identity_ok(
            identity_value.get(
                "motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_generation_identity"
            )
        )
        airgap_stress_harmonic_torque_identity_ok = (
            _airgap_stress_harmonic_torque_identity_ok(
                identity_value.get(
                    "airgap_stress_harmonic_sector_periodicity_origin_sampling_alias_radius_torque_generation_identity"
                )
            )
        )
        laminated_core_loss_identity_ok = _laminated_core_loss_identity_ok(
            identity_value.get(
                "laminated_core_hysteresis_eddy_excess_frequency_flux_lamination_volume_result_generation_identity"
            )
        )
        magnet_demag_volume_fraction_identity_ok = _magnet_demag_volume_fraction_identity_ok(
            identity_value.get(
                "magnet_demag_recoil_knee_field_volume_generation_identity"
            )
        )
        linear_motor_wave_end_effect_identity_ok = _linear_motor_wave_end_effect_identity_ok(
            identity_value.get(
                "linear_motor_end_phase_wave_pitch_force_generation_identity"
            )
        )
        maglev_stiffness_identity_ok = _maglev_stiffness_identity_ok(
            identity_value.get(
                "maglev_force_stiffness_position_current_derivative_frame_mesh_result_identity"
            )
        )
        cogging_torque_sampling_identity_ok = _cogging_torque_sampling_identity_ok(
            identity_value.get(
                "cogging_torque_slot_pole_period_origin_sampling_harmonic_phase_mesh_result_identity"
            )
        )
        bem_near_singular_force_identity_ok = _bem_near_singular_force_identity_ok(
            identity_value.get(
                "bem_near_singular_gap_quadrature_normal_order_force_reciprocity_geometry_result_identity"
            )
        )
        hysteresis_minor_loop_identity_ok = _hysteresis_minor_loop_identity_ok(
            identity_value.get(
                "hysteresis_minor_loop_state_reversal_return_memory_remanence_energy_time_material_identity"
            )
        )
        maglev_equilibrium_closure_identity_ok = _maglev_equilibrium_closure_identity_ok(
            identity_value.get(
                "maglev_equilibrium_force_displacement_derivative_stiffness_gravity_mesh_result_identity"
            )
        )
        bem_surface_charge_closure_identity_ok = _bem_surface_charge_closure_identity_ok(
            identity_value.get(
                "bem_surface_charge_gauge_normal_energy_reciprocity_geometry_owner_result_identity"
            )
        )
        halbach_harmonic_closure_identity_ok = _halbach_harmonic_closure_identity_ok(
            identity_value.get(
                "halbach_harmonic_magnetization_order_pitch_phase_grid_field_energy_force_geometry_owner_result_identity"
            )
        )
        magnetic_bearing_linearization_identity_ok = (
            _magnetic_bearing_linearization_identity_ok(
                identity_value.get(
                    "magnetic_bearing_force_current_displacement_stiffness_reciprocity_bias_frame_mesh_result_identity"
                )
            )
        )
        magnetic_bearing_dynamic_identity_ok = _magnetic_bearing_dynamic_identity_ok(
            identity_value.get(
                "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
            )
        )
        moving_conductor_drag_identity_ok = _moving_conductor_drag_identity_ok(
            identity_value.get(
                "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
            )
        )
        magnetic_gear_identity_ok = _magnetic_gear_identity_ok(
            identity_value.get(
                "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
            )
        )
        demag_bem_charge_identity_ok = _demag_bem_charge_identity_ok(
            identity_value.get(
                "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
            )
        )
        magnetic_bearing_bias_sweep_identity_ok = (
            _magnetic_bearing_bias_sweep_identity_ok(
                identity_value.get(
                    "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
                )
            )
        )
        pm_demag_recoil_identity_ok = _pm_demag_recoil_identity_ok(
            identity_value.get(
                "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
            )
        )
        maglev_dynamic_stiffness_identity_ok = _maglev_dynamic_stiffness_identity_ok(
            identity_value.get(
                "maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"
            )
        )
        bem_demag_reciprocity_identity_ok = _bem_demag_reciprocity_identity_ok(
            identity_value.get(
                "bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"
            )
        )
        eddy_maglev_power_identity_ok = _eddy_maglev_power_identity_ok(
            identity_value.get(
                "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
            )
        )
        pm_coupling_energy_identity_ok = _pm_coupling_energy_identity_ok(
            identity_value.get(
                "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
            )
        )
        thin_conductor_surface_impedance_identity_ok = _thin_conductor_surface_impedance_identity_ok(
            identity_value.get(
                "thin_conductor_surface_impedance_skin_sheetcurrent_fieldjump_complexpower_surface_owner_result_identity"
            )
        )
        magnetic_gear_action_reaction_identity_ok = _magnetic_gear_action_reaction_identity_ok(
            identity_value.get(
                "magnetic_gear_harmonic_polepair_modulation_phase_ratio_torque_actionreaction_power_owner_result_identity"
            )
        )
        multilayer_shield_closure_identity_ok = _multilayer_shield_closure_identity_ok(
            identity_value.get(
                "multilayer_magnetic_shield_permeability_thickness_radius_interface_flux_attenuation_leakage_energy_geometry_result_identity"
            )
        )
        transformer_energy_force_identity_ok = _transformer_energy_force_identity_ok(
            identity_value.get(
                "transformer_leakage_mutual_inductance_fluxlinkage_reciprocity_psd_coenergy_force_winding_result_identity"
            )
        )
        maglev_equilibrium_energy_identity_ok = _maglev_equilibrium_energy_identity_ok(
            identity_value.get(
                "maglev_equilibrium_airgap_position_force_gradient_stiffness_potential_energy_stability_geometry_result_identity"
            )
        )
        eddy_shield_frequency_identity_ok = _eddy_shield_frequency_identity_ok(
            identity_value.get(
                "eddy_shield_frequency_conductivity_permeability_skin_depth_thickness_phase_attenuation_loss_energy_geometry_result_identity"
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
        "bem_demag_surface_shares_orientation_magnetization_volume_material_and_frame": (
            bem_demag_surface_material_frame_generation_identity_ok
        ),
        "linear_motor_thrust_ripple_shares_period_position_phase_frame_and_order": (
            linear_motor_thrust_ripple_generation_identity_ok
        ),
        "levitation_force_stiffness_and_energy_derivative_share_displacement_frame_sign_and_generation": (
            levitation_gradient_energy_identity_ok
        ),
        "cogging_torque_uses_current_position_periodicity_mesh_interpolation_reference_and_result": (
            cogging_periodic_interpolation_identity_ok
        ),
        "bem_demag_force_uses_current_panel_orientation_magnetization_self_term_energy_mesh_and_result": (
            bem_panel_self_term_energy_force_identity_ok
        ),
        "motor_reduced_torque_uses_current_basis_snapshots_operating_point_weights_residual_and_result": (
            motor_reduced_basis_torque_identity_ok
        ),
        "maglev_stiffness_uses_current_displacement_coordinate_force_geometry_mesh_and_solution": (
            maglev_force_stiffness_identity_ok
        ),
        "motor_torque_uses_current_winding_harmonics_currents_phases_angles_coenergy_mesh_and_result": (
            motor_coenergy_torque_identity_ok
        ),
        "airgap_torque_uses_current_sector_sampling_alias_harmonics_geometry_mesh_and_result": (
            airgap_stress_harmonic_torque_identity_ok
        ),
        "laminated_core_loss_uses_current_frequency_flux_lamination_volume_components_and_result": (
            laminated_core_loss_identity_ok
        ),
        "magnet_demag_uses_current_recoil_knee_temperature_local_field_mask_volume_and_result": (
            magnet_demag_volume_fraction_identity_ok
        ),
        "linear_motor_force_uses_current_end_effect_phase_wave_pitch_positions_ripple_and_result": (
            linear_motor_wave_end_effect_identity_ok
        ),
        "maglev_stiffness_uses_fixed_current_symmetric_positions_force_derivative_frame_mesh_and_result": (
            maglev_stiffness_identity_ok
        ),
        "cogging_torque_uses_slot_pole_period_origin_sampling_harmonics_phase_mesh_and_result": (
            cogging_torque_sampling_identity_ok
        ),
        "bem_near_contact_force_uses_gap_adaptive_quadrature_normals_order_reciprocity_geometry_and_result": (
            bem_near_singular_force_identity_ok
        ),
        "hysteresis_minor_loop_uses_initial_state_reversals_return_memory_remanence_energy_time_and_material": (
            hysteresis_minor_loop_identity_ok
        ),
        "maglev_equilibrium_uses_upward_force_global_displacement_central_stiffness_gravity_mesh_and_result": (
            maglev_equilibrium_closure_identity_ok
        ),
        "bem_surface_charge_uses_neutral_charge_mean_zero_gauge_opposed_normals_energy_reciprocity_geometry_and_result": (
            bem_surface_charge_closure_identity_ok
        ),
        "halbach_harmonics_use_current_magnetization_pitch_phase_grid_field_energy_force_geometry_owner_and_result": (
            halbach_harmonic_closure_identity_ok
        ),
        "magnetic_bearing_uses_current_force_jacobians_bias_frame_reciprocal_positive_stiffness_mesh_and_result": (
            magnetic_bearing_linearization_identity_ok
        ),
        "magnetic_bearing_dynamics_use_current_force_perturbations_stiffness_damping_coordinates_stability_operating_point_owner_and_result": (
            magnetic_bearing_dynamic_identity_ok
        ),
        "moving_conductor_uses_current_velocity_frame_drag_lift_power_skin_depth_slip_mesh_owner_field_and_result": (
            moving_conductor_drag_identity_ok
        ),
        "magnetic_gear_uses_current_poles_harmonic_torque_phase_power_frame_mesh_owner_and_result": (
            magnetic_gear_identity_ok
        ),
        "demag_bem_uses_neutral_surface_charge_outward_normals_jump_farfield_energy_mesh_owner_and_solution": (
            demag_bem_charge_identity_ok
        ),
        "magnetic_bearing_bias_sweep_uses_current_bias_symmetric_force_derivative_crosscoupling_frame_owner_and_result": (
            magnetic_bearing_bias_sweep_identity_ok
        ),
        "pm_demag_uses_temperature_adjusted_recoil_knee_loadline_irreversible_loss_orientation_mesh_owner_and_result": (
            pm_demag_recoil_identity_ok
        ),
        "maglev_dynamics_close_bias_equilibrium_frequency_stiffness_damping_phase_frame_owner_and_result": (
            maglev_dynamic_stiffness_identity_ok
        ),
        "bem_demag_closes_reciprocal_energy_field_magnetization_surface_volume_mesh_solution_and_result": (
            bem_demag_reciprocity_identity_ok
        ),
        "eddy_maglev_closes_velocity_frequency_skin_depth_lift_drag_joule_power_mesh_owner_and_result": (
            eddy_maglev_power_identity_ok
        ),
        "pm_coupling_closes_pole_periodic_energy_derivative_action_reaction_frame_mesh_owner_and_result": (
            pm_coupling_energy_identity_ok
        ),
        "thin_conductors_close_surface_impedance_skin_current_field_jump_complex_power_owner_and_result": (
            thin_conductor_surface_impedance_identity_ok
        ),
        "magnetic_gears_close_harmonics_poles_ratio_torque_reaction_power_owner_and_result": (
            magnetic_gear_action_reaction_identity_ok
        ),
        "multilayer_shields_close_material_thickness_flux_attenuation_leakage_energy_geometry_and_result": (
            multilayer_shield_closure_identity_ok
        ),
        "transformers_close_inductance_leakage_flux_reciprocity_psd_coenergy_force_winding_and_result": (
            transformer_energy_force_identity_ok
        ),
        "maglev_equilibrium_closes_gap_force_gradient_stiffness_energy_stability_owner_and_result": (
            maglev_equilibrium_energy_identity_ok
        ),
        "eddy_shields_close_frequency_skin_depth_phase_attenuation_loss_energy_owner_and_result": (
            eddy_shield_frequency_identity_ok
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
