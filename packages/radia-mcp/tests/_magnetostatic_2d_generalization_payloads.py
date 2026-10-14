"""Payload builders for ``test_magnetostatic_2d_generalization.py``.

Not collected by pytest (leading underscore).  The ``_identity_v19`` ..
``_identity_v28`` chain feeds ``force_coenergy_displacement_gate``; every
later builder (v45..v56) feeds its own ``validate_public_identity`` gate.
"""

from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v49 import (
    HARMONIC as HARMONIC_V49,
    NONLINEAR as NONLINEAR_V49,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v50 import CIRCUIT as CIRCUIT_V50, IRON as IRON_V50
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v51 import (
    INCREMENTAL as INCREMENTAL_V51,
    WEIGHTED_FORCE as WEIGHTED_FORCE_V51,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v52 import (
    FROZEN_INDUCTANCE as FROZEN_INDUCTANCE_V52,
    MAGNETIC_PRESSURE as MAGNETIC_PRESSURE_V52,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v53 import (
    HYSTERESIS as HYSTERESIS_V53,
    VIRTUAL_WORK as VIRTUAL_WORK_V53,
)
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v54 import FORCE as FORCE_V54, POWER as POWER_V54
from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v55 import (
    CAPACITANCE as CAPACITANCE_V55,
    INDUCTION as INDUCTION_V55,
)
from radia_mcp.radia_ngsolve.electromagnetic_force_heat_identity_v56 import FORCE as FORCE_V56, HEAT as HEAT_V56
from radia_mcp.radia_ngsolve.electromagnetic_semantic_identity_v48 import (
    ELECTROSTATIC as ELECTROSTATIC_V48,
    INCREMENTAL as INCREMENTAL_V48,
)
from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from radia_mcp.radia_ngsolve.magnetic_artifact_lineage_v47 import BH as BH_V47, FORCE as FORCE_V47
from test_force_coenergy_gate import _artifact_identity, _quadratic_case


# --- force/coenergy displacement gate chain (v19 .. v28) -------------------


def _identity_v19(sample_count):
    identity = _artifact_identity(sample_count)
    identity["weighted_stress_air_mask_nodal_weight_mesh_generation_identity"] = {
        "field_mesh_generation": "mesh-21",
        "air_mask_mesh_generation": "mesh-21",
        "nodal_weight_mesh_generation": "mesh-21",
        "force_integral_mesh_generation": "mesh-21",
        "air_region_ids": [1, 2],
        "mask_air_region_ids": [1, 2],
        "nodal_weight_node_ids": [101, 102, 103],
        "force_weight_node_ids": [101, 102, 103],
        "air_mask_sha256": "a" * 64,
        "force_air_mask_sha256": "a" * 64,
        "nodal_weight_sha256": "b" * 64,
        "force_nodal_weight_sha256": "b" * 64,
    }
    identity["sliding_band_periodic_angle_rotor_position_generation_identity"] = {
        "rotor_position_generation": "rotor-position-21",
        "sliding_band_rotor_position_generation": "rotor-position-21",
        "torque_rotor_position_generation": "rotor-position-21",
        "periodic_angle_generation": "periodic-angle-21",
        "sliding_band_periodic_angle_generation": "periodic-angle-21",
        "torque_periodic_angle_generation": "periodic-angle-21",
        "rotor_angle_deg": 15.0,
        "sliding_band_rotor_angle_deg": 15.0,
        "periodic_angle_pairs_deg": [[0.0, 30.0], [30.0, 60.0]],
        "torque_periodic_angle_pairs_deg": [[0.0, 30.0], [30.0, 60.0]],
        "sliding_band_map_sha256": "c" * 64,
        "torque_sliding_band_map_sha256": "c" * 64,
    }
    return identity


def _gate(identity):
    positions, coenergy, forces = _quadratic_case()
    return force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )


def _identity_v20(sample_count):
    identity = _identity_v19(sample_count)
    solve_generations = ["solve-22-a", "solve-22-b", "solve-22-c"]
    identity["coenergy_torque_angle_difference_remesh_state_generation_identity"] = {
        "torque_generation": "torque-22",
        "derivative_torque_generation": "torque-22",
        "state_table_generation": "state-table-22",
        "derivative_state_table_generation": "state-table-22",
        "coenergy_solve_generations": solve_generations,
        "mesh_remap_solve_generations": solve_generations,
        "excitation_solve_generations": solve_generations,
        "angle_state_solve_generations": solve_generations,
        "angles_deg": [14.0, 15.0, 16.0],
        "derivative_angles_deg": [14.0, 15.0, 16.0],
        "angle_spacing_generation": "angle-spacing-22",
        "derivative_angle_spacing_generation": "angle-spacing-22",
        "coenergy_state_table_sha256": "1" * 64,
        "derivative_state_table_sha256": "1" * 64,
    }
    identity[
        "axisymmetric_henrotte_hodge_radius_weight_coordinate_generation_identity"
    ] = {
        "mesh_geometry_generation": "axi-mesh-22",
        "field_mesh_geometry_generation": "axi-mesh-22",
        "radius_weight_mesh_geometry_generation": "axi-mesh-22",
        "cylindrical_coordinate_mesh_geometry_generation": "axi-mesh-22",
        "node_ids": [101, 102, 103],
        "field_node_ids": [101, 102, 103],
        "radius_m": [0.01, 0.02, 0.03],
        "hodge_radius_weight_m": [0.01, 0.02, 0.03],
        "cylindrical_r_coordinate_m": [0.01, 0.02, 0.03],
        "radius_weight_table_sha256": "2" * 64,
        "force_radius_weight_table_sha256": "2" * 64,
        "coordinate_table_sha256": "3" * 64,
        "field_coordinate_table_sha256": "3" * 64,
    }
    return identity


def _identity_v21(sample_count):
    identity = _identity_v20(sample_count)
    identity["weighted_stress_force_mask_material_mesh_generation_identity"] = {
        "solve_generation": "magnetostatic-solve-31",
        "weighted_stress_mask_solve_generation": "magnetostatic-solve-31",
        "material_label_solve_generation": "magnetostatic-solve-31",
        "mesh_field_solve_generation": "magnetostatic-solve-31",
        "force_integral_solve_generation": "magnetostatic-solve-31",
        "body_group_ids": [4, 5],
        "mask_body_group_ids": [4, 5],
        "material_labels": ["steel", "magnet"],
        "resolved_material_labels": ["steel", "magnet"],
        "weighted_force_n": [12.5, -3.2],
        "reported_force_n": [12.5, -3.2],
        "mask_table_sha256": "1" * 64,
        "force_mask_table_sha256": "1" * 64,
        "mesh_field_sha256": "2" * 64,
        "force_mesh_field_sha256": "2" * 64,
    }
    identity["harmonic_loss_phase_frequency_lamination_generation_identity"] = {
        "analysis_generation": "harmonic-loss-31",
        "phase_convention_analysis_generation": "harmonic-loss-31",
        "frequency_analysis_generation": "harmonic-loss-31",
        "lamination_analysis_generation": "harmonic-loss-31",
        "material_loss_analysis_generation": "harmonic-loss-31",
        "loss_result_analysis_generation": "harmonic-loss-31",
        "phase_convention": "exp(+jwt)",
        "loss_phase_convention": "exp(+jwt)",
        "frequency_hz": 400.0,
        "loss_frequency_hz": 400.0,
        "lamination_orientations": ["in_plane", "stacking_z"],
        "loss_lamination_orientations": ["in_plane", "stacking_z"],
        "material_loss_coefficients": [[0.02, 1.6], [0.01, 2.0]],
        "loss_material_coefficients": [[0.02, 1.6], [0.01, 2.0]],
        "material_loss_table_sha256": "3" * 64,
        "loss_material_table_sha256": "3" * 64,
    }
    return identity


def _identity_v22(sample_count):
    identity = _identity_v21(sample_count)
    identity["axisymmetric_force_energy_measure_depth_coordinate_generation_identity"] = {
        "solve_generation": "axisym-solve-41",
        "measure_solve_generation": "axisym-solve-41",
        "depth_solve_generation": "axisym-solve-41",
        "coordinate_solve_generation": "axisym-solve-41",
        "result_solve_generation": "axisym-solve-41",
        "problem_type": "axisymmetric",
        "result_problem_type": "axisymmetric",
        "measure_convention": "2*pi*r",
        "result_measure_convention": "2*pi*r",
        "planar_depth_m": 1.0,
        "result_planar_depth_m": 1.0,
        "coordinate_convention": "r_z",
        "result_coordinate_convention": "r_z",
        "force_energy_values": [12.5, 0.031],
        "reported_force_energy_values": [12.5, 0.031],
        "normalization_table_sha256": "1" * 64,
        "result_normalization_table_sha256": "1" * 64,
    }
    identity["nonlinear_incremental_mu_force_branch_perturbation_generation_identity"] = {
        "operating_point_generation": "nonlinear-op-41",
        "branch_operating_point_generation": "nonlinear-op-41",
        "differential_mu_operating_point_generation": "nonlinear-op-41",
        "perturbation_operating_point_generation": "nonlinear-op-41",
        "force_operating_point_generation": "nonlinear-op-41",
        "branch_id": "up-sweep:17",
        "force_branch_id": "up-sweep:17",
        "perturbation_current_a": 0.01,
        "force_perturbation_current_a": 0.01,
        "differential_mu_sha256": "2" * 64,
        "force_differential_mu_sha256": "2" * 64,
        "incremental_force_n": [0.12, -0.03],
        "reported_incremental_force_n": [0.12, -0.03],
        "incremental_state_sha256": "3" * 64,
        "force_incremental_state_sha256": "3" * 64,
    }
    return identity


def _identity_v23(sample_count):
    identity = _identity_v22(sample_count)
    identity["nonlinear_bh_branch_operating_point_force_mesh_generation_identity"] = {
        "solve_generation": "nonlinear-force-51",
        "branch_solve_generation": "nonlinear-force-51",
        "operating_point_solve_generation": "nonlinear-force-51",
        "permeability_solve_generation": "nonlinear-force-51",
        "force_mesh_solve_generation": "nonlinear-force-51",
        "force_result_solve_generation": "nonlinear-force-51",
        "branch_id": "ascending:23",
        "force_branch_id": "ascending:23",
        "operating_point_current_a": 7.5,
        "force_operating_point_current_a": 7.5,
        "operating_point_flux_density_t": [1.35, 1.42],
        "force_operating_point_flux_density_t": [1.35, 1.42],
        "permeability_state_sha256": "1" * 64,
        "force_permeability_state_sha256": "1" * 64,
        "force_mesh_sha256": "2" * 64,
        "integrated_force_mesh_sha256": "2" * 64,
        "force_n": [18.2, -0.4],
        "reported_force_n": [18.2, -0.4],
    }
    identity["sliding_band_angle_mesh_harmonic_torque_generation_identity"] = {
        "sweep_generation": "sliding-sweep-51",
        "angle_sweep_generation": "sliding-sweep-51",
        "airgap_mesh_sweep_generation": "sliding-sweep-51",
        "phase_current_sweep_generation": "sliding-sweep-51",
        "torque_sample_sweep_generation": "sliding-sweep-51",
        "harmonic_sweep_generation": "sliding-sweep-51",
        "rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "torque_rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "airgap_mesh_sha256": "3" * 64,
        "torque_airgap_mesh_sha256": "3" * 64,
        "phase_current_table_sha256": "4" * 64,
        "torque_phase_current_table_sha256": "4" * 64,
        "torque_samples_nm": [1.0, 1.2, 0.9, 1.1],
        "harmonic_torque_samples_nm": [1.0, 1.2, 0.9, 1.1],
        "harmonic_orders": [0, 1, 2],
        "reported_harmonic_orders": [0, 1, 2],
        "harmonic_amplitudes_nm": [1.05, 0.12, 0.04],
        "reported_harmonic_amplitudes_nm": [1.05, 0.12, 0.04],
        "torque_sample_table_sha256": "5" * 64,
        "harmonic_sample_table_sha256": "5" * 64,
    }
    return identity


def _identity_v24(sample_count):
    identity = _identity_v23(sample_count)
    identity[
        "weighted_stress_energy_derivative_force_mesh_frame_unit_generation_identity"
    ] = {
        "force_generation": "force-101",
        "weighted_stress_force_generation": "force-101",
        "energy_derivative_force_generation": "force-101",
        "mesh_force_generation": "force-101",
        "displacement_frame_force_generation": "force-101",
        "unit_force_generation": "force-101",
        "result_force_generation": "force-101",
        "weighted_stress_mesh_sha256": "1" * 64,
        "energy_derivative_mesh_sha256": "1" * 64,
        "displacement_frame_id": "world:x",
        "energy_derivative_displacement_frame_id": "world:x",
        "displacement_unit": "m",
        "energy_derivative_displacement_unit": "m",
        "force_unit": "N",
        "energy_derivative_force_unit": "N",
        "weighted_stress_force_n": [12.4, -0.2],
        "energy_derivative_force_n": [12.4, -0.2],
        "force_result_sha256": "2" * 64,
        "energy_derivative_result_sha256": "2" * 64,
    }
    identity[
        "axisymmetric_revolved_energy_force_2pir_jacobian_derham_generation_identity"
    ] = {
        "axisymmetric_generation": "axisym-101",
        "jacobian_axisymmetric_generation": "axisym-101",
        "field_axisymmetric_generation": "axisym-101",
        "material_axisymmetric_generation": "axisym-101",
        "mesh_axisymmetric_generation": "axisym-101",
        "revolved_result_axisymmetric_generation": "axisym-101",
        "jacobian_measure": "2*pi*r",
        "revolved_jacobian_measure": "2*pi*r",
        "field_state_sha256": "3" * 64,
        "revolved_field_state_sha256": "3" * 64,
        "material_map_sha256": "4" * 64,
        "revolved_material_map_sha256": "4" * 64,
        "axisymmetric_mesh_sha256": "5" * 64,
        "revolved_source_mesh_sha256": "5" * 64,
        "revolution_angle_deg": 360.0,
        "axisymmetric_energy_j": 0.125,
        "revolved_energy_j": 0.125,
        "axisymmetric_force_n": [4.2, 0.0],
        "revolved_force_n": [4.2, 0.0],
        "derham_sequence_id": "H1-axisym:Hodge-2pir",
        "revolved_derham_sequence_id": "H1-axisym:Hodge-2pir",
        "result_sha256": "6" * 64,
        "revolved_result_sha256": "6" * 64,
    }
    return identity


def _identity_v25(sample_count):
    identity = _identity_v24(sample_count)
    identity[
        "nonlinear_bh_incremental_energy_coenergy_force_branch_mesh_generation_identity"
    ] = {
        "nonlinear_generation": "nonlinear-111",
        "bh_curve_nonlinear_generation": "nonlinear-111",
        "material_nonlinear_generation": "nonlinear-111",
        "branch_nonlinear_generation": "nonlinear-111",
        "incremental_state_nonlinear_generation": "nonlinear-111",
        "mesh_nonlinear_generation": "nonlinear-111",
        "energy_nonlinear_generation": "nonlinear-111",
        "coenergy_nonlinear_generation": "nonlinear-111",
        "force_nonlinear_generation": "nonlinear-111",
        "result_nonlinear_generation": "nonlinear-111",
        "nonlinear_material_ids": ["core"],
        "result_nonlinear_material_ids": ["core"],
        "bh_curve_sha256": "1" * 64,
        "result_bh_curve_sha256": "1" * 64,
        "material_map_sha256": "2" * 64,
        "result_material_map_sha256": "2" * 64,
        "branch_id": "ascending:step-12",
        "result_branch_id": "ascending:step-12",
        "load_current_a": 8.0,
        "result_load_current_a": 8.0,
        "incremental_state_sha256": "3" * 64,
        "result_incremental_state_sha256": "3" * 64,
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "magnetic_energy_j": 0.42,
        "result_magnetic_energy_j": 0.42,
        "magnetic_coenergy_j": 0.58,
        "result_magnetic_coenergy_j": 0.58,
        "incremental_force_n": [12.2, -0.1],
        "result_incremental_force_n": [12.2, -0.1],
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    identity[
        "open_boundary_domain_decay_multipole_moment_material_generation_identity"
    ] = {
        "boundary_generation": "open-boundary-111",
        "domain_boundary_generation": "open-boundary-111",
        "mesh_boundary_generation": "open-boundary-111",
        "material_boundary_generation": "open-boundary-111",
        "multipole_boundary_generation": "open-boundary-111",
        "decay_boundary_generation": "open-boundary-111",
        "result_boundary_generation": "open-boundary-111",
        "boundary_type": "asymptotic_multipole",
        "result_boundary_type": "asymptotic_multipole",
        "source_radius_m": 0.05,
        "result_source_radius_m": 0.05,
        "outer_radius_m": 0.4,
        "result_outer_radius_m": 0.4,
        "multipole_order": 3,
        "result_multipole_order": 3,
        "material_map_sha256": "6" * 64,
        "result_material_map_sha256": "6" * 64,
        "mesh_sha256": "7" * 64,
        "result_mesh_sha256": "7" * 64,
        "multipole_moment_sha256": "8" * 64,
        "result_multipole_moment_sha256": "8" * 64,
        "decay_sample_radii_m": [0.2, 0.3, 0.4],
        "result_decay_sample_radii_m": [0.2, 0.3, 0.4],
        "decay_flux_density_t": [1.0e-3, 3.0e-4, 1.2e-4],
        "result_decay_flux_density_t": [1.0e-3, 3.0e-4, 1.2e-4],
        "result_sha256": "9" * 64,
        "accepted_result_sha256": "9" * 64,
    }
    return identity


def _identity_v26(sample_count):
    identity = _identity_v25(sample_count)
    identity["weighted_stress_tensor_mask_mesh_region_force_torque_energy_generation_identity"] = {
        "force_generation": "force-131", "mask_force_generation": "force-131",
        "mesh_force_generation": "force-131", "region_force_generation": "force-131",
        "energy_force_generation": "force-131", "torque_force_generation": "force-131",
        "result_force_generation": "force-131", "body_group_ids": [3, 4],
        "result_body_group_ids": [3, 4], "weighted_mask_node_ids": [101, 102, 103, 104],
        "result_weighted_mask_node_ids": [101, 102, 103, 104],
        "integration_region": "surrounding-air", "result_integration_region": "surrounding-air",
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "mask_sha256": "2" * 64, "result_mask_sha256": "2" * 64,
        "integration_region_sha256": "3" * 64, "result_integration_region_sha256": "3" * 64,
        "weighted_force_n": [12.5, -0.2, 0.0], "result_weighted_force_n": [12.5, -0.2, 0.0],
        "weighted_torque_nm": [0.0, 0.0, 1.8], "result_weighted_torque_nm": [0.0, 0.0, 1.8],
        "magnetic_energy_j": 0.44, "result_magnetic_energy_j": 0.44,
        "magnetic_coenergy_j": 0.61, "result_magnetic_coenergy_j": 0.61,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    total_force = 1.2
    identity["axisymmetric_planar_depth_two_pi_r_force_normalization_coordinate_unit_generation_identity"] = {
        "normalization_generation": "normalization-131",
        "planar_depth_normalization_generation": "normalization-131",
        "radius_normalization_generation": "normalization-131",
        "coordinate_normalization_generation": "normalization-131",
        "unit_normalization_generation": "normalization-131",
        "mesh_normalization_generation": "normalization-131",
        "result_normalization_generation": "normalization-131",
        "planar_depth_m": 0.05, "result_planar_depth_m": 0.05,
        "planar_force_n_per_m": 24.0, "result_planar_force_n_per_m": 24.0,
        "planar_total_force_n": total_force, "result_planar_total_force_n": total_force,
        "axisymmetric_radius_m": 0.03, "result_axisymmetric_radius_m": 0.03,
        "axisymmetric_meridian_force_n_per_rad": total_force / (2.0 * math.pi),
        "result_axisymmetric_meridian_force_n_per_rad": total_force / (2.0 * math.pi),
        "axisymmetric_total_force_n": total_force, "result_axisymmetric_total_force_n": total_force,
        "radius_measure_convention": "2*pi*r", "result_radius_measure_convention": "2*pi*r",
        "coordinate_convention": "r_z_right_handed", "result_coordinate_convention": "r_z_right_handed",
        "force_unit": "N_total_3d", "result_force_unit": "N_total_3d",
        "mesh_sha256": "5" * 64, "result_mesh_sha256": "5" * 64,
    }
    return identity


def _identity_v27(sample_count):
    identity = _identity_v26(sample_count)
    identity["nonlinear_bh_minor_loop_branch_interpolation_state_coenergy_force_generation_identity"] = {
        "nonlinear_generation": "bh-minor-loop-141",
        "branch_nonlinear_generation": "bh-minor-loop-141",
        "interpolation_nonlinear_generation": "bh-minor-loop-141",
        "state_nonlinear_generation": "bh-minor-loop-141",
        "coenergy_nonlinear_generation": "bh-minor-loop-141",
        "mesh_nonlinear_generation": "bh-minor-loop-141",
        "force_nonlinear_generation": "bh-minor-loop-141",
        "result_nonlinear_generation": "bh-minor-loop-141",
        "bh_branch": "ascending-minor-loop",
        "result_bh_branch": "ascending-minor-loop",
        "interpolation_rule": "monotone-cubic-h",
        "result_interpolation_rule": "monotone-cubic-h",
        "state_point_am": [120.0, 800.0],
        "result_state_point_am": [120.0, 800.0],
        "magnetic_coenergy_j": 0.37,
        "result_magnetic_coenergy_j": 0.37,
        "force_n": [18.0, -0.4],
        "result_force_n": [18.0, -0.4],
        "bh_table_sha256": "1" * 64,
        "result_bh_table_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "solution_sha256": "3" * 64,
        "accepted_solution_sha256": "3" * 64,
    }
    mu0 = 4.0e-7 * math.pi
    sigma = 5.8e7
    frequency = 2000.0
    skin_depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * mu0 * sigma))
    identity["harmonic_eddy_phasor_conductivity_skin_depth_frequency_loss_mesh_generation_identity"] = {
        "eddy_generation": "harmonic-eddy-141",
        "phasor_eddy_generation": "harmonic-eddy-141",
        "conductivity_eddy_generation": "harmonic-eddy-141",
        "skin_eddy_generation": "harmonic-eddy-141",
        "frequency_eddy_generation": "harmonic-eddy-141",
        "loss_eddy_generation": "harmonic-eddy-141",
        "mesh_eddy_generation": "harmonic-eddy-141",
        "result_eddy_generation": "harmonic-eddy-141",
        "phasor_convention": "exp(+jwt)",
        "result_phasor_convention": "exp(+jwt)",
        "conductivity_s_m": sigma,
        "result_conductivity_s_m": sigma,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "minimum_elements_per_skin_depth": 4.0,
        "result_minimum_elements_per_skin_depth": 4.0,
        "joule_loss_w": 42.0,
        "result_joule_loss_w": 42.0,
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "loss_result_sha256": "5" * 64,
        "accepted_loss_result_sha256": "5" * 64,
    }
    return identity


def _identity_v28(sample_count):
    identity = _identity_v27(sample_count)
    identity["axisymmetric_aphi_radial_weight_region_energy_force_mesh_solution_generation_identity"] = {
        "axisymmetric_generation": "axisym-force-151", "aphi_axisymmetric_generation": "axisym-force-151",
        "weight_axisymmetric_generation": "axisym-force-151", "region_axisymmetric_generation": "axisym-force-151",
        "energy_axisymmetric_generation": "axisym-force-151", "force_axisymmetric_generation": "axisym-force-151",
        "mesh_axisymmetric_generation": "axisym-force-151", "solution_axisymmetric_generation": "axisym-force-151",
        "result_axisymmetric_generation": "axisym-force-151", "formulation": "Aphi", "result_formulation": "Aphi",
        "radial_weighting": "2*pi*r", "result_radial_weighting": "2*pi*r",
        "region_ids": [4, 5], "result_region_ids": [4, 5],
        "displacement_m": [-0.0001, 0.0, 0.0001],
        "magnetic_coenergy_j": [0.4982, 0.5, 0.5018],
        "result_magnetic_coenergy_j": [0.4982, 0.5, 0.5018],
        "force_axis": "z", "result_force_axis": "z",
        "force_from_energy_n": 18.0, "result_force_from_energy_n": 18.0,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "solution_sha256": "2" * 64, "accepted_solution_sha256": "2" * 64,
    }
    identity["permanent_magnet_recoil_temperature_operating_point_demag_force_generation_identity"] = {
        "magnet_generation": "pm-operating-point-151", "recoil_magnet_generation": "pm-operating-point-151",
        "temperature_magnet_generation": "pm-operating-point-151", "operating_point_magnet_generation": "pm-operating-point-151",
        "frame_magnet_generation": "pm-operating-point-151", "demag_magnet_generation": "pm-operating-point-151",
        "force_magnet_generation": "pm-operating-point-151", "mesh_magnet_generation": "pm-operating-point-151",
        "result_magnet_generation": "pm-operating-point-151",
        "recoil_relative_permeability": 1.05, "result_recoil_relative_permeability": 1.05,
        "remanence_t": 1.2, "result_remanence_t": 1.2,
        "magnet_temperature_c": 80.0, "result_magnet_temperature_c": 80.0,
        "operating_point_bh": [0.82, -302000.0], "result_operating_point_bh": [0.82, -302000.0],
        "demag_margin_a_m": 95000.0, "result_demag_margin_a_m": 95000.0,
        "magnetization_frame_sha256": "3" * 64, "result_magnetization_frame_sha256": "3" * 64,
        "force_n": [12.0, 0.2], "result_force_n": [12.0, 0.2],
        "mesh_sha256": "4" * 64, "result_mesh_sha256": "4" * 64,
        "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
    }
    return identity


# --- axisymmetric_v44_identity (v45) ----------------------------------------


def _identity_v45():
    generation = "test-845"
    return {
        "v45_public_axisymmetric_force_torque_coenergy_stress_contour_owner_mismatch": {
            "generation": generation,
            **{key: generation for key in ("force_generation", "torque_generation", "coenergy_generation", "stress_contour_generation", "axis_factor_generation", "mesh_generation", "result_generation")},
            "force_method": "weighted_stress_tensor", "result_force_method": "weighted_stress_tensor", "torque_method": "airgap_contour", "result_torque_method": "airgap_contour",
            "force_n": [12.0, -3.0], "result_force_n": [12.0, -3.0], "torque_nm": 0.8, "result_torque_nm": 0.8, "coenergy_j": 2.1, "result_coenergy_j": 2.1,
            "axisymmetric_factor": 2.0 * math.pi, "result_axisymmetric_factor": 2.0 * math.pi, "contour_owner": "contour:test-845", "result_contour_owner": "contour:test-845", "mesh_owner": "mesh:test-845", "result_mesh_owner": "mesh:test-845",
            "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v45_public_electrostatic_fringe_charge_energy_capacitance_interface_flux_axisfactor_mismatch": {
            "generation": generation,
            **{key: generation for key in ("charge_generation", "energy_generation", "capacitance_generation", "interface_generation", "axis_factor_generation", "mesh_generation", "result_generation")},
            "voltage_v": 100.0, "result_voltage_v": 100.0, "capacitance_f": 1.0e-9, "result_capacitance_f": 1.0e-9, "charge_c": 1.0e-7, "result_charge_c": 1.0e-7,
            "stored_energy_j": 5.0e-6, "result_stored_energy_j": 5.0e-6, "interface_flux_c": 1.0e-7, "result_interface_flux_c": 1.0e-7, "axisymmetric_factor": 2.0 * math.pi, "result_axisymmetric_factor": 2.0 * math.pi,
            "mesh_owner": "mesh:test-845", "result_mesh_owner": "mesh:test-845", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


# --- axisymmetric_v46_identity (v46) ----------------------------------------

FORCE_V46 = "v46_public_axisymmetric_force_weighted_stress_contour_open_boundary_partial_mismatch"
BH_V46 = "v46_public_nonlinear_bh_curve_unit_scale_convergence_status_nan_mismatch"


def _identity_v46() -> dict[str, object]:
    generation = "test-force-open-v46"
    force = {
        "generation": generation,
        **{key: generation for key in ("force_generation", "contour_generation", "open_boundary_generation", "solve_generation", "result_generation")},
        "force_method": "weighted_stress_tensor", "result_force_method": "weighted_stress_tensor",
        "open_boundary": "outer_kelvin", "result_open_boundary": "outer_kelvin",
        "solve_state": "complete", "result_solve_state": "complete",
        "partial_solve_status": "none", "result_partial_solve_status": "none",
        "force_n": [12.0, -3.0], "result_force_n": [12.0, -3.0],
        "axisymmetric_factor": 2.0 * math.pi, "result_axisymmetric_factor": 2.0 * math.pi,
        "contour_owner": "contour:test-force-open-v46", "result_contour_owner": "contour:test-force-open-v46",
        "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
    }
    generation = "test-bh-unit-v46"
    bh = {
        "generation": generation,
        **{key: generation for key in ("bh_curve_generation", "unit_generation", "convergence_generation", "result_generation")},
        "bh_unit": "tesla_ampere_per_meter", "result_bh_unit": "tesla_ampere_per_meter",
        "convergence_status": "converged", "result_convergence_status": "converged",
        "nonfinite_status": "none", "result_nonfinite_status": "none",
        "bh_points": [[0.0, 0.0], [0.8, 500.0]], "result_bh_points": [[0.0, 0.0], [0.8, 500.0]],
        "material_owner": "material:test-bh-unit-v46", "result_material_owner": "material:test-bh-unit-v46",
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    return {FORCE_V46: force, BH_V46: bh}


# --- magnetic_artifact_lineage_v47 (v47) ------------------------------------


def _identity_v47() -> dict[str, object]:
    force_generation = "force-v47"
    bh_generation = "bh-v47"
    keys = ["current=1", "current=2", "current=3"]
    branches = ["ascending", "ascending", "descending"]
    return {
        FORCE_V47: {
            "generation": force_generation,
            "force_method_generation": force_generation,
            "body_owner_generation": force_generation,
            "displacement_pair_generation": force_generation,
            "result_generation": force_generation,
            "force_method": "weighted_stress_tensor",
            "result_force_method": "weighted_stress_tensor",
            "body_owner": "group:moving",
            "result_body_owner": "group:moving",
            "force_sign_convention": "positive_displacement_direction",
            "result_force_sign_convention": "positive_displacement_direction",
            "displacement_pair_m": [0.0, 0.001],
            "result_displacement_pair_m": [0.0, 0.001],
            "coenergy_pair_j": [1.0, 1.01],
            "result_coenergy_pair_j": [1.0, 1.01],
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        BH_V47: {
            "generation": bh_generation,
            "operating_point_generation": bh_generation,
            "branch_generation": bh_generation,
            "history_generation": bh_generation,
            "result_generation": bh_generation,
            "operating_point_row_keys": keys,
            "result_operating_point_row_keys": keys,
            "hysteresis_branches": branches,
            "result_hysteresis_branches": branches,
            "excitation_history_sha256": "2" * 64,
            "result_excitation_history_sha256": "2" * 64,
            "material_owner": "material:steel",
            "result_material_owner": "material:steel",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


# --- electromagnetic_semantic_identity_v48 (v48) ----------------------------


def _identity_v48() -> dict[str, object]:
    incremental_generation = "incremental-bias-v48-901"
    electrostatic_generation = "capacitance-sweep-v48-901"
    tangent = {"material:core": [[0.0020, 0.0001], [0.0001, 0.0018]]}
    voltage = [0.0, 1.0, 2.0, 3.0]
    capacitance = 2.0e-9
    charge = [capacitance * value for value in voltage]
    energy = [0.5 * capacitance * value * value for value in voltage]
    owners = ["conductor:electrode-a" for _ in voltage]
    return {
        INCREMENTAL_V48: {
            "generation": incremental_generation,
            "bias_generation": incremental_generation,
            "phasor_generation": incremental_generation,
            "tangent_generation": incremental_generation,
            "operating_point_generation": incremental_generation,
            "result_generation": incremental_generation,
            "frozen_bias_sha256": "6" * 64,
            "result_frozen_bias_sha256": "6" * 64,
            "harmonic_phasor_a": [1.5, -0.25],
            "result_harmonic_phasor_a": [1.5, -0.25],
            "material_tangent_h_per_m": tangent,
            "result_material_tangent_h_per_m": tangent,
            "operating_point_owner": "operating-point:incremental-v48-901",
            "result_operating_point_owner": "operating-point:incremental-v48-901",
            "result_sha256": "7" * 64,
            "accepted_result_sha256": "7" * 64,
        },
        ELECTROSTATIC_V48: {
            "generation": electrostatic_generation,
            "voltage_generation": electrostatic_generation,
            "charge_generation": electrostatic_generation,
            "energy_generation": electrostatic_generation,
            "conductor_generation": electrostatic_generation,
            "result_generation": electrostatic_generation,
            "voltage_v": voltage,
            "result_voltage_v": voltage,
            "charge_c": charge,
            "result_charge_c": charge,
            "field_energy_j": energy,
            "result_field_energy_j": energy,
            "capacitance_f": capacitance,
            "result_capacitance_f": capacitance,
            "conductor_owner_rows": owners,
            "result_conductor_owner_rows": owners,
            "result_sha256": "8" * 64,
            "accepted_result_sha256": "8" * 64,
        },
    }


# --- electromagnetic_artifact_identity_v49 (v49) ----------------------------


def _identity_v49() -> dict[str, object]:
    nonlinear_generation = "nonlinear-material-v49-901"
    harmonic_generation = "harmonic-model-v49-901"
    bh_rows = [[0.0, 0.0], [0.7, 120.0], [1.35, 620.0], [1.62, 2400.0]]
    incremental = [[0.0023, 0.0001], [0.0001, 0.0020]]
    lamination = {"fill_factor": 0.95, "direction": "in-plane", "sheet_thickness_m": 0.00035}
    circuit = {"name": "coil-a", "current_a": [3.0, -1.0], "turns": 120}
    losses = {"joule_w": 4.25, "core_w": 1.75}
    return {
        NONLINEAR_V49: {
            "generation": nonlinear_generation,
            "branch_generation": nonlinear_generation,
            "incremental_generation": nonlinear_generation,
            "temperature_generation": nonlinear_generation,
            "lamination_generation": nonlinear_generation,
            "result_generation": nonlinear_generation,
            "bh_branch": "ascending:first-quadrant",
            "result_bh_branch": "ascending:first-quadrant",
            "bh_rows_t_a_per_m": bh_rows,
            "result_bh_rows_t_a_per_m": bh_rows,
            "incremental_permeability_h_per_m": incremental,
            "result_incremental_permeability_h_per_m": incremental,
            "temperature_c": 80.0,
            "result_temperature_c": 80.0,
            "lamination": lamination,
            "result_lamination": lamination,
            "material_owner": "material:core-v49-901",
            "result_material_owner": "material:core-v49-901",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        HARMONIC_V49: {
            "generation": harmonic_generation,
            "geometry_generation": harmonic_generation,
            "frequency_generation": harmonic_generation,
            "circuit_generation": harmonic_generation,
            "loss_generation": harmonic_generation,
            "result_generation": harmonic_generation,
            "geometry_type": "planar",
            "result_geometry_type": "planar",
            "depth_m": 0.04,
            "result_depth_m": 0.04,
            "frequency_hz": 400.0,
            "result_frequency_hz": 400.0,
            "phase_convention": "exp(+jwt)",
            "result_phase_convention": "exp(+jwt)",
            "circuit": circuit,
            "result_circuit": circuit,
            "losses_w": losses,
            "result_losses_w": losses,
            "result_owner": "harmonic-result:v49-901",
            "accepted_result_owner": "harmonic-result:v49-901",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


# --- electromagnetic_artifact_identity_v50 (v50) ----------------------------


def _identity_v50() -> dict[str, object]:
    generation = "magnetostatic-v50"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    waveform = [0.0, 0.7, 1.2, 0.7, 0.0, -0.7, -1.2, -0.7, 0.0]
    components = {"hysteresis_w": 1.8, "eddy_w": 0.9, "excess_w": 0.3, "total_w": 3.0}
    return {
        CIRCUIT_V50: {
            "generation": generation, "current_generation": generation, "connection_generation": generation,
            "turn_generation": generation, "depth_generation": generation, "owner_generation": generation,
            "result_generation": generation, "current_phasor_a": [3.0, -1.0], "result_current_phasor_a": [3.0, -1.0],
            "winding_connection": "series", "result_winding_connection": "series", "parallel_paths": 1,
            "result_parallel_paths": 1, "turns": 120, "result_turns": 120, "depth_m": 0.04, "result_depth_m": 0.04,
            "circuit_owner": "circuit:coil-v50", "result_circuit_owner": "circuit:coil-v50", **result,
        },
        IRON_V50: {
            "generation": generation, "coefficient_generation": generation, "frequency_generation": generation,
            "waveform_generation": generation, "component_generation": generation, "owner_generation": generation,
            "result_generation": generation, "steinmetz_coefficients": {"kh": 0.021, "ke": 0.00018, "alpha": 1.62},
            "result_steinmetz_coefficients": {"kh": 0.021, "ke": 0.00018, "alpha": 1.62},
            "frequency_hz": 400.0, "result_frequency_hz": 400.0, "flux_density_waveform_t": waveform,
            "result_flux_density_waveform_t": waveform, "loss_components": components, "result_loss_components": components,
            "material_owner": "material:lamination-v50", "result_material_owner": "material:lamination-v50", **result,
        },
    }


# --- electromagnetic_artifact_identity_v51 (v51) ----------------------------


def _identity_v51() -> dict[str, object]:
    generation = "magnetostatic-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    tangent = {"radial": 410.0, "tangential": 370.0}
    force = [12.5, -0.4]
    return {
        INCREMENTAL_V51: {
            "generation": generation, "bias_generation": generation, "harmonic_generation": generation,
            "tangent_generation": generation, "branch_generation": generation, "owner_generation": generation,
            "result_generation": generation, "analysis_mode": "incremental_permeability",
            "result_analysis_mode": "incremental_permeability", "frozen_bias_solution_sha256": "1" * 64,
            "result_frozen_bias_solution_sha256": "1" * 64, "harmonic_frequency_hz": 1000.0,
            "result_harmonic_frequency_hz": 1000.0, "tangent_permeability_relative": tangent,
            "result_tangent_permeability_relative": tangent, "branch_state": "ascending_major_loop",
            "result_branch_state": "ascending_major_loop", "operating_point_owner": "operating-point:v51",
            "result_operating_point_owner": "operating-point:v51", **result,
        },
        WEIGHTED_FORCE_V51: {
            "generation": generation, "mask_generation": generation, "air_generation": generation,
            "axisym_generation": generation, "force_generation": generation, "frame_generation": generation,
            "owner_generation": generation, "result_generation": generation, "weighted_stress_mask_sha256": "2" * 64,
            "result_weighted_stress_mask_sha256": "2" * 64, "air_element_ids": [101, 102, 103],
            "result_air_element_ids": [101, 102, 103], "axisymmetric_radius_m": 0.025,
            "result_axisymmetric_radius_m": 0.025, "axisymmetric_factor_m": 2.0 * 3.141592653589793 * 0.025,
            "result_axisymmetric_factor_m": 2.0 * 3.141592653589793 * 0.025, "force_n": force,
            "result_force_n": force, "force_frame": "global_rz", "result_force_frame": "global_rz",
            "force_owner": "force:weighted-v51", "result_force_owner": "force:weighted-v51", **result,
        },
    }


# --- electromagnetic_artifact_identity_v52 (v52) ----------------------------


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity_v52() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    field_jump = {"normal_h_a_per_m": 1200.0, "tangential_b_t": 0.42}
    normal = [0.6, 0.8, 0.0]
    traction = [3200.0, 4200.0, 0.0]
    bias = {"current_a": 8.0, "solution_sha256": "b" * 64}
    perturbation = {"delta_current_a": 0.02, "frequency_hz": 400.0}
    return {
        MAGNETIC_PRESSURE_V52: {
            **_generations("magnetic-pressure-v52", ("field_generation", "normal_generation", "traction_generation", "owner_generation", "result_generation")),
            "field_jump": field_jump,
            "result_field_jump": field_jump,
            "boundary_normal": normal,
            "result_boundary_normal": normal,
            "traction_n_per_m2": traction,
            "result_traction_n_per_m2": traction,
            "integration_measure": "surface_area_m2",
            "result_integration_measure": "surface_area_m2",
            "field_owner": "field:magnetic-pressure-v52",
            "result_field_owner": "field:magnetic-pressure-v52",
            **result,
        },
        FROZEN_INDUCTANCE_V52: {
            **_generations("frozen-inductance-v52", ("permeability_generation", "bias_generation", "perturbation_generation", "inductance_generation", "owner_generation", "result_generation")),
            "permeability_mode": "frozen_at_bias",
            "result_permeability_mode": "frozen_at_bias",
            "bias_point": bias,
            "result_bias_point": bias,
            "perturbation": perturbation,
            "result_perturbation": perturbation,
            "incremental_inductance_h": 0.014,
            "result_incremental_inductance_h": 0.014,
            "solution_owner": "solution:frozen-inductance-v52",
            "result_solution_owner": "solution:frozen-inductance-v52",
            **result,
        },
    }


# --- electromagnetic_artifact_identity_v53 (v53) ----------------------------


def _identity_v53():
    frequency = 400.0; mu = {"real": 220.0, "imag": -18.0}; h_rms = 120.0
    loss = 2.0 * math.pi * frequency * 4.0e-7 * math.pi * (-mu["imag"]) * h_rms**2
    hysteresis = {**_generations("hyst-v53", ("permeability_generation", "phasor_generation", "loss_generation", "material_generation", "owner_generation", "result_generation")), "relative_permeability": mu, "result_relative_permeability": mu, "phasor_convention": "exp(+j_omega_t)", "result_phasor_convention": "exp(+j_omega_t)", "frequency_hz": frequency, "result_frequency_hz": frequency, "h_rms_a_per_m": h_rms, "result_h_rms_a_per_m": h_rms, "loss_density_w_m3": loss, "result_loss_density_w_m3": loss, "material_id": "material:steel", "result_material_id": "material:steel", "material_owner": "material-owner:v53", "result_material_owner": "material-owner:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    virtual_work = {**_generations("vw-v53", ("path_generation", "constraint_generation", "energy_generation", "force_generation", "owner_generation", "result_generation")), "virtual_work_path": "constant_voltage_coenergy", "result_virtual_work_path": "constant_voltage_coenergy", "constraint_mode": "fixed_voltage", "result_constraint_mode": "fixed_voltage", "voltage_v": 800.0, "result_voltage_v": 800.0, "charge_c": 2.0e-8, "result_charge_c": 2.0e-8, "virtual_displacement_m": 1.0e-6, "result_virtual_displacement_m": 1.0e-6, "coenergy_before_j": 0.01, "result_coenergy_before_j": 0.01, "coenergy_after_j": 0.010004, "result_coenergy_after_j": 0.010004, "force_n": 4.0, "result_force_n": 4.0, "force_owner": "force:v53", "result_force_owner": "force:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {HYSTERESIS_V53: hysteresis, VIRTUAL_WORK_V53: virtual_work}


# --- electromagnetic_artifact_identity_v54 (v54) ----------------------------


def _generation(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payload_v54():
    power = {**_generation("power-v54", ("power_generation", "loss_generation", "frequency_generation", "circuit_generation", "owner_generation", "result_generation")), "complex_power_va": {"real": 120.0, "imag": 45.0}, "result_complex_power_va": {"real": 120.0, "imag": 45.0}, "active_power_w": 120.0, "result_active_power_w": 120.0, "reactive_power_var": 45.0, "result_reactive_power_var": 45.0, "loss_components_w": {"copper": 80.0, "core": 40.0}, "result_loss_components_w": {"copper": 80.0, "core": 40.0}, "frequency_hz": 400.0, "result_frequency_hz": 400.0, "circuit_id": "circuit:a", "result_circuit_id": "circuit:a", "circuit_owner": "circuit-owner:v54", "result_circuit_owner": "circuit-owner:v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    selection = ["block:armature", "interface:airgap"]
    force = {**_generation("force-v54", ("radius_generation", "stress_generation", "selection_generation", "direction_generation", "owner_generation", "result_generation")), "radius_weighting": "2*pi*r", "result_radius_weighting": "2*pi*r", "stress_measure": "weighted_stress_tensor", "result_stress_measure": "weighted_stress_tensor", "integration_selection": selection, "result_integration_selection": selection, "force_direction_rz": [1.0, 0.0], "result_force_direction_rz": [1.0, 0.0], "force_n": 12.5, "result_force_n": 12.5, "mesh_owner": "mesh:v54", "result_mesh_owner": "mesh:v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {POWER_V54: power, FORCE_V54: force}


# --- electromagnetic_artifact_identity_v55 (v55) ----------------------------


def _payload_v55():
    frequency = 1000.0; conductivity = 5.8e7; permeability = 4.0e-7 * math.pi
    depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * permeability * conductivity))
    induction = {**_generation("induction-v55", ("skin_generation", "field_generation", "loss_generation", "frequency_generation", "material_generation", "owner_generation", "result_generation")), "frequency_hz": frequency, "result_frequency_hz": frequency, "conductivity_s_m": conductivity, "result_conductivity_s_m": conductivity, "permeability_h_m": permeability, "result_permeability_h_m": permeability, "skin_depth_m": depth, "result_skin_depth_m": depth, "complex_magnetic_field_a_m": {"real": 1200.0, "imag": -350.0}, "result_complex_magnetic_field_a_m": {"real": 1200.0, "imag": -350.0}, "joule_loss_w": 12.5, "result_joule_loss_w": 12.5, "conductor_owner": "conductor:v55", "result_conductor_owner": "conductor:v55", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    matrix = [[2.0e-12, -1.0e-12], [-1.0e-12, 2.0e-12]]
    capacitance = {**_generation("capacitance-v55", ("capacitance_generation", "charge_generation", "voltage_generation", "energy_generation", "symmetry_generation", "owner_generation", "result_generation")), "conductor_order": ["conductor:1", "conductor:2"], "result_conductor_order": ["conductor:1", "conductor:2"], "capacitance_matrix_f": matrix, "result_capacitance_matrix_f": matrix, "voltage_v": [1.0, 0.0], "result_voltage_v": [1.0, 0.0], "charge_c": [2.0e-12, -1.0e-12], "result_charge_c": [2.0e-12, -1.0e-12], "stored_energy_j": 1.0e-12, "result_stored_energy_j": 1.0e-12, "solution_owner": "solution:v55", "result_solution_owner": "solution:v55", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {INDUCTION_V55: induction, CAPACITANCE_V55: capacitance}


# --- electromagnetic_force_heat_identity_v56 (v56) --------------------------


def _identity_v56() -> dict[str, object]:
    generation = "magnetostatic-public-v56-test"
    generations = lambda fields: {field: generation for field in fields}
    samples = [{"displacement_m": -1.0e-4, "coenergy_j": 0.9998}, {"displacement_m": 1.0e-4, "coenergy_j": 1.0002}]
    temperatures = {"minimum_k": 293.15, "maximum_k": 331.5, "mean_k": 307.2}
    regions = {"region:coil": 32.0, "region:core": 18.0}
    return {
        FORCE_V56: {"generation": generation, **generations(("energy_generation", "coenergy_generation", "displacement_generation", "derivative_generation", "force_generation", "owner_generation", "result_generation")), "magnetic_energy_j": 0.98, "result_magnetic_energy_j": 0.98, "coenergy_samples": samples, "result_coenergy_samples": samples, "displacement_step_m": 1.0e-4, "result_displacement_step_m": 1.0e-4, "coenergy_derivative_n": 2.0, "result_coenergy_derivative_n": 2.0, "force_n": 2.0, "result_force_n": 2.0, "solution_owner": "solution:force-v56", "result_solution_owner": "solution:force-v56", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64},
        HEAT_V56: {"generation": generation, **generations(("source_generation", "temperature_generation", "flux_generation", "balance_generation", "region_generation", "owner_generation", "result_generation")), "joule_source_w": 50.0, "result_joule_source_w": 50.0, "temperature_field_k": temperatures, "result_temperature_field_k": temperatures, "outward_boundary_flux_w": 50.0, "result_outward_boundary_flux_w": 50.0, "region_joule_balance_w": regions, "result_region_joule_balance_w": regions, "balance_residual_w": 0.0, "result_balance_residual_w": 0.0, "solution_owner": "solution:heat-v56", "result_solution_owner": "solution:heat-v56", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64},
    }
