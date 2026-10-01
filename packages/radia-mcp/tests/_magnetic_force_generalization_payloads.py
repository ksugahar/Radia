"""Payload builders for test_magnetic_force_generalization.py.

Each ``_summary_vN`` / ``_identity_vN`` extends the previous cumulative payload
of its chain; standalone gates keep their own builders. Data only, not collected.
"""
from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.bem_hysteresis_identity_v48 import (
    BEM as V48_BEM,
    HYSTERESIS,
)
from radia_mcp.radia_ngsolve.bem_motion_identity_v50 import (
    BEM as V50_BEM,
    MOTION,
)
from radia_mcp.radia_ngsolve.demag_virtual_work_identity_v49 import (
    DEMAG as V49_DEMAG,
    VIRTUAL_WORK,
)
from radia_mcp.radia_ngsolve.energy_derivative_identity_v52 import (
    MAGNET_TORQUE,
    VIRTUAL_FORCE,
)
from radia_mcp.radia_ngsolve.energy_derivative_identity_v53 import (
    MAGLEV,
    QUADRATURE,
)
from radia_mcp.radia_ngsolve.energy_derivative_identity_v54 import (
    CHARGE,
    STIFFNESS,
)
from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from radia_mcp.radia_ngsolve.magnetic_force_artifact_lineage_v47 import (
    FORCE,
    MOTOR,
    MOTOR_LANES,
)
from radia_mcp.radia_ngsolve.magnetostatic_energy_identity_v55 import (
    BEARING,
    DEMAG as V55_DEMAG,
)
from radia_mcp.radia_ngsolve.potential_bem_identity_v51 import (
    BEM as V51_BEM,
    POTENTIAL,
)

from test_force_coenergy_gate import (
    _artifact_identity,
    _quadratic_case,
)
from test_magnetic_force_method_profile_gate import _summary_v22


def _summary_v23():
    summary = _summary_v22()
    identity = summary["artifact_identity"]
    identity["bem_panel_normal_material_region_demag_force_generation_identity"] = {
        "solve_generation": "bem-force-51",
        "panel_mesh_solve_generation": "bem-force-51",
        "outward_normal_solve_generation": "bem-force-51",
        "material_region_solve_generation": "bem-force-51",
        "demag_result_solve_generation": "bem-force-51",
        "force_result_solve_generation": "bem-force-51",
        "panel_ids": [101, 102, 103],
        "result_panel_ids": [101, 102, 103],
        "outward_normals": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "result_outward_normals": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        "material_region_ids": [1, 1, 2],
        "result_material_region_ids": [1, 1, 2],
        "demag_field_a_per_m": [-125000.0, -118000.0, -62000.0],
        "result_demag_field_a_per_m": [-125000.0, -118000.0, -62000.0],
        "force_vectors_n": [[12.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 3.0]],
        "result_force_vectors_n": [[12.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 3.0]],
        "force_coordinate_frame": "global_xyz",
        "result_force_coordinate_frame": "global_xyz",
        "panel_force_table_sha256": "1" * 64,
        "result_panel_force_table_sha256": "1" * 64,
    }
    identity["motor_harmonic_rotor_angle_current_phase_force_frame_generation_identity"] = {
        "sweep_generation": "motor-harmonic-51",
        "rotor_angle_sweep_generation": "motor-harmonic-51",
        "current_phase_sweep_generation": "motor-harmonic-51",
        "harmonic_bin_sweep_generation": "motor-harmonic-51",
        "force_frame_sweep_generation": "motor-harmonic-51",
        "force_result_sweep_generation": "motor-harmonic-51",
        "rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "result_rotor_angles_deg": [0.0, 5.0, 10.0, 15.0],
        "current_phase_deg": [0.0, -120.0, 120.0],
        "result_current_phase_deg": [0.0, -120.0, 120.0],
        "harmonic_bins": [0, 1, 2, 3],
        "result_harmonic_bins": [0, 1, 2, 3],
        "force_coordinate_frame": "rotor_dq",
        "result_force_coordinate_frame": "rotor_dq",
        "force_harmonics_n": [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        "result_force_harmonics_n": [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        "harmonic_force_table_sha256": "2" * 64,
        "result_harmonic_force_table_sha256": "2" * 64,
    }
    return summary


def _summary_v24():
    summary = _summary_v23()
    identity = summary["artifact_identity"]
    identity[
        "maglev_force_stiffness_equilibrium_energy_finite_difference_generation_identity"
    ] = {
        "maglev_generation": "maglev-101",
        "equilibrium_maglev_generation": "maglev-101",
        "displacement_maglev_generation": "maglev-101",
        "energy_maglev_generation": "maglev-101",
        "force_maglev_generation": "maglev-101",
        "stiffness_maglev_generation": "maglev-101",
        "coordinate_frame_maglev_generation": "maglev-101",
        "result_maglev_generation": "maglev-101",
        "equilibrium_displacement_m": 0.0,
        "result_equilibrium_displacement_m": 0.0,
        "equilibrium_force_n": 0.0,
        "result_equilibrium_force_n": 0.0,
        "displacement_samples_m": [-0.001, 0.0, 0.001],
        "energy_samples_j": [0.005, 0.0, 0.005],
        "force_samples_n": [10.0, 0.0, -10.0],
        "energy_finite_difference_force_n": [10.0, 0.0, -10.0],
        "stiffness_n_m": 10000.0,
        "reported_stiffness_n_m": 10000.0,
        "force_energy_sign_convention": "force=-dW/dx",
        "result_force_energy_sign_convention": "force=-dW/dx",
        "coordinate_frame": "global_z",
        "result_coordinate_frame": "global_z",
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "maglev_result_sha256": "2" * 64,
        "reported_maglev_result_sha256": "2" * 64,
    }
    identity[
        "motor_dual_lane_geometry_excitation_force_frame_harmonic_alignment_generation_identity"
    ] = {
        "comparison_generation": "dual-motor-101",
        "geometry_comparison_generation": "dual-motor-101",
        "excitation_comparison_generation": "dual-motor-101",
        "force_frame_comparison_generation": "dual-motor-101",
        "harmonic_comparison_generation": "dual-motor-101",
        "rotor_angle_comparison_generation": "dual-motor-101",
        "result_comparison_generation": "dual-motor-101",
        "lane_ids": ["ngsolve-age", "hdiv-mmm-hcurl-eddy-bubble"],
        "result_lane_ids": ["ngsolve-age", "hdiv-mmm-hcurl-eddy-bubble"],
        "geometry_revision_sha256": ["3" * 64, "3" * 64],
        "result_geometry_revision_sha256": ["3" * 64, "3" * 64],
        "excitation_table_sha256": ["4" * 64, "4" * 64],
        "result_excitation_table_sha256": ["4" * 64, "4" * 64],
        "force_coordinate_frames": ["rotor_dq", "rotor_dq"],
        "result_force_coordinate_frames": ["rotor_dq", "rotor_dq"],
        "harmonic_bins": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "result_harmonic_bins": [[0, 1, 2, 3], [0, 1, 2, 3]],
        "rotor_angles_deg": [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]],
        "result_rotor_angles_deg": [[0.0, 5.0, 10.0], [0.0, 5.0, 10.0]],
        "force_harmonics_n": [
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        ],
        "result_force_harmonics_n": [
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
            [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
        ],
        "comparison_result_sha256": "5" * 64,
        "reported_comparison_result_sha256": "5" * 64,
    }
    return summary


def _summary_v25():
    summary = _summary_v24()
    identity = summary["artifact_identity"]
    normals = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
    identity[
        "bem_demag_surface_orientation_magnetization_volume_material_frame_generation_identity"
    ] = {
        "solve_generation": "bem-demag-201",
        "surface_mesh_solve_generation": "bem-demag-201",
        "surface_orientation_solve_generation": "bem-demag-201",
        "magnetization_solve_generation": "bem-demag-201",
        "body_volume_solve_generation": "bem-demag-201",
        "material_region_solve_generation": "bem-demag-201",
        "coordinate_frame_solve_generation": "bem-demag-201",
        "result_solve_generation": "bem-demag-201",
        "surface_ids": [101, 102, 103, 104, 105, 106],
        "result_surface_ids": [101, 102, 103, 104, 105, 106],
        "surface_orientation": "outward_from_magnet",
        "result_surface_orientation": "outward_from_magnet",
        "outward_normals": normals,
        "result_outward_normals": normals,
        "magnetization_vector_a_per_m": [0.0, 0.0, 900000.0],
        "result_magnetization_vector_a_per_m": [0.0, 0.0, 900000.0],
        "body_volume_m3": 1.0e-6,
        "result_body_volume_m3": 1.0e-6,
        "material_region_id": 7,
        "result_material_region_id": 7,
        "coordinate_frame": "global_xyz",
        "result_coordinate_frame": "global_xyz",
        "surface_mesh_sha256": "1" * 64,
        "result_surface_mesh_sha256": "1" * 64,
        "demag_result_sha256": "2" * 64,
        "reported_demag_result_sha256": "2" * 64,
    }
    positions = [0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
    phases = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0, 360.0]
    thrust = [100.0, 102.0, 101.0, 99.0, 98.0, 100.0, 100.0]
    identity[
        "linear_motor_thrust_ripple_period_position_phase_frame_generation_identity"
    ] = {
        "sweep_generation": "linear-thrust-201",
        "period_sweep_generation": "linear-thrust-201",
        "position_sweep_generation": "linear-thrust-201",
        "phase_sweep_generation": "linear-thrust-201",
        "force_frame_sweep_generation": "linear-thrust-201",
        "sample_order_sweep_generation": "linear-thrust-201",
        "result_sweep_generation": "linear-thrust-201",
        "mechanical_period_m": 0.03,
        "result_mechanical_period_m": 0.03,
        "mover_positions_m": positions,
        "result_mover_positions_m": positions,
        "excitation_phase_deg": phases,
        "result_excitation_phase_deg": phases,
        "sample_order": list(range(7)),
        "result_sample_order": list(range(7)),
        "force_coordinate_frame": "global_x",
        "result_force_coordinate_frame": "global_x",
        "thrust_samples_n": thrust,
        "result_thrust_samples_n": thrust,
        "thrust_ripple_peak_to_peak_n": 4.0,
        "reported_thrust_ripple_peak_to_peak_n": 4.0,
        "thrust_table_sha256": "3" * 64,
        "result_thrust_table_sha256": "3" * 64,
    }
    return summary


def _summary_v26():
    summary = _summary_v25()
    identity = summary["artifact_identity"]
    identity["levitation_force_displacement_gradient_stiffness_energy_derivative_frame_generation_identity"] = {
        "levitation_generation": "levitation-301", "displacement_levitation_generation": "levitation-301",
        "force_levitation_generation": "levitation-301", "energy_levitation_generation": "levitation-301",
        "gradient_levitation_generation": "levitation-301", "frame_levitation_generation": "levitation-301",
        "result_levitation_generation": "levitation-301", "displacement_m": [-0.001, 0.0, 0.001],
        "result_displacement_m": [-0.001, 0.0, 0.001], "force_n": [10.0, 0.0, -10.0],
        "result_force_n": [10.0, 0.0, -10.0], "magnetic_energy_j": [0.005, 0.0, 0.005],
        "result_magnetic_energy_j": [0.005, 0.0, 0.005],
        "negative_energy_derivative_force_n": [10.0, 0.0, -10.0],
        "result_negative_energy_derivative_force_n": [10.0, 0.0, -10.0],
        "restoring_stiffness_n_m": 10000.0, "result_restoring_stiffness_n_m": 10000.0,
        "coordinate_frame": "global_z_up", "result_coordinate_frame": "global_z_up",
        "force_sign_convention": "restoring_negative_gradient",
        "result_force_sign_convention": "restoring_negative_gradient",
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    positions = [0.0, 7.5, 15.0, 22.5, 30.0, 37.5, 45.0]
    torque = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0]
    identity["cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_identity"] = {
        "cogging_generation": "cogging-301", "position_cogging_generation": "cogging-301",
        "periodicity_cogging_generation": "cogging-301", "mesh_cogging_generation": "cogging-301",
        "interpolation_cogging_generation": "cogging-301", "reference_cogging_generation": "cogging-301",
        "result_cogging_generation": "cogging-301", "mechanical_positions_deg": positions,
        "result_mechanical_positions_deg": positions, "cogging_torque_nm": torque,
        "result_cogging_torque_nm": torque, "periodicity": 8, "result_periodicity": 8,
        "mechanical_period_deg": 45.0, "result_mechanical_period_deg": 45.0,
        "reference_angle_deg": 0.0, "result_reference_angle_deg": 0.0,
        "interpolation_method": "periodic_cubic", "result_interpolation_method": "periodic_cubic",
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


def _summary_v27():
    summary = _summary_v26()
    identity = summary["artifact_identity"]
    normals = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ]
    identity["bem_panel_orientation_magnetization_frame_self_term_energy_force_generation_identity"] = {
        "panel_generation": "bem-panel-311",
        "orientation_panel_generation": "bem-panel-311",
        "magnetization_panel_generation": "bem-panel-311",
        "self_term_panel_generation": "bem-panel-311",
        "energy_panel_generation": "bem-panel-311",
        "force_panel_generation": "bem-panel-311",
        "mesh_panel_generation": "bem-panel-311",
        "result_panel_generation": "bem-panel-311",
        "panel_ids": [11, 12, 13, 14],
        "result_panel_ids": [11, 12, 13, 14],
        "outward_unit_normals": normals,
        "result_outward_unit_normals": normals,
        "panel_area_m2": [0.5, 0.5, 0.5, 0.5],
        "result_panel_area_m2": [0.5, 0.5, 0.5, 0.5],
        "magnetization_a_m": [[0.0, 0.0, 800000.0]] * 4,
        "result_magnetization_a_m": [[0.0, 0.0, 800000.0]] * 4,
        "magnetization_frame": "global-cartesian",
        "result_magnetization_frame": "global-cartesian",
        "singular_self_term": "analytic-solid-angle",
        "result_singular_self_term": "analytic-solid-angle",
        "displacement_m": [-0.001, 0.0, 0.001],
        "result_displacement_m": [-0.001, 0.0, 0.001],
        "magnetic_energy_j": [0.005, 0.0, 0.005],
        "result_magnetic_energy_j": [0.005, 0.0, 0.005],
        "negative_energy_derivative_force_n": [10.0, 0.0, -10.0],
        "result_force_n": [10.0, 0.0, -10.0],
        "panel_mesh_sha256": "1" * 64,
        "result_panel_mesh_sha256": "1" * 64,
        "force_result_sha256": "2" * 64,
        "accepted_force_result_sha256": "2" * 64,
    }
    identity["motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_generation_identity"] = {
        "reduced_generation": "motor-rom-311",
        "basis_reduced_generation": "motor-rom-311",
        "snapshot_reduced_generation": "motor-rom-311",
        "operating_point_reduced_generation": "motor-rom-311",
        "weight_reduced_generation": "motor-rom-311",
        "torque_reduced_generation": "motor-rom-311",
        "residual_reduced_generation": "motor-rom-311",
        "result_reduced_generation": "motor-rom-311",
        "basis_dimension": 3,
        "result_basis_dimension": 3,
        "snapshot_ids": ["snap-a", "snap-b", "snap-c"],
        "result_snapshot_ids": ["snap-a", "snap-b", "snap-c"],
        "snapshot_operating_points": [[1000.0, 10.0, 0.0], [2000.0, 20.0, 30.0], [3000.0, 30.0, 60.0]],
        "result_snapshot_operating_points": [[1000.0, 10.0, 0.0], [2000.0, 20.0, 30.0], [3000.0, 30.0, 60.0]],
        "query_operating_point": [2200.0, 22.0, 36.0],
        "result_query_operating_point": [2200.0, 22.0, 36.0],
        "interpolation_weights": [0.2, 0.6, 0.2],
        "result_interpolation_weights": [0.2, 0.6, 0.2],
        "snapshot_torque_nm": [1.0, 2.0, 3.0],
        "result_snapshot_torque_nm": [1.0, 2.0, 3.0],
        "reduced_torque_nm": 2.0,
        "result_reduced_torque_nm": 2.0,
        "relative_residual": 1.0e-6,
        "accepted_relative_residual": 1.0e-4,
        "result_relative_residual": 1.0e-6,
        "basis_sha256": "3" * 64,
        "loaded_basis_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _summary_v28():
    summary = _summary_v27()
    identity = summary["artifact_identity"]
    identity[
        "maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_generation_identity"
    ] = {
        "stiffness_generation": "maglev-stiffness-321",
        "displacement_stiffness_generation": "maglev-stiffness-321",
        "coordinate_stiffness_generation": "maglev-stiffness-321",
        "geometry_stiffness_generation": "maglev-stiffness-321",
        "mesh_stiffness_generation": "maglev-stiffness-321",
        "force_stiffness_generation": "maglev-stiffness-321",
        "derivative_stiffness_generation": "maglev-stiffness-321",
        "solution_stiffness_generation": "maglev-stiffness-321",
        "result_stiffness_generation": "maglev-stiffness-321",
        "displacement_m": [-0.001, 0.0, 0.001],
        "result_displacement_m": [-0.001, 0.0, 0.001],
        "displacement_step_m": 0.001,
        "result_displacement_step_m": 0.001,
        "coordinate_direction": "global-z-positive",
        "result_coordinate_direction": "global-z-positive",
        "force_n": [12.0, 10.0, 8.0],
        "result_force_n": [12.0, 10.0, 8.0],
        "derivative_convention": "stiffness-equals-negative-force-derivative",
        "result_derivative_convention": "stiffness-equals-negative-force-derivative",
        "stiffness_n_m": 2000.0,
        "result_stiffness_n_m": 2000.0,
        "geometry_sha256": "1" * 64,
        "result_geometry_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "solution_sha256": "3" * 64,
        "accepted_solution_sha256": "3" * 64,
    }
    identity[
        "motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_generation_identity"
    ] = {
        "torque_generation": "motor-coenergy-321",
        "winding_torque_generation": "motor-coenergy-321",
        "harmonic_torque_generation": "motor-coenergy-321",
        "current_torque_generation": "motor-coenergy-321",
        "phase_torque_generation": "motor-coenergy-321",
        "angle_torque_generation": "motor-coenergy-321",
        "coenergy_torque_generation": "motor-coenergy-321",
        "mesh_torque_generation": "motor-coenergy-321",
        "result_torque_generation": "motor-coenergy-321",
        "phase_order": ["U", "V", "W"],
        "result_phase_order": ["U", "V", "W"],
        "harmonic_orders": [1, 5, 7],
        "result_harmonic_orders": [1, 5, 7],
        "phase_current_harmonic_a": [
            [100.0, -50.0, -50.0],
            [8.0, -4.0, -4.0],
            [5.0, -2.5, -2.5],
        ],
        "result_phase_current_harmonic_a": [
            [100.0, -50.0, -50.0],
            [8.0, -4.0, -4.0],
            [5.0, -2.5, -2.5],
        ],
        "current_phase_deg": [0.0, -120.0, 120.0],
        "result_current_phase_deg": [0.0, -120.0, 120.0],
        "rotor_mechanical_angle_deg": [0.0, 1.0, 2.0],
        "result_rotor_mechanical_angle_deg": [0.0, 1.0, 2.0],
        "coenergy_j": [0.100, 0.102, 0.104],
        "result_coenergy_j": [0.100, 0.102, 0.104],
        "torque_convention": "positive-coenergy-angle-derivative",
        "result_torque_convention": "positive-coenergy-angle-derivative",
        "torque_nm": [0.11459155902616465, 0.11459155902616465],
        "result_torque_nm": [0.11459155902616465, 0.11459155902616465],
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return summary


def _summary_v29():
    summary = _summary_v28()
    identity = summary["artifact_identity"]
    identity["magnet_demag_recoil_knee_field_volume_generation_identity"] = {
        "demag_generation": "demag-state-331",
        "recoil_demag_generation": "demag-state-331",
        "knee_demag_generation": "demag-state-331",
        "temperature_demag_generation": "demag-state-331",
        "field_demag_generation": "demag-state-331",
        "mask_demag_generation": "demag-state-331",
        "volume_demag_generation": "demag-state-331",
        "mesh_demag_generation": "demag-state-331",
        "result_demag_generation": "demag-state-331",
        "recoil_relative_permeability": 1.05,
        "result_recoil_relative_permeability": 1.05,
        "knee_field_a_m": -700000.0,
        "result_knee_field_a_m": -700000.0,
        "temperature_c": 120.0,
        "result_temperature_c": 120.0,
        "element_ids": [1, 2, 3],
        "result_element_ids": [1, 2, 3],
        "local_recoil_axis_field_a_m": [-600000.0, -800000.0, -500000.0],
        "result_local_recoil_axis_field_a_m": [-600000.0, -800000.0, -500000.0],
        "irreversible_mask": [False, True, False],
        "result_irreversible_mask": [False, True, False],
        "element_volumes_m3": [1.0e-6, 2.0e-6, 1.0e-6],
        "result_element_volumes_m3": [1.0e-6, 2.0e-6, 1.0e-6],
        "magnet_volume_m3": 4.0e-6,
        "result_magnet_volume_m3": 4.0e-6,
        "irreversible_volume_fraction": 0.5,
        "result_irreversible_volume_fraction": 0.5,
        "material_state_sha256": "1" * 64,
        "result_material_state_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    identity["linear_motor_end_phase_wave_pitch_force_generation_identity"] = {
        "linear_motor_generation": "linear-force-331",
        "end_effect_linear_motor_generation": "linear-force-331",
        "phase_linear_motor_generation": "linear-force-331",
        "wave_linear_motor_generation": "linear-force-331",
        "pitch_linear_motor_generation": "linear-force-331",
        "position_linear_motor_generation": "linear-force-331",
        "force_linear_motor_generation": "linear-force-331",
        "ripple_linear_motor_generation": "linear-force-331",
        "result_linear_motor_generation": "linear-force-331",
        "phase_sequence": ["U", "V", "W"],
        "result_phase_sequence": ["U", "V", "W"],
        "traveling_wave_direction": "global-x-positive",
        "result_traveling_wave_direction": "global-x-positive",
        "pole_pitch_m": 0.06,
        "result_pole_pitch_m": 0.06,
        "position_m": [0.0, 0.015, 0.03, 0.045, 0.06],
        "result_position_m": [0.0, 0.015, 0.03, 0.045, 0.06],
        "end_effect_factor": [0.9, 1.0, 1.0, 1.0, 0.9],
        "result_end_effect_factor": [0.9, 1.0, 1.0, 1.0, 0.9],
        "force_n": [100.0, 110.0, 100.0, 90.0, 100.0],
        "result_force_n": [100.0, 110.0, 100.0, 90.0, 100.0],
        "mean_force_n": 100.0,
        "result_mean_force_n": 100.0,
        "force_ripple_peak_to_peak_n": 20.0,
        "result_force_ripple_peak_to_peak_n": 20.0,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _force_summary_v29():
    positions, _, _ = _quadratic_case()
    identity = _artifact_identity(len(positions))
    generation = "airgap-stress-161"
    identity["airgap_stress_harmonic_sector_periodicity_origin_sampling_alias_radius_torque_generation_identity"] = {
        "stress_generation": generation, "harmonic_stress_generation": generation,
        "sector_stress_generation": generation, "sampling_stress_generation": generation,
        "alias_stress_generation": generation, "geometry_stress_generation": generation,
        "mesh_stress_generation": generation, "result_stress_generation": generation,
        "sector_pitch_deg": 30.0, "result_sector_pitch_deg": 30.0,
        "sector_count": 12, "result_sector_count": 12,
        "angular_origin_deg": 0.0, "result_angular_origin_deg": 0.0,
        "angular_sample_count": 720, "result_angular_sample_count": 720,
        "sector_sample_count": 60, "result_sector_sample_count": 60,
        "harmonic_orders": [0, 6, 12, 18], "result_harmonic_orders": [0, 6, 12, 18],
        "torque_harmonics_nm": [5.0, 0.2, 0.05, 0.01],
        "result_torque_harmonics_nm": [5.0, 0.2, 0.05, 0.01],
        "alias_filter": "truncate_below_nyquist", "result_alias_filter": "truncate_below_nyquist",
        "alias_cutoff_order": 24, "result_alias_cutoff_order": 24,
        "airgap_radius_m": 0.05, "result_airgap_radius_m": 0.05,
        "axial_length_m": 0.1, "result_axial_length_m": 0.1,
        "torque_nm": 5.26, "result_torque_nm": 5.26,
        "airgap_mesh_sha256": "1" * 64, "result_airgap_mesh_sha256": "1" * 64,
        "torque_result_sha256": "2" * 64, "accepted_torque_result_sha256": "2" * 64,
    }
    frequency, flux, thickness, volume = 400.0, 1.2, 0.00035, 0.001
    kh, alpha, ke, kex = 50.0, 1.6, 0.02, 0.5
    hysteresis = kh * frequency * flux**alpha * volume
    eddy = ke * frequency**2 * flux**2 * thickness**2 * volume
    excess = kex * frequency**1.5 * flux**1.5 * volume
    generation = "laminated-loss-161"
    identity["laminated_core_hysteresis_eddy_excess_frequency_flux_lamination_volume_result_generation_identity"] = {
        "loss_generation": generation, "hysteresis_loss_generation": generation,
        "eddy_loss_generation": generation, "excess_loss_generation": generation,
        "frequency_loss_generation": generation, "flux_loss_generation": generation,
        "lamination_loss_generation": generation, "volume_loss_generation": generation,
        "result_loss_generation": generation,
        "frequency_hz": frequency, "result_frequency_hz": frequency,
        "peak_flux_density_t": flux, "result_peak_flux_density_t": flux,
        "lamination_thickness_m": thickness, "result_lamination_thickness_m": thickness,
        "magnetic_volume_m3": volume, "result_magnetic_volume_m3": volume,
        "hysteresis_coefficient": kh, "result_hysteresis_coefficient": kh,
        "hysteresis_exponent": alpha, "result_hysteresis_exponent": alpha,
        "eddy_coefficient": ke, "result_eddy_coefficient": ke,
        "excess_coefficient": kex, "result_excess_coefficient": kex,
        "hysteresis_loss_w": hysteresis, "result_hysteresis_loss_w": hysteresis,
        "eddy_loss_w": eddy, "result_eddy_loss_w": eddy,
        "excess_loss_w": excess, "result_excess_loss_w": excess,
        "total_core_loss_w": hysteresis + eddy + excess,
        "result_total_core_loss_w": hysteresis + eddy + excess,
        "material_sha256": "3" * 64, "result_material_sha256": "3" * 64,
        "loss_result_sha256": "4" * 64, "accepted_loss_result_sha256": "4" * 64,
    }
    return identity


def _force_gate(identity):
    positions, coenergy, forces = _quadratic_case()
    return force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )


def _identity_v30():
    identity = _force_summary_v29()
    generation = "nonlinear-coenergy-171"
    dx = 1.0e-4
    coenergy = [2.005, 2.0, 1.995]
    force = -(coenergy[2] - coenergy[0]) / (2.0 * dx)
    identity[
        "nonlinear_coenergy_force_current_perturbation_remesh_central_difference_frame_result_identity"
    ] = {
        "force_generation": generation,
        "current_force_generation": generation,
        "coenergy_force_generation": generation,
        "remesh_force_generation": generation,
        "difference_force_generation": generation,
        "frame_force_generation": generation,
        "result_force_generation": generation,
        "nonlinear_material": True,
        "result_nonlinear_material": True,
        "current_constraint": "fixed_current",
        "result_current_constraint": "fixed_current",
        "nominal_current_a": 10.0,
        "branch_currents_a": [10.0, 10.0, 10.0],
        "result_branch_currents_a": [10.0, 10.0, 10.0],
        "displacements_m": [-dx, 0.0, dx],
        "result_displacements_m": [-dx, 0.0, dx],
        "coenergy_j": coenergy,
        "result_coenergy_j": coenergy,
        "difference_rule": "negative_central_difference",
        "result_difference_rule": "negative_central_difference",
        "branch_mesh_generations": ["mesh-minus-171", "mesh-center-171", "mesh-plus-171"],
        "result_branch_mesh_generations": ["mesh-minus-171", "mesh-center-171", "mesh-plus-171"],
        "displacement_frame": "global_x",
        "result_displacement_frame": "global_x",
        "force_n": force,
        "result_force_n": force,
        "branch_result_sha256": "1" * 64,
        "accepted_branch_result_sha256": "1" * 64,
    }
    generation = "axisym-stress-171"
    radius, jacobian, stress = 0.025, 0.001, 1200.0
    radial_weight = 2.0 * math.pi * radius
    force = stress * radial_weight * jacobian
    identity[
        "axisymmetric_force_radial_weight_jacobian_coordinate_stress_contour_material_mesh_result_identity"
    ] = {
        "axisymmetric_generation": generation,
        "radial_weight_generation": generation,
        "jacobian_generation": generation,
        "coordinate_generation": generation,
        "stress_contour_generation": generation,
        "material_side_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "coordinate_convention": "r_z_axisymmetric",
        "result_coordinate_convention": "r_z_axisymmetric",
        "radius_m": radius,
        "result_radius_m": radius,
        "radial_weight": radial_weight,
        "result_radial_weight": radial_weight,
        "line_jacobian_m": jacobian,
        "result_line_jacobian_m": jacobian,
        "stress_normal_pa": stress,
        "result_stress_normal_pa": stress,
        "stress_contour_closed": True,
        "result_stress_contour_closed": True,
        "stress_contour_orientation": "counterclockwise",
        "result_stress_contour_orientation": "counterclockwise",
        "stress_contour_material_side": "air",
        "result_stress_contour_material_side": "air",
        "force_n": force,
        "result_force_n": force,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "force_result_sha256": "3" * 64,
        "accepted_force_result_sha256": "3" * 64,
    }
    return identity


def _gate(identity):
    positions, coenergy, forces = _quadratic_case()
    return force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )


def _summary_v30():
    summary = _summary_v29(); identity = summary["artifact_identity"]
    generation = "maglev-stiffness-341"
    identity["maglev_force_stiffness_position_current_derivative_frame_mesh_result_identity"] = {
        "maglev_generation": generation, "position_maglev_generation": generation,
        "current_maglev_generation": generation, "force_maglev_generation": generation,
        "derivative_maglev_generation": generation, "frame_maglev_generation": generation,
        "mesh_maglev_generation": generation, "result_maglev_generation": generation,
        "position_m": [-0.001, 0.0, 0.001], "result_position_m": [-0.001, 0.0, 0.001],
        "current_a": [10.0, 10.0, 10.0], "result_current_a": [10.0, 10.0, 10.0],
        "force_z_n": [12.0, 10.0, 8.0], "result_force_z_n": [12.0, 10.0, 8.0],
        "stiffness_n_per_m": -2000.0, "result_stiffness_n_per_m": -2000.0,
        "coordinate_frame": "global_z_positive_up", "result_coordinate_frame": "global_z_positive_up",
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "cogging-341"
    identity["cogging_torque_slot_pole_period_origin_sampling_harmonic_phase_mesh_result_identity"] = {
        "cogging_generation": generation, "slot_cogging_generation": generation,
        "pole_cogging_generation": generation, "period_cogging_generation": generation,
        "sampling_cogging_generation": generation, "harmonic_cogging_generation": generation,
        "mesh_cogging_generation": generation, "result_cogging_generation": generation,
        "slot_count": 12, "result_slot_count": 12, "pole_count": 10, "result_pole_count": 10,
        "cogging_period_mechanical_deg": 6.0, "result_cogging_period_mechanical_deg": 6.0,
        "angular_origin_deg": 0.0, "result_angular_origin_deg": 0.0,
        "sample_angles_deg": [0.0, 1.5, 3.0, 4.5, 6.0],
        "result_sample_angles_deg": [0.0, 1.5, 3.0, 4.5, 6.0],
        "harmonic_orders": [1, 2], "result_harmonic_orders": [1, 2],
        "harmonic_phase_deg": [0.0, 90.0], "result_harmonic_phase_deg": [0.0, 90.0],
        "torque_nm": [0.0, 0.1, 0.0, -0.1, 0.0], "result_torque_nm": [0.0, 0.1, 0.0, -0.1, 0.0],
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


def _identity_v31():
    identity = _identity_v30()
    generation = "incremental-inductance-181"
    identity[
        "nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mesh_result_identity"
    ] = {
        "analysis_generation": generation,
        **{key: generation for key in (
            "current_path_generation", "magnetic_state_generation",
            "perturbation_generation", "derivative_generation",
            "circuit_generation", "mesh_generation", "result_generation",
        )},
        "nonlinear_material": True,
        "result_nonlinear_material": True,
        "current_path_a": [8.0, 9.0, 10.0, 11.0, 12.0],
        "result_current_path_a": [8.0, 9.0, 10.0, 11.0, 12.0],
        "current_state_index": 2,
        "result_current_state_index": 2,
        "nominal_current_a": 10.0,
        "result_nominal_current_a": 10.0,
        "perturbation_a": 0.1,
        "result_perturbation_a": 0.1,
        "current_samples_a": [9.9, 10.0, 10.1],
        "result_current_samples_a": [9.9, 10.0, 10.1],
        "flux_linkage_wb_turn": [0.792, 0.8, 0.808],
        "result_flux_linkage_wb_turn": [0.792, 0.8, 0.808],
        "derivative_rule": "symmetric_central_difference",
        "result_derivative_rule": "symmetric_central_difference",
        "incremental_inductance_h": 0.08,
        "result_incremental_inductance_h": 0.08,
        "circuit_name": "coil_a",
        "result_circuit_name": "coil_a",
        "magnetic_state_sha256": "1" * 64,
        "result_magnetic_state_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    generation = "lamination-loss-181"
    components = {"hysteresis_w": 3.2, "classical_eddy_w": 1.1, "excess_w": 0.4}
    coefficients = {"kh": 0.021, "kc": 0.00031, "ke": 0.0012}
    identity[
        "lamination_anisotropy_fill_orientation_frequency_loss_volume_balance_result_identity"
    ] = {
        "analysis_generation": generation,
        **{key: generation for key in (
            "anisotropy_generation", "fill_generation", "orientation_generation",
            "frequency_generation", "loss_generation", "volume_generation",
            "result_generation",
        )},
        "mu_axis_order": ["rolling", "transverse"],
        "result_mu_axis_order": ["rolling", "transverse"],
        "relative_permeability_axes": [1500.0, 40.0],
        "result_relative_permeability_axes": [1500.0, 40.0],
        "lamination_fill_factor": 0.95,
        "result_lamination_fill_factor": 0.95,
        "stacking_direction": "global_z",
        "result_stacking_direction": "global_z",
        "material_frame": "local_xy",
        "result_material_frame": "local_xy",
        "frequency_hz": 400.0,
        "result_frequency_hz": 400.0,
        "loss_coefficient_basis": "bertotti_three_term",
        "result_loss_coefficient_basis": "bertotti_three_term",
        "loss_coefficients": coefficients,
        "result_loss_coefficients": dict(coefficients),
        "gross_volume_m3": 0.001,
        "result_gross_volume_m3": 0.001,
        "active_iron_volume_m3": 0.00095,
        "result_active_iron_volume_m3": 0.00095,
        "loss_components_w": components,
        "result_loss_components_w": dict(components),
        "total_core_loss_w": 4.7,
        "result_total_core_loss_w": 4.7,
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return identity


def _summary_v31():
    summary = _summary_v30(); identity = summary["artifact_identity"]
    g = "near-gap-bem-351"
    identity["bem_near_singular_gap_quadrature_normal_order_force_reciprocity_geometry_result_identity"] = {
        "bem_generation": g, **{key: g for key in ("gap_bem_generation", "quadrature_bem_generation", "normal_bem_generation", "order_bem_generation", "force_bem_generation", "reciprocity_bem_generation", "geometry_bem_generation", "result_bem_generation")},
        "gap_m": 2e-5, "result_gap_m": 2e-5, "panel_size_m": 1e-3, "result_panel_size_m": 1e-3,
        "quadrature_policy": "gap_adaptive_duffy", "result_quadrature_policy": "gap_adaptive_duffy", "quadrature_order": 12, "result_quadrature_order": 12,
        "source_normal": [0., 0., 1.], "result_source_normal": [0., 0., 1.], "target_normal": [0., 0., -1.], "result_target_normal": [0., 0., -1.],
        "source_target_order": ["body_a", "body_b"], "result_source_target_order": ["body_a", "body_b"],
        "force_on_source_n": [0., 0., 5.], "result_force_on_source_n": [0., 0., 5.], "force_on_target_n": [0., 0., -5.], "result_force_on_target_n": [0., 0., -5.],
        "action_reaction_residual_n": 0., "result_action_reaction_residual_n": 0., "geometry_sha256": "1" * 64, "result_geometry_sha256": "1" * 64, "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    g = "minor-loop-351"
    identity["hysteresis_minor_loop_state_reversal_return_memory_remanence_energy_time_material_identity"] = {
        "loop_generation": g, **{key: g for key in ("state_loop_generation", "reversal_loop_generation", "memory_loop_generation", "remanence_loop_generation", "energy_loop_generation", "time_loop_generation", "material_loop_generation", "result_loop_generation")},
        "initial_state_sha256": "3" * 64, "result_initial_state_sha256": "3" * 64, "time_s": [0., .1, .2, .3, .4], "result_time_s": [0., .1, .2, .3, .4],
        "drive_h_a_per_m": [0., 100., 20., 100., 0.], "result_drive_h_a_per_m": [0., 100., 20., 100., 0.], "reversal_indices": [1, 2, 3], "result_reversal_indices": [1, 2, 3],
        "return_point_memory_closed": True, "result_return_point_memory_closed": True, "remanence_t": .35, "result_remanence_t": .35, "loop_energy_j_per_m3": 42., "result_loop_energy_j_per_m3": 42.,
        "material_owner": "steel-a", "result_material_owner": "steel-a", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


def _identity_v32():
    identity = _identity_v31()
    generation = "axisymmetric-force-closure-191"
    identity[
        "axisymmetric_weighted_stress_coenergy_contour_displacement_radius_weight_material_mesh_owner_result_identity"
    ] = {
        "force_generation": generation,
        **{
            key: generation
            for key in (
                "stress_generation",
                "coenergy_generation",
                "contour_generation",
                "displacement_generation",
                "weight_generation",
                "material_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "stress_method": "weighted_stress_tensor",
        "result_stress_method": "weighted_stress_tensor",
        "coenergy_method": "symmetric_virtual_displacement",
        "result_coenergy_method": "symmetric_virtual_displacement",
        "contour_radius_m": 0.025,
        "result_contour_radius_m": 0.025,
        "virtual_displacement_m": 1.0e-5,
        "result_virtual_displacement_m": 1.0e-5,
        "axisymmetric_weight": "2*pi*r",
        "result_axisymmetric_weight": "2*pi*r",
        "material_side": "air_gap",
        "result_material_side": "air_gap",
        "weighted_stress_force_n": 12.5,
        "coenergy_force_n": 12.5,
        "force_relative_tolerance": 1.0e-8,
        "force_mesh_sha256": "1" * 64,
        "result_force_mesh_sha256": "1" * 64,
        "force_owner": "axisymmetric/plunger/group1",
        "result_force_owner": "axisymmetric/plunger/group1",
        "force_result_sha256": "2" * 64,
        "accepted_force_result_sha256": "2" * 64,
    }
    generation = "laminated-diffusion-closure-191"
    conductivity = [[2.0e6, 0.0], [0.0, 1.0e3]]
    identity[
        "laminated_diffusion_conductivity_skin_depth_frequency_phasor_power_volume_mesh_loss_result_identity"
    ] = {
        "diffusion_generation": generation,
        **{
            key: generation
            for key in (
                "conductivity_generation",
                "skin_depth_generation",
                "frequency_generation",
                "phasor_generation",
                "power_generation",
                "volume_generation",
                "mesh_generation",
                "loss_generation",
                "result_generation",
            )
        },
        "conductivity_tensor_s_per_m": conductivity,
        "result_conductivity_tensor_s_per_m": [row[:] for row in conductivity],
        "skin_depth_m": 5.0e-4,
        "result_skin_depth_m": 5.0e-4,
        "frequency_hz": 400.0,
        "result_frequency_hz": 400.0,
        "phasor_convention": "exp(+jwt)_rms",
        "result_phasor_convention": "exp(+jwt)_rms",
        "complex_power_va_ri": [4.7, 1.2],
        "result_complex_power_va_ri": [4.7, 1.2],
        "active_volume_m3": 9.5e-4,
        "result_active_volume_m3": 9.5e-4,
        "laminated_mesh_sha256": "3" * 64,
        "result_laminated_mesh_sha256": "3" * 64,
        "laminated_loss_w": 4.7,
        "result_laminated_loss_w": 4.7,
        "laminated_result_sha256": "4" * 64,
        "accepted_laminated_result_sha256": "4" * 64,
    }
    return identity


def _summary_v32():
    summary = _summary_v31()
    identity = summary["artifact_identity"]
    generation = "maglev-equilibrium-361"
    identity[
        "maglev_equilibrium_force_displacement_derivative_stiffness_gravity_mesh_result_identity"
    ] = {
        "equilibrium_generation": generation,
        **{
            key: generation
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
        },
        "force_sign_convention": "positive_up",
        "result_force_sign_convention": "positive_up",
        "displacement_frame": "global_z_up",
        "result_displacement_frame": "global_z_up",
        "displacement_samples_m": [-1.0e-4, 0.0, 1.0e-4],
        "result_displacement_samples_m": [-1.0e-4, 0.0, 1.0e-4],
        "magnetic_force_samples_n": [10.001, 9.81, 9.619],
        "result_magnetic_force_samples_n": [10.001, 9.81, 9.619],
        "derivative_stencil": "symmetric_central_difference",
        "result_derivative_stencil": "symmetric_central_difference",
        "force_derivative_n_per_m": -1910.0,
        "result_force_derivative_n_per_m": -1910.0,
        "vertical_stiffness_n_per_m": 1910.0,
        "result_vertical_stiffness_n_per_m": 1910.0,
        "gravity_force_n": -9.81,
        "result_gravity_force_n": -9.81,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_owner": "maglev/case-361/equilibrium-z",
        "accepted_result_owner": "maglev/case-361/equilibrium-z",
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "bem-charge-361"
    identity[
        "bem_surface_charge_gauge_normal_energy_reciprocity_geometry_owner_result_identity"
    ] = {
        "bem_generation": generation,
        **{
            key: generation
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
        },
        "net_surface_charge": 0.0,
        "result_net_surface_charge": 0.0,
        "charge_balance_tolerance": 1.0e-12,
        "gauge_reference": "mean_zero_scalar_potential",
        "result_gauge_reference": "mean_zero_scalar_potential",
        "source_normal": [0.0, 0.0, 1.0],
        "result_source_normal": [0.0, 0.0, 1.0],
        "target_normal": [0.0, 0.0, -1.0],
        "result_target_normal": [0.0, 0.0, -1.0],
        "field_energy_j": 0.25,
        "result_field_energy_j": 0.25,
        "reciprocity_residual": 1.0e-12,
        "result_reciprocity_residual": 1.0e-12,
        "reciprocity_tolerance": 1.0e-9,
        "geometry_sha256": "3" * 64,
        "result_geometry_sha256": "3" * 64,
        "result_owner": "bem/case-361/scalar-potential",
        "accepted_result_owner": "bem/case-361/scalar-potential",
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _identity_v33():
    identity = _identity_v32()
    generation = "electrostatic-capacitance-201"
    matrix = [
        [2e-12, -1e-12, -1e-12],
        [-1e-12, 2e-12, -1e-12],
        [-1e-12, -1e-12, 2e-12],
    ]
    identity[
        "electrostatic_capacitance_conductor_reciprocity_neutrality_voltage_charge_energy_unit_mesh_owner_result_identity"
    ] = {
        "electrostatic_generation": generation,
        **{
            key: generation
            for key in (
                "conductor_generation",
                "reciprocity_generation",
                "neutrality_generation",
                "voltage_generation",
                "charge_generation",
                "energy_generation",
                "unit_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "conductor_order": ["left", "shield", "right"],
        "result_conductor_order": ["left", "shield", "right"],
        "capacitance_matrix_f": matrix,
        "result_capacitance_matrix_f": [row[:] for row in matrix],
        "voltages_v": [1.0, 0.0, -1.0],
        "result_voltages_v": [1.0, 0.0, -1.0],
        "charges_c": [3e-12, 0.0, -3e-12],
        "result_charges_c": [3e-12, 0.0, -3e-12],
        "electrostatic_energy_j": 3e-12,
        "result_electrostatic_energy_j": 3e-12,
        "capacitance_unit": "F",
        "result_capacitance_unit": "F",
        "charge_unit": "C",
        "result_charge_unit": "C",
        "voltage_unit": "V",
        "result_voltage_unit": "V",
        "energy_unit": "J",
        "result_energy_unit": "J",
        "electrostatic_mesh_sha256": "1" * 64,
        "result_electrostatic_mesh_sha256": "1" * 64,
        "result_owner": "electrostatic/capacitance-matrix-201",
        "accepted_result_owner": "electrostatic/capacitance-matrix-201",
        "electrostatic_result_sha256": "2" * 64,
        "accepted_electrostatic_result_sha256": "2" * 64,
    }
    generation = "axisymmetric-heat-balance-201"
    identity[
        "axisymmetric_heat_conduction_convection_source_boundary_flux_weight_temperature_mesh_owner_result_identity"
    ] = {
        "heat_generation": generation,
        **{
            key: generation
            for key in (
                "conduction_generation",
                "convection_generation",
                "source_generation",
                "boundary_generation",
                "weight_generation",
                "temperature_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "conduction_w": 40.0,
        "result_conduction_w": 40.0,
        "convection_w": 30.0,
        "result_convection_w": 30.0,
        "volume_source_w": 100.0,
        "result_volume_source_w": 100.0,
        "boundary_flux_w": 30.0,
        "result_boundary_flux_w": 30.0,
        "axisymmetric_weight": "2*pi*r",
        "result_axisymmetric_weight": "2*pi*r",
        "temperature_reference_k": 293.15,
        "result_temperature_reference_k": 293.15,
        "heat_balance_tolerance_w": 1e-9,
        "heat_mesh_sha256": "3" * 64,
        "result_heat_mesh_sha256": "3" * 64,
        "heat_result_owner": "axisymmetric/thermal-body-201",
        "result_heat_result_owner": "axisymmetric/thermal-body-201",
        "heat_result_sha256": "4" * 64,
        "accepted_heat_result_sha256": "4" * 64,
    }
    return identity


def _summary_v33():
    summary = _summary_v32()
    identity = summary["artifact_identity"]
    generation = "halbach-harmonic-371"
    identity[
        "halbach_harmonic_magnetization_order_pitch_phase_grid_field_energy_force_geometry_owner_result_identity"
    ] = {
        "halbach_generation": generation,
        **{
            key: generation
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
        },
        "magnetization_angles_deg": [0.0, 90.0, 180.0, 270.0],
        "result_magnetization_angles_deg": [0.0, 90.0, 180.0, 270.0],
        "pole_pitch_m": 0.02,
        "result_pole_pitch_m": 0.02,
        "harmonic_orders": [1, 3, 5],
        "result_harmonic_orders": [1, 3, 5],
        "harmonic_phase_deg": [0.0, 30.0, -15.0],
        "result_harmonic_phase_deg": [0.0, 30.0, -15.0],
        "sampling_grid_m": [0.0, 0.005, 0.01, 0.015, 0.02],
        "result_sampling_grid_m": [0.0, 0.005, 0.01, 0.015, 0.02],
        "field_harmonic_amplitude_t": [1.0, 0.1, 0.03],
        "result_field_harmonic_amplitude_t": [1.0, 0.1, 0.03],
        "magnetic_energy_j": 0.5,
        "result_magnetic_energy_j": 0.5,
        "force_direction": "+x",
        "result_force_direction": "+x",
        "force_n": 10.0,
        "result_force_n": 10.0,
        "geometry_sha256": "1" * 64,
        "result_geometry_sha256": "1" * 64,
        "result_owner": "halbach/case-371/harmonic-line",
        "accepted_result_owner": "halbach/case-371/harmonic-line",
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "magnetic-bearing-371"
    stiffness = [[1000.0, 50.0], [50.0, 900.0]]
    eigenvalues = [879.2893218813452, 1020.7106781186548]
    identity[
        "magnetic_bearing_force_current_displacement_stiffness_reciprocity_bias_frame_mesh_result_identity"
    ] = {
        "bearing_generation": generation,
        **{
            key: generation
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
        },
        "coordinate_frame": "global_xy_right_handed",
        "result_coordinate_frame": "global_xy_right_handed",
        "bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "result_bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "bias_displacement_m": [0.0, 0.0],
        "result_bias_displacement_m": [0.0, 0.0],
        "force_current_jacobian_n_per_a": [[10.0, -10.0, 0.0, 0.0], [0.0, 0.0, 10.0, -10.0]],
        "result_force_current_jacobian_n_per_a": [[10.0, -10.0, 0.0, 0.0], [0.0, 0.0, 10.0, -10.0]],
        "force_displacement_jacobian_n_per_m": [[-1000.0, -50.0], [-50.0, -900.0]],
        "result_force_displacement_jacobian_n_per_m": [[-1000.0, -50.0], [-50.0, -900.0]],
        "stiffness_matrix_n_per_m": stiffness,
        "result_stiffness_matrix_n_per_m": [list(row) for row in stiffness],
        "stiffness_eigenvalues_n_per_m": eigenvalues,
        "result_stiffness_eigenvalues_n_per_m": list(eigenvalues),
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _identity_v34():
    identity = _identity_v33()
    generation = "kelvin-open-boundary-211"
    identity[
        "kelvin_transform_radius_permeability_jacobian_interface_energy_flux_far_field_mesh_owner_result_identity"
    ] = {
        "kelvin_generation": generation,
        **{key: generation for key in (
            "radius_generation", "mapping_generation", "jacobian_generation",
            "interface_generation", "energy_generation", "flux_generation",
            "far_field_generation", "mesh_generation", "owner_generation", "result_generation")},
        "kelvin_radius_m": 1.0, "result_kelvin_radius_m": 1.0,
        "mapped_permeability_relative": [1.0, 4.0, 16.0],
        "result_mapped_permeability_relative": [1.0, 4.0, 16.0],
        "mapping_jacobian_determinants": [1.0, 0.25, 0.0625],
        "result_mapping_jacobian_determinants": [1.0, 0.25, 0.0625],
        "interface_potential_jump": 0.0, "result_interface_potential_jump": 0.0,
        "interface_normal_flux_jump_wb": 0.0, "result_interface_normal_flux_jump_wb": 0.0,
        "magnetic_energy_j": 0.012, "result_magnetic_energy_j": 0.012,
        "inner_flux_wb": 0.02, "result_inner_flux_wb": 0.02,
        "outer_flux_wb": -0.02, "result_outer_flux_wb": -0.02,
        "far_field_radius_m": [2.0, 4.0, 8.0], "result_far_field_radius_m": [2.0, 4.0, 8.0],
        "far_field_potential": [0.5, 0.25, 0.125], "result_far_field_potential": [0.5, 0.25, 0.125],
        "kelvin_mesh_sha256": "1" * 64, "result_kelvin_mesh_sha256": "1" * 64,
        "kelvin_result_owner": "open-boundary/kelvin-211",
        "accepted_kelvin_result_owner": "open-boundary/kelvin-211",
        "kelvin_result_sha256": "2" * 64, "accepted_kelvin_result_sha256": "2" * 64,
    }
    generation = "airgap-force-211"
    identity[
        "force_airgap_contour_stress_fourier_virtual_work_symmetry_torque_origin_field_result_identity"
    ] = {
        "force_generation": generation,
        **{key: generation for key in (
            "contour_generation", "stress_generation", "fourier_generation",
            "virtual_work_generation", "symmetry_generation", "torque_generation",
            "field_generation", "owner_generation", "result_generation")},
        "airgap_contour_ids": [1, 2], "result_airgap_contour_ids": [1, 2],
        "contour_forces_n": [[10.0, 0.0], [10.00000001, 0.0]],
        "result_contour_forces_n": [[10.0, 0.0], [10.00000001, 0.0]],
        "fourier_stress_harmonics_n": [[0, 10.000000005], [1, 0.0], [2, 0.0]],
        "result_fourier_stress_harmonics_n": [[0, 10.000000005], [1, 0.0], [2, 0.0]],
        "virtual_work_force_n": 10.000000005, "result_virtual_work_force_n": 10.000000005,
        "virtual_work_displacement_m": 1.0e-5, "result_virtual_work_displacement_m": 1.0e-5,
        "torque_origin_m": [0.0, 0.0], "result_torque_origin_m": [0.0, 0.0],
        "torque_nm": 0.0, "result_torque_nm": 0.0,
        "force_symmetry": "mirror_y", "result_force_symmetry": "mirror_y",
        "force_field_sha256": "3" * 64, "result_force_field_sha256": "3" * 64,
        "force_result_owner": "force/airgap-211", "accepted_force_result_owner": "force/airgap-211",
        "force_result_sha256": "4" * 64, "accepted_force_result_sha256": "4" * 64,
    }
    return identity


def _summary_v34():
    summary = _summary_v33()
    identity = summary["artifact_identity"]

    generation = "bearing-dynamic-381"
    identity[
        "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
    ] = {
        "bearing_dynamic_generation": generation,
        **{
            key: generation
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
        },
        "coordinate_order": ["x", "y"],
        "result_coordinate_order": ["x", "y"],
        "displacement_perturbations_m": [[1.0e-4, 0.0], [-1.0e-4, 0.0], [0.0, 1.0e-4], [0.0, -1.0e-4]],
        "result_displacement_perturbations_m": [[1.0e-4, 0.0], [-1.0e-4, 0.0], [0.0, 1.0e-4], [0.0, -1.0e-4]],
        "force_perturbations_n": [[-0.1, -0.005], [0.1, 0.005], [-0.005, -0.09], [0.005, 0.09]],
        "result_force_perturbations_n": [[-0.1, -0.005], [0.1, 0.005], [-0.005, -0.09], [0.005, 0.09]],
        "stiffness_matrix_n_per_m": [[1000.0, 50.0], [50.0, 900.0]],
        "result_stiffness_matrix_n_per_m": [[1000.0, 50.0], [50.0, 900.0]],
        "damping_matrix_n_s_per_m": [[10.0, 2.0], [-2.0, 12.0]],
        "result_damping_matrix_n_s_per_m": [[10.0, 2.0], [-2.0, 12.0]],
        "state_eigenvalues_per_s": [[-5.0, 30.0], [-5.0, -30.0], [-6.0, 28.0], [-6.0, -28.0]],
        "result_state_eigenvalues_per_s": [[-5.0, 30.0], [-5.0, -30.0], [-6.0, 28.0], [-6.0, -28.0]],
        "operating_displacement_m": [0.0, 0.0],
        "result_operating_displacement_m": [0.0, 0.0],
        "operating_velocity_m_s": [0.0, 0.0],
        "result_operating_velocity_m_s": [0.0, 0.0],
        "bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "result_bias_currents_a": [5.0, 5.0, 5.0, 5.0],
        "bearing_mesh_sha256": "1" * 64,
        "result_bearing_mesh_sha256": "1" * 64,
        "bearing_result_owner": "bearing/case-381/linearization",
        "accepted_bearing_result_owner": "bearing/case-381/linearization",
        "bearing_result_sha256": "2" * 64,
        "accepted_bearing_result_sha256": "2" * 64,
    }

    generation = "moving-conductor-381"
    conductivity = 3.5e7
    frequency = 50.0
    skin_depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * 4.0e-7 * math.pi * conductivity))
    identity[
        "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
    ] = {
        "moving_conductor_generation": generation,
        **{
            key: generation
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
        },
        "coordinate_frame": "global_xyz_right_handed",
        "result_coordinate_frame": "global_xyz_right_handed",
        "velocity_m_s": [10.0, 0.0, 0.0],
        "result_velocity_m_s": [10.0, 0.0, 0.0],
        "drag_force_n": [-100.0, 0.0, 0.0],
        "result_drag_force_n": [-100.0, 0.0, 0.0],
        "lift_force_n": [0.0, 20.0, 0.0],
        "result_lift_force_n": [0.0, 20.0, 0.0],
        "joule_power_w": 1000.0,
        "result_joule_power_w": 1000.0,
        "mechanical_drag_power_w": 1000.0,
        "result_mechanical_drag_power_w": 1000.0,
        "conductivity_s_m": conductivity,
        "result_conductivity_s_m": conductivity,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "excitation_frequency_hz": frequency,
        "result_excitation_frequency_hz": frequency,
        "spatial_period_m": 0.2,
        "result_spatial_period_m": 0.2,
        "slip_frequency_hz": frequency,
        "result_slip_frequency_hz": frequency,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "conductor_mesh_sha256": "3" * 64,
        "result_conductor_mesh_sha256": "3" * 64,
        "field_sha256": "4" * 64,
        "accepted_field_sha256": "4" * 64,
        "conductor_result_owner": "moving-conductor/case-381",
        "accepted_conductor_result_owner": "moving-conductor/case-381",
        "conductor_result_sha256": "5" * 64,
        "accepted_conductor_result_sha256": "5" * 64,
    }
    return summary


def _identity_v35():
    identity = _identity_v34()
    generation = "inductance-matrix-235"
    identity[
        "inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_coil_mesh_owner_result_identity"
    ] = {
        "inductance_generation": generation,
        **{key: generation for key in (
            "reciprocity_generation", "psd_generation", "flux_generation",
            "current_generation", "energy_generation", "coil_generation",
            "mesh_generation", "owner_generation", "result_generation")},
        "coil_order": ["phase_a", "phase_b"], "result_coil_order": ["phase_a", "phase_b"],
        "currents_a": [2.0, -1.0], "result_currents_a": [2.0, -1.0],
        "inductance_matrix_h": [[0.008, 0.002], [0.002, 0.006]],
        "result_inductance_matrix_h": [[0.008, 0.002], [0.002, 0.006]],
        "flux_linkages_wb_turn": [0.014, -0.002],
        "result_flux_linkages_wb_turn": [0.014, -0.002],
        "stored_energy_j": 0.015, "result_stored_energy_j": 0.015,
        "minimum_eigenvalue_h": 0.00476393202250021,
        "result_minimum_eigenvalue_h": 0.00476393202250021,
        "coil_mesh_sha256": "1" * 64, "result_coil_mesh_sha256": "1" * 64,
        "inductance_result_owner": "magnetic/coils-235",
        "accepted_inductance_result_owner": "magnetic/coils-235",
        "inductance_result_sha256": "2" * 64,
        "accepted_inductance_result_sha256": "2" * 64,
    }
    generation = "nonlinear-bh-235"
    identity[
        "nonlinear_bh_interpolation_differential_permeability_branch_energy_coenergy_operating_material_solution_identity"
    ] = {
        "bh_generation": generation,
        **{key: generation for key in (
            "interpolation_generation", "differential_generation", "branch_generation",
            "energy_generation", "coenergy_generation", "operating_generation",
            "material_generation", "solution_generation", "result_generation")},
        "b_samples_t": [0.0, 0.5, 1.0, 1.4], "result_b_samples_t": [0.0, 0.5, 1.0, 1.4],
        "h_samples_a_m": [0.0, 100.0, 400.0, 1200.0],
        "result_h_samples_a_m": [0.0, 100.0, 400.0, 1200.0],
        "differential_permeability_h_m": [0.005, 1.0 / 600.0, 0.0005],
        "result_differential_permeability_h_m": [0.005, 1.0 / 600.0, 0.0005],
        "branch": "ascending", "result_branch": "ascending",
        "operating_point_b_t": 1.4, "result_operating_point_b_t": 1.4,
        "operating_point_h_a_m": 1200.0, "result_operating_point_h_a_m": 1200.0,
        "magnetic_energy_density_j_m3": 470.0, "result_magnetic_energy_density_j_m3": 470.0,
        "magnetic_coenergy_density_j_m3": 1210.0,
        "result_magnetic_coenergy_density_j_m3": 1210.0,
        "material_owner": "materials/nonlinear-core-235",
        "accepted_material_owner": "materials/nonlinear-core-235",
        "nonlinear_solution_sha256": "3" * 64,
        "accepted_nonlinear_solution_sha256": "3" * 64,
    }
    return identity


def _summary_v35():
    summary = _summary_v34()
    identity = summary["artifact_identity"]
    generation = "magnetic-gear-401"
    identity[
        "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
    ] = {
        "magnetic_gear_generation": generation,
        **{
            key: generation
            for key in (
                "pole_generation", "harmonic_generation", "torque_generation",
                "phase_generation", "power_generation", "frame_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "high_speed_pole_pairs": 4, "result_high_speed_pole_pairs": 4,
        "low_speed_pole_pairs": 22, "result_low_speed_pole_pairs": 22,
        "modulator_pole_count": 26, "result_modulator_pole_count": 26,
        "transmitted_harmonic_order": 22, "result_transmitted_harmonic_order": 22,
        "high_speed_torque_nm": -10.0, "result_high_speed_torque_nm": -10.0,
        "low_speed_torque_nm": 55.0, "result_low_speed_torque_nm": 55.0,
        "high_speed_angular_velocity_rad_s": 110.0,
        "result_high_speed_angular_velocity_rad_s": 110.0,
        "low_speed_angular_velocity_rad_s": 20.0,
        "result_low_speed_angular_velocity_rad_s": 20.0,
        "high_speed_harmonic_phase_rad": 0.5,
        "result_high_speed_harmonic_phase_rad": 0.5,
        "low_speed_harmonic_phase_rad": -0.1,
        "result_low_speed_harmonic_phase_rad": -0.1,
        "modulator_phase_rad": 0.2, "result_modulator_phase_rad": 0.2,
        "transmitted_phase_rad": 0.8, "result_transmitted_phase_rad": 0.8,
        "coordinate_frame": "global_xyz_right_handed",
        "result_coordinate_frame": "global_xyz_right_handed",
        "gear_mesh_sha256": "1" * 64, "result_gear_mesh_sha256": "1" * 64,
        "gear_result_owner": "magnetic-gear/case-401",
        "accepted_gear_result_owner": "magnetic-gear/case-401",
        "gear_result_sha256": "2" * 64,
        "accepted_gear_result_sha256": "2" * 64,
    }
    generation = "demag-bem-401"
    identity[
        "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
    ] = {
        "demag_bem_generation": generation,
        **{
            key: generation
            for key in (
                "charge_generation", "normal_generation", "jump_generation",
                "farfield_generation", "energy_generation", "mesh_generation",
                "owner_generation", "solution_generation",
            )
        },
        "panel_areas_m2": [1.0, 1.0, 2.0, 2.0],
        "result_panel_areas_m2": [1.0, 1.0, 2.0, 2.0],
        "surface_charge_density_a_m": [2.0, -2.0, 1.0, -1.0],
        "result_surface_charge_density_a_m": [2.0, -2.0, 1.0, -1.0],
        "surface_charge_integral_a_m": 0.0, "result_surface_charge_integral_a_m": 0.0,
        "outward_normals": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        "result_outward_normals": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        "outward_orientation_verified": True, "result_outward_orientation_verified": True,
        "normal_field_jump_a_m": [2.0, -2.0, 1.0, -1.0],
        "result_normal_field_jump_a_m": [2.0, -2.0, 1.0, -1.0],
        "farfield_radius_m": [2.0, 4.0, 8.0], "result_farfield_radius_m": [2.0, 4.0, 8.0],
        "farfield_potential_a": [0.25, 0.0625, 0.015625],
        "result_farfield_potential_a": [0.25, 0.0625, 0.015625],
        "farfield_field_a_m": [0.125, 0.015625, 0.001953125],
        "result_farfield_field_a_m": [0.125, 0.015625, 0.001953125],
        "magnetic_energy_j": 0.75, "result_magnetic_energy_j": 0.75,
        "boundary_mesh_sha256": "3" * 64, "result_boundary_mesh_sha256": "3" * 64,
        "demag_solution_owner": "demag-bem/case-401",
        "accepted_demag_solution_owner": "demag-bem/case-401",
        "demag_solution_sha256": "4" * 64,
        "accepted_demag_solution_sha256": "4" * 64,
    }
    return summary


def _identity_v36():
    identity = _identity_v35()
    generation = "nonlinear-incremental-236"
    identity[
        "nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_current_mesh_owner_solution_identity"
    ] = {
        "nonlinear_generation": generation,
        **{key: generation for key in (
            "branch_generation", "incremental_generation", "energy_generation",
            "coenergy_generation", "inductance_generation", "current_generation",
            "mesh_generation", "owner_generation", "solution_generation", "result_generation")},
        "bh_branch": "ascending", "result_bh_branch": "ascending",
        "current_points_a": [1.0, 2.0, 3.0], "result_current_points_a": [1.0, 2.0, 3.0],
        "flux_linkages_wb_turn": [0.01, 0.018, 0.024],
        "result_flux_linkages_wb_turn": [0.01, 0.018, 0.024],
        "incremental_permeability_h_m": 0.0012, "result_incremental_permeability_h_m": 0.0012,
        "differential_inductance_h": 0.006, "result_differential_inductance_h": 0.006,
        "magnetic_energy_j": 0.03, "result_magnetic_energy_j": 0.03,
        "magnetic_coenergy_j": 0.042, "result_magnetic_coenergy_j": 0.042,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "field_owner": "nonlinear/core-236", "accepted_field_owner": "nonlinear/core-236",
        "solution_sha256": "2" * 64, "accepted_solution_sha256": "2" * 64,
    }
    generation = "weighted-stress-force-236"
    identity[
        "weighted_stress_force_region_air_contour_mesh_convergence_direction_owner_result_identity"
    ] = {
        "force_generation": generation,
        **{key: generation for key in (
            "weighting_generation", "air_generation", "contour_generation",
            "convergence_generation", "direction_generation", "field_generation",
            "owner_generation", "result_generation")},
        "weighting_region_id": "air:force-mask", "result_weighting_region_id": "air:force-mask",
        "air_enclosure_id": "air:enclosure", "result_air_enclosure_id": "air:enclosure",
        "weighted_force_n": [12.0, -0.2], "result_weighted_force_n": [12.0, -0.2],
        "contour_force_samples_n": [[12.0, -0.2], [12.01, -0.19], [11.99, -0.21]],
        "result_contour_force_samples_n": [[12.0, -0.2], [12.01, -0.19], [11.99, -0.21]],
        "mesh_sizes_m": [0.004, 0.002, 0.001], "result_mesh_sizes_m": [0.004, 0.002, 0.001],
        "mesh_force_sequence_n": [[11.5, -0.3], [11.9, -0.22], [12.0, -0.2]],
        "result_mesh_force_sequence_n": [[11.5, -0.3], [11.9, -0.22], [12.0, -0.2]],
        "force_direction_unit": [1.0, 0.0], "result_force_direction_unit": [1.0, 0.0],
        "field_owner": "force/field-236", "accepted_field_owner": "force/field-236",
        "force_result_sha256": "3" * 64, "accepted_force_result_sha256": "3" * 64,
    }
    return identity


def _summary_v36():
    summary = _summary_v35()
    identity = summary["artifact_identity"]
    generation = "bearing-linearization-402"
    identity[
        "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
    ] = {
        "bearing_generation": generation,
        **{
            key: generation
            for key in (
                "bias_generation", "displacement_generation", "force_generation",
                "stiffness_generation", "crosscoupling_generation", "frame_generation",
                "owner_generation", "result_generation",
            )
        },
        "bias_current_a": 5.0, "result_bias_current_a": 5.0,
        "displacement_samples_m": [-0.001, 0.0, 0.001],
        "result_displacement_samples_m": [-0.001, 0.0, 0.001],
        "force_x_samples_n": [10.0, 0.0, -10.0],
        "result_force_x_samples_n": [10.0, 0.0, -10.0],
        "stiffness_matrix_n_m": [[10000.0, 100.0], [100.0, 9000.0]],
        "result_stiffness_matrix_n_m": [[10000.0, 100.0], [100.0, 9000.0]],
        "coordinate_frame": "global_xyz_right_handed",
        "result_coordinate_frame": "global_xyz_right_handed",
        "bearing_owner": "bearing/case-402",
        "accepted_bearing_owner": "bearing/case-402",
        "bearing_result_sha256": "1" * 64,
        "accepted_bearing_result_sha256": "1" * 64,
    }
    generation = "pm-demag-402"
    reference_temperature = 20.0
    operating_temperature = 100.0
    remanence_reference = 1.2
    coefficient = -0.001
    remanence_temperature = remanence_reference * (
        1.0 + coefficient * (operating_temperature - reference_temperature)
    )
    recoil_mu = 1.05
    h_points = [-900000.0, -800000.0, -600000.0]
    b_points = [
        remanence_temperature + 4.0e-7 * math.pi * recoil_mu * field
        for field in h_points
    ]
    identity[
        "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
    ] = {
        "demag_generation": generation,
        **{
            key: generation
            for key in (
                "recoil_generation", "knee_generation", "loadline_generation",
                "temperature_generation", "irreversible_generation", "orientation_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "reference_temperature_c": reference_temperature,
        "result_reference_temperature_c": reference_temperature,
        "operating_temperature_c": operating_temperature,
        "result_operating_temperature_c": operating_temperature,
        "remanence_reference_t": remanence_reference,
        "result_remanence_reference_t": remanence_reference,
        "remanence_temperature_coefficient_per_c": coefficient,
        "result_remanence_temperature_coefficient_per_c": coefficient,
        "temperature_adjusted_remanence_t": remanence_temperature,
        "result_temperature_adjusted_remanence_t": remanence_temperature,
        "recoil_relative_permeability": recoil_mu,
        "result_recoil_relative_permeability": recoil_mu,
        "knee_field_a_m": -800000.0, "result_knee_field_a_m": -800000.0,
        "loadline_h_a_m": h_points, "result_loadline_h_a_m": h_points,
        "loadline_b_t": b_points, "result_loadline_b_t": b_points,
        "knee_crossed": True, "result_knee_crossed": True,
        "irreversible_flux_loss_fraction": 0.05,
        "result_irreversible_flux_loss_fraction": 0.05,
        "field_orientation": "magnetization_antiparallel_h",
        "result_field_orientation": "magnetization_antiparallel_h",
        "mesh_sha256": "2" * 64, "result_mesh_sha256": "2" * 64,
        "demag_owner": "pm/demag-402", "accepted_demag_owner": "pm/demag-402",
        "demag_result_sha256": "3" * 64,
        "accepted_demag_result_sha256": "3" * 64,
    }
    return summary


def _identity_v37():
    identity = _identity_v36()
    generation = "axisym-revolution-246"
    factor = 2.0 * math.pi
    area = 0.01
    radius = 0.05
    identity[
        "axisymmetric_3d_force_revolution_volume_energy_coenergy_direction_displacement_field_mesh_result_identity"
    ] = {
        "revolution_generation": generation,
        **{
            key: generation
            for key in (
                "factor_generation",
                "volume_generation",
                "energy_generation",
                "coenergy_generation",
                "direction_generation",
                "displacement_generation",
                "field_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "revolution_factor": factor,
        "result_revolution_factor": factor,
        "meridional_area_m2": area,
        "result_meridional_area_m2": area,
        "centroid_radius_m": radius,
        "result_centroid_radius_m": radius,
        "swept_volume_m3": factor * radius * area,
        "result_swept_volume_m3": factor * radius * area,
        "axisymmetric_energy_j_per_rad": 2.0,
        "result_axisymmetric_energy_j_per_rad": 2.0,
        "revolved_energy_j": factor * 2.0,
        "result_revolved_energy_j": factor * 2.0,
        "axisymmetric_coenergy_j_per_rad": 2.1,
        "result_axisymmetric_coenergy_j_per_rad": 2.1,
        "revolved_coenergy_j": factor * 2.1,
        "result_revolved_coenergy_j": factor * 2.1,
        "virtual_displacement_m": [0.0, 0.0, 1.0e-4],
        "result_virtual_displacement_m": [0.0, 0.0, 1.0e-4],
        "force_direction_unit": [0.0, 0.0, 1.0],
        "result_force_direction_unit": [0.0, 0.0, 1.0],
        "force_n": [0.0, 0.0, 12.0],
        "result_force_n": [0.0, 0.0, 12.0],
        "field_owner": "axisym/field-246",
        "accepted_field_owner": "axisym/field-246",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "force_result_sha256": "2" * 64,
        "accepted_force_result_sha256": "2" * 64,
    }
    generation = "harmonic-circuit-246"
    voltage = [10.0, 2.0]
    current = [2.0, -1.0]
    denominator = current[0] ** 2 + current[1] ** 2
    impedance = [
        (voltage[0] * current[0] + voltage[1] * current[1]) / denominator,
        (voltage[1] * current[0] - voltage[0] * current[1]) / denominator,
    ]
    power = [
        0.5 * (voltage[0] * current[0] + voltage[1] * current[1]),
        0.5 * (voltage[1] * current[0] - voltage[0] * current[1]),
    ]
    identity[
        "harmonic_circuit_voltage_current_impedance_complex_power_copper_field_loss_rms_owner_result_identity"
    ] = {
        "circuit_generation": generation,
        **{
            key: generation
            for key in (
                "voltage_generation",
                "current_generation",
                "impedance_generation",
                "power_generation",
                "loss_generation",
                "rms_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "voltage_peak_phasor_v": voltage,
        "result_voltage_peak_phasor_v": voltage,
        "current_peak_phasor_a": current,
        "result_current_peak_phasor_a": current,
        "impedance_ohm": impedance,
        "result_impedance_ohm": impedance,
        "complex_power_va": power,
        "result_complex_power_va": power,
        "copper_loss_w": 7.0,
        "result_copper_loss_w": 7.0,
        "field_loss_w": 2.0,
        "result_field_loss_w": 2.0,
        "phasor_convention": "peak_cosine",
        "result_phasor_convention": "peak_cosine",
        "circuit_owner": "circuit:winding-246",
        "accepted_circuit_owner": "circuit:winding-246",
        "circuit_result_sha256": "3" * 64,
        "accepted_circuit_result_sha256": "3" * 64,
    }
    return identity


def _summary_v37():
    summary = _summary_v36()
    identity = summary["artifact_identity"]
    generation = "maglev-dynamic-246"
    frequency, stiffness, damping = 100.0, 10000.0, 10.0
    omega = 2.0 * math.pi * frequency
    displacement = [1.0e-5, 0.0]
    dynamic_stiffness = [stiffness, omega * damping]
    force = [dynamic_stiffness[0] * displacement[0], dynamic_stiffness[1] * displacement[0]]
    identity["maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"] = {
        "maglev_generation": generation,
        **{key: generation for key in (
            "bias_generation", "equilibrium_generation", "frequency_generation",
            "stiffness_generation", "damping_generation", "force_generation",
            "displacement_generation", "phase_generation", "frame_generation",
            "owner_generation", "result_generation")},
        "bias_current_a": 5.0, "result_bias_current_a": 5.0,
        "equilibrium_gap_m": 0.005, "result_equilibrium_gap_m": 0.005,
        "equilibrium_force_n": 100.0, "result_equilibrium_force_n": 100.0,
        "supported_load_n": 100.0, "result_supported_load_n": 100.0,
        "excitation_frequency_hz": frequency, "result_excitation_frequency_hz": frequency,
        "complex_stiffness_n_m": dynamic_stiffness, "result_complex_stiffness_n_m": dynamic_stiffness,
        "viscous_damping_n_s_m": damping, "result_viscous_damping_n_s_m": damping,
        "displacement_phasor_m": displacement, "result_displacement_phasor_m": displacement,
        "force_phasor_n": force, "result_force_phasor_n": force,
        "force_displacement_phase_rad": math.atan2(force[1], force[0]),
        "result_force_displacement_phase_rad": math.atan2(force[1], force[0]),
        "coordinate_frame": "global_z_up_force_positive",
        "result_coordinate_frame": "global_z_up_force_positive",
        "maglev_owner": "maglev/dynamic-246", "accepted_maglev_owner": "maglev/dynamic-246",
        "maglev_result_sha256": "1" * 64, "accepted_maglev_result_sha256": "1" * 64,
    }
    generation = "bem-demag-246"
    identity["bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"] = {
        "demag_generation": generation,
        **{key: generation for key in (
            "reciprocity_generation", "energy_generation", "field_generation",
            "magnetization_generation", "surface_generation", "volume_generation",
            "mesh_generation", "solution_generation", "result_generation")},
        "interaction_energy_12_j": -0.01, "result_interaction_energy_12_j": -0.01,
        "interaction_energy_21_j": -0.01, "result_interaction_energy_21_j": -0.01,
        "field_1_due_2_a_m": [-1000.0, 0.0, 0.0], "result_field_1_due_2_a_m": [-1000.0, 0.0, 0.0],
        "field_2_due_1_a_m": [0.0, -1000.0, 0.0], "result_field_2_due_1_a_m": [0.0, -1000.0, 0.0],
        "magnetization_1_a_m": [800000.0, 0.0, 0.0], "result_magnetization_1_a_m": [800000.0, 0.0, 0.0],
        "magnetization_2_a_m": [0.0, 800000.0, 0.0], "result_magnetization_2_a_m": [0.0, 800000.0, 0.0],
        "surface_orientation": "outward_right_handed", "result_surface_orientation": "outward_right_handed",
        "region_volumes_m3": [1.0e-5, 1.0e-5], "result_region_volumes_m3": [1.0e-5, 1.0e-5],
        "mesh_owner": "mesh/bem-demag-246", "accepted_mesh_owner": "mesh/bem-demag-246",
        "solution_owner": "solution/bem-demag-246", "accepted_solution_owner": "solution/bem-demag-246",
        "demag_result_sha256": "2" * 64, "accepted_demag_result_sha256": "2" * 64,
    }
    return summary


def _identity_v38():
    identity = _identity_v37()
    generation = "harmonic-conductor-258"
    frequency = 10_000.0
    conductivity = 5.8e7
    skin_depth = math.sqrt(
        2.0
        / (
            2.0
            * math.pi
            * frequency
            * (4.0e-7 * math.pi)
            * conductivity
        )
    )
    voltage = [1.0, 0.5]
    current = [10.0, -2.0]
    denominator = current[0] ** 2 + current[1] ** 2
    impedance = [
        (voltage[0] * current[0] + voltage[1] * current[1]) / denominator,
        (voltage[1] * current[0] - voltage[0] * current[1]) / denominator,
    ]
    loss = 0.5 * (
        voltage[0] * current[0] + voltage[1] * current[1]
    )
    identity[
        "harmonic_conductor_skin_proximity_impedance_current_voltage_loss_poynting_frequency_mesh_owner_result_identity"
    ] = {
        "harmonic_generation": generation,
        **{
            key: generation
            for key in (
                "skin_generation", "proximity_generation",
                "impedance_generation", "current_generation",
                "voltage_generation", "loss_generation",
                "poynting_generation", "frequency_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "conductivity_s_m": conductivity,
        "result_conductivity_s_m": conductivity,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "proximity_current_density_peak_a_m2": [1.0e6, 1.2e6, 0.9e6],
        "result_proximity_current_density_peak_a_m2": [1.0e6, 1.2e6, 0.9e6],
        "terminal_voltage_peak_phasor_v": voltage,
        "result_terminal_voltage_peak_phasor_v": voltage,
        "terminal_current_peak_phasor_a": current,
        "result_terminal_current_peak_phasor_a": current,
        "complex_impedance_ohm": impedance,
        "result_complex_impedance_ohm": impedance,
        "copper_loss_w": loss,
        "result_copper_loss_w": loss,
        "inward_poynting_power_w": loss,
        "result_inward_poynting_power_w": loss,
        "phasor_convention": "peak_cosine",
        "result_phasor_convention": "peak_cosine",
        "mesh_levels": [1, 2, 3],
        "result_mesh_levels": [1, 2, 3],
        "impedance_relative_changes": [0.10, 0.02, 0.004],
        "result_impedance_relative_changes": [0.10, 0.02, 0.004],
        "maximum_final_relative_change": 0.01,
        "result_maximum_final_relative_change": 0.01,
        "field_owner": "harmonic:conductor-258",
        "accepted_field_owner": "harmonic:conductor-258",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "harmonic_result_sha256": "2" * 64,
        "accepted_harmonic_result_sha256": "2" * 64,
    }

    generation = "heat-radiation-258"
    sigma = 5.670374419e-8
    coefficient = 10.0
    emissivity = 0.8
    ambient = 300.0
    boundary = 400.0
    convection = coefficient * (boundary - ambient)
    radiation = emissivity * sigma * (boundary**4 - ambient**4)
    conductive = convection + radiation
    area = 0.2 * 0.1
    loss = conductive * area
    identity[
        "heat_convection_radiation_emissivity_ambient_flux_temperature_geometry_mesh_energy_result_identity"
    ] = {
        "heat_generation": generation,
        **{
            key: generation
            for key in (
                "convection_generation", "radiation_generation",
                "emissivity_generation", "ambient_generation",
                "flux_generation", "temperature_generation",
                "geometry_generation", "mesh_generation",
                "energy_generation", "result_generation",
            )
        },
        "convection_coefficient_w_m2_k": coefficient,
        "result_convection_coefficient_w_m2_k": coefficient,
        "emissivity": emissivity,
        "result_emissivity": emissivity,
        "stefan_boltzmann_w_m2_k4": sigma,
        "result_stefan_boltzmann_w_m2_k4": sigma,
        "ambient_temperature_k": ambient,
        "result_ambient_temperature_k": ambient,
        "boundary_temperature_k": boundary,
        "result_boundary_temperature_k": boundary,
        "convection_flux_w_m2": convection,
        "result_convection_flux_w_m2": convection,
        "radiation_flux_w_m2": radiation,
        "result_radiation_flux_w_m2": radiation,
        "conductive_outward_flux_w_m2": conductive,
        "result_conductive_outward_flux_w_m2": conductive,
        "geometry_weighting": "planar_depth",
        "result_geometry_weighting": "planar_depth",
        "boundary_length_m": 0.2,
        "result_boundary_length_m": 0.2,
        "planar_depth_m": 0.1,
        "result_planar_depth_m": 0.1,
        "effective_boundary_area_m2": area,
        "result_effective_boundary_area_m2": area,
        "boundary_heat_loss_w": loss,
        "result_boundary_heat_loss_w": loss,
        "domain_heat_generation_w": loss,
        "result_domain_heat_generation_w": loss,
        "energy_balance_residual_w": 0.0,
        "result_energy_balance_residual_w": 0.0,
        "energy_tolerance_w": 1.0e-9,
        "result_energy_tolerance_w": 1.0e-9,
        "mesh_owner": "heat:mesh-258",
        "accepted_mesh_owner": "heat:mesh-258",
        "heat_result_sha256": "3" * 64,
        "accepted_heat_result_sha256": "3" * 64,
    }
    return identity


def _summary_v38():
    summary = _summary_v37()
    identity = summary["artifact_identity"]
    generation = "eddy-maglev-258"
    velocity = 20.0
    pole_pitch = 0.1
    frequency = velocity / pole_pitch
    conductivity = 3.5e7
    relative_permeability = 1.0
    skin_depth = math.sqrt(
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
    lift = 100.0
    drag = 20.0
    drag_power = drag * velocity
    identity[
        "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
    ] = {
        "maglev_generation": generation,
        **{
            key: generation
            for key in (
                "velocity_generation", "frequency_generation",
                "conductivity_generation", "skin_generation", "force_generation",
                "loss_generation", "power_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )
        },
        "plate_velocity_m_s": velocity,
        "result_plate_velocity_m_s": velocity,
        "pole_pitch_m": pole_pitch,
        "result_pole_pitch_m": pole_pitch,
        "excitation_frequency_hz": frequency,
        "result_excitation_frequency_hz": frequency,
        "plate_conductivity_s_m": conductivity,
        "result_plate_conductivity_s_m": conductivity,
        "relative_permeability": relative_permeability,
        "result_relative_permeability": relative_permeability,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "lift_force_n": lift,
        "result_lift_force_n": lift,
        "drag_force_n": drag,
        "result_drag_force_n": drag,
        "joule_loss_w": drag_power,
        "result_joule_loss_w": drag_power,
        "mechanical_drag_power_w": drag_power,
        "result_mechanical_drag_power_w": drag_power,
        "power_balance_residual_w": 0.0,
        "result_power_balance_residual_w": 0.0,
        "power_tolerance_w": 1.0e-9,
        "result_power_tolerance_w": 1.0e-9,
        "mesh_owner": "mesh:eddy-maglev-258",
        "accepted_mesh_owner": "mesh:eddy-maglev-258",
        "maglev_result_sha256": "1" * 64,
        "accepted_maglev_result_sha256": "1" * 64,
    }

    generation = "pm-coupling-258"
    pole_pairs = 4
    period = 2.0 * math.pi / pole_pairs
    angle = math.pi / 16.0
    delta = 1.0e-4

    def energy(theta: float) -> float:
        return -0.5 * math.cos(pole_pairs * theta)

    energy_minus = energy(angle - delta)
    energy_center = energy(angle)
    energy_plus = energy(angle + delta)
    derivative_torque = -(energy_plus - energy_minus) / (2.0 * delta)
    identity[
        "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
    ] = {
        "coupling_generation": generation,
        **{
            key: generation
            for key in (
                "angle_generation", "periodicity_generation", "energy_generation",
                "derivative_generation", "torque_generation", "reaction_generation",
                "frame_generation", "mesh_generation", "owner_generation",
                "result_generation",
            )
        },
        "relative_angle_rad": angle,
        "result_relative_angle_rad": angle,
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "pole_period_rad": period,
        "result_pole_period_rad": period,
        "angle_perturbation_rad": delta,
        "result_angle_perturbation_rad": delta,
        "energy_minus_j": energy_minus,
        "result_energy_minus_j": energy_minus,
        "energy_center_j": energy_center,
        "result_energy_center_j": energy_center,
        "energy_plus_j": energy_plus,
        "result_energy_plus_j": energy_plus,
        "periodic_energy_j": energy(angle + period),
        "result_periodic_energy_j": energy(angle + period),
        "energy_derivative_torque_nm": derivative_torque,
        "result_energy_derivative_torque_nm": derivative_torque,
        "driver_torque_nm": derivative_torque,
        "result_driver_torque_nm": derivative_torque,
        "driven_torque_nm": -derivative_torque,
        "result_driven_torque_nm": -derivative_torque,
        "torque_frame": "relative_angle_driver_positive",
        "result_torque_frame": "relative_angle_driver_positive",
        "mesh_owner": "mesh:pm-coupling-258",
        "accepted_mesh_owner": "mesh:pm-coupling-258",
        "coupling_result_sha256": "2" * 64,
        "accepted_coupling_result_sha256": "2" * 64,
    }
    return summary


_NONLINEAR = "nonlinear_magnetic_circuit_bh_flux_linkage_coenergy_incremental_inductance_force_power_mesh_result_generation_identity"


_ELECTROSTATIC = "electrostatic_capacitance_matrix_charge_energy_reciprocity_gauge_conductor_mesh_result_generation_identity"


def _identity_v39():
    identity = _identity_v38()
    generation = "nonlinear-magnetic-circuit-271"
    identity[_NONLINEAR] = {
        "magnetic_generation": generation,
        **{key: generation for key in ("bh_generation", "flux_generation", "coenergy_generation", "inductance_generation", "force_generation", "power_generation", "mesh_generation", "result_generation")},
        "bh_curve_h_a_per_m": [0.0, 100.0, 500.0, 1000.0],
        "result_bh_curve_h_a_per_m": [0.0, 100.0, 500.0, 1000.0],
        "bh_curve_b_t": [0.0, 0.5, 1.2, 1.45],
        "result_bh_curve_b_t": [0.0, 0.5, 1.2, 1.45],
        "operating_h_a_per_m": 500.0,
        "result_operating_h_a_per_m": 500.0,
        "operating_b_t": 1.2,
        "result_operating_b_t": 1.2,
        "terminal_current_a": 10.0,
        "result_terminal_current_a": 10.0,
        "flux_linkage_wb_turn": 2.0e-2,
        "result_flux_linkage_wb_turn": 2.0e-2,
        "magnetic_energy_j": 8.0e-2,
        "result_magnetic_energy_j": 8.0e-2,
        "coenergy_j": 1.2e-1,
        "result_coenergy_j": 1.2e-1,
        "current_increment_a": 5.0e-1,
        "result_current_increment_a": 5.0e-1,
        "flux_linkage_increment_wb_turn": 1.0e-3,
        "result_flux_linkage_increment_wb_turn": 1.0e-3,
        "incremental_inductance_h": 2.0e-3,
        "result_incremental_inductance_h": 2.0e-3,
        "virtual_displacement_m": 1.0e-3,
        "result_virtual_displacement_m": 1.0e-3,
        "coenergy_increment_j": 2.0e-2,
        "result_coenergy_increment_j": 2.0e-2,
        "virtual_work_force_n": 20.0,
        "result_virtual_work_force_n": 20.0,
        "terminal_voltage_v": 5.0,
        "result_terminal_voltage_v": 5.0,
        "terminal_power_w": 50.0,
        "result_terminal_power_w": 50.0,
        "mesh_owner": "magnetic:mesh-271",
        "accepted_mesh_owner": "magnetic:mesh-271",
        "nonlinear_result_sha256": "1" * 64,
        "accepted_nonlinear_result_sha256": "1" * 64,
    }
    generation = "electrostatic-capacitance-271"
    identity[_ELECTROSTATIC] = {
        "electrostatic_generation": generation,
        **{key: generation for key in ("matrix_generation", "charge_generation", "energy_generation", "reciprocity_generation", "gauge_generation", "conductor_generation", "mesh_generation", "result_generation")},
        "capacitance_matrix_f": [[2.0e-11, -2.0e-11], [-2.0e-11, 2.0e-11]],
        "result_capacitance_matrix_f": [[2.0e-11, -2.0e-11], [-2.0e-11, 2.0e-11]],
        "terminal_voltage_v": [100.0, 0.0],
        "result_terminal_voltage_v": [100.0, 0.0],
        "terminal_charge_c": [2.0e-9, -2.0e-9],
        "result_terminal_charge_c": [2.0e-9, -2.0e-9],
        "field_energy_j": 1.0e-7,
        "result_field_energy_j": 1.0e-7,
        "reciprocity_residual_f": 0.0,
        "result_reciprocity_residual_f": 0.0,
        "reference_gauge": "conductor:2=0V",
        "result_reference_gauge": "conductor:2=0V",
        "conductor_owner": "electrostatic:conductors-271",
        "accepted_conductor_owner": "electrostatic:conductors-271",
        "mesh_owner": "electrostatic:mesh-271",
        "accepted_mesh_owner": "electrostatic:mesh-271",
        "electrostatic_result_sha256": "2" * 64,
        "accepted_electrostatic_result_sha256": "2" * 64,
    }
    return identity


_THIN_KEY = (
    "thin_conductor_surface_impedance_skin_sheetcurrent_fieldjump_complexpower_"
    "surface_owner_result_identity"
)


_GEAR_KEY = (
    "magnetic_gear_harmonic_polepair_modulation_phase_ratio_torque_"
    "actionreaction_power_owner_result_identity"
)


def _summary_v39() -> dict:
    summary = _summary_v38()
    identity = summary["artifact_identity"]

    generation = "thin-conductor-271"
    frequency = 100_000.0
    conductivity = 5.8e7
    permeability = 4.0e-7 * math.pi
    skin_depth = math.sqrt(
        2.0 / (2.0 * math.pi * frequency * permeability * conductivity)
    )
    surface_resistance = math.sqrt(
        math.pi * frequency * permeability / conductivity
    )
    sheet_current = [100.0, 0.0]
    area = 2.0e-2
    loss = 0.5 * surface_resistance * sheet_current[0] ** 2 * area
    identity[_THIN_KEY] = {
        "thin_generation": generation,
        **{
            key: generation
            for key in (
                "impedance_generation",
                "skin_generation",
                "current_generation",
                "field_generation",
                "power_generation",
                "surface_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "conductivity_s_m": conductivity,
        "result_conductivity_s_m": conductivity,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "surface_impedance_ohm": [surface_resistance, surface_resistance],
        "result_surface_impedance_ohm": [surface_resistance, surface_resistance],
        "sheet_current_peak_a_m": sheet_current,
        "result_sheet_current_peak_a_m": sheet_current,
        "tangential_field_jump_peak_a_m": sheet_current,
        "result_tangential_field_jump_peak_a_m": sheet_current,
        "surface_area_m2": area,
        "result_surface_area_m2": area,
        "joule_loss_w": loss,
        "result_joule_loss_w": loss,
        "reactive_power_var": loss,
        "result_reactive_power_var": loss,
        "surface_owner": "surface:thin-conductor-271",
        "accepted_surface_owner": "surface:thin-conductor-271",
        "thin_result_sha256": "1" * 64,
        "accepted_thin_result_sha256": "1" * 64,
    }

    generation = "magnetic-gear-271"
    high_speed = 10.0
    low_speed = -20.0
    high_torque = 20.0
    low_torque = 10.0
    identity[_GEAR_KEY] = {
        "gear_generation": generation,
        **{
            key: generation
            for key in (
                "harmonic_generation",
                "pole_generation",
                "phase_generation",
                "ratio_generation",
                "torque_generation",
                "reaction_generation",
                "power_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "high_speed_pole_pairs": 4,
        "result_high_speed_pole_pairs": 4,
        "low_speed_pole_pairs": 2,
        "result_low_speed_pole_pairs": 2,
        "modulator_segment_count": 6,
        "result_modulator_segment_count": 6,
        "working_harmonic_order": 6,
        "result_working_harmonic_order": 6,
        "mechanical_phase_rad": math.pi / 12.0,
        "result_mechanical_phase_rad": math.pi / 12.0,
        "high_speed_rad_s": high_speed,
        "result_high_speed_rad_s": high_speed,
        "low_speed_rad_s": low_speed,
        "result_low_speed_rad_s": low_speed,
        "gear_ratio": -2.0,
        "result_gear_ratio": -2.0,
        "high_speed_torque_nm": high_torque,
        "result_high_speed_torque_nm": high_torque,
        "low_speed_torque_nm": low_torque,
        "result_low_speed_torque_nm": low_torque,
        "modulator_reaction_torque_nm": -(high_torque + low_torque),
        "result_modulator_reaction_torque_nm": -(high_torque + low_torque),
        "power_balance_residual_w": high_torque * high_speed + low_torque * low_speed,
        "result_power_balance_residual_w": high_torque * high_speed + low_torque * low_speed,
        "model_owner": "gear:magnetic-271",
        "accepted_model_owner": "gear:magnetic-271",
        "gear_result_sha256": "2" * 64,
        "accepted_gear_result_sha256": "2" * 64,
    }
    return summary


_HEAT = "heat_conduction_contact_convection_flux_temperature_energy_mesh_result_generation_identity"


_CURRENT = "current_flow_electrode_voltage_current_resistance_joule_power_reciprocity_conductor_mesh_result_generation_identity"


def _identity_v40():
    identity = _identity_v39()
    generation = "heat-contact-311"
    area = 1.0e-2
    length = 2.0e-2
    conductivity = 15.0
    contact_resistance_area = 1.0e-4
    convection = 25.0
    hot = 373.15
    ambient = 293.15
    contact_resistance = contact_resistance_area / area
    convection_resistance = 1.0 / (convection * area)
    heat_rate = (hot - ambient) / (
        length / (conductivity * area) + contact_resistance + convection_resistance
    )
    identity[_HEAT] = {
        "heat_generation": generation,
        **{
            key: generation
            for key in (
                "conductivity_generation", "contact_generation", "convection_generation",
                "flux_generation", "temperature_generation", "energy_generation",
                "mesh_generation", "result_generation",
            )
        },
        "conductivity_w_per_m_k": conductivity,
        "result_conductivity_w_per_m_k": conductivity,
        "conduction_length_m": length,
        "result_conduction_length_m": length,
        "boundary_area_m2": area,
        "result_boundary_area_m2": area,
        "contact_resistance_m2_k_per_w": contact_resistance_area,
        "result_contact_resistance_m2_k_per_w": contact_resistance_area,
        "convection_coefficient_w_per_m2_k": convection,
        "result_convection_coefficient_w_per_m2_k": convection,
        "hot_temperature_k": hot,
        "result_hot_temperature_k": hot,
        "ambient_temperature_k": ambient,
        "result_ambient_temperature_k": ambient,
        "boundary_heat_flux_w_per_m2": heat_rate / area,
        "result_boundary_heat_flux_w_per_m2": heat_rate / area,
        "interface_temperature_jump_k": heat_rate * contact_resistance,
        "result_interface_temperature_jump_k": heat_rate * contact_resistance,
        "convection_surface_temperature_k": ambient + heat_rate * convection_resistance,
        "result_convection_surface_temperature_k": ambient + heat_rate * convection_resistance,
        "total_heat_rate_w": heat_rate,
        "result_total_heat_rate_w": heat_rate,
        "energy_balance_residual_w": 0.0,
        "result_energy_balance_residual_w": 0.0,
        "mesh_owner": "heat:mesh-311",
        "accepted_mesh_owner": "heat:mesh-311",
        "heat_result_sha256": "5" * 64,
        "accepted_heat_result_sha256": "5" * 64,
    }
    generation = "current-flow-311"
    identity[_CURRENT] = {
        "current_generation": generation,
        **{
            key: generation
            for key in (
                "electrode_generation", "voltage_generation", "current_generation_id",
                "resistance_generation", "joule_generation", "power_generation",
                "reciprocity_generation", "conductor_generation", "mesh_generation",
                "result_generation",
            )
        },
        "electrode_voltage_v": [10.0, 0.0],
        "result_electrode_voltage_v": [10.0, 0.0],
        "terminal_current_a": [2.0, -2.0],
        "result_terminal_current_a": [2.0, -2.0],
        "effective_resistance_ohm": 5.0,
        "result_effective_resistance_ohm": 5.0,
        "joule_loss_w": 20.0,
        "result_joule_loss_w": 20.0,
        "terminal_power_w": 20.0,
        "result_terminal_power_w": 20.0,
        "reciprocity_residual_ohm": 0.0,
        "result_reciprocity_residual_ohm": 0.0,
        "conductor_owner": "current:conductors-311",
        "accepted_conductor_owner": "current:conductors-311",
        "mesh_owner": "current:mesh-311",
        "accepted_mesh_owner": "current:mesh-311",
        "current_result_sha256": "6" * 64,
        "accepted_current_result_sha256": "6" * 64,
    }
    return identity


_SHIELD = "multilayer_magnetic_shield_permeability_thickness_radius_interface_flux_attenuation_leakage_energy_geometry_result_identity"


_TRANSFORMER = "transformer_leakage_mutual_inductance_fluxlinkage_reciprocity_psd_coenergy_force_winding_result_identity"


def _summary_v40() -> dict:
    summary = _summary_v39()
    identity = summary["artifact_identity"]
    generation = "multilayer-shield-280"
    permeability = [2000.0, 5000.0]
    thickness = [1.0e-3, 0.8e-3]
    radii = [0.12, 0.10]
    factors = [1.0 + (mu_r - 1.0) * layer / (2.0 * radius) for mu_r, layer, radius in zip(permeability, thickness, radii)]
    attenuation = math.prod(factors)
    external = 1.0e-3
    leakage = external / attenuation
    volume = 4.0e-3
    energy = leakage**2 * volume / (2.0 * 4.0e-7 * math.pi)
    values = {
        "relative_permeability": permeability, "layer_thickness_m": thickness,
        "layer_mean_radius_m": radii, "layer_shielding_factor": factors,
        "interface_normal_flux_t": [leakage, leakage, leakage],
        "external_field_t": external, "attenuation_factor": attenuation,
        "leakage_field_t": leakage, "cavity_volume_m3": volume,
        "stored_energy_j": energy,
    }
    identity[_SHIELD] = {
        "shield_generation": generation,
        **{key: generation for key in ("material_generation", "thickness_generation", "geometry_generation", "flux_generation", "attenuation_generation", "field_generation", "energy_generation", "owner_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "geometry_owner": "geometry:multilayer-shield-280",
        "accepted_geometry_owner": "geometry:multilayer-shield-280",
        "shield_result_sha256": "1" * 64,
        "accepted_shield_result_sha256": "1" * 64,
    }

    generation = "transformer-coupling-280"
    primary = 12.0e-3
    secondary = 3.0e-3
    mutual = 5.4e-3
    currents = [4.0, -6.0]
    linkages = [primary * currents[0] + mutual * currents[1], mutual * currents[0] + secondary * currents[1]]
    coenergy = 0.5 * sum(current * linkage for current, linkage in zip(currents, linkages))
    gradient = -2.0e-2
    values = {
        "inductance_matrix_h": [[primary, mutual], [mutual, secondary]],
        "primary_leakage_inductance_h": primary - mutual**2 / secondary,
        "winding_currents_a": currents, "flux_linkages_wb_turn": linkages,
        "reciprocity_residual_h": 0.0, "coenergy_j": coenergy,
        "mutual_inductance_gradient_h_per_m": gradient,
        "force_n": gradient * currents[0] * currents[1],
    }
    identity[_TRANSFORMER] = {
        "transformer_generation": generation,
        **{key: generation for key in ("inductance_generation", "leakage_generation", "flux_generation", "reciprocity_generation", "energy_generation", "force_generation", "winding_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "winding_owner": "winding:transformer-280",
        "accepted_winding_owner": "winding:transformer-280",
        "transformer_result_sha256": "2" * 64,
        "accepted_transformer_result_sha256": "2" * 64,
    }
    return summary


_MAGNET = (
    "permanent_magnet_recoil_loadline_operating_demag_energy_virtualwork_"
    "force_mesh_result_generation_identity"
)


_CAPACITANCE = (
    "electrostatic_capacitance_matrix_symmetry_psd_charge_energy_reciprocity_"
    "mesh_result_generation_identity"
)


def _identity_v41():
    identity = _identity_v40()
    generation = "pm-loadline-724"
    mu0 = 4.0e-7 * math.pi
    recoil_mu_r = 1.05
    remanence = 1.2
    coercive_field = remanence / (mu0 * recoil_mu_r)
    permeance = 2.0
    operating_h = -remanence / (mu0 * (recoil_mu_r + permeance))
    operating_b = remanence + mu0 * recoil_mu_r * operating_h
    demag_knee_h = -8.0e5
    demag_margin = abs(demag_knee_h) - abs(operating_h)
    volume = 1.0e-6
    energy = 0.5 * operating_b * operating_b * volume / mu0
    displacement = 1.0e-4
    energy_minus = energy - 1.0e-3
    energy_plus = energy + 1.0e-3
    force = -(energy_plus - energy_minus) / (2.0 * displacement)
    identity[_MAGNET] = {
        "magnet_generation": generation,
        **{
            key: generation
            for key in (
                "recoil_generation",
                "loadline_generation",
                "operating_generation",
                "demag_generation",
                "energy_generation",
                "force_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "recoil_relative_permeability": recoil_mu_r,
        "result_recoil_relative_permeability": recoil_mu_r,
        "remanence_t": remanence,
        "result_remanence_t": remanence,
        "coercive_field_a_per_m": coercive_field,
        "result_coercive_field_a_per_m": coercive_field,
        "loadline_permeance_coefficient": permeance,
        "result_loadline_permeance_coefficient": permeance,
        "operating_h_a_per_m": operating_h,
        "result_operating_h_a_per_m": operating_h,
        "operating_b_t": operating_b,
        "result_operating_b_t": operating_b,
        "demag_knee_h_a_per_m": demag_knee_h,
        "result_demag_knee_h_a_per_m": demag_knee_h,
        "demag_margin_a_per_m": demag_margin,
        "result_demag_margin_a_per_m": demag_margin,
        "magnet_volume_m3": volume,
        "result_magnet_volume_m3": volume,
        "field_energy_j": energy,
        "result_field_energy_j": energy,
        "virtual_work_displacement_m": displacement,
        "result_virtual_work_displacement_m": displacement,
        "energy_minus_j": energy_minus,
        "result_energy_minus_j": energy_minus,
        "energy_plus_j": energy_plus,
        "result_energy_plus_j": energy_plus,
        "virtual_work_force_n": force,
        "result_virtual_work_force_n": force,
        "mesh_owner": "magnetics:pm-mesh-724",
        "accepted_mesh_owner": "magnetics:pm-mesh-724",
        "magnet_result_sha256": "5" * 64,
        "accepted_magnet_result_sha256": "5" * 64,
    }

    generation = "capacitance-matrix-724"
    matrix = [[2.0e-9, -0.5e-9], [-0.5e-9, 1.5e-9]]
    voltages = [10.0, 0.0]
    charges = [
        sum(matrix[row][column] * voltages[column] for column in range(2))
        for row in range(2)
    ]
    energy = 0.5 * sum(voltages[index] * charges[index] for index in range(2))
    identity[_CAPACITANCE] = {
        "capacitance_generation": generation,
        **{
            key: generation
            for key in (
                "matrix_generation",
                "symmetry_generation",
                "psd_generation",
                "charge_generation",
                "energy_generation",
                "reciprocity_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "conductor_names": ["electrode_1", "electrode_2"],
        "result_conductor_names": ["electrode_1", "electrode_2"],
        "capacitance_matrix_f": matrix,
        "result_capacitance_matrix_f": matrix,
        "drive_voltage_v": voltages,
        "result_drive_voltage_v": voltages,
        "conductor_charge_c": charges,
        "result_conductor_charge_c": charges,
        "stored_energy_j": energy,
        "result_stored_energy_j": energy,
        "symmetry_residual_f": 0.0,
        "result_symmetry_residual_f": 0.0,
        "reciprocity_residual_f": 0.0,
        "result_reciprocity_residual_f": 0.0,
        "mesh_owner": "electrostatics:cap-mesh-724",
        "accepted_mesh_owner": "electrostatics:cap-mesh-724",
        "capacitance_result_sha256": "6" * 64,
        "accepted_capacitance_result_sha256": "6" * 64,
    }
    return identity


_MAGLEV = "maglev_equilibrium_airgap_position_force_gradient_stiffness_potential_energy_stability_geometry_result_identity"


_EDDY = "eddy_shield_frequency_conductivity_permeability_skin_depth_thickness_phase_attenuation_loss_energy_geometry_result_identity"


def _summary_v41() -> dict:
    summary = _summary_v40()
    identity = summary["artifact_identity"]
    generation = "maglev-equilibrium-724"
    positions = [-0.01, 0.0, 0.01]
    stiffness = 200.0
    values = {
        "air_gap_m": 0.02,
        "sample_position_m": positions,
        "force_n": [-stiffness * position for position in positions],
        "equilibrium_position_m": 0.0,
        "equilibrium_force_n": 0.0,
        "force_gradient_n_per_m": -stiffness,
        "stiffness_n_per_m": stiffness,
        "potential_energy_j": [0.5 * stiffness * position**2 for position in positions],
        "energy_curvature_j_per_m2": stiffness,
        "stability": "stable",
    }
    identity[_MAGLEV] = {
        "maglev_generation": generation,
        **{key: generation for key in ("gap_generation", "position_generation", "force_generation", "gradient_generation", "stiffness_generation", "energy_generation", "stability_generation", "geometry_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "geometry_owner": "geometry:maglev-equilibrium-724",
        "accepted_geometry_owner": "geometry:maglev-equilibrium-724",
        "maglev_result_sha256": "5" * 64,
        "accepted_maglev_result_sha256": "5" * 64,
    }

    generation = "eddy-shield-724"
    mu0 = 4.0e-7 * math.pi
    frequency = 1000.0
    omega = 2.0 * math.pi * frequency
    mu_r = 1.0
    conductivity = 5.8e7
    permeability = mu0 * mu_r
    skin_depth = math.sqrt(2.0 / (omega * permeability * conductivity))
    thickness = 2.0e-3
    attenuation = math.exp(-thickness / skin_depth)
    incident = 1.0e-3
    area = 0.05
    surface_resistance = 1.0 / (conductivity * skin_depth)
    transmitted = incident * attenuation
    volume = area * thickness
    values = {
        "frequency_hz": frequency,
        "angular_frequency_rad_s": omega,
        "relative_permeability": mu_r,
        "conductivity_s_per_m": conductivity,
        "skin_depth_m": skin_depth,
        "shield_thickness_m": thickness,
        "attenuation_factor": attenuation,
        "phase_lag_rad": -thickness / skin_depth,
        "incident_field_t": incident,
        "transmitted_field_t": transmitted,
        "surface_resistance_ohm": surface_resistance,
        "shield_area_m2": area,
        "eddy_loss_w": 0.5 * surface_resistance * (incident / permeability) ** 2 * area,
        "shield_volume_m3": volume,
        "stored_energy_j": transmitted**2 * volume / (2.0 * permeability),
    }
    identity[_EDDY] = {
        "eddy_shield_generation": generation,
        **{key: generation for key in ("frequency_generation", "material_generation", "skin_depth_generation", "thickness_generation", "phase_generation", "field_generation", "loss_generation", "energy_generation", "geometry_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "geometry_owner": "geometry:eddy-shield-724",
        "accepted_geometry_owner": "geometry:eddy-shield-724",
        "eddy_shield_result_sha256": "6" * 64,
        "accepted_eddy_shield_result_sha256": "6" * 64,
    }
    return summary


_FORCE = (
    "weightedstress_airgapcontour_virtualwork_force_direction_mesh_"
    "fieldowner_result_generation_identity"
)


_CONDUCTOR = (
    "harmonicconductor_skin_depth_complexcurrent_jouleloss_impedance_"
    "power_mesh_result_generation_identity"
)


def _identity_v42():
    identity = _identity_v41()
    generation = "force-closure-842"
    force = [12.0, -3.0]
    contours = [[12.02, -2.99], [11.98, -3.01], [12.01, -3.0]]
    force_norm = math.hypot(*force)
    spread = max(
        math.hypot(sample[0] - force[0], sample[1] - force[1])
        for sample in contours
    ) / force_norm
    displacement = 1.0e-4
    energy_minus = 1.0012
    energy_plus = 0.9988
    virtual_force = -(energy_plus - energy_minus) / (2.0 * displacement)
    mesh_forces = [[11.9, -2.95], force]
    identity[_FORCE] = {
        "force_generation": generation,
        **{
            key: generation
            for key in (
                "weighted_stress_generation",
                "airgap_contour_generation",
                "virtual_work_generation",
                "direction_generation",
                "mesh_generation",
                "field_generation",
                "result_generation",
            )
        },
        "weighted_stress_force_n": force,
        "result_weighted_stress_force_n": force,
        "airgap_contour_forces_n": contours,
        "result_airgap_contour_forces_n": contours,
        "contour_independence_relative_spread": spread,
        "result_contour_independence_relative_spread": spread,
        "virtual_work_displacement_m": displacement,
        "result_virtual_work_displacement_m": displacement,
        "energy_minus_j": energy_minus,
        "result_energy_minus_j": energy_minus,
        "energy_plus_j": energy_plus,
        "result_energy_plus_j": energy_plus,
        "virtual_work_direction": [1.0, 0.0],
        "result_virtual_work_direction": [1.0, 0.0],
        "virtual_work_force_n": virtual_force,
        "result_virtual_work_force_n": virtual_force,
        "mesh_refinement_force_samples_n": mesh_forces,
        "result_mesh_refinement_force_samples_n": mesh_forces,
        "field_owner": "field:force-closure-842",
        "accepted_field_owner": "field:force-closure-842",
        "mesh_owner": "mesh:force-closure-842",
        "accepted_mesh_owner": "mesh:force-closure-842",
        "force_result_sha256": "1" * 64,
        "accepted_force_result_sha256": "1" * 64,
    }

    generation = "harmonic-conductor-842"
    frequency = 1000.0
    conductivity = 5.8e7
    permeability = 4.0e-7 * math.pi
    skin_depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * permeability * conductivity))
    current_density = [[1.0e6, -2.0e5], [0.8e6, -1.5e5]]
    integrated_j_squared = 2.0 * conductivity * 20.8
    identity[_CONDUCTOR] = {
        "conductor_generation": generation,
        **{
            key: generation
            for key in (
                "skin_depth_generation",
                "current_density_generation",
                "joule_generation",
                "impedance_generation",
                "power_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "conductivity_s_per_m": conductivity,
        "result_conductivity_s_per_m": conductivity,
        "absolute_permeability_h_per_m": permeability,
        "result_absolute_permeability_h_per_m": permeability,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "complex_current_density_a_per_m2": current_density,
        "result_complex_current_density_a_per_m2": current_density,
        "integrated_abs_current_density_sq_a2_per_m": integrated_j_squared,
        "result_integrated_abs_current_density_sq_a2_per_m": integrated_j_squared,
        "joule_loss_w": 20.8,
        "result_joule_loss_w": 20.8,
        "terminal_current_a": [10.0, -2.0],
        "result_terminal_current_a": [10.0, -2.0],
        "terminal_voltage_v": [4.4, 1.2],
        "result_terminal_voltage_v": [4.4, 1.2],
        "terminal_impedance_ohm": [0.4, 0.2],
        "result_terminal_impedance_ohm": [0.4, 0.2],
        "complex_power_va": [20.8, 10.4],
        "result_complex_power_va": [20.8, 10.4],
        "mesh_owner": "mesh:harmonic-conductor-842",
        "accepted_mesh_owner": "mesh:harmonic-conductor-842",
        "conductor_result_sha256": "2" * 64,
        "accepted_conductor_result_sha256": "2" * 64,
    }
    return identity


_DEMAG = (
    "demagnetizing_tensor_symmetry_trace_eigenvalue_reciprocity_energy_mesh_"
    "magnetization_result_identity"
)


_LINEAR = (
    "linear_motor_cogging_force_position_periodicity_work_coenergy_phase_"
    "thrust_mesh_result_identity"
)


def _matvec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(left * right for left, right in zip(row, vector)) for row in matrix]


def _dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _summary_v42() -> dict:
    summary = _summary_v41()
    identity = summary["artifact_identity"]

    generation = "demag-tensor-842"
    tensor = [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.5]]
    magnetization = [1.0e5, 2.0e5, 3.0e5]
    probe = [-2.0e5, 1.0e5, 0.5e5]
    tensor_m = _matvec(tensor, magnetization)
    tensor_probe = _matvec(tensor, probe)
    volume = 1.0e-6
    values = {
        "demag_tensor": tensor,
        "tensor_trace": 1.0,
        "tensor_eigenvalues": [0.2, 0.3, 0.5],
        "magnetization_a_per_m": magnetization,
        "probe_magnetization_a_per_m": probe,
        "demag_field_a_per_m": [-value for value in tensor_m],
        "reciprocity_left": _dot(magnetization, tensor_probe),
        "reciprocity_right": _dot(probe, tensor_m),
        "magnet_volume_m3": volume,
        "demag_energy_j": (
            0.5
            * 4.0e-7
            * math.pi
            * volume
            * _dot(magnetization, tensor_m)
        ),
    }
    identity[_DEMAG] = {
        "demag_tensor_generation": generation,
        **{
            key: generation
            for key in (
                "symmetry_generation", "trace_generation", "eigenvalue_generation",
                "reciprocity_generation", "energy_generation", "mesh_generation",
                "magnetization_generation", "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:demag-tensor-842",
        "accepted_mesh_owner": "mesh:demag-tensor-842",
        "magnetization_owner": "magnetization:demag-tensor-842",
        "accepted_magnetization_owner": "magnetization:demag-tensor-842",
        "demag_result_sha256": "1" * 64,
        "accepted_demag_result_sha256": "1" * 64,
    }

    generation = "linear-motor-period-842"
    period = 0.04
    positions = [0.0, 0.01, 0.02, 0.03, 0.04]
    wave_number = 2.0 * math.pi / period
    amplitude = 10.0 / wave_number
    current_peak = 5.0
    force = [
        amplitude * wave_number * math.sin(wave_number * position)
        for position in positions
    ]
    phase_currents = [
        [
            current_peak * math.sin(wave_number * position),
            current_peak * math.sin(wave_number * position - 2.0 * math.pi / 3.0),
            current_peak * math.sin(wave_number * position + 2.0 * math.pi / 3.0),
        ]
        for position in positions
    ]
    values = {
        "position_m": positions,
        "period_m": period,
        "phase_order": ["U", "V", "W"],
        "phase_currents_a": phase_currents,
        "current_peak_a": current_peak,
        "coenergy_amplitude_j": amplitude,
        "coenergy_j": [-amplitude * math.cos(wave_number * position) for position in positions],
        "cogging_force_n": force,
        "base_thrust_n": 100.0,
        "thrust_n": [100.0 + value for value in force],
        "periodic_work_j": 0.0,
    }
    identity[_LINEAR] = {
        "linear_motor_generation": generation,
        **{
            key: generation
            for key in (
                "position_generation", "periodicity_generation", "phase_generation",
                "force_generation", "thrust_generation", "work_generation",
                "coenergy_generation", "mesh_generation", "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:linear-motor-period-842",
        "accepted_mesh_owner": "mesh:linear-motor-period-842",
        "linear_motor_result_sha256": "2" * 64,
        "accepted_linear_motor_result_sha256": "2" * 64,
    }
    return summary


_SOLENOID = "axisymmetric_solenoid_flux_inductance_force_coenergy_axisfactor_mesh_generation_identity"


_DIELECTRIC = "dielectric_interface_capacitance_charge_flux_energy_reciprocity_mesh_generation_identity"


def _identity_v43():
    identity = _identity_v42()
    generation = "axisym-solenoid-843"
    identity[_SOLENOID] = {
        "solenoid_generation": generation,
        **{key: generation for key in ("flux_generation", "inductance_generation", "force_generation", "coenergy_generation", "axis_factor_generation", "mesh_generation", "result_generation")},
        "current_a": 3.0, "result_current_a": 3.0, "mean_radius_m": 0.02, "result_mean_radius_m": 0.02,
        "flux_linkage_wb_turn": 0.12, "result_flux_linkage_wb_turn": 0.12, "inductance_h": 0.04, "result_inductance_h": 0.04,
        "axial_force_n": 12.0, "result_axial_force_n": 12.0, "coenergy_force_derivative_n": 12.0, "result_coenergy_force_derivative_n": 12.0,
        "two_pi_r_factor_m": 2.0 * math.pi * 0.02, "result_two_pi_r_factor_m": 2.0 * math.pi * 0.02,
        "mesh_owner": "mesh:axisym-solenoid-843", "result_mesh_owner": "mesh:axisym-solenoid-843",
        "solenoid_result_sha256": "5" * 64, "accepted_solenoid_result_sha256": "5" * 64,
    }
    generation = "dielectric-interface-843"
    identity[_DIELECTRIC] = {
        "dielectric_generation": generation,
        **{key: generation for key in ("capacitance_generation", "charge_generation", "flux_generation", "energy_generation", "interface_generation", "mesh_generation", "result_generation")},
        "voltage_v": 100.0, "result_voltage_v": 100.0, "capacitance_f": 1.0e-9, "result_capacitance_f": 1.0e-9,
        "conductor_charge_c": 1.0e-7, "result_conductor_charge_c": 1.0e-7, "normal_displacement_flux_c": 1.0e-7, "result_normal_displacement_flux_c": 1.0e-7,
        "stored_energy_j": 5.0e-6, "result_stored_energy_j": 5.0e-6, "interface_continuity_residual_c": 1.0e-12, "result_interface_continuity_residual_c": 1.0e-12,
        "reciprocity_residual_f": 1.0e-12, "result_reciprocity_residual_f": 1.0e-12,
        "mesh_owner": "mesh:dielectric-interface-843", "result_mesh_owner": "mesh:dielectric-interface-843",
        "dielectric_result_sha256": "6" * 64, "accepted_dielectric_result_sha256": "6" * 64,
    }
    return identity


_BEARING = "magneticbearing_stiffnessmatrix_symmetry_crosscoupling_force_energy_stability_mesh_result_identity"


_HYSTERESIS = "hysteresis_minorloop_remanence_coercivity_loss_path_energy_material_mesh_result_identity"


def _summary_v43() -> dict:
    summary = _summary_v42()
    generation = "bearing-stiffness-843"
    summary["artifact_identity"][_BEARING] = {
        "bearing_generation": generation,
        **{key: generation for key in ("equilibrium_generation", "stiffness_generation", "symmetry_generation", "energy_generation", "stability_generation", "gap_generation", "mesh_generation", "result_generation")},
        "force_equilibrium_n": [0.0, 0.0], "result_force_equilibrium_n": [0.0, 0.0],
        "stiffness_matrix_n_per_m": [[100.0, -5.0], [-5.0, 80.0]],
        "result_stiffness_matrix_n_per_m": [[100.0, -5.0], [-5.0, 80.0]],
        "energy_curvature_n_per_m": 79.0, "result_energy_curvature_n_per_m": 79.0,
        "stability_sign": "stable", "result_stability_sign": "stable",
        "gap_m": 1.0e-3, "result_gap_m": 1.0e-3,
        "mesh_owner": "mesh:bearing-stiffness-843", "result_mesh_owner": "mesh:bearing-stiffness-843",
        "bearing_result_sha256": "5" * 64, "accepted_bearing_result_sha256": "5" * 64,
    }
    generation = "hyst-minor-843"
    summary["artifact_identity"][_HYSTERESIS] = {
        "hysteresis_generation": generation,
        **{key: generation for key in ("path_generation", "branch_generation", "remanence_generation", "coercivity_generation", "loss_generation", "energy_generation", "material_generation", "mesh_generation", "result_generation")},
        "field_path_a_per_m": [-1.0, 1.0, -1.0], "result_field_path_a_per_m": [-1.0, 1.0, -1.0],
        "magnetization_a_per_m": [-0.8, 1.0, -0.8], "result_magnetization_a_per_m": [-0.8, 1.0, -0.8],
        "remanence_a_per_m": 0.8, "result_remanence_a_per_m": 0.8,
        "coercivity_a_per_m": 0.4, "result_coercivity_a_per_m": 0.4,
        "loop_area_loss_j_per_m3": 0.16, "result_loop_area_loss_j_per_m3": 0.16,
        "cycle_energy_j": 0.16, "result_cycle_energy_j": 0.16,
        "material_owner": "material:hyst-minor-843", "result_material_owner": "material:hyst-minor-843",
        "mesh_owner": "mesh:hyst-minor-843", "result_mesh_owner": "mesh:hyst-minor-843",
        "hysteresis_result_sha256": "6" * 64, "accepted_hysteresis_result_sha256": "6" * 64,
    }
    return summary


_DYNAMIC = "magneticbearing_dynamicstiffness_phase_damping_force_power_stability_mesh_result_identity"


_DEMAG_V44 = "demag_minorloop_fieldpath_remanence_loss_energy_temperature_material_mesh_result_identity"


def _identity_v44() -> dict:
    return {
        _DYNAMIC: {
            "bearing_dynamic_generation": "bearing-dynamic-844",
            **{key: "bearing-dynamic-844" for key in ("dynamic_stiffness_generation", "phase_generation", "damping_generation", "force_generation", "power_generation", "stability_generation", "mesh_generation", "result_generation")},
            "frequency_hz": [100.0, 200.0, 300.0], "result_frequency_hz": [100.0, 200.0, 300.0],
            "dynamic_stiffness_n_per_m": [100.0, 110.0, 120.0], "result_dynamic_stiffness_n_per_m": [100.0, 110.0, 120.0],
            "phase_deg": [0.0, 10.0, 20.0], "result_phase_deg": [0.0, 10.0, 20.0],
            "damping_n_s_per_m": [2.0, 2.5, 3.0], "result_damping_n_s_per_m": [2.0, 2.5, 3.0],
            "force_n": [10.0, 11.0, 12.0], "result_force_n": [10.0, 11.0, 12.0],
            "power_w": [1.0, 1.2, 1.4], "result_power_w": [1.0, 1.2, 1.4],
            "stability_sign": "stable", "result_stability_sign": "stable",
            "mesh_owner": "mesh:bearing-dynamic-844", "result_mesh_owner": "mesh:bearing-dynamic-844",
            "bearing_dynamic_result_sha256": "9" * 64, "accepted_bearing_dynamic_result_sha256": "9" * 64,
        },
        _DEMAG_V44: {
            "demag_minor_generation": "demag-minor-844",
            **{key: "demag-minor-844" for key in ("fieldpath_generation", "remanence_generation", "branch_generation", "loss_generation", "energy_generation", "temperature_generation", "material_generation", "mesh_generation", "result_generation")},
            "field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0], "result_field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0],
            "magnetization_a_per_m": [-0.8, 0.4, -0.3, 0.9], "result_magnetization_a_per_m": [-0.8, 0.4, -0.3, 0.9],
            "remanence_a_per_m": 0.8, "result_remanence_a_per_m": 0.8,
            "coercivity_a_per_m": 0.4, "result_coercivity_a_per_m": 0.4,
            "loss_energy_j": 0.16, "result_loss_energy_j": 0.16,
            "temperature_k": 293.15, "result_temperature_k": 293.15,
            "material_owner": "material:demag-minor-844", "result_material_owner": "material:demag-minor-844",
            "mesh_owner": "mesh:demag-minor-844", "result_mesh_owner": "mesh:demag-minor-844",
            "demag_minor_result_sha256": "a" * 64, "accepted_demag_minor_result_sha256": "a" * 64,
        },
    }


def _identity_v47() -> dict[str, object]:
    motor_generation = "motor-v47"
    force_generation = "force-v47"
    excitation = {"phase_order": ["A", "B", "C"], "current_a": [10.0, -5.0, -5.0]}
    components = {"core": 80.0, "magnet": 20.0}
    return {
        MOTOR: {
            "generation": motor_generation,
            **{
                key: motor_generation
                for key in (
                    "geometry_generation",
                    "material_generation",
                    "excitation_generation",
                    "operating_point_generation",
                    "lane_a_generation",
                    "lane_b_generation",
                    "result_generation",
                )
            },
            "lane_ids": MOTOR_LANES,
            "result_lane_ids": MOTOR_LANES,
            "geometry_identity_sha256": "1" * 64,
            "lane_a_geometry_identity_sha256": "1" * 64,
            "lane_b_geometry_identity_sha256": "1" * 64,
            "material_identity": "material:motor-v47",
            "lane_a_material_identity": "material:motor-v47",
            "lane_b_material_identity": "material:motor-v47",
            "excitation_identity": excitation,
            "lane_a_excitation_identity": excitation,
            "lane_b_excitation_identity": excitation,
            "operating_point_key": "speed=3000rpm,current=10A",
            "lane_a_operating_point_key": "speed=3000rpm,current=10A",
            "lane_b_operating_point_key": "speed=3000rpm,current=10A",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        FORCE: {
            "generation": force_generation,
            **{
                key: force_generation
                for key in (
                    "coenergy_generation",
                    "displacement_generation",
                    "body_owner_generation",
                    "aggregation_generation",
                    "result_generation",
                )
            },
            "displacement_pair_m": [0.0, 0.001],
            "result_displacement_pair_m": [0.0, 0.001],
            "coenergy_pair_j": [1.0, 1.1],
            "result_coenergy_pair_j": [1.0, 1.1],
            "body_owner": "body:moving-assembly",
            "result_body_owner": "body:moving-assembly",
            "component_force_n": components,
            "result_component_force_n": components,
            "aggregated_force_n": 100.0,
            "result_aggregated_force_n": 100.0,
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def _identity_v48() -> dict[str, object]:
    bem_generation = "bem-panel-v48-901"
    hysteresis_generation = "minor-loop-v48-901"
    panels = ["panel:1", "panel:2", "panel:3", "panel:4"]
    normals = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]]
    angles = [math.pi, math.pi, math.pi, math.pi]
    points = [[0.2, 120.0], [0.8, 480.0], [0.3, 160.0]]
    state = {"branch": "descending", "last_reversal": 0.8, "memory_depth": 2}
    return {
        V48_BEM: {
            "generation": bem_generation,
            "quadrature_generation": bem_generation,
            "normal_generation": bem_generation,
            "solid_angle_generation": bem_generation,
            "mesh_generation": bem_generation,
            "result_generation": bem_generation,
            "panel_ids": panels,
            "result_panel_ids": panels,
            "near_quadrature_order": [8, 8, 8, 8],
            "result_near_quadrature_order": [8, 8, 8, 8],
            "far_quadrature_order": [4, 4, 4, 4],
            "result_far_quadrature_order": [4, 4, 4, 4],
            "panel_normals": normals,
            "result_panel_normals": normals,
            "panel_solid_angles_sr": angles,
            "result_panel_solid_angles_sr": angles,
            "mesh_revision": "mesh:bem-v48-901",
            "result_mesh_revision": "mesh:bem-v48-901",
            "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
        },
        HYSTERESIS: {
            "generation": hysteresis_generation,
            "return_point_generation": hysteresis_generation,
            "state_generation": hysteresis_generation,
            "environment_generation": hysteresis_generation,
            "material_generation": hysteresis_generation,
            "result_generation": hysteresis_generation,
            "return_points": points,
            "result_return_points": points,
            "internal_state": state,
            "result_internal_state": state,
            "temperature_k": 353.15,
            "result_temperature_k": 353.15,
            "frequency_hz": 50.0,
            "result_frequency_hz": 50.0,
            "material_owner": "material:hysteresis-v48-901",
            "result_material_owner": "material:hysteresis-v48-901",
            "result_sha256": "7" * 64,
            "accepted_result_sha256": "7" * 64,
        },
    }


def _identity_v49() -> dict[str, object]:
    demag_generation = "demag-recoil-v49-901"
    virtual_work_generation = "virtual-work-v49-901"
    recoil = [[0.9, -120000.0], [1.0, -80000.0], [1.1, -40000.0]]
    displacement = [-0.0001, 0.0, 0.0001]
    energy = [1.00012, 1.00000, 0.99988]
    return {
        V49_DEMAG: {
            "generation": demag_generation,
            "branch_generation": demag_generation,
            "temperature_generation": demag_generation,
            "loadstep_generation": demag_generation,
            "result_generation": demag_generation,
            "recoil_branch_t_a_per_m": recoil,
            "result_recoil_branch_t_a_per_m": recoil,
            "temperature_c": 140.0,
            "result_temperature_c": 140.0,
            "load_step": 12,
            "result_load_step": 12,
            "magnet_owner": "magnet:rotor-v49-901",
            "result_magnet_owner": "magnet:rotor-v49-901",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        VIRTUAL_WORK: {
            "generation": virtual_work_generation,
            "displacement_generation": virtual_work_generation,
            "frame_generation": virtual_work_generation,
            "mesh_generation": virtual_work_generation,
            "energy_generation": virtual_work_generation,
            "result_generation": virtual_work_generation,
            "displacement_m": displacement,
            "result_displacement_m": displacement,
            "displacement_frame": "frame:global-x",
            "result_displacement_frame": "frame:global-x",
            "mesh_state_sha256": "2" * 64,
            "result_mesh_state_sha256": "2" * 64,
            "energy_j": energy,
            "result_energy_j": energy,
            "force_owner": "force:body-v49-901",
            "result_force_owner": "force:body-v49-901",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def _identity_v50() -> dict[str, object]:
    bem_generation = "bem-quadrature-v50-901"
    motion_generation = "motion-emf-v50-901"
    nearfield = {"distance_ratio": 0.15, "regularization": "adaptive-subdivision", "max_depth": 6}
    velocity = [12.0, 0.0, 0.0]
    path = [[0.0, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.2, 0.0]]
    return {
        V50_BEM: {
            "generation": bem_generation, "quadrature_generation": bem_generation, "self_panel_generation": bem_generation,
            "nearfield_generation": bem_generation, "mesh_generation": bem_generation, "result_generation": bem_generation,
            "singular_quadrature": "duffy-order-8", "result_singular_quadrature": "duffy-order-8",
            "self_panel_treatment": "analytic-solid-angle", "result_self_panel_treatment": "analytic-solid-angle",
            "nearfield_regularization": nearfield, "result_nearfield_regularization": nearfield,
            "mesh_owner": "mesh:bem-v50-901", "result_mesh_owner": "mesh:bem-v50-901",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        MOTION: {
            "generation": motion_generation, "velocity_generation": motion_generation, "frame_generation": motion_generation,
            "path_generation": motion_generation, "direction_generation": motion_generation, "result_generation": motion_generation,
            "velocity_m_s": velocity, "result_velocity_m_s": velocity,
            "velocity_frame": "frame:global", "result_velocity_frame": "frame:global",
            "conductor_path_m": path, "result_conductor_path_m": path,
            "integration_direction": "path-forward", "result_integration_direction": "path-forward",
            "emf_result_owner": "emf:conductor-v50-901", "result_emf_owner": "emf:conductor-v50-901",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def _identity_v51() -> dict[str, object]:
    generation = "magnetic-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    scalar = ["air", "iron"]
    vector = ["coil", "magnet"]
    traces = [{"scalar_domain": "air", "vector_domain": "magnet", "trace": "tangential_continuity"}]
    return {
        POTENTIAL: {
            "generation": generation, "gauge_generation": generation, "domain_generation": generation,
            "interface_generation": generation, "trace_generation": generation, "owner_generation": generation,
            "result_generation": generation, "gauge": "coulomb", "result_gauge": "coulomb",
            "scalar_potential_domains": scalar, "result_scalar_potential_domains": scalar,
            "vector_potential_domains": vector, "result_vector_potential_domains": vector,
            "interface_traces": traces, "result_interface_traces": traces, "interface_trace_sha256": "1" * 64,
            "result_interface_trace_sha256": "1" * 64, "solution_owner": "solution:coupled-v51",
            "result_solution_owner": "solution:coupled-v51", **result,
        },
        V51_BEM: {
            "generation": generation, "reciprocity_generation": generation, "symmetry_generation": generation,
            "orientation_generation": generation, "cache_generation": generation, "owner_generation": generation,
            "result_generation": generation, "matrix_shape": [128, 128], "result_matrix_shape": [128, 128],
            "reciprocity_relative_error": 2.0e-13, "result_reciprocity_relative_error": 2.0e-13,
            "symmetry_class": "symmetric", "result_symmetry_class": "symmetric", "panel_orientation": "outward",
            "result_panel_orientation": "outward", "panel_orientation_sha256": "2" * 64,
            "result_panel_orientation_sha256": "2" * 64, "cache_revision": "cache:bem-v51",
            "result_cache_revision": "cache:bem-v51", "matrix_owner": "matrix:bem-v51",
            "result_matrix_owner": "matrix:bem-v51", **result,
        },
    }


def _generations_v52(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity_v52() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    wrapped = [350.0, 355.0, 0.0, 5.0, 10.0]
    unwrapped = [350.0, 355.0, 360.0, 365.0, 370.0]
    energy = [1.00, 0.95, 0.90, 0.85, 0.80]
    torque = -((energy[3] - energy[1]) / math.radians(unwrapped[3] - unwrapped[1]))
    return {
        VIRTUAL_FORCE: {
            **_generations_v52("virtual-force-v52", ("mesh_generation", "energy_generation", "displacement_generation", "force_generation", "owner_generation", "result_generation")),
            "minus_mesh_sha256": "a" * 64, "result_minus_mesh_sha256": "a" * 64,
            "plus_mesh_sha256": "b" * 64, "result_plus_mesh_sha256": "b" * 64,
            "displacement_axis": [1.0, 0.0, 0.0], "result_displacement_axis": [1.0, 0.0, 0.0],
            "displacement_step_m": 1.0e-4, "result_displacement_step_m": 1.0e-4,
            "energy_minus_j": 1.002, "result_energy_minus_j": 1.002,
            "energy_plus_j": 0.998, "result_energy_plus_j": 0.998,
            "force_n": [20.0, 0.0, 0.0], "result_force_n": [20.0, 0.0, 0.0],
            "force_sign_convention": "negative_energy_gradient", "result_force_sign_convention": "negative_energy_gradient",
            "solution_owner": "solution:virtual-force-v52", "result_solution_owner": "solution:virtual-force-v52", **result,
        },
        MAGNET_TORQUE: {
            **_generations_v52("magnet-torque-v52", ("angle_generation", "unwrap_generation", "energy_generation", "derivative_generation", "owner_generation", "result_generation")),
            "angles_wrapped_deg": wrapped, "result_angles_wrapped_deg": wrapped,
            "angles_unwrapped_deg": unwrapped, "result_angles_unwrapped_deg": unwrapped,
            "angular_energy_j": energy, "result_angular_energy_j": energy,
            "torque_at_center_nm": torque, "result_torque_at_center_nm": torque,
            "derivative_sign_convention": "negative_energy_gradient", "result_derivative_sign_convention": "negative_energy_gradient",
            "magnet_owner": "magnet:torque-v52", "result_magnet_owner": "magnet:torque-v52", **result,
        },
    }


def _generations_v53(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity_v53():
    interactions = [{"source_panel": 11, "target_panel": 11, "separation_over_size": 0.0, "classification": "self", "rule": "duffy_singular"}, {"source_panel": 11, "target_panel": 12, "separation_over_size": 0.4, "classification": "near", "rule": "adaptive_near"}, {"source_panel": 11, "target_panel": 91, "separation_over_size": 4.5, "classification": "far", "rule": "gauss_far"}]
    quadrature = {**_generations_v53("quad-v53", ("classification_generation", "quadrature_generation", "panel_generation", "owner_generation", "result_generation")), "panel_interactions": interactions, "result_panel_interactions": interactions, "panel_owner": "panel-set:v53", "result_panel_owner": "panel-set:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    displacement = [-0.001, 0.0, 0.001]; force = [1.2, 0.0, -1.2]
    maglev = {**_generations_v53("maglev-v53", ("equilibrium_generation", "force_generation", "stiffness_generation", "displacement_generation", "owner_generation", "result_generation")), "displacement_path_m": displacement, "result_displacement_path_m": displacement, "force_path_n": force, "result_force_path_n": force, "equilibrium_index": 1, "result_equilibrium_index": 1, "equilibrium_force_n": 0.0, "result_equilibrium_force_n": 0.0, "stiffness_n_per_m": -1200.0, "result_stiffness_n_per_m": -1200.0, "body_owner": "body:v53", "result_body_owner": "body:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {QUADRATURE: quadrature, MAGLEV: maglev}


def _generations_v54(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload_v54():
    charges = [{"panel": 11, "magnetic_charge_a_m": 0.25}, {"panel": 12, "magnetic_charge_a_m": -0.25}]; normals = {"11": [0.0, 0.0, 1.0], "12": [0.0, 0.0, -1.0]}; regions = {"11": "region:magnet", "12": "region:air"}; orientations = {"11": 1, "12": -1}
    charge = {**_generations_v54("charge-v54", ("charge_generation", "normal_generation", "region_generation", "orientation_generation", "owner_generation", "result_generation")), "surface_charges": charges, "result_surface_charges": charges, "surface_normals": normals, "result_surface_normals": normals, "material_region_map": regions, "result_material_region_map": regions, "boundary_orientation": orientations, "result_boundary_orientation": orientations, "solution_owner": "solution:v54", "result_solution_owner": "solution:v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    direction = [0.0, 0.0, 1.0]; load_point = [0.0, 0.0, 0.025]
    stiffness = {**_generations_v54("stiffness-v54", ("gradient_generation", "coordinate_generation", "loadpoint_generation", "increment_generation", "owner_generation", "result_generation")), "force_gradient_n_per_m": -1200.0, "result_force_gradient_n_per_m": -1200.0, "stiffness_n_per_m": 1200.0, "result_stiffness_n_per_m": 1200.0, "coordinate_direction": direction, "result_coordinate_direction": direction, "load_point_m": load_point, "result_load_point_m": load_point, "displacement_increment_m": 1.0e-5, "result_displacement_increment_m": 1.0e-5, "body_owner": "body:v54", "result_body_owner": "body:v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {CHARGE: charge, STIFFNESS: stiffness}


def _payload_v55():
    gen = lambda name, fields: {"generation": name, **{field: name for field in fields}}
    samples = [{"h_a_per_m": 0.0, "b_t": 1.2}, {"h_a_per_m": -1.0e5, "b_t": 0.8}, {"h_a_per_m": -3.0e5, "b_t": 0.4}]
    demag = {**gen("demag-v55", ("hb_generation", "branch_generation", "volume_generation", "energy_generation", "material_generation", "owner_generation", "result_generation")), "hb_samples": samples, "result_hb_samples": samples, "curve_branch": "descending_demag", "result_curve_branch": "descending_demag", "material_volume_m3": 1.0e-5, "result_material_volume_m3": 1.0e-5, "demag_energy_j": 1.0, "result_demag_energy_j": 1.0, "material_revision": "magnet-v55-r4", "result_material_revision": "magnet-v55-r4", "material_owner": "material:magnet-v55", "result_material_owner": "material:magnet-v55", "solution_owner": "solution:demag-v55", "result_solution_owner": "solution:demag-v55", "result_sha256": "9" * 64, "accepted_result_sha256": "9" * 64}
    matrix = [[1200.0, 50.0], [50.0, 900.0]]; basis = {"x": [1.0, 0.0, 0.0], "y": [0.0, 1.0, 0.0]}; load = [0.0, 0.0, 0.025]
    bearing = {**gen("bearing-v55", ("matrix_generation", "basis_generation", "reciprocity_generation", "loadpoint_generation", "owner_generation", "result_generation")), "cross_stiffness_n_per_m": matrix, "result_cross_stiffness_n_per_m": matrix, "coordinate_basis": basis, "result_coordinate_basis": basis, "reciprocity_tolerance": 1.0e-10, "result_reciprocity_tolerance": 1.0e-10, "load_point_m": load, "result_load_point_m": load, "body_owner": "body:bearing-v55", "result_body_owner": "body:bearing-v55", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64}
    return {V55_DEMAG: demag, BEARING: bearing}
