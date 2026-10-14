"""Magnetic-force generalization gate contracts (v23-v55).

One positive closure per gate / chain head; every distinct failure signal is kept.
"""
from __future__ import annotations

import math
from copy import deepcopy

from radia_mcp.radia_ngsolve.bem_hysteresis_identity_v48 import (
    BEM as V48_BEM,
    HYSTERESIS,
    validate_public_identity as validate_v48_identity,
)
from radia_mcp.radia_ngsolve.bem_motion_identity_v50 import (
    BEM as V50_BEM,
    MOTION,
    validate_public_identity as validate_v50_identity,
)
from radia_mcp.radia_ngsolve.demag_virtual_work_identity_v49 import (
    DEMAG as V49_DEMAG,
    VIRTUAL_WORK,
    validate_public_identity as validate_v49_identity,
)
from radia_mcp.radia_ngsolve.energy_derivative_identity_v52 import (
    MAGNET_TORQUE,
    VIRTUAL_FORCE,
    validate_public_identity as validate_v52_identity,
)
from radia_mcp.radia_ngsolve.energy_derivative_identity_v53 import (
    MAGLEV,
    QUADRATURE,
    validate_public_identity as validate_v53_identity,
)
from radia_mcp.radia_ngsolve.energy_derivative_identity_v54 import (
    CHARGE,
    STIFFNESS,
    validate_public_identity as validate_v54_identity,
)
from radia_mcp.radia_ngsolve.magnetic_force_artifact_lineage_v47 import (
    FORCE,
    MOTOR,
    validate_public_v47_identity,
)
from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from radia_mcp.radia_ngsolve.magnetic_force_v44_identity import (
    validate_public_identity as validate_v44_identity,
)
from radia_mcp.radia_ngsolve.magnetostatic_energy_identity_v55 import (
    BEARING,
    DEMAG as V55_DEMAG,
    validate_public_identity as validate_v55_identity,
)
from radia_mcp.radia_ngsolve.potential_bem_identity_v51 import (
    BEM as V51_BEM,
    POTENTIAL,
    validate_public_identity as validate_v51_identity,
)

from _magnetic_force_generalization_payloads import (
    _BEARING,
    _CAPACITANCE,
    _CONDUCTOR,
    _CURRENT,
    _DEMAG,
    _DEMAG_V44,
    _DIELECTRIC,
    _DYNAMIC,
    _EDDY,
    _ELECTROSTATIC,
    _FORCE,
    _GEAR_KEY,
    _HEAT,
    _HYSTERESIS,
    _LINEAR,
    _MAGLEV,
    _MAGNET,
    _NONLINEAR,
    _SHIELD,
    _SOLENOID,
    _THIN_KEY,
    _TRANSFORMER,
    _force_gate,
    _force_summary_v29,
    _gate,
    _identity_v30,
    _identity_v31,
    _identity_v32,
    _identity_v33,
    _identity_v34,
    _identity_v35,
    _identity_v36,
    _identity_v37,
    _identity_v38,
    _identity_v39,
    _identity_v40,
    _identity_v41,
    _identity_v42,
    _identity_v43,
    _identity_v44,
    _identity_v47,
    _identity_v48,
    _identity_v49,
    _identity_v50,
    _identity_v51,
    _identity_v52,
    _identity_v53,
    _payload_v54,
    _payload_v55,
    _summary_v23,
    _summary_v24,
    _summary_v25,
    _summary_v26,
    _summary_v27,
    _summary_v28,
    _summary_v29,
    _summary_v30,
    _summary_v31,
    _summary_v32,
    _summary_v33,
    _summary_v34,
    _summary_v35,
    _summary_v36,
    _summary_v37,
    _summary_v38,
    _summary_v39,
    _summary_v40,
    _summary_v41,
    _summary_v42,
    _summary_v43,
)


def test_v23_public_bem_panel_normal_material_region_demag_force_generation_mismatch() -> None:
    summary = _summary_v23()
    summary["artifact_identity"][
        "bem_panel_normal_material_region_demag_force_generation_identity"
    ].update(
        {
            "panel_mesh_solve_generation": "bem-force-50",
            "outward_normal_solve_generation": "bem-force-49",
            "material_region_solve_generation": "bem-force-48",
            "demag_result_solve_generation": "bem-force-47",
            "force_result_solve_generation": "bem-force-46",
            "result_panel_ids": [101, 103, 104],
            "result_outward_normals": [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
            "result_material_region_ids": [1, 2, 3],
            "result_demag_field_a_per_m": [-90000.0, -62000.0, -30000.0],
            "result_force_vectors_n": [[-8.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
            "result_force_coordinate_frame": "local_xyz",
            "result_panel_force_table_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_demag_force_uses_current_panels_normals_materials_and_results"
    ]


def test_v23_public_motor_harmonic_rotor_angle_current_phase_force_frame_generation_mismatch() -> None:
    summary = _summary_v23()
    summary["artifact_identity"][
        "motor_harmonic_rotor_angle_current_phase_force_frame_generation_identity"
    ].update(
        {
            "rotor_angle_sweep_generation": "motor-harmonic-50",
            "current_phase_sweep_generation": "motor-harmonic-49",
            "harmonic_bin_sweep_generation": "motor-harmonic-48",
            "force_frame_sweep_generation": "motor-harmonic-47",
            "force_result_sweep_generation": "motor-harmonic-46",
            "result_rotor_angles_deg": [15.0, 10.0, 5.0, 0.0],
            "result_current_phase_deg": [0.0, 120.0, -120.0],
            "result_harmonic_bins": [0, 2, 1, 4],
            "result_force_coordinate_frame": "global_xyz",
            "result_force_harmonics_n": [[50.0, 0.0], [1.0, 0.0], [4.0, 2.0], [0.8, 0.2]],
            "result_harmonic_force_table_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_force_harmonics_use_current_angles_phases_bins_and_frame"
    ]


def test_v24_public_maglev_force_stiffness_energy_generation_mismatch() -> None:
    summary = _summary_v24()
    summary["artifact_identity"][
        "maglev_force_stiffness_equilibrium_energy_finite_difference_generation_identity"
    ].update(
        {
            "equilibrium_maglev_generation": "maglev-100",
            "displacement_maglev_generation": "maglev-99",
            "energy_maglev_generation": "maglev-98",
            "force_maglev_generation": "maglev-97",
            "stiffness_maglev_generation": "maglev-96",
            "coordinate_frame_maglev_generation": "maglev-95",
            "result_equilibrium_displacement_m": 0.0002,
            "result_equilibrium_force_n": 2.0,
            "energy_finite_difference_force_n": [-10.0, 0.0, 10.0],
            "reported_stiffness_n_m": -10000.0,
            "result_force_energy_sign_convention": "force=+dW/dx",
            "result_coordinate_frame": "local_z",
            "result_mesh_sha256": "8" * 64,
            "reported_maglev_result_sha256": "9" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "maglev_force_stiffness_energy_and_equilibrium_share_current_state"
    ]


def test_v24_public_motor_dual_lane_alignment_generation_mismatch() -> None:
    summary = _summary_v24()
    summary["artifact_identity"][
        "motor_dual_lane_geometry_excitation_force_frame_harmonic_alignment_generation_identity"
    ].update(
        {
            "geometry_comparison_generation": "dual-motor-100",
            "excitation_comparison_generation": "dual-motor-99",
            "force_frame_comparison_generation": "dual-motor-98",
            "harmonic_comparison_generation": "dual-motor-97",
            "rotor_angle_comparison_generation": "dual-motor-96",
            "result_geometry_revision_sha256": ["3" * 64, "a" * 64],
            "result_excitation_table_sha256": ["4" * 64, "b" * 64],
            "result_force_coordinate_frames": ["rotor_dq", "global_xyz"],
            "result_harmonic_bins": [[0, 1, 2, 3], [0, 2, 1, 4]],
            "result_rotor_angles_deg": [[0.0, 5.0, 10.0], [10.0, 5.0, 0.0]],
            "result_force_harmonics_n": [
                [[100.0, 0.0], [3.0, -1.0], [1.0, 0.5], [0.3, -0.1]],
                [[50.0, 0.0], [1.0, 0.0], [4.0, 2.0], [0.8, 0.2]],
            ],
            "reported_comparison_result_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_dual_lanes_share_geometry_excitation_frame_harmonics_and_angles"
    ]


def test_v25_public_bem_demag_surface_material_frame_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["artifact_identity"][
        "bem_demag_surface_orientation_magnetization_volume_material_frame_generation_identity"
    ]
    identity.update(
        {
            "surface_orientation_solve_generation": "bem-demag-200",
            "magnetization_solve_generation": "bem-demag-199",
            "body_volume_solve_generation": "bem-demag-198",
            "material_region_solve_generation": "bem-demag-197",
            "coordinate_frame_solve_generation": "bem-demag-196",
            "result_surface_orientation": "inward_to_magnet",
            "result_outward_normals": [[-1.0, 0.0, 0.0]],
            "result_magnetization_vector_a_per_m": [0.0, 900000.0, 0.0],
            "result_body_volume_m3": 2.0e-6,
            "result_material_region_id": 8,
            "result_coordinate_frame": "local_xyz",
            "result_surface_mesh_sha256": "7" * 64,
            "reported_demag_result_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_demag_surface_shares_orientation_magnetization_volume_material_and_frame"
    ]


def test_v25_public_linear_motor_thrust_ripple_generation_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["artifact_identity"][
        "linear_motor_thrust_ripple_period_position_phase_frame_generation_identity"
    ]
    identity.update(
        {
            "period_sweep_generation": "linear-thrust-200",
            "position_sweep_generation": "linear-thrust-199",
            "phase_sweep_generation": "linear-thrust-198",
            "force_frame_sweep_generation": "linear-thrust-197",
            "sample_order_sweep_generation": "linear-thrust-196",
            "result_mechanical_period_m": 0.04,
            "result_mover_positions_m": [0.03, 0.0, 0.01],
            "result_excitation_phase_deg": [0.0, 120.0, 60.0],
            "result_sample_order": [6, 0, 2],
            "result_force_coordinate_frame": "mover_local_x",
            "result_thrust_samples_n": [100.0, 98.0, 102.0],
            "reported_thrust_ripple_peak_to_peak_n": 2.0,
            "result_thrust_table_sha256": "9" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "linear_motor_thrust_ripple_shares_period_position_phase_frame_and_order"
    ]


def test_v26_public_levitation_force_displacement_gradient_stiffness_energy_derivative_frame_mismatch():
    summary = _summary_v26()
    summary["artifact_identity"]["levitation_force_displacement_gradient_stiffness_energy_derivative_frame_generation_identity"].update({
        "displacement_levitation_generation": "levitation-300", "result_force_n": [-10.0, 0.0, 10.0],
        "result_magnetic_energy_j": [0.0, 0.005, 0.0],
        "result_negative_energy_derivative_force_n": [-10.0, 0.0, 10.0],
        "result_restoring_stiffness_n_m": -10000.0, "result_coordinate_frame": "local_z_down",
        "result_force_sign_convention": "positive_gradient"})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["levitation_force_stiffness_and_energy_derivative_share_displacement_frame_sign_and_generation"]


def test_v26_public_cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_mismatch():
    summary = _summary_v26()
    summary["artifact_identity"]["cogging_torque_position_periodicity_mesh_interpolation_reference_angle_generation_identity"].update({
        "position_cogging_generation": "cogging-300",
        "result_mechanical_positions_deg": [5.0, 12.5, 20.0, 27.5, 35.0, 42.5],
        "result_cogging_torque_nm": [0.5, 1.0, 0.0, -1.0, 0.0, -0.5],
        "result_periodicity": 6, "result_mechanical_period_deg": 60.0,
        "result_reference_angle_deg": 5.0, "result_interpolation_method": "linear_open",
        "result_mesh_sha256": "9" * 64})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["cogging_torque_uses_current_position_periodicity_mesh_interpolation_reference_and_result"]


def test_v27_public_demag_bem_panel_orientation_magnetization_frame_self_term_energy_force_mismatch():
    summary = _summary_v27()
    summary["artifact_identity"][
        "bem_panel_orientation_magnetization_frame_self_term_energy_force_generation_identity"
    ].update({
        "orientation_panel_generation": "bem-panel-310",
        "self_term_panel_generation": "bem-panel-309",
        "mesh_panel_generation": "bem-panel-308",
        "result_outward_unit_normals": [[-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        "result_panel_area_m2": [0.5, -0.5],
        "result_magnetization_a_m": [[800000.0, 0.0, 0.0]],
        "result_magnetization_frame": "panel-local",
        "result_singular_self_term": "omitted",
        "result_magnetic_energy_j": [0.0, 0.005, 0.0],
        "result_force_n": [-10.0, 0.0, 10.0],
        "result_panel_mesh_sha256": "8" * 64,
        "accepted_force_result_sha256": "9" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_demag_force_uses_current_panel_orientation_magnetization_self_term_energy_mesh_and_result"
    ]


def test_v27_public_motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_mismatch():
    summary = _summary_v27()
    summary["artifact_identity"][
        "motor_reduced_basis_snapshot_operating_point_interpolation_torque_residual_generation_identity"
    ].update({
        "basis_reduced_generation": "motor-rom-310",
        "snapshot_reduced_generation": "motor-rom-309",
        "weight_reduced_generation": "motor-rom-308",
        "result_basis_dimension": 2,
        "result_snapshot_ids": ["snap-a", "snap-old"],
        "result_snapshot_operating_points": [[1000.0, 10.0, 0.0]],
        "result_query_operating_point": [4000.0, 50.0, 90.0],
        "result_interpolation_weights": [0.5, 0.5, 0.5],
        "result_snapshot_torque_nm": [1.0, 2.0],
        "result_reduced_torque_nm": 4.0,
        "result_relative_residual": 0.2,
        "loaded_basis_sha256": "a" * 64,
        "accepted_result_sha256": "b" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_reduced_torque_uses_current_basis_snapshots_operating_point_weights_residual_and_result"
    ]


def test_v28_public_maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_mismatch():
    summary = _summary_v28()
    identity = summary["artifact_identity"][
        "maglev_force_stiffness_displacement_step_coordinate_mesh_solution_derivative_generation_identity"
    ]
    identity.update(
        {
            "displacement_stiffness_generation": "maglev-stiffness-320",
            "mesh_stiffness_generation": "maglev-stiffness-319",
            "result_displacement_m": [0.0, 0.002, 0.004],
            "result_displacement_step_m": 0.002,
            "result_coordinate_direction": "local-x-negative",
            "result_force_n": [8.0, 10.0, 12.0],
            "result_derivative_convention": "stiffness-equals-force-derivative",
            "result_stiffness_n_m": -1000.0,
            "result_geometry_sha256": "a" * 64,
            "result_mesh_sha256": "b" * 64,
            "accepted_solution_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "maglev_stiffness_uses_current_displacement_coordinate_force_geometry_mesh_and_solution"
    ]


def test_v28_public_motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_mismatch():
    summary = _summary_v28()
    identity = summary["artifact_identity"][
        "motor_winding_harmonic_current_phase_rotor_angle_coenergy_torque_result_generation_identity"
    ]
    identity.update(
        {
            "winding_torque_generation": "motor-coenergy-320",
            "angle_torque_generation": "motor-coenergy-319",
            "result_phase_order": ["U", "W", "V"],
            "result_harmonic_orders": [1, 3, 5],
            "result_phase_current_harmonic_a": [[100.0, -40.0, -50.0]],
            "result_current_phase_deg": [0.0, 120.0, -120.0],
            "result_rotor_mechanical_angle_deg": [0.0, 2.0, 4.0],
            "result_coenergy_j": [0.100, 0.101, 0.099],
            "result_torque_convention": "negative-coenergy-angle-derivative",
            "result_torque_nm": [-0.03, 0.06],
            "result_mesh_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "motor_torque_uses_current_winding_harmonics_currents_phases_angles_coenergy_mesh_and_result"
    ]


def test_v29_public_magnet_demag_recoil_knee_temperature_local_field_volume_fraction_mismatch():
    summary = _summary_v29()
    identity = summary["artifact_identity"][
        "magnet_demag_recoil_knee_field_volume_generation_identity"
    ]
    identity.update({
        "recoil_demag_generation": "demag-state-330",
        "field_demag_generation": "demag-state-329",
        "result_recoil_relative_permeability": 1.2,
        "result_knee_field_a_m": -500000.0,
        "result_temperature_c": 20.0,
        "result_element_ids": [3, 2, 1],
        "result_local_recoil_axis_field_a_m": [-400000.0, -450000.0, -300000.0],
        "result_irreversible_mask": [False, False, False],
        "result_element_volumes_m3": [1.0e-6, 1.0e-6, 1.0e-6],
        "result_magnet_volume_m3": 3.0e-6,
        "result_irreversible_volume_fraction": 0.0,
        "result_material_state_sha256": "9" * 64,
        "result_mesh_sha256": "a" * 64,
        "accepted_result_sha256": "b" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnet_demag_uses_current_recoil_knee_temperature_local_field_mask_volume_and_result"
    ]


def test_v29_public_linear_motor_end_effect_phase_sequence_travel_wave_force_ripple_pitch_mismatch():
    summary = _summary_v29()
    identity = summary["artifact_identity"][
        "linear_motor_end_phase_wave_pitch_force_generation_identity"
    ]
    identity.update({
        "phase_linear_motor_generation": "linear-force-330",
        "end_effect_linear_motor_generation": "linear-force-329",
        "result_phase_sequence": ["U", "W", "V"],
        "result_traveling_wave_direction": "global-x-negative",
        "result_pole_pitch_m": 0.08,
        "result_position_m": [0.0, 0.02, 0.04, 0.08],
        "result_end_effect_factor": [1.0, 1.0, 1.0, 1.0],
        "result_force_n": [80.0, 120.0, 70.0, 130.0],
        "result_mean_force_n": 90.0,
        "result_force_ripple_peak_to_peak_n": 60.0,
        "accepted_result_sha256": "c" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "linear_motor_force_uses_current_end_effect_phase_wave_pitch_positions_ripple_and_result"
    ]


def test_v29_public_airgap_stress_harmonic_sector_periodicity_sampling_alias_torque_mismatch():
    summary = _force_summary_v29()
    summary["airgap_stress_harmonic_sector_periodicity_origin_sampling_alias_radius_torque_generation_identity"].update({
        "harmonic_stress_generation": "airgap-stress-160", "sampling_stress_generation": "airgap-stress-159",
        "result_stress_generation": "airgap-stress-158", "result_sector_pitch_deg": 45.0,
        "result_sector_count": 10, "result_angular_origin_deg": 15.0,
        "result_angular_sample_count": 30, "result_sector_sample_count": 3,
        "result_harmonic_orders": [0, 7, 13, 25], "result_torque_harmonics_nm": [5.0, -1.0],
        "result_alias_filter": "none", "result_alias_cutoff_order": 400,
        "result_airgap_radius_m": 0.06, "result_axial_length_m": 0.08,
        "result_torque_nm": 8.0, "result_airgap_mesh_sha256": "b" * 64,
        "accepted_torque_result_sha256": "c" * 64,
    })
    result = _force_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["airgap_torque_uses_current_sector_sampling_alias_harmonics_geometry_mesh_and_result"]


def test_v29_public_laminated_core_hysteresis_eddy_excess_loss_frequency_flux_volume_mismatch():
    summary = _force_summary_v29()
    summary["laminated_core_hysteresis_eddy_excess_frequency_flux_lamination_volume_result_generation_identity"].update({
        "hysteresis_loss_generation": "laminated-loss-160", "frequency_loss_generation": "laminated-loss-159",
        "result_loss_generation": "laminated-loss-158", "result_frequency_hz": 50.0,
        "result_peak_flux_density_t": 0.8, "result_lamination_thickness_m": 0.0005,
        "result_magnetic_volume_m3": 0.01, "result_hysteresis_coefficient": 30.0,
        "result_hysteresis_exponent": 2.0, "result_eddy_coefficient": 0.2,
        "result_excess_coefficient": 2.0, "result_hysteresis_loss_w": 1.0,
        "result_eddy_loss_w": 2.0, "result_excess_loss_w": -1.0,
        "result_total_core_loss_w": 99.0, "result_material_sha256": "d" * 64,
        "accepted_loss_result_sha256": "e" * 64,
    })
    result = _force_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["laminated_core_loss_uses_current_frequency_flux_lamination_volume_components_and_result"]


def test_v30_public_nonlinear_coenergy_force_current_perturbation_remesh_central_difference_mismatch():
    identity = _identity_v30()
    identity[
        "nonlinear_coenergy_force_current_perturbation_remesh_central_difference_frame_result_identity"
    ].update({
        "current_force_generation": "nonlinear-coenergy-170",
        "remesh_force_generation": "nonlinear-coenergy-169",
        "result_force_generation": "nonlinear-coenergy-168",
        "result_current_constraint": "fixed_flux",
        "result_branch_currents_a": [9.0, 10.0, 11.0],
        "result_displacements_m": [-2.0e-4, 0.0, 1.0e-4],
        "result_coenergy_j": [2.005, 2.0, 2.02],
        "result_difference_rule": "forward_difference",
        "result_branch_mesh_generations": ["mesh-old", "mesh-center-171", "mesh-plus-171"],
        "result_displacement_frame": "local_r",
        "result_force_n": -150.0,
        "accepted_branch_result_sha256": "8" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_coenergy_force_uses_fixed_current_symmetric_displacement_remesh_frame_and_result"
    ]


def test_v30_public_axisymmetric_force_r_weight_jacobian_coordinate_stress_contour_mismatch():
    identity = _identity_v30()
    identity[
        "axisymmetric_force_radial_weight_jacobian_coordinate_stress_contour_material_mesh_result_identity"
    ].update({
        "radial_weight_generation": "axisym-stress-170",
        "stress_contour_generation": "axisym-stress-169",
        "result_generation": "axisym-stress-168",
        "result_coordinate_convention": "x_y_planar",
        "result_radius_m": 0.05,
        "result_radial_weight": 1.0,
        "result_line_jacobian_m": 0.002,
        "result_stress_normal_pa": -1200.0,
        "result_stress_contour_closed": False,
        "result_stress_contour_orientation": "clockwise",
        "result_stress_contour_material_side": "iron",
        "result_force_n": 12.0,
        "result_mesh_sha256": "9" * 64,
        "accepted_force_result_sha256": "a" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_stress_force_uses_two_pi_r_jacobian_air_contour_mesh_and_result"
    ]


def test_v30_public_maglev_force_stiffness_position_current_perturbation_derivative_mesh_mismatch():
    summary = _summary_v30()
    summary["artifact_identity"]["maglev_force_stiffness_position_current_derivative_frame_mesh_result_identity"].update({
        "position_maglev_generation": "maglev-stiffness-340", "mesh_maglev_generation": "maglev-stiffness-339",
        "result_position_m": [0.0, 0.001, 0.002], "result_current_a": [9.0, 10.0, 11.0],
        "result_force_z_n": [8.0, 10.0, 12.0], "result_stiffness_n_per_m": 2000.0,
        "result_coordinate_frame": "local_r_positive_down", "result_mesh_sha256": "7" * 64,
        "accepted_result_sha256": "8" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["maglev_stiffness_uses_fixed_current_symmetric_positions_force_derivative_frame_mesh_and_result"]


def test_v30_public_cogging_torque_slot_pole_period_angular_sampling_harmonic_phase_mismatch():
    summary = _summary_v30()
    summary["artifact_identity"]["cogging_torque_slot_pole_period_origin_sampling_harmonic_phase_mesh_result_identity"].update({
        "slot_cogging_generation": "cogging-340", "sampling_cogging_generation": "cogging-339",
        "result_slot_count": 9, "result_pole_count": 8, "result_cogging_period_mechanical_deg": 15.0,
        "result_angular_origin_deg": 2.0, "result_sample_angles_deg": [2.0, 4.0, 6.0],
        "result_harmonic_orders": [3, 7], "result_harmonic_phase_deg": [45.0, -20.0],
        "result_torque_nm": [0.2, 0.1, 0.3], "result_mesh_sha256": "9" * 64,
        "accepted_result_sha256": "a" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["cogging_torque_uses_slot_pole_period_origin_sampling_harmonics_phase_mesh_and_result"]


def test_v31_public_nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mismatch():
    identity = _identity_v31()
    record = identity[
        "nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mesh_result_identity"
    ]
    record.update({
        "current_path_generation": "incremental-inductance-180",
        "result_current_state_index": 3,
        "result_current_samples_a": [9.9, 10.0, 10.2],
        "result_derivative_rule": "forward_difference",
        "result_incremental_inductance_h": 0.19,
        "result_magnetic_state_sha256": "a" * 64,
        "accepted_result_sha256": "b" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_incremental_inductance_uses_current_path_state_symmetric_derivative_circuit_mesh_and_result"
    ]


def test_v31_public_lamination_anisotropy_fill_factor_orientation_frequency_loss_balance_mismatch():
    identity = _identity_v31()
    record = identity[
        "lamination_anisotropy_fill_orientation_frequency_loss_volume_balance_result_identity"
    ]
    record.update({
        "anisotropy_generation": "lamination-loss-180",
        "result_mu_axis_order": ["transverse", "rolling"],
        "result_lamination_fill_factor": 1.0,
        "result_frequency_hz": 50.0,
        "result_active_iron_volume_m3": 0.001,
        "result_total_core_loss_w": 9.0,
        "accepted_result_sha256": "c" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "laminated_loss_uses_current_anisotropy_fill_orientation_frequency_volume_balance_and_result"
    ]


def test_v31_public_bem_near_singular_gap_quadrature_normal_force_reciprocity_mismatch():
    summary = _summary_v31(); record = summary["artifact_identity"]["bem_near_singular_gap_quadrature_normal_order_force_reciprocity_geometry_result_identity"]
    record.update({"gap_bem_generation": "old", "result_quadrature_policy": "fixed_gauss", "result_source_normal": [0., 0., -1.], "result_force_on_target_n": [0., 0., 4.], "result_action_reaction_residual_n": 9., "accepted_result_sha256": "8" * 64})
    result = magnetic_force_method_profile_gate(summary); assert result["status"] == "needs_attention"; assert not result["checks"]["bem_near_contact_force_uses_gap_adaptive_quadrature_normals_order_reciprocity_geometry_and_result"]


def test_v31_public_hysteresis_minor_loop_state_remanence_return_point_energy_dissipation_mismatch():
    summary = _summary_v31(); record = summary["artifact_identity"]["hysteresis_minor_loop_state_reversal_return_memory_remanence_energy_time_material_identity"]
    record.update({"state_loop_generation": "old", "result_time_s": [0., .2, .1, .3, .4], "result_return_point_memory_closed": False, "result_remanence_t": -.2, "result_loop_energy_j_per_m3": -42., "result_material_owner": "old", "accepted_result_sha256": "a" * 64})
    result = magnetic_force_method_profile_gate(summary); assert result["status"] == "needs_attention"; assert not result["checks"]["hysteresis_minor_loop_uses_initial_state_reversals_return_memory_remanence_energy_time_and_material"]


def test_v32_public_axisymmetric_force_weighted_stress_coenergy_contour_radius_mesh_mismatch():
    identity = _identity_v32()
    record = identity[
        "axisymmetric_weighted_stress_coenergy_contour_displacement_radius_weight_material_mesh_owner_result_identity"
    ]
    record.update(
        {
            "stress_generation": "axisymmetric-force-190",
            "mesh_generation": "axisymmetric-force-189",
            "result_generation": "axisymmetric-force-188",
            "result_stress_method": "contour_maxwell_stress",
            "result_coenergy_method": "forward_displacement",
            "result_contour_radius_m": 0.03,
            "result_virtual_displacement_m": 1.0e-3,
            "result_axisymmetric_weight": "1",
            "result_material_side": "steel",
            "coenergy_force_n": -10.0,
            "result_force_mesh_sha256": "9" * 64,
            "result_force_owner": "planar/group2",
            "accepted_force_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_force_closes_weighted_stress_coenergy_contour_displacement_weight_material_mesh_owner_and_result"
    ]


def test_v32_public_laminated_diffusion_skin_depth_complex_power_frequency_mesh_mismatch():
    identity = _identity_v32()
    record = identity[
        "laminated_diffusion_conductivity_skin_depth_frequency_phasor_power_volume_mesh_loss_result_identity"
    ]
    record.update(
        {
            "conductivity_generation": "laminated-diffusion-190",
            "frequency_generation": "laminated-diffusion-189",
            "result_generation": "laminated-diffusion-188",
            "result_conductivity_tensor_s_per_m": [[1.0e3, 0.0], [0.0, 2.0e6]],
            "result_skin_depth_m": 5.0e-3,
            "result_frequency_hz": 50.0,
            "result_phasor_convention": "exp(-jwt)_peak",
            "result_complex_power_va_ri": [1.2, -4.7],
            "result_active_volume_m3": 1.0e-3,
            "result_laminated_mesh_sha256": "b" * 64,
            "result_laminated_loss_w": 9.0,
            "accepted_laminated_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "laminated_diffusion_uses_current_conductivity_skin_depth_frequency_phasor_power_volume_mesh_loss_and_result"
    ]


def test_v32_public_maglev_equilibrium_force_displacement_stiffness_derivative_sign_mesh_mismatch():
    summary = _summary_v32()
    record = summary["artifact_identity"][
        "maglev_equilibrium_force_displacement_derivative_stiffness_gravity_mesh_result_identity"
    ]
    record.update(
        {
            "force_equilibrium_generation": "maglev-equilibrium-360",
            "mesh_equilibrium_generation": "maglev-equilibrium-359",
            "result_equilibrium_generation": "maglev-equilibrium-358",
            "result_force_sign_convention": "positive_down",
            "result_displacement_frame": "local_y_down",
            "result_displacement_samples_m": [0.0, 1.0e-4, 2.0e-4],
            "result_magnetic_force_samples_n": [-9.81, -9.7, -9.5],
            "result_derivative_stencil": "forward_difference",
            "result_force_derivative_n_per_m": 1100.0,
            "result_vertical_stiffness_n_per_m": -1100.0,
            "result_gravity_force_n": 9.81,
            "result_mesh_sha256": "7" * 64,
            "accepted_result_owner": "maglev/old-case",
            "accepted_result_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "maglev_equilibrium_uses_upward_force_global_displacement_central_stiffness_gravity_mesh_and_result"
    ]


def test_v32_public_bem_surface_charge_net_zero_gauge_reference_energy_reciprocity_mismatch():
    summary = _summary_v32()
    record = summary["artifact_identity"][
        "bem_surface_charge_gauge_normal_energy_reciprocity_geometry_owner_result_identity"
    ]
    record.update(
        {
            "charge_bem_generation": "bem-charge-360",
            "geometry_bem_generation": "bem-charge-359",
            "result_bem_generation": "bem-charge-358",
            "result_net_surface_charge": 0.1,
            "result_gauge_reference": "pin_first_node",
            "result_source_normal": [0.0, 0.0, -1.0],
            "result_target_normal": [0.0, 0.0, -1.0],
            "result_field_energy_j": -0.25,
            "result_reciprocity_residual": 0.2,
            "result_geometry_sha256": "9" * 64,
            "accepted_result_owner": "bem/old-case",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_surface_charge_uses_neutral_charge_mean_zero_gauge_opposed_normals_energy_reciprocity_geometry_and_result"
    ]


def test_v33_public_electrostatic_capacitance_matrix_reciprocity_charge_neutrality_energy_mismatch():
    identity = _identity_v33()
    identity[
        "electrostatic_capacitance_conductor_reciprocity_neutrality_voltage_charge_energy_unit_mesh_owner_result_identity"
    ].update(
        {
            "reciprocity_generation": "electrostatic-capacitance-200",
            "energy_generation": "electrostatic-capacitance-199",
            "result_generation": "electrostatic-capacitance-198",
            "result_conductor_order": ["right", "shield", "left"],
            "result_capacitance_matrix_f": [[2e-12, 1e-12, -1e-12], [-2e-12, 2e-12, 0.0], [-1e-12, -1e-12, 3e-12]],
            "result_voltages_v": [-1.0, 0.0, 1.0],
            "result_charges_c": [3e-12, 1e-12, -2e-12],
            "result_electrostatic_energy_j": -3e-12,
            "result_capacitance_unit": "pF",
            "result_charge_unit": "nC",
            "result_electrostatic_mesh_sha256": "9" * 64,
            "accepted_result_owner": "stale/result",
            "accepted_electrostatic_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrostatic_capacitance_closes_conductors_reciprocity_neutrality_charge_energy_units_mesh_owner_and_result"
    ]


def test_v33_public_axisymmetric_heat_flux_conduction_convection_source_2pir_balance_mismatch():
    identity = _identity_v33()
    identity[
        "axisymmetric_heat_conduction_convection_source_boundary_flux_weight_temperature_mesh_owner_result_identity"
    ].update(
        {
            "source_generation": "axisymmetric-heat-balance-200",
            "weight_generation": "axisymmetric-heat-balance-199",
            "result_generation": "axisymmetric-heat-balance-198",
            "result_conduction_w": 20.0,
            "result_convection_w": 10.0,
            "result_boundary_flux_w": 5.0,
            "result_axisymmetric_weight": "1",
            "result_temperature_reference_k": 20.0,
            "result_heat_mesh_sha256": "b" * 64,
            "result_heat_result_owner": "planar/old",
            "accepted_heat_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_heat_closes_conduction_convection_source_flux_weight_temperature_mesh_owner_and_result"
    ]


def test_v33_public_halbach_harmonic_magnetization_order_pole_pitch_phase_field_energy_force_mismatch():
    summary = _summary_v33()
    record = summary["artifact_identity"][
        "halbach_harmonic_magnetization_order_pitch_phase_grid_field_energy_force_geometry_owner_result_identity"
    ]
    record.update(
        {
            "magnetization_generation": "halbach-harmonic-370",
            "field_generation": "halbach-harmonic-369",
            "result_generation": "halbach-harmonic-368",
            "result_magnetization_angles_deg": [270.0, 180.0, 90.0, 0.0],
            "result_pole_pitch_m": 0.04,
            "result_harmonic_orders": [1, 2, 4],
            "result_harmonic_phase_deg": [180.0, 0.0, 0.0],
            "result_sampling_grid_m": [0.0, 0.01, 0.02],
            "result_field_harmonic_amplitude_t": [0.2, 0.4, 0.8],
            "result_magnetic_energy_j": -0.5,
            "result_force_direction": "-z",
            "result_force_n": -10.0,
            "result_geometry_sha256": "7" * 64,
            "accepted_result_owner": "halbach/old-case",
            "accepted_result_sha256": "8" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "halbach_harmonics_use_current_magnetization_pitch_phase_grid_field_energy_force_geometry_owner_and_result"
    ]


def test_v33_public_magnetic_bearing_force_current_displacement_stiffness_matrix_reciprocity_stability_mismatch():
    summary = _summary_v33()
    record = summary["artifact_identity"][
        "magnetic_bearing_force_current_displacement_stiffness_reciprocity_bias_frame_mesh_result_identity"
    ]
    record.update(
        {
            "jacobian_generation": "magnetic-bearing-370",
            "stiffness_generation": "magnetic-bearing-369",
            "result_generation": "magnetic-bearing-368",
            "result_coordinate_frame": "local_yz_left_handed",
            "result_bias_currents_a": [2.0, 2.0],
            "result_bias_displacement_m": [0.001, -0.001],
            "result_force_current_jacobian_n_per_a": [[1.0, 2.0]],
            "result_force_displacement_jacobian_n_per_m": [[100.0, 200.0], [-50.0, 100.0]],
            "result_stiffness_matrix_n_per_m": [[-100.0, -200.0], [50.0, -100.0]],
            "result_stiffness_eigenvalues_n_per_m": [-100.0, -100.0],
            "result_mesh_sha256": "9" * 64,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_bearing_uses_current_force_jacobians_bias_frame_reciprocal_positive_stiffness_mesh_and_result"
    ]


def test_v34_public_kelvin_transform_open_boundary_radius_energy_flux_far_field_mismatch():
    identity = _identity_v34()
    identity["kelvin_transform_radius_permeability_jacobian_interface_energy_flux_far_field_mesh_owner_result_identity"].update({
        "radius_generation": "kelvin-open-boundary-210", "flux_generation": "kelvin-open-boundary-209", "result_generation": "kelvin-open-boundary-208",
        "result_kelvin_radius_m": 2.0, "result_mapped_permeability_relative": [1.0, 1.0, 1.0],
        "result_mapping_jacobian_determinants": [1.0, -0.25, 0.0],
        "result_interface_potential_jump": 0.1, "result_interface_normal_flux_jump_wb": 0.02,
        "result_magnetic_energy_j": -0.012, "result_outer_flux_wb": 0.01,
        "result_far_field_potential": [0.5, 0.6, 0.7], "result_kelvin_mesh_sha256": "9" * 64,
        "accepted_kelvin_result_owner": "stale/owner", "accepted_kelvin_result_sha256": "a" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "kelvin_open_boundary_closes_mapping_interface_energy_flux_far_field_mesh_owner_and_result"
    ]


def test_v34_public_force_contour_path_airgap_stress_fourier_virtual_work_symmetry_mismatch():
    identity = _identity_v34()
    identity["force_airgap_contour_stress_fourier_virtual_work_symmetry_torque_origin_field_result_identity"].update({
        "contour_generation": "airgap-force-210", "virtual_work_generation": "airgap-force-209", "result_generation": "airgap-force-208",
        "result_airgap_contour_ids": [2, 1], "result_contour_forces_n": [[10.0, 0.0], [-5.0, 2.0]],
        "result_fourier_stress_harmonics_n": [[0, 10.0], [1, 4.0], [2, -3.0]],
        "result_virtual_work_force_n": -10.0, "result_virtual_work_displacement_m": -1.0e-5,
        "result_torque_origin_m": [1.0, 0.0], "result_torque_nm": 5.0,
        "result_force_symmetry": "none", "result_force_field_sha256": "b" * 64,
        "accepted_force_result_owner": "stale/force", "accepted_force_result_sha256": "c" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "airgap_force_closes_contours_fourier_virtual_work_symmetry_torque_field_owner_and_result"
    ]


def test_v34_rejects_self_consistent_nondecaying_kelvin_far_field():
    identity = _identity_v34()
    row = identity["kelvin_transform_radius_permeability_jacobian_interface_energy_flux_far_field_mesh_owner_result_identity"]
    row["far_field_potential"] = [0.5, 0.5, 0.5]
    row["result_far_field_potential"] = [0.5, 0.5, 0.5]
    assert _gate(identity)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_force_contour_disagreement():
    identity = _identity_v34()
    row = identity["force_airgap_contour_stress_fourier_virtual_work_symmetry_torque_origin_field_result_identity"]
    row["contour_forces_n"] = [[10.0, 0.0], [5.0, 0.0]]
    row["result_contour_forces_n"] = [[10.0, 0.0], [5.0, 0.0]]
    assert _gate(identity)["status"] == "needs_attention"


def test_v34_public_magnetic_bearing_force_matrix_cross_coupled_stiffness_damping_stability_mismatch():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
    ]
    record.update(
        {
            "force_generation": "bearing-dynamic-380",
            "stability_generation": "bearing-dynamic-379",
            "result_generation": "bearing-dynamic-378",
            "result_coordinate_order": ["y", "x"],
            "result_displacement_perturbations_m": [[0.0, 1.0e-3]],
            "result_force_perturbations_n": [[1.0, 1.0]],
            "result_stiffness_matrix_n_per_m": [[-1000.0, 500.0], [-50.0, -900.0]],
            "result_damping_matrix_n_s_per_m": [[-10.0, 20.0], [2.0, -12.0]],
            "result_state_eigenvalues_per_s": [[5.0, 30.0], [6.0, -28.0]],
            "result_operating_displacement_m": [0.001, -0.001],
            "result_operating_velocity_m_s": [1.0, 0.0],
            "result_bias_currents_a": [1.0, 2.0],
            "result_bearing_mesh_sha256": "a" * 64,
            "accepted_bearing_result_owner": "bearing/old",
            "accepted_bearing_result_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_bearing_dynamics_use_current_force_perturbations_stiffness_damping_coordinates_stability_operating_point_owner_and_result"
    ]


def test_v34_public_moving_conductor_eddy_drag_lift_power_velocity_skin_depth_sign_mismatch():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
    ]
    record.update(
        {
            "velocity_generation": "moving-conductor-380",
            "power_generation": "moving-conductor-379",
            "result_generation": "moving-conductor-378",
            "result_coordinate_frame": "body_left_handed",
            "result_velocity_m_s": [-10.0, 0.0, 0.0],
            "result_drag_force_n": [100.0, 0.0, 0.0],
            "result_lift_force_n": [100.0, 20.0, 0.0],
            "result_joule_power_w": -1000.0,
            "result_mechanical_drag_power_w": 500.0,
            "result_conductivity_s_m": -3.5e7,
            "result_relative_permeability": -1.0,
            "result_excitation_frequency_hz": 25.0,
            "result_spatial_period_m": 0.1,
            "result_slip_frequency_hz": -50.0,
            "result_skin_depth_m": -0.01,
            "result_conductor_mesh_sha256": "c" * 64,
            "accepted_field_sha256": "d" * 64,
            "accepted_conductor_result_owner": "moving-conductor/old",
            "accepted_conductor_result_sha256": "e" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "moving_conductor_uses_current_velocity_frame_drag_lift_power_skin_depth_slip_mesh_owner_field_and_result"
    ]


def test_v34_public_rejects_self_consistent_negative_bearing_damping():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "magnetic_bearing_perturbation_cross_coupled_stiffness_damping_coordinate_stability_operating_owner_result_identity"
    ]
    damping = [[-10.0, 2.0], [-2.0, -12.0]]
    record["damping_matrix_n_s_per_m"] = damping
    record["result_damping_matrix_n_s_per_m"] = damping
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_drag_aligned_with_velocity():
    summary = _summary_v34()
    record = summary["artifact_identity"][
        "moving_conductor_velocity_frame_drag_lift_joule_work_skin_depth_frequency_slip_mesh_owner_field_result_identity"
    ]
    record["drag_force_n"] = [100.0, 0.0, 0.0]
    record["result_drag_force_n"] = [100.0, 0.0, 0.0]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v35_public_inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_mismatch():
    identity = _identity_v35()
    identity["inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_coil_mesh_owner_result_identity"].update({
        "reciprocity_generation": "inductance-matrix-234", "energy_generation": "inductance-matrix-233",
        "result_generation": "inductance-matrix-232", "result_coil_order": ["phase_b", "phase_a"],
        "result_currents_a": [-1.0, 2.0], "result_inductance_matrix_h": [[0.008, 0.004], [-0.003, -0.006]],
        "result_flux_linkages_wb_turn": [0.2, 0.1], "result_stored_energy_j": -0.015,
        "result_minimum_eigenvalue_h": -0.01, "result_coil_mesh_sha256": "9" * 64,
        "accepted_inductance_result_owner": "stale/coils", "accepted_inductance_result_sha256": "a" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["inductance_matrix_closes_reciprocity_psd_flux_current_energy_coils_mesh_owner_and_result"]


def test_v35_public_nonlinear_bh_interpolation_differential_permeability_energy_coenergy_branch_mismatch():
    identity = _identity_v35()
    identity["nonlinear_bh_interpolation_differential_permeability_branch_energy_coenergy_operating_material_solution_identity"].update({
        "interpolation_generation": "nonlinear-bh-234", "energy_generation": "nonlinear-bh-233",
        "result_generation": "nonlinear-bh-232", "result_b_samples_t": [0.0, 0.5, 0.4, 1.4],
        "result_h_samples_a_m": [0.0, 100.0, 80.0, 1200.0],
        "result_differential_permeability_h_m": [0.005, -0.005, 0.001],
        "result_branch": "descending", "result_operating_point_b_t": 1.0,
        "result_operating_point_h_a_m": 400.0, "result_magnetic_energy_density_j_m3": -470.0,
        "result_magnetic_coenergy_density_j_m3": 470.0,
        "accepted_material_owner": "stale/material", "accepted_nonlinear_solution_sha256": "b" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["nonlinear_bh_closes_interpolation_differential_mu_branch_energy_coenergy_operating_material_and_solution"]


def test_v35_rejects_self_consistent_inductance_energy_mismatch():
    identity = _identity_v35()
    row = identity["inductance_matrix_reciprocity_psd_fluxlinkage_current_energy_coil_mesh_owner_result_identity"]
    row["stored_energy_j"] = row["result_stored_energy_j"] = 0.03
    assert _gate(identity)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_nonmonotone_bh_branch():
    identity = _identity_v35()
    row = identity["nonlinear_bh_interpolation_differential_permeability_branch_energy_coenergy_operating_material_solution_identity"]
    row["h_samples_a_m"] = row["result_h_samples_a_m"] = [0.0, 100.0, 80.0, 1200.0]
    assert _gate(identity)["status"] == "needs_attention"


def test_v35_public_magnetic_gear_polepair_harmonic_torque_phase_power_balance_mismatch():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
    ]
    record.update(
        {
            "pole_generation": "magnetic-gear-400",
            "power_generation": "magnetic-gear-399",
            "result_generation": "magnetic-gear-398",
            "result_high_speed_pole_pairs": 5,
            "result_low_speed_pole_pairs": 20,
            "result_modulator_pole_count": 24,
            "result_transmitted_harmonic_order": 21,
            "result_high_speed_torque_nm": 10.0,
            "result_low_speed_torque_nm": 20.0,
            "result_high_speed_angular_velocity_rad_s": 100.0,
            "result_low_speed_angular_velocity_rad_s": 40.0,
            "result_high_speed_harmonic_phase_rad": -0.5,
            "result_low_speed_harmonic_phase_rad": 0.7,
            "result_modulator_phase_rad": -0.2,
            "result_transmitted_phase_rad": 2.5,
            "result_coordinate_frame": "rotor_left_handed",
            "result_gear_mesh_sha256": "9" * 64,
            "accepted_gear_result_owner": "magnetic-gear/old",
            "accepted_gear_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_gear_uses_current_poles_harmonic_torque_phase_power_frame_mesh_owner_and_result"
    ]


def test_v35_public_demag_bem_surface_charge_neutrality_normal_farfield_energy_mesh_mismatch():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
    ]
    record.update(
        {
            "charge_generation": "demag-bem-400",
            "farfield_generation": "demag-bem-399",
            "solution_generation": "demag-bem-398",
            "result_panel_areas_m2": [1.0, -1.0],
            "result_surface_charge_density_a_m": [2.0, 2.0],
            "result_surface_charge_integral_a_m": 4.0,
            "result_outward_normals": [[0.0, 0.0, 0.0]],
            "result_outward_orientation_verified": False,
            "result_normal_field_jump_a_m": [-2.0, 2.0, -1.0, 1.0],
            "result_farfield_radius_m": [8.0, 4.0, 2.0],
            "result_farfield_potential_a": [0.25, 0.25, 0.25],
            "result_farfield_field_a_m": [0.125, 0.125, 0.125],
            "result_magnetic_energy_j": -0.75,
            "result_boundary_mesh_sha256": "b" * 64,
            "accepted_demag_solution_owner": "demag-bem/old",
            "accepted_demag_solution_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "demag_bem_uses_neutral_surface_charge_outward_normals_jump_farfield_energy_mesh_owner_and_solution"
    ]


def test_v35_public_rejects_self_consistent_non_neutral_demag_charge():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
    ]
    charges = [2.0, 2.0, 1.0, 1.0]
    record["surface_charge_density_a_m"] = charges
    record["result_surface_charge_density_a_m"] = charges
    record["normal_field_jump_a_m"] = charges
    record["result_normal_field_jump_a_m"] = charges
    record["surface_charge_integral_a_m"] = 8.0
    record["result_surface_charge_integral_a_m"] = 8.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v35_public_rejects_self_consistent_magnetic_gear_power_imbalance():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
    ]
    record["low_speed_torque_nm"] = 50.0
    record["result_low_speed_torque_nm"] = 50.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v36_public_nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_mismatch():
    identity = _identity_v36()
    identity["nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_current_mesh_owner_solution_identity"].update({
        "branch_generation": "nonlinear-incremental-235", "mesh_generation": "nonlinear-incremental-234",
        "result_generation": "nonlinear-incremental-233", "result_bh_branch": "descending",
        "result_current_points_a": [3.0, 2.0, 1.0], "result_flux_linkages_wb_turn": [0.024, 0.018, 0.01],
        "result_incremental_permeability_h_m": -0.0012, "result_differential_inductance_h": -0.006,
        "result_magnetic_energy_j": -0.03, "result_magnetic_coenergy_j": 0.01,
        "result_mesh_sha256": "a" * 64, "accepted_field_owner": "stale/core",
        "accepted_solution_sha256": "b" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["nonlinear_incremental_solution_closes_branch_mu_energy_coenergy_differential_inductance_current_mesh_owner_and_result"]


def test_v36_public_weighted_stress_force_region_weighting_contour_independence_mesh_convergence_mismatch():
    identity = _identity_v36()
    identity["weighted_stress_force_region_air_contour_mesh_convergence_direction_owner_result_identity"].update({
        "weighting_generation": "weighted-stress-force-235", "direction_generation": "weighted-stress-force-234",
        "result_generation": "weighted-stress-force-233", "result_weighting_region_id": "iron:body",
        "result_air_enclosure_id": "air:old", "result_weighted_force_n": [-12.0, 2.0],
        "result_contour_force_samples_n": [[3.0, 4.0], [-5.0, 1.0]],
        "result_mesh_sizes_m": [0.001, 0.004], "result_mesh_force_sequence_n": [[3.0, 0.0], [20.0, 0.0]],
        "result_force_direction_unit": [-1.0, 2.0], "accepted_field_owner": "stale/field",
        "accepted_force_result_sha256": "c" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["weighted_stress_force_closes_region_air_contours_mesh_convergence_direction_owner_and_result"]


def test_v36_rejects_self_consistent_incremental_energy_partition_error():
    identity = _identity_v36()
    row = identity["nonlinear_bh_incremental_permeability_energy_coenergy_differential_inductance_current_mesh_owner_solution_identity"]
    row["magnetic_coenergy_j"] = row["result_magnetic_coenergy_j"] = 0.02
    assert _gate(identity)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_nonconvergent_weighted_force_sequence():
    identity = _identity_v36()
    row = identity["weighted_stress_force_region_air_contour_mesh_convergence_direction_owner_result_identity"]
    row["mesh_force_sequence_n"] = row["result_mesh_force_sequence_n"] = [[11.5, -0.3], [10.0, -0.2], [12.0, -0.2]]
    assert _gate(identity)["status"] == "needs_attention"


def test_v36_public_magnetic_bearing_stiffness_force_displacement_bias_current_linearization_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
    ]
    record.update(
        {
            "bias_generation": "bearing-linearization-401",
            "stiffness_generation": "bearing-linearization-400",
            "result_generation": "bearing-linearization-399",
            "result_bias_current_a": -5.0,
            "result_displacement_samples_m": [0.0, 0.001],
            "result_force_x_samples_n": [10.0, 10.0],
            "result_stiffness_matrix_n_m": [[-10000.0, 5000.0], [-100.0, -9000.0]],
            "result_coordinate_frame": "rotor_left_handed",
            "accepted_bearing_owner": "bearing/old",
            "accepted_bearing_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_bearing_bias_sweep_uses_current_bias_symmetric_force_derivative_crosscoupling_frame_owner_and_result"
    ]


def test_v36_public_pm_demag_recoil_knee_loadline_temperature_irreversible_loss_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
    ]
    record.update(
        {
            "recoil_generation": "pm-demag-401",
            "temperature_generation": "pm-demag-400",
            "result_generation": "pm-demag-399",
            "result_temperature_adjusted_remanence_t": -1.0,
            "result_recoil_relative_permeability": -1.05,
            "result_knee_crossed": False,
            "result_irreversible_flux_loss_fraction": -0.2,
            "result_field_orientation": "parallel_h",
            "result_mesh_sha256": "b" * 64,
            "accepted_demag_owner": "pm/old",
            "accepted_demag_result_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pm_demag_uses_temperature_adjusted_recoil_knee_loadline_irreversible_loss_orientation_mesh_owner_and_result"
    ]


def test_v36_public_rejects_self_consistent_bearing_derivative_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
    ]
    stiffness = [[8000.0, 100.0], [100.0, 9000.0]]
    record["stiffness_matrix_n_m"] = stiffness
    record["result_stiffness_matrix_n_m"] = stiffness
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_demag_loadline_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
    ]
    loadline = [0.5, 0.5, 0.5]
    record["loadline_b_t"] = loadline
    record["result_loadline_b_t"] = loadline
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v37_public_axisymmetric_to_3d_force_revolution_volume_energy_direction_owner_mismatch():
    identity = _identity_v37()
    row = identity[
        "axisymmetric_3d_force_revolution_volume_energy_coenergy_direction_displacement_field_mesh_result_identity"
    ]
    row.update(
        {
            "factor_generation": "axisym-revolution-245",
            "direction_generation": "axisym-revolution-244",
            "result_generation": "axisym-revolution-243",
            "result_revolution_factor": math.pi,
            "result_swept_volume_m3": -1.0,
            "result_revolved_energy_j": 2.0,
            "result_revolved_coenergy_j": 2.1,
            "result_virtual_displacement_m": [1.0e-4, 0.0, 0.0],
            "result_force_direction_unit": [-1.0, 0.0, 0.0],
            "result_force_n": [-12.0, 0.0, 0.0],
            "accepted_field_owner": "stale:field",
            "accepted_mesh_sha256": "a" * 64,
            "accepted_force_result_sha256": "b" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_3d_forces_use_current_revolution_volume_energy_coenergy_direction_displacement_field_mesh_and_result"
    ]


def test_v37_public_harmonic_circuit_impedance_voltage_current_complex_power_loss_owner_mismatch():
    identity = _identity_v37()
    row = identity[
        "harmonic_circuit_voltage_current_impedance_complex_power_copper_field_loss_rms_owner_result_identity"
    ]
    row.update(
        {
            "voltage_generation": "harmonic-circuit-245",
            "power_generation": "harmonic-circuit-244",
            "result_generation": "harmonic-circuit-243",
            "result_voltage_peak_phasor_v": [10.0, -2.0],
            "result_current_peak_phasor_a": [-2.0, 1.0],
            "result_impedance_ohm": [-4.0, 1.0],
            "result_complex_power_va": [-9.0, -7.0],
            "result_copper_loss_w": -7.0,
            "result_field_loss_w": 20.0,
            "result_phasor_convention": "rms_sine",
            "accepted_circuit_owner": "stale:circuit",
            "accepted_circuit_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_circuits_use_current_voltage_current_impedance_power_losses_rms_owner_and_result"
    ]


def test_v37_public_rejects_self_consistent_wrong_revolution_factor():
    identity = _identity_v37()
    row = identity[
        "axisymmetric_3d_force_revolution_volume_energy_coenergy_direction_displacement_field_mesh_result_identity"
    ]
    row["revolution_factor"] = row["result_revolution_factor"] = math.pi
    row["swept_volume_m3"] = row["result_swept_volume_m3"] = (
        math.pi * row["centroid_radius_m"] * row["meridional_area_m2"]
    )
    row["revolved_energy_j"] = row["result_revolved_energy_j"] = math.pi * 2.0
    row["revolved_coenergy_j"] = row["result_revolved_coenergy_j"] = math.pi * 2.1
    assert _gate(identity)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_harmonic_power_loss_imbalance():
    identity = _identity_v37()
    row = identity[
        "harmonic_circuit_voltage_current_impedance_complex_power_copper_field_loss_rms_owner_result_identity"
    ]
    row["copper_loss_w"] = row["result_copper_loss_w"] = 1.0
    row["field_loss_w"] = row["result_field_loss_w"] = 1.0
    assert _gate(identity)["status"] == "needs_attention"


def test_v37_public_maglev_dynamic_stiffness_frequency_damping_phase_bias_equilibrium_owner_mismatch():
    summary = _summary_v37()
    row = summary["artifact_identity"]["maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"]
    row.update({"bias_generation": "maglev-dynamic-245", "phase_generation": "maglev-dynamic-244",
                "result_generation": "maglev-dynamic-243", "result_bias_current_a": -5.0,
                "result_equilibrium_gap_m": -0.005, "result_equilibrium_force_n": -100.0,
                "result_supported_load_n": 50.0, "result_excitation_frequency_hz": -100.0,
                "result_complex_stiffness_n_m": [-10000.0, -1.0], "result_viscous_damping_n_s_m": -10.0,
                "result_displacement_phasor_m": [0.0, 1.0e-5], "result_force_phasor_n": [-0.1, -0.1],
                "result_force_displacement_phase_rad": -2.0, "result_coordinate_frame": "rotor_down_left_handed",
                "accepted_maglev_owner": "stale/maglev", "accepted_maglev_result_sha256": "a" * 64})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["maglev_dynamics_close_bias_equilibrium_frequency_stiffness_damping_phase_frame_owner_and_result"]


def test_v37_public_bem_demag_reciprocity_energy_field_magnetization_surface_owner_mismatch():
    summary = _summary_v37()
    row = summary["artifact_identity"]["bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"]
    row.update({"reciprocity_generation": "bem-demag-245", "surface_generation": "bem-demag-244",
                "result_generation": "bem-demag-243", "result_interaction_energy_21_j": 0.02,
                "result_field_1_due_2_a_m": [1000.0, 0.0, 0.0], "result_field_2_due_1_a_m": [0.0, 1000.0, 0.0],
                "result_surface_orientation": "inward_left_handed", "result_region_volumes_m3": [1.0e-5, -1.0e-5],
                "accepted_mesh_owner": "stale/mesh", "accepted_solution_owner": "stale/solution",
                "accepted_demag_result_sha256": "b" * 64})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["bem_demag_closes_reciprocal_energy_field_magnetization_surface_volume_mesh_solution_and_result"]


def test_v37_public_rejects_self_consistent_negative_dynamic_damping():
    summary = _summary_v37()
    row = summary["artifact_identity"]["maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"]
    row["viscous_damping_n_s_m"] = row["result_viscous_damping_n_s_m"] = -10.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_bem_nonreciprocity():
    summary = _summary_v37()
    row = summary["artifact_identity"]["bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"]
    row["interaction_energy_21_j"] = row["result_interaction_energy_21_j"] = 0.02
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v38_public_harmonic_conductor_skin_proximity_impedance_current_loss_poynting_mesh_mismatch():
    identity = _identity_v38()
    row = identity[
        "harmonic_conductor_skin_proximity_impedance_current_voltage_loss_poynting_frequency_mesh_owner_result_identity"
    ]
    row.update(
        {
            "skin_generation": "harmonic-conductor-257",
            "result_frequency_hz": -10_000.0,
            "result_skin_depth_m": -1.0,
            "result_complex_impedance_ohm": [-1.0, 2.0],
            "result_copper_loss_w": -4.5,
            "result_inward_poynting_power_w": 45.0,
            "result_mesh_levels": [3, 2, 1],
            "accepted_field_owner": "stale:harmonic",
        }
    )
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_heat_radiation_convection_emissivity_ambient_flux_temperature_energy_mismatch():
    identity = _identity_v38()
    row = identity[
        "heat_convection_radiation_emissivity_ambient_flux_temperature_geometry_mesh_energy_result_identity"
    ]
    row.update(
        {
            "radiation_generation": "heat-radiation-257",
            "result_emissivity": 1.5,
            "result_ambient_temperature_k": -300.0,
            "result_radiation_flux_w_m2": -1.0,
            "result_geometry_weighting": "axisymmetric",
            "result_boundary_heat_loss_w": -1.0,
            "result_energy_balance_residual_w": 11.0,
            "accepted_mesh_owner": "stale:heat",
        }
    )
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_but_wrong_skin_depth():
    identity = _identity_v38()
    row = identity[
        "harmonic_conductor_skin_proximity_impedance_current_voltage_loss_poynting_frequency_mesh_owner_result_identity"
    ]
    row["skin_depth_m"] *= 2.0
    row["result_skin_depth_m"] = row["skin_depth_m"]
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_but_wrong_radiation_flux():
    identity = _identity_v38()
    row = identity[
        "heat_convection_radiation_emissivity_ambient_flux_temperature_geometry_mesh_energy_result_identity"
    ]
    row["radiation_flux_w_m2"] *= 0.5
    row["result_radiation_flux_w_m2"] = row["radiation_flux_w_m2"]
    assert _gate(identity)["status"] == "needs_attention"


def test_v38_public_eddy_current_maglev_plate_velocity_skin_depth_lift_drag_loss_power_mismatch():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
    ]
    row.update(
        {
            "skin_generation": "eddy-maglev-257",
            "power_generation": "eddy-maglev-256",
            "result_generation": "eddy-maglev-255",
            "result_plate_velocity_m_s": -20.0,
            "result_excitation_frequency_hz": -200.0,
            "result_plate_conductivity_s_m": -3.5e7,
            "result_skin_depth_m": -1.0,
            "result_lift_force_n": -100.0,
            "result_drag_force_n": -20.0,
            "result_joule_loss_w": -400.0,
            "result_mechanical_drag_power_w": 40.0,
            "result_power_balance_residual_w": 440.0,
            "accepted_mesh_owner": "stale:maglev",
            "accepted_maglev_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "eddy_maglev_closes_velocity_frequency_skin_depth_lift_drag_joule_power_mesh_owner_and_result"
    ]


def test_v38_public_pm_coupling_torque_angle_periodicity_energy_derivative_action_reaction_mismatch():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
    ]
    row.update(
        {
            "periodicity_generation": "pm-coupling-257",
            "reaction_generation": "pm-coupling-256",
            "result_generation": "pm-coupling-255",
            "result_relative_angle_rad": -1.0,
            "result_pole_pairs": 0,
            "result_pole_period_rad": -1.0,
            "result_periodic_energy_j": 9.0,
            "result_energy_derivative_torque_nm": -9.0,
            "result_driver_torque_nm": 5.0,
            "result_driven_torque_nm": 5.0,
            "result_torque_frame": "left_handed_local",
            "accepted_mesh_owner": "stale:coupling",
            "accepted_coupling_result_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pm_coupling_closes_pole_periodic_energy_derivative_action_reaction_frame_mesh_owner_and_result"
    ]


def test_v38_public_rejects_self_consistent_wrong_maglev_skin_depth():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "eddy_current_maglev_plate_velocity_frequency_conductivity_skin_depth_lift_drag_loss_power_mesh_owner_result_identity"
    ]
    row["skin_depth_m"] *= 2.0
    row["result_skin_depth_m"] = row["skin_depth_m"]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_same_sign_coupling_torque():
    summary = _summary_v38()
    row = summary["artifact_identity"][
        "pm_coupling_angle_pole_periodicity_energy_derivative_driver_driven_torque_action_reaction_frame_mesh_owner_result_identity"
    ]
    row["driven_torque_nm"] = row["driver_torque_nm"]
    row["result_driven_torque_nm"] = row["driven_torque_nm"]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v39_public_nonlinear_magnetic_circuit_bh_coenergy_incremental_inductance_force_mismatch():
    identity = _identity_v39()
    identity[_NONLINEAR].update({"bh_generation": "nonlinear-magnetic-circuit-270", "force_generation": "nonlinear-magnetic-circuit-269", "result_generation": "nonlinear-magnetic-circuit-268", "result_operating_h_a_per_m": 1000.0, "result_operating_b_t": -1.2, "result_flux_linkage_wb_turn": -2.0e-2, "result_magnetic_energy_j": -8.0e-2, "result_coenergy_j": -1.2e-1, "result_current_increment_a": -5.0e-1, "result_flux_linkage_increment_wb_turn": -1.0e-3, "result_incremental_inductance_h": -2.0e-3, "result_coenergy_increment_j": -2.0e-2, "result_virtual_work_force_n": -20.0, "result_terminal_power_w": -50.0, "accepted_mesh_owner": "stale:mesh", "accepted_nonlinear_result_sha256": "a" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["nonlinear_magnetic_circuits_close_bh_flux_coenergy_incremental_inductance_force_power_mesh_and_result"]


def test_v39_public_electrostatic_capacitance_matrix_charge_energy_reciprocity_gauge_mismatch():
    identity = _identity_v39()
    identity[_ELECTROSTATIC].update({"matrix_generation": "electrostatic-capacitance-270", "gauge_generation": "electrostatic-capacitance-269", "result_generation": "electrostatic-capacitance-268", "result_capacitance_matrix_f": [[2.0e-11, 1.0e-11], [-2.0e-11, 2.0e-11]], "result_terminal_charge_c": [1.0e-9, 1.0e-9], "result_field_energy_j": -1.0e-7, "result_reciprocity_residual_f": 1.0e-11, "result_reference_gauge": "floating", "accepted_conductor_owner": "stale:conductors", "accepted_mesh_owner": "stale:mesh", "accepted_electrostatic_result_sha256": "b" * 64})
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["electrostatic_capacitance_closes_matrix_charge_energy_reciprocity_gauge_owners_and_result"]


def test_v39_public_rejects_self_consistent_wrong_coenergy():
    identity = _identity_v39()
    identity[_NONLINEAR]["coenergy_j"] = 0.2
    identity[_NONLINEAR]["result_coenergy_j"] = 0.2
    assert _gate(identity)["status"] == "needs_attention"


def test_v39_public_accepts_nonlinear_energy_coenergy_partition():
    identity = _identity_v39()
    row = identity[_NONLINEAR]
    assert row["magnetic_energy_j"] != row["coenergy_j"]
    assert row["coenergy_j"] != 0.5 * row["terminal_current_a"] * row["flux_linkage_wb_turn"]
    assert _gate(identity)["status"] == "ok"


def test_v39_public_rejects_self_consistent_nonsymmetric_capacitance():
    identity = _identity_v39()
    matrix = [[2.0e-11, -1.0e-11], [-2.0e-11, 2.0e-11]]
    identity[_ELECTROSTATIC]["capacitance_matrix_f"] = matrix
    identity[_ELECTROSTATIC]["result_capacitance_matrix_f"] = matrix
    assert _gate(identity)["status"] == "needs_attention"


def test_v39_public_thin_conductor_eddy_surface_impedance_skin_current_complex_power_mismatch() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_THIN_KEY]
    row.update(
        {
            "impedance_generation": "thin-conductor-270",
            "power_generation": "thin-conductor-269",
            "result_generation": "thin-conductor-268",
            "result_skin_depth_m": -1.0,
            "result_surface_impedance_ohm": [-1.0, 1.0],
            "result_sheet_current_peak_a_m": [-100.0, 0.0],
            "result_tangential_field_jump_peak_a_m": [0.0, 0.0],
            "result_joule_loss_w": -1.0,
            "result_reactive_power_var": -1.0,
            "accepted_surface_owner": "stale:surface",
            "accepted_thin_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thin_conductors_close_surface_impedance_skin_current_field_jump_complex_power_owner_and_result"
    ]


def test_v39_public_magnetic_gear_harmonic_polepair_phase_torque_actionreaction_mismatch() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_GEAR_KEY]
    row.update(
        {
            "pole_generation": "magnetic-gear-270",
            "reaction_generation": "magnetic-gear-269",
            "result_generation": "magnetic-gear-268",
            "result_modulator_segment_count": 5,
            "result_working_harmonic_order": 3,
            "result_mechanical_phase_rad": -1.0,
            "result_gear_ratio": 2.0,
            "result_low_speed_rad_s": 20.0,
            "result_low_speed_torque_nm": -10.0,
            "result_modulator_reaction_torque_nm": 0.0,
            "result_power_balance_residual_w": 400.0,
            "accepted_model_owner": "stale:gear",
            "accepted_gear_result_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_gears_close_harmonics_poles_ratio_torque_reaction_power_owner_and_result"
    ]


def test_v39_public_rejects_self_consistent_wrong_surface_impedance() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_THIN_KEY]
    row["surface_impedance_ohm"] = [2.0 * value for value in row["surface_impedance_ohm"]]
    row["result_surface_impedance_ohm"] = list(row["surface_impedance_ohm"])
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_wrong_gear_direction() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_GEAR_KEY]
    row["low_speed_rad_s"] = 20.0
    row["result_low_speed_rad_s"] = 20.0
    row["gear_ratio"] = 2.0
    row["result_gear_ratio"] = 2.0
    row["power_balance_residual_w"] = 400.0
    row["result_power_balance_residual_w"] = 400.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_heat_contact_mismatch():
    identity = _identity_v40()
    identity[_HEAT].update(
        {
            "contact_generation": "heat-contact-310",
            "energy_generation": "heat-contact-309",
            "result_generation": "heat-contact-308",
            "result_conductivity_w_per_m_k": -15.0,
            "result_contact_resistance_m2_k_per_w": -1.0e-4,
            "result_convection_coefficient_w_per_m2_k": 0.0,
            "result_boundary_heat_flux_w_per_m2": -1.0,
            "result_interface_temperature_jump_k": -1.0,
            "result_convection_surface_temperature_k": 400.0,
            "result_total_heat_rate_w": -20.0,
            "result_energy_balance_residual_w": 5.0,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_heat_result_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "heat_contact_models_close_conduction_contact_convection_flux_temperature_energy_mesh_and_result"
    ]


def test_v40_public_current_flow_mismatch():
    identity = _identity_v40()
    identity[_CURRENT].update(
        {
            "electrode_generation": "current-flow-310",
            "power_generation": "current-flow-309",
            "result_generation": "current-flow-308",
            "result_electrode_voltage_v": [0.0, 10.0],
            "result_terminal_current_a": [2.0, 2.0],
            "result_effective_resistance_ohm": -5.0,
            "result_joule_loss_w": -20.0,
            "result_terminal_power_w": 0.0,
            "result_reciprocity_residual_ohm": 1.0,
            "accepted_conductor_owner": "stale:conductor",
            "accepted_mesh_owner": "stale:mesh",
            "accepted_current_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "current_flow_models_close_electrodes_currents_resistance_joule_terminal_power_reciprocity_owners_and_result"
    ]


def test_v40_public_rejects_self_consistent_wrong_contact_jump():
    identity = _identity_v40()
    identity[_HEAT]["interface_temperature_jump_k"] = 1.0
    identity[_HEAT]["result_interface_temperature_jump_k"] = 1.0
    assert _gate(identity)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_current_nonconservation():
    identity = _identity_v40()
    identity[_CURRENT]["terminal_current_a"] = [2.0, -1.5]
    identity[_CURRENT]["result_terminal_current_a"] = [2.0, -1.5]
    assert _gate(identity)["status"] == "needs_attention"


def test_v40_public_multilayer_magnetic_shield_permeability_thickness_attenuation_flux_energy_mismatch() -> None:
    summary = _summary_v40()
    summary["artifact_identity"][_SHIELD].update({"material_generation": "multilayer-shield-279", "result_attenuation_factor": 0.5, "result_interface_normal_flux_t": [1.0, 2.0], "result_stored_energy_j": -1.0, "accepted_geometry_owner": "stale:geometry"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_transformer_leakage_mutual_inductance_fluxlinkage_reciprocity_coenergy_force_mismatch() -> None:
    summary = _summary_v40()
    summary["artifact_identity"][_TRANSFORMER].update({"inductance_generation": "transformer-coupling-279", "result_inductance_matrix_h": [[-1.0, 2.0], [3.0, -1.0]], "result_coenergy_j": -1.0, "result_force_n": -1.0, "accepted_winding_owner": "stale:winding"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_wrong_shield_attenuation() -> None:
    summary = _summary_v40()
    row = summary["artifact_identity"][_SHIELD]
    row["attenuation_factor"] = row["result_attenuation_factor"] = 2.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_nonreciprocal_transformer() -> None:
    summary = _summary_v40()
    row = summary["artifact_identity"][_TRANSFORMER]
    row["inductance_matrix_h"] = row["result_inductance_matrix_h"] = [[12.0e-3, 5.4e-3], [4.0e-3, 3.0e-3]]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_permanent_magnet_mismatch():
    identity = _identity_v41()
    identity[_MAGNET].update(
        {
            "recoil_generation": "pm-loadline-723",
            "force_generation": "pm-loadline-722",
            "result_generation": "pm-loadline-721",
            "result_recoil_relative_permeability": -1.05,
            "result_loadline_permeance_coefficient": -2.0,
            "result_operating_h_a_per_m": 1.0e6,
            "result_operating_b_t": -1.0,
            "result_demag_margin_a_per_m": -1.0,
            "result_field_energy_j": -1.0,
            "result_virtual_work_force_n": 10.0,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_magnet_result_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "permanent_magnets_close_recoil_loadline_operating_demag_energy_virtualwork_mesh_and_result"
    ]


def test_v41_public_capacitance_matrix_mismatch():
    identity = _identity_v41()
    identity[_CAPACITANCE].update(
        {
            "matrix_generation": "capacitance-matrix-723",
            "energy_generation": "capacitance-matrix-722",
            "result_generation": "capacitance-matrix-721",
            "result_conductor_names": ["electrode_2", "electrode_1"],
            "result_capacitance_matrix_f": [[-2.0e-9, 1.0e-9], [0.0, 1.5e-9]],
            "result_conductor_charge_c": [0.0, 0.0],
            "result_stored_energy_j": -1.0,
            "result_symmetry_residual_f": 1.0,
            "result_reciprocity_residual_f": 1.0,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_capacitance_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrostatic_capacitance_matrices_close_symmetry_psd_charge_energy_reciprocity_mesh_and_result"
    ]


def test_v41_public_rejects_self_consistent_wrong_loadline_point():
    identity = _identity_v41()
    identity[_MAGNET]["operating_b_t"] = 0.6
    identity[_MAGNET]["result_operating_b_t"] = 0.6
    assert _gate(identity)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_indefinite_capacitance_matrix():
    identity = _identity_v41()
    matrix = [[1.0e-9, 2.0e-9], [2.0e-9, 1.0e-9]]
    charges = [10.0e-9, 20.0e-9]
    identity[_CAPACITANCE]["capacitance_matrix_f"] = matrix
    identity[_CAPACITANCE]["result_capacitance_matrix_f"] = matrix
    identity[_CAPACITANCE]["conductor_charge_c"] = charges
    identity[_CAPACITANCE]["result_conductor_charge_c"] = charges
    identity[_CAPACITANCE]["stored_energy_j"] = 0.5e-7
    identity[_CAPACITANCE]["result_stored_energy_j"] = 0.5e-7
    assert _gate(identity)["status"] == "needs_attention"


def test_v41_public_maglev_equilibrium_forcegradient_stiffness_energy_stability_gap_mismatch() -> None:
    summary = _summary_v41()
    summary["artifact_identity"][_MAGLEV].update({"force_generation": "maglev-equilibrium-723", "result_air_gap_m": -0.01, "result_stiffness_n_per_m": -200.0, "result_stability": "unstable", "accepted_geometry_owner": "stale:geometry"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_eddyshield_skin_depth_phase_lag_loss_frequency_energy_mismatch() -> None:
    summary = _summary_v41()
    summary["artifact_identity"][_EDDY].update({"frequency_generation": "eddy-shield-723", "result_skin_depth_m": -1.0, "result_attenuation_factor": 2.0, "result_eddy_loss_w": -1.0, "accepted_geometry_owner": "stale:geometry"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_maglev_energy() -> None:
    summary = _summary_v41()
    row = summary["artifact_identity"][_MAGLEV]
    row["potential_energy_j"] = row["result_potential_energy_j"] = [0.0, 1.0, 0.0]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_eddy_skin_depth() -> None:
    summary = _summary_v41()
    row = summary["artifact_identity"][_EDDY]
    row["skin_depth_m"] = row["result_skin_depth_m"] = 1.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_weighted_stress_force_mismatch():
    identity = _identity_v42()
    identity[_FORCE].update(
        {
            "airgap_contour_generation": "force-closure-841",
            "mesh_generation": "force-closure-840",
            "result_generation": "force-closure-839",
            "result_weighted_stress_force_n": [-12.0, 3.0],
            "result_airgap_contour_forces_n": [[30.0, 0.0]],
            "result_contour_independence_relative_spread": 2.0,
            "result_energy_plus_j": 1.0012,
            "result_virtual_work_direction": [-1.0, 0.0],
            "result_virtual_work_force_n": -12.0,
            "result_mesh_refinement_force_samples_n": [[1.0, 1.0], [-12.0, 3.0]],
            "accepted_field_owner": "stale:field",
            "accepted_mesh_owner": "stale:mesh",
            "accepted_force_result_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "weighted_stress_forces_close_airgap_contours_virtual_work_direction_mesh_field_owner_and_result"
    ]


def test_v42_public_harmonic_conductor_mismatch():
    identity = _identity_v42()
    identity[_CONDUCTOR].update(
        {
            "skin_depth_generation": "harmonic-conductor-841",
            "power_generation": "harmonic-conductor-840",
            "result_generation": "harmonic-conductor-839",
            "result_skin_depth_m": 0.1,
            "result_complex_current_density_a_per_m2": [[-1.0e6, 2.0e5]],
            "result_integrated_abs_current_density_sq_a2_per_m": -1.0,
            "result_joule_loss_w": -20.8,
            "result_terminal_impedance_ohm": [-0.4, -0.2],
            "result_complex_power_va": [-20.8, 10.4],
            "accepted_mesh_owner": "stale:mesh",
            "accepted_conductor_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_conductors_close_skin_depth_complex_current_joule_impedance_power_mesh_and_result"
    ]


def test_v42_public_rejects_self_consistent_wrong_contour_spread():
    identity = _identity_v42()
    record = identity[_FORCE]
    record["airgap_contour_forces_n"] = [[30.0, 0.0]] * 3
    record["result_airgap_contour_forces_n"] = [[30.0, 0.0]] * 3
    spread = math.hypot(18.0, 3.0) / math.hypot(12.0, -3.0)
    record["contour_independence_relative_spread"] = spread
    record["result_contour_independence_relative_spread"] = spread
    assert _gate(identity)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_complex_power():
    identity = _identity_v42()
    record = identity[_CONDUCTOR]
    record["complex_power_va"] = [20.8, -10.4]
    record["result_complex_power_va"] = [20.8, -10.4]
    assert _gate(identity)["status"] == "needs_attention"


def test_v42_public_demag_tensor_mismatch() -> None:
    summary = _summary_v42()
    summary["artifact_identity"][_DEMAG].update(
        {
            "symmetry_generation": "demag-tensor-841",
            "result_tensor_trace": 1.2,
            "result_tensor_eigenvalues": [-0.1, 0.3, 1.0],
            "result_reciprocity_right": -1.0,
            "result_demag_energy_j": -1.0,
            "accepted_mesh_owner": "stale:mesh",
        }
    )
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_linear_motor_mismatch() -> None:
    summary = _summary_v42()
    summary["artifact_identity"][_LINEAR].update(
        {
            "position_generation": "linear-motor-period-841",
            "result_position_m": [0.04, 0.03, 0.02, 0.01, 0.0],
            "result_phase_currents_a": [[5.0, 5.0, 5.0]],
            "result_periodic_work_j": 1.0,
            "accepted_mesh_owner": "stale:mesh",
        }
    )
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_demag_trace() -> None:
    summary = _summary_v42()
    row = summary["artifact_identity"][_DEMAG]
    tensor = [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.6]]
    row["demag_tensor"] = row["result_demag_tensor"] = tensor
    row["tensor_trace"] = row["result_tensor_trace"] = 1.1
    row["tensor_eigenvalues"] = row["result_tensor_eigenvalues"] = [0.2, 0.3, 0.6]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_cogging_force() -> None:
    summary = _summary_v42()
    row = summary["artifact_identity"][_LINEAR]
    row["cogging_force_n"] = row["result_cogging_force_n"] = [0.0] * 5
    row["thrust_n"] = row["result_thrust_n"] = [100.0] * 5
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v43_public_positive_contracts():
    assert _gate(_identity_v43())["status"] == "ok"


def test_v43_public_solenoid_mismatch():
    identity = _identity_v43()
    identity[_SOLENOID]["result_inductance_h"] = 0.08
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["axisymmetric_solenoids_close_flux_inductance_coenergy_force_two_pi_r_mesh_and_result"]


def test_v43_public_dielectric_mismatch():
    identity = _identity_v43()
    identity[_DIELECTRIC]["result_stored_energy_j"] = -1.0
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["dielectric_interfaces_close_capacitance_charge_flux_energy_reciprocity_mesh_and_result"]


def test_v43_public_rejects_self_consistent_wrong_two_pi_r_factor():
    identity = _identity_v43()
    record = identity[_SOLENOID]
    record["two_pi_r_factor_m"] = record["result_two_pi_r_factor_m"] = 1.0
    assert _gate(identity)["status"] == "needs_attention"


def test_v43_public_positive_bearing_and_hysteresis_closure() -> None:
    assert magnetic_force_method_profile_gate(_summary_v43())["status"] == "ok"


def test_v43_public_bearing_mismatch() -> None:
    summary = _summary_v43()
    summary["artifact_identity"][_BEARING]["result_stability_sign"] = "unstable"
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["magnetic_bearings_close_stiffness_symmetry_crosscoupling_force_energy_stability_gap_mesh_and_result"]


def test_v43_public_hysteresis_mismatch() -> None:
    summary = _summary_v43()
    summary["artifact_identity"][_HYSTERESIS]["result_remanence_a_per_m"] = -0.8
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["hysteresis_minorloops_close_field_path_branch_remanence_coercivity_loss_energy_material_mesh_and_result"]


def test_v44_public_dynamic_and_demag_identity_positive() -> None:
    result = validate_v44_identity(_identity_v44())
    assert result == {
        "magnetic_force_v44_dynamic_bearing_identity": True,
        "magnetic_force_v44_demag_minorloop_identity": True,
    }


def test_v44_public_identity_rejects_phase_and_temperature_mutations() -> None:
    identity = _identity_v44()
    identity[_DYNAMIC]["result_phase_deg"] = [0.0, -10.0, -20.0]
    identity[_DEMAG_V44]["result_temperature_k"] = 350.0
    result = validate_v44_identity(identity)
    assert result["magnetic_force_v44_dynamic_bearing_identity"] is False
    assert result["magnetic_force_v44_demag_minorloop_identity"] is False


def test_v47_positive_dual_lane_and_force_artifacts_are_accepted() -> None:
    assert all(validate_public_v47_identity(_identity_v47()).values())


def test_v47_dual_lane_shared_physics_mutation_is_rejected() -> None:
    identity = _identity_v47()
    identity[MOTOR]["lane_b_geometry_identity_sha256"] = "a" * 64
    identity[MOTOR]["lane_b_material_identity"] = "material:other"
    identity[MOTOR]["lane_b_excitation_identity"] = {
        "phase_order": ["A", "C", "B"],
        "current_a": [10.0, -5.0, -5.0],
    }
    identity[MOTOR]["lane_b_operating_point_key"] = "speed=6000rpm,current=5A"
    assert not all(validate_public_v47_identity(identity).values())


def test_v47_force_pair_body_aggregation_mutation_is_rejected() -> None:
    identity = _identity_v47()
    identity[FORCE]["result_displacement_pair_m"] = [0.001, 0.0]
    identity[FORCE]["result_body_owner"] = "body:fixed"
    identity[FORCE]["result_component_force_n"] = {"core": 80.0}
    identity[FORCE]["result_aggregated_force_n"] = 80.0
    assert not all(validate_public_v47_identity(identity).values())


def test_v48_positive_bem_and_hysteresis_artifacts_are_accepted() -> None:
    assert all(validate_v48_identity(_identity_v48()).values())


def test_v48_bem_mesh_and_quadrature_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v48())
    identity[V48_BEM]["result_near_quadrature_order"] = [8, 4, 8, 8]
    identity[V48_BEM]["result_mesh_revision"] = "mesh:old"
    assert validate_v48_identity(identity)["bem_v48_panel_quadrature_normal_solid_angle_mesh"] is False


def test_v48_hysteresis_history_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v48())
    identity[HYSTERESIS]["result_internal_state"] = {"branch": "ascending", "last_reversal": 0.3, "memory_depth": 1}
    identity[HYSTERESIS]["result_material_owner"] = "material:old"
    assert validate_v48_identity(identity)["hysteresis_v48_return_state_environment_material_owner"] is False


def test_v49_positive_demag_and_virtual_work_artifacts_are_accepted() -> None:
    assert all(validate_v49_identity(_identity_v49()).values())


def test_v49_demag_branch_temperature_step_digest_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v49())
    identity[V49_DEMAG]["result_recoil_branch_t_a_per_m"] = list(reversed(identity[V49_DEMAG]["recoil_branch_t_a_per_m"]))
    identity[V49_DEMAG]["result_temperature_c"] = 20.0
    identity[V49_DEMAG]["accepted_result_sha256"] = "8" * 64
    identity[V49_DEMAG]["result_magnet_owner"] = "magnet:old"
    assert validate_v49_identity(identity)["demag_v49_recoil_temperature_loadstep_digest_owner"] is False


def test_v49_virtual_work_frame_mesh_energy_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity_v49())
    identity[VIRTUAL_WORK]["result_displacement_frame"] = "frame:local-y"
    identity[VIRTUAL_WORK]["result_mesh_state_sha256"] = "9" * 64
    identity[VIRTUAL_WORK]["result_energy_j"] = list(reversed(identity[VIRTUAL_WORK]["energy_j"]))
    identity[VIRTUAL_WORK]["result_force_owner"] = "force:old"
    assert validate_v49_identity(identity)["force_v49_virtual_work_displacement_frame_mesh_energy_owner"] is False


def test_v50_positive_bem_and_motional_emf_artifacts_are_accepted() -> None:
    assert all(validate_v50_identity(_identity_v50()).values())


def test_v50_bem_quadrature_self_panel_nearfield_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity_v50())
    identity[V50_BEM]["result_singular_quadrature"] = "gauss-order-2"
    identity[V50_BEM]["result_self_panel_treatment"] = "centroid-sample"
    identity[V50_BEM]["result_nearfield_regularization"] = {"distance_ratio": 1.0, "regularization": "none", "max_depth": 0}
    identity[V50_BEM]["result_mesh_owner"] = "mesh:foreign"
    assert validate_v50_identity(identity)["magnetic_force_v50_bem_singular_self_nearfield_mesh_owner"] is False


def test_v50_motion_emf_velocity_frame_path_direction_and_owner_drift_is_rejected() -> None:
    identity = deepcopy(_identity_v50())
    identity[MOTION]["result_velocity_m_s"] = [-12.0, 0.0, 0.0]
    identity[MOTION]["result_velocity_frame"] = "frame:body"
    identity[MOTION]["result_conductor_path_m"] = list(reversed(identity[MOTION]["conductor_path_m"]))
    identity[MOTION]["result_integration_direction"] = "path-reverse"
    identity[MOTION]["result_emf_owner"] = "emf:foreign"
    assert validate_v50_identity(identity)["magnetic_force_v50_motion_emf_velocity_frame_path_direction_owner"] is False


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_v51_identity(_identity_v51()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity_v51())
    identity[POTENTIAL].update({"result_gauge": "tree_cotree", "result_solution_owner": "solution:stale"})
    identity[V51_BEM].update({"result_panel_orientation": "inward", "result_matrix_owner": "matrix:stale"})
    assert not all(validate_v51_identity(identity).values())


def test_v51_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity_v51())
    identity[POTENTIAL]["gauge"] = identity[POTENTIAL]["result_gauge"] = "tree_cotree"
    identity[V51_BEM]["panel_orientation"] = identity[V51_BEM]["result_panel_orientation"] = "inward"
    assert not all(validate_v51_identity(identity).values())


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_v52_identity(_identity_v52()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity_v52())
    identity[VIRTUAL_FORCE]["result_force_n"] = [-20.0, 0.0, 0.0]
    identity[MAGNET_TORQUE]["result_angles_unwrapped_deg"] = [350.0, 355.0, 0.0, 5.0, 10.0]
    assert not all(validate_v52_identity(identity).values())


def test_v52_self_consistent_wrong_derivative_sign_is_rejected() -> None:
    identity = deepcopy(_identity_v52())
    identity[VIRTUAL_FORCE]["force_sign_convention"] = identity[VIRTUAL_FORCE]["result_force_sign_convention"] = "positive_energy_gradient"
    identity[MAGNET_TORQUE]["torque_at_center_nm"] = identity[MAGNET_TORQUE]["result_torque_at_center_nm"] = -identity[MAGNET_TORQUE]["torque_at_center_nm"]
    assert not all(validate_v52_identity(identity).values())


def test_v53_positive_public_artifacts_are_accepted():
    assert all(validate_v53_identity(_identity_v53()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    identity = deepcopy(_identity_v53())
    identity[QUADRATURE]["result_panel_owner"] = "panel-set:stale"
    identity[MAGLEV]["result_stiffness_n_per_m"] = 1200.0
    assert not all(validate_v53_identity(identity).values())


def test_v53_self_consistent_wrong_physics_is_rejected():
    identity = deepcopy(_identity_v53())
    identity[QUADRATURE]["panel_interactions"][0]["classification"] = identity[QUADRATURE]["result_panel_interactions"][0]["classification"] = "far"
    identity[MAGLEV]["stiffness_n_per_m"] = identity[MAGLEV]["result_stiffness_n_per_m"] = 1200.0
    assert not all(validate_v53_identity(identity).values())


def test_v54_positive_identities_are_accepted():
    assert all(validate_v54_identity(_payload_v54()).values())


def test_v54_frozen_mutations_are_rejected():
    payload = deepcopy(_payload_v54())
    payload[CHARGE]["result_surface_charges"] = [{"panel": 11, "magnetic_charge_a_m": 0.25}]
    payload[STIFFNESS]["result_coordinate_direction"] = [1.0, 0.0, 0.0]
    assert not all(validate_v54_identity(payload).values())


def test_v54_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload_v54())
    bad_charges = [{"panel": 11, "magnetic_charge_a_m": 0.25}, {"panel": 12, "magnetic_charge_a_m": 0.1}]
    payload[CHARGE]["surface_charges"] = payload[CHARGE]["result_surface_charges"] = bad_charges
    payload[STIFFNESS]["force_gradient_n_per_m"] = payload[STIFFNESS]["result_force_gradient_n_per_m"] = 1200.0
    assert not all(validate_v54_identity(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload_v54())
    payload[CHARGE]["boundary_orientation"] = {"11": [1], "12": -1}
    payload[STIFFNESS]["coordinate_direction"] = [[0.0], 0.0, 1.0]
    assert not all(validate_v54_identity(payload).values())


def test_v55_positive_identities_are_accepted():
    assert all(validate_v55_identity(_payload_v55()).values())


def test_v55_frozen_mutations_are_rejected():
    payload = deepcopy(_payload_v55()); payload[V55_DEMAG]["result_solution_owner"] = "solution:stale"; payload[BEARING]["result_body_owner"] = "body:stale"
    assert not all(validate_v55_identity(payload).values())


def test_v55_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload_v55()); payload[V55_DEMAG]["demag_energy_j"] = payload[V55_DEMAG]["result_demag_energy_j"] = -1.0
    payload[BEARING]["cross_stiffness_n_per_m"] = payload[BEARING]["result_cross_stiffness_n_per_m"] = [[1200.0, 500.0], [-50.0, 900.0]]
    assert not all(validate_v55_identity(payload).values())


def test_v55_malformed_values_reject_without_raising():
    payload = deepcopy(_payload_v55()); payload[V55_DEMAG]["hb_samples"] = [{"h_a_per_m": [0.0], "b_t": 1.2}]; payload[BEARING]["coordinate_basis"] = {"x": [[1.0], 0.0, 0.0]}
    assert not all(validate_v55_identity(payload).values())


def test_v55_numeric_sha256_values_are_rejected():
    payload = _payload_v55()
    numeric_digest = int("9" * 64)
    for row in payload.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert not all(validate_v55_identity(payload).values())
