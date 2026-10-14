"""Payload builders for the multiphysics generalization identity tests.

Builders are data, not tests: each ``_with_v*`` adds identity blocks to the
rotational eddy-brake summary chain, and each ``_records_v*`` returns the
replay records for one public identity validator.
"""
from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.adjoint_weakform_identity_v52 import (
    ADJOINT,
    WEAK_FORM,
)
from radia_mcp.radia_ngsolve.conservation_identity_v54 import (
    PIEZO,
    SPECIES,
)
from radia_mcp.radia_ngsolve.coupled_periodic_identity_v53 import (
    ELECTROTHERMAL,
    FLOQUET,
)
from radia_mcp.radia_ngsolve.dissipation_reaction_identity_v55 import (
    ELECTROCHEM,
    THERMO,
)
from radia_mcp.radia_ngsolve.multiphysics_identity_v56 import (
    INDUCTION,
    MODAL,
)
from test_rotational_eddy_brake_energy_gate import (
    _with_v28_thermoelastic_and_field_circuit_identity,
)


def _with_v29_sliding_and_radiation_identity(summary: dict) -> dict:
    summary = _with_v28_thermoelastic_and_field_circuit_identity(summary)
    generation = "sliding-interface-161"
    summary[
        "rotating_sliding_interface_sector_pitch_azimuth_interpolation_frame_periodicity_mesh_torque_generation_identity"
    ] = {
        "sliding_generation": generation,
        "sector_sliding_generation": generation,
        "azimuth_sliding_generation": generation,
        "interpolation_sliding_generation": generation,
        "frame_sliding_generation": generation,
        "periodicity_sliding_generation": generation,
        "mesh_sliding_generation": generation,
        "result_sliding_generation": generation,
        "sector_pitch_deg": 30.0,
        "result_sector_pitch_deg": 30.0,
        "sector_count": 12,
        "result_sector_count": 12,
        "azimuth_origin_deg": 0.0,
        "result_azimuth_origin_deg": 0.0,
        "source_interface_tag": "rotor_if",
        "result_source_interface_tag": "rotor_if",
        "target_interface_tag": "stator_if",
        "result_target_interface_tag": "stator_if",
        "interpolation": "conservative_mortar_azimuth",
        "result_interpolation": "conservative_mortar_azimuth",
        "rotor_frame": "rotor_cylindrical",
        "result_rotor_frame": "rotor_cylindrical",
        "periodic_phase_deg": 0.0,
        "result_periodic_phase_deg": 0.0,
        "azimuth_samples_deg": [0.0, 7.5, 15.0, 22.5, 30.0],
        "result_azimuth_samples_deg": [0.0, 7.5, 15.0, 22.5, 30.0],
        "torque_nm": [10.0, 10.5, 10.0, 9.5, 10.0],
        "result_torque_nm": [10.0, 10.5, 10.0, 9.5, 10.0],
        "sliding_mesh_sha256": "1" * 64,
        "result_sliding_mesh_sha256": "1" * 64,
        "torque_result_sha256": "2" * 64,
        "accepted_torque_result_sha256": "2" * 64,
    }
    generation = "acoustic-radiation-161"
    summary[
        "acoustic_radiation_impedance_modal_trace_reference_area_pressure_velocity_power_frequency_result_generation_identity"
    ] = {
        "radiation_generation": generation,
        "modal_radiation_generation": generation,
        "trace_radiation_generation": generation,
        "area_radiation_generation": generation,
        "convention_radiation_generation": generation,
        "power_radiation_generation": generation,
        "frequency_radiation_generation": generation,
        "result_radiation_generation": generation,
        "modal_basis_id": "p1_interface_modes_mass_normalized",
        "result_modal_basis_id": "p1_interface_modes_mass_normalized",
        "mode_indices": [1, 2, 3],
        "result_mode_indices": [1, 2, 3],
        "trace_projection": "l2_p1_boundary",
        "result_trace_projection": "l2_p1_boundary",
        "reference_area_m2": 0.1,
        "result_reference_area_m2": 0.1,
        "pressure_velocity_convention": "outward_positive_velocity",
        "result_pressure_velocity_convention": "outward_positive_velocity",
        "frequency_grid_hz": [100.0, 200.0, 500.0],
        "result_frequency_grid_hz": [100.0, 200.0, 500.0],
        "radiation_impedance_ri": [[20.0, 5.0], [30.0, 8.0], [50.0, 12.0]],
        "result_radiation_impedance_ri": [[20.0, 5.0], [30.0, 8.0], [50.0, 12.0]],
        "outward_power_flux_w": [1.0, 1.5, 2.0],
        "result_outward_power_flux_w": [1.0, 1.5, 2.0],
        "radiation_mesh_sha256": "3" * 64,
        "result_radiation_mesh_sha256": "3" * 64,
        "radiation_result_sha256": "4" * 64,
        "accepted_radiation_result_sha256": "4" * 64,
    }
    return summary


def _with_v30_joule_and_eigenmode_identity(summary: dict) -> dict:
    summary = _with_v29_sliding_and_radiation_identity(summary)
    generation = "joule-heat-closure-171"
    summary[
        "joule_heat_source_current_density_resistivity_temperature_frame_time_average_energy_mesh_result_generation_identity"
    ] = {
        "joule_generation": generation,
        "mapping_joule_generation": generation,
        "resistivity_joule_generation": generation,
        "temperature_joule_generation": generation,
        "frame_joule_generation": generation,
        "averaging_joule_generation": generation,
        "energy_joule_generation": generation,
        "mesh_joule_generation": generation,
        "result_joule_generation": generation,
        "current_density_field_id": "ec.J-current-171",
        "result_current_density_field_id": "ec.J-current-171",
        "temperature_field_id": "ht.T-current-171",
        "result_temperature_field_id": "ht.T-current-171",
        "resistivity_model_id": "rho(T)-copper-171",
        "result_resistivity_model_id": "rho(T)-copper-171",
        "source_frame": "material_spatial",
        "result_source_frame": "material_spatial",
        "averaging_window_s": [0.02, 0.04],
        "result_averaging_window_s": [0.02, 0.04],
        "time_average_method": "trapezoidal_period_average",
        "result_time_average_method": "trapezoidal_period_average",
        "electric_loss_w": 12.5,
        "result_electric_loss_w": 12.5,
        "heat_source_integral_w": 12.5,
        "result_heat_source_integral_w": 12.5,
        "energy_balance_relative_tolerance": 1.0e-9,
        "coupled_mesh_sha256": "1" * 64,
        "result_coupled_mesh_sha256": "1" * 64,
        "joule_heat_result_sha256": "2" * 64,
        "accepted_joule_heat_result_sha256": "2" * 64,
    }
    generation = "nonlinear-eigenmode-171"
    summary[
        "nonlinear_eigenmode_continuation_parameter_normalization_phase_mac_branch_eigenvalue_mesh_result_generation_identity"
    ] = {
        "eigenmode_generation": generation,
        "continuation_eigenmode_generation": generation,
        "normalization_eigenmode_generation": generation,
        "phase_eigenmode_generation": generation,
        "mac_eigenmode_generation": generation,
        "branch_eigenmode_generation": generation,
        "eigenvalue_eigenmode_generation": generation,
        "mesh_eigenmode_generation": generation,
        "result_eigenmode_generation": generation,
        "continuation_parameter_name": "prestress_scale",
        "result_continuation_parameter_name": "prestress_scale",
        "continuation_parameter_values": [0.0, 0.5, 1.0],
        "result_continuation_parameter_values": [0.0, 0.5, 1.0],
        "mode_normalization": "unit_mass",
        "result_mode_normalization": "unit_mass",
        "phase_anchor_dof": "tip-z",
        "result_phase_anchor_dof": "tip-z",
        "mac_reference_branch_ids": [1, 2],
        "result_mac_reference_branch_ids": [1, 2],
        "mode_branch_ids": [[1, 2], [1, 2], [1, 2]],
        "result_mode_branch_ids": [[1, 2], [1, 2], [1, 2]],
        "eigenvalues_ri": [
            [[100.0, 0.0], [150.0, 0.0]],
            [[102.0, 0.0], [148.0, 0.0]],
            [[105.0, 0.0], [145.0, 0.0]],
        ],
        "result_eigenvalues_ri": [
            [[100.0, 0.0], [150.0, 0.0]],
            [[102.0, 0.0], [148.0, 0.0]],
            [[105.0, 0.0], [145.0, 0.0]],
        ],
        "mac_assignment_sha256": "3" * 64,
        "result_mac_assignment_sha256": "3" * 64,
        "eigenmode_mesh_sha256": "4" * 64,
        "result_eigenmode_mesh_sha256": "4" * 64,
        "eigenmode_result_sha256": "5" * 64,
        "accepted_eigenmode_result_sha256": "5" * 64,
    }
    return summary


def _with_v31_transform_and_force_identity(summary: dict) -> dict:
    summary = _with_v30_joule_and_eigenmode_identity(summary)
    generation = "frequency-time-closure-181"
    summary[
        "frequency_time_hermitian_spacing_window_group_delay_parseval_mesh_result_generation_identity"
    ] = {
        "transform_generation": generation,
        "spectrum_transform_generation": generation,
        "window_transform_generation": generation,
        "delay_transform_generation": generation,
        "energy_transform_generation": generation,
        "mesh_transform_generation": generation,
        "result_transform_generation": generation,
        "frequencies_hz": [0.0, 100.0, 200.0],
        "result_frequencies_hz": [0.0, 100.0, 200.0],
        "frequency_spacing_hz": 100.0,
        "result_frequency_spacing_hz": 100.0,
        "spectrum_ri": [[1.0, 0.0], [0.5, -0.25], [0.2, 0.0]],
        "result_spectrum_ri": [[1.0, 0.0], [0.5, -0.25], [0.2, 0.0]],
        "hermitian_completion": "conjugate_negative_frequencies",
        "result_hermitian_completion": "conjugate_negative_frequencies",
        "window_name": "hann_periodic",
        "result_window_name": "hann_periodic",
        "window_coherent_gain": 0.5,
        "result_window_coherent_gain": 0.5,
        "group_delay_s": 0.001,
        "result_group_delay_s": 0.001,
        "time_origin_s": 0.0,
        "result_time_origin_s": 0.0,
        "frequency_energy": 2.0,
        "time_energy": 2.0,
        "parseval_relative_tolerance": 1.0e-9,
        "transform_mesh_sha256": "1" * 64,
        "result_transform_mesh_sha256": "1" * 64,
        "time_trace_sha256": "2" * 64,
        "accepted_time_trace_sha256": "2" * 64,
    }
    generation = "rotating-force-balance-181"
    summary[
        "rotating_force_virtual_work_stress_phase_frame_lever_angle_power_mesh_result_generation_identity"
    ] = {
        "force_generation": generation,
        "virtual_work_force_generation": generation,
        "stress_force_generation": generation,
        "phase_force_generation": generation,
        "frame_force_generation": generation,
        "angle_force_generation": generation,
        "power_force_generation": generation,
        "mesh_force_generation": generation,
        "result_force_generation": generation,
        "phasor_convention": "exp_positive_jwt_rms",
        "result_phasor_convention": "exp_positive_jwt_rms",
        "coordinate_frame": "rotor_material",
        "result_coordinate_frame": "rotor_material",
        "lever_arm_m": [0.05, 0.0, 0.0],
        "result_lever_arm_m": [0.05, 0.0, 0.0],
        "mechanical_angles_rad": [0.0, 0.1, 0.2],
        "result_mechanical_angles_rad": [0.0, 0.1, 0.2],
        "virtual_work_torque_nm": [1.0, 1.2, 1.1],
        "stress_tensor_torque_nm": [1.0, 1.2, 1.1],
        "torque_relative_tolerance": 1.0e-9,
        "mechanical_power_w": 110.0,
        "airgap_power_w": 110.0,
        "power_relative_tolerance": 1.0e-9,
        "force_mesh_sha256": "3" * 64,
        "result_force_mesh_sha256": "3" * 64,
        "force_result_sha256": "4" * 64,
        "accepted_force_result_sha256": "4" * 64,
    }
    return summary


def _with_v32_nonlinear_and_eigenmode_identity(summary: dict) -> dict:
    summary = _with_v31_transform_and_force_identity(summary)
    generation = "nonlinear-segregated-closure-191"
    summary[
        "nonlinear_segregated_group_relaxation_residual_jacobian_continuation_mesh_result_generation_identity"
    ] = {
        "nonlinear_generation": generation,
        "segregated_group_generation": generation,
        "relaxation_generation": generation,
        "residual_generation": generation,
        "jacobian_generation": generation,
        "continuation_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "segregated_group_order": ["magnetic", "thermal", "structural"],
        "result_segregated_group_order": ["magnetic", "thermal", "structural"],
        "relaxation_schedule": [0.5, 0.7, 1.0],
        "result_relaxation_schedule": [0.5, 0.7, 1.0],
        "residual_norms": [1.0e-2, 2.0e-5, 4.0e-9],
        "accepted_residual_norms": [1.0e-2, 2.0e-5, 4.0e-9],
        "residual_relative_tolerance": 1.0e-8,
        "continuation_parameter": 1.0,
        "accepted_continuation_parameter": 1.0,
        "continuation_unit": "1",
        "accepted_continuation_unit": "1",
        "jacobian_sha256": "1" * 64,
        "accepted_jacobian_sha256": "1" * 64,
        "nonlinear_mesh_sha256": "2" * 64,
        "accepted_nonlinear_mesh_sha256": "2" * 64,
        "nonlinear_solution_sha256": "3" * 64,
        "accepted_nonlinear_solution_sha256": "3" * 64,
    }
    generation = "degenerate-eigenmode-closure-191"
    summary[
        "degenerate_eigenmode_subspace_phase_normalization_participation_mass_mesh_owner_result_generation_identity"
    ] = {
        "mode_generation": generation,
        "subspace_generation": generation,
        "phase_generation": generation,
        "normalization_generation": generation,
        "participation_generation": generation,
        "mass_generation": generation,
        "mesh_generation": generation,
        "result_generation": generation,
        "eigenvalues": [100.0, 100.0],
        "result_eigenvalues": [100.0, 100.0],
        "degenerate_subspace_sha256": "4" * 64,
        "result_degenerate_subspace_sha256": "4" * 64,
        "phase_anchor_dofs": [12, 37],
        "result_phase_anchor_dofs": [12, 37],
        "normalization": "mass_orthonormal",
        "result_normalization": "mass_orthonormal",
        "participation_factors": [0.8, 0.6],
        "result_participation_factors": [0.8, 0.6],
        "effective_masses_kg": [0.64, 0.36],
        "result_effective_masses_kg": [0.64, 0.36],
        "mode_owner": "component1/solid/eig1",
        "result_mode_owner": "component1/solid/eig1",
        "eigenmode_mesh_sha256": "5" * 64,
        "result_eigenmode_mesh_sha256": "5" * 64,
        "eigenmode_result_sha256": "6" * 64,
        "accepted_eigenmode_result_sha256": "6" * 64,
    }
    return summary


def _with_v33_contact_and_dae_identity(summary: dict) -> dict:
    summary = _with_v32_nonlinear_and_eigenmode_identity(summary)
    generation = "contact-complementarity-201"
    summary[
        "contact_gap_pressure_active_set_friction_dissipation_normal_mesh_result_generation_identity"
    ] = {
        "contact_generation": generation,
        **{
            key: generation
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
        },
        "contact_pair": "contact1",
        "result_contact_pair": "contact1",
        "active_contact_ids": [3, 7],
        "result_active_contact_ids": [3, 7],
        "normal_gap_m": [0.0, 0.0],
        "result_normal_gap_m": [0.0, 0.0],
        "normal_pressure_pa": [1.0e6, 2.0e6],
        "result_normal_pressure_pa": [1.0e6, 2.0e6],
        "tangential_slip_m": [0.0, 1.0e-4],
        "result_tangential_slip_m": [0.0, 1.0e-4],
        "friction_traction_pa": [0.0, 2.0e5],
        "result_friction_traction_pa": [0.0, 2.0e5],
        "contact_area_m2": [1.0e-4, 1.0e-4],
        "result_contact_area_m2": [1.0e-4, 1.0e-4],
        "friction_coefficient": 0.2,
        "result_friction_coefficient": 0.2,
        "friction_dissipation_j": 0.002,
        "result_friction_dissipation_j": 0.002,
        "normal_orientation": "outward_slave_to_master",
        "result_normal_orientation": "outward_slave_to_master",
        "contact_mesh_sha256": "1" * 64,
        "result_contact_mesh_sha256": "1" * 64,
        "contact_result_sha256": "2" * 64,
        "accepted_contact_result_sha256": "2" * 64,
    }
    generation = "field-circuit-dae-201"
    summary[
        "field_circuit_dae_charge_current_event_energy_time_dataset_result_generation_identity"
    ] = {
        "dae_generation": generation,
        **{
            key: generation
            for key in (
                "charge_generation",
                "current_generation",
                "event_generation",
                "energy_generation",
                "time_generation",
                "dataset_generation",
                "result_generation",
            )
        },
        "time_s": [0.0, 0.5e-3, 1.0e-3],
        "result_time_s": [0.0, 0.5e-3, 1.0e-3],
        "switch_event_time_s": 0.5e-3,
        "result_switch_event_time_s": 0.5e-3,
        "event_side": "right_limit_after_event",
        "result_event_side": "right_limit_after_event",
        "charge_c": [0.0, 1.0e-6, 1.5e-6],
        "result_charge_c": [0.0, 1.0e-6, 1.5e-6],
        "integrated_current_c": [0.0, 1.0e-6, 1.5e-6],
        "result_integrated_current_c": [0.0, 1.0e-6, 1.5e-6],
        "algebraic_residual_c": [0.0, 1.0e-14, 0.0],
        "accepted_algebraic_residual_c": [0.0, 1.0e-14, 0.0],
        "algebraic_tolerance_c": 1.0e-12,
        "current_sign_convention": "positive_into_field_device",
        "result_current_sign_convention": "positive_into_field_device",
        "stored_energy_before_j": 0.002,
        "result_stored_energy_before_j": 0.002,
        "stored_energy_after_j": 0.0018,
        "result_stored_energy_after_j": 0.0018,
        "switch_dissipation_j": 0.0002,
        "result_switch_dissipation_j": 0.0002,
        "dataset_owner": "dset1/sol2",
        "result_dataset_owner": "dset1/sol2",
        "dae_dataset_sha256": "3" * 64,
        "result_dae_dataset_sha256": "3" * 64,
        "dae_result_sha256": "4" * 64,
        "accepted_dae_result_sha256": "4" * 64,
    }
    return summary


def _with_v34_arclength_and_electrochemical_identity(summary: dict) -> dict:
    summary = _with_v33_contact_and_dae_identity(summary)
    generation = "arclength-211"
    summary[
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    ] = {
        "continuation_generation": generation,
        **{
            key: generation
            for key in (
                "arclength_generation",
                "tangent_generation",
                "branch_generation",
                "turning_generation",
                "residual_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "previous_augmented_state": [1.0, 0.0, 1.0],
        "result_previous_augmented_state": [1.0, 0.0, 1.0],
        "predictor_tangent": [0.6, 0.0, 0.8],
        "result_predictor_tangent": [0.6, 0.0, 0.8],
        "arclength_step": 0.05,
        "result_arclength_step": 0.05,
        "predictor_augmented_state": [1.03, 0.0, 1.04],
        "result_predictor_augmented_state": [1.03, 0.0, 1.04],
        "corrected_augmented_state": [1.03, 0.0, 1.04],
        "result_corrected_augmented_state": [1.03, 0.0, 1.04],
        "branch_id": "upper_branch",
        "result_branch_id": "upper_branch",
        "turning_point_side": "pre_turn_positive_parameter_tangent",
        "result_turning_point_side": "pre_turn_positive_parameter_tangent",
        "corrected_residual_norm": 1.0e-10,
        "result_corrected_residual_norm": 1.0e-10,
        "residual_tolerance": 1.0e-8,
        "result_residual_tolerance": 1.0e-8,
        "continuation_mesh_sha256": "1" * 64,
        "result_continuation_mesh_sha256": "1" * 64,
        "continuation_result_sha256": "2" * 64,
        "accepted_continuation_result_sha256": "2" * 64,
    }
    generation = "electrochemical-211"
    summary[
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    ] = {
        "electrochemical_generation": generation,
        **{
            key: generation
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
        },
        "species_order": ["A_plus", "B_minus", "C_neutral"],
        "result_species_order": ["A_plus", "B_minus", "C_neutral"],
        "charge_numbers": [1, -1, 0],
        "result_charge_numbers": [1, -1, 0],
        "molar_mass_basis": [1.0, 1.0, 1.0],
        "result_molar_mass_basis": [1.0, 1.0, 1.0],
        "reaction_stoichiometry": [-1.0, -1.0, 2.0],
        "result_reaction_stoichiometry": [-1.0, -1.0, 2.0],
        "reaction_extent_mol": 0.5,
        "result_reaction_extent_mol": 0.5,
        "initial_inventory_mol": [1.0, 1.0, 0.0],
        "result_initial_inventory_mol": [1.0, 1.0, 0.0],
        "final_inventory_mol": [0.5, 0.5, 1.0],
        "result_final_inventory_mol": [0.5, 0.5, 1.0],
        "integrated_boundary_flux_mol": [0.0, 0.0, 0.0],
        "result_integrated_boundary_flux_mol": [0.0, 0.0, 0.0],
        "integrated_electric_current_c": 0.0,
        "result_integrated_electric_current_c": 0.0,
        "initial_free_energy_j": 2.0,
        "result_initial_free_energy_j": 2.0,
        "final_free_energy_j": 1.8,
        "result_final_free_energy_j": 1.8,
        "dissipated_free_energy_j": 0.2,
        "result_dissipated_free_energy_j": 0.2,
        "time_s": [0.0, 1.0],
        "result_time_s": [0.0, 1.0],
        "electrochemical_mesh_sha256": "3" * 64,
        "result_electrochemical_mesh_sha256": "3" * 64,
        "electrochemical_result_sha256": "4" * 64,
        "accepted_electrochemical_result_sha256": "4" * 64,
    }
    return summary


def _with_v35_multirate_and_adjoint_identity(summary: dict) -> dict:
    summary = _with_v34_arclength_and_electrochemical_identity(summary)
    generation = "multirate-coupling-221"
    summary[
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    ] = {
        "coupling_generation": generation,
        **{
            key: generation
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
        },
        "electrical_time_s": [0.0, 0.0005, 0.001],
        "result_electrical_time_s": [0.0, 0.0005, 0.001],
        "mechanical_time_s": [0.0, 0.00025, 0.0005, 0.00075, 0.001],
        "result_mechanical_time_s": [0.0, 0.00025, 0.0005, 0.00075, 0.001],
        "event_time_s": 0.0005,
        "result_event_time_s": 0.0005,
        "event_interpolation_side": "right_continuous_after_event",
        "result_event_interpolation_side": "right_continuous_after_event",
        "substep_owner": "coupler:electrical2_mechanical4",
        "result_substep_owner": "coupler:electrical2_mechanical4",
        "coordinate_frame": "stationary_xyz",
        "result_coordinate_frame": "stationary_xyz",
        "electrical_input_energy_j": 0.012,
        "result_electrical_input_energy_j": 0.012,
        "mechanical_output_work_j": 0.009,
        "result_mechanical_output_work_j": 0.009,
        "dissipated_energy_j": 0.003,
        "result_dissipated_energy_j": 0.003,
        "energy_balance_tolerance_j": 1.0e-10,
        "result_energy_balance_tolerance_j": 1.0e-10,
        "coupling_mesh_sha256": "1" * 64,
        "result_coupling_mesh_sha256": "1" * 64,
        "coupling_result_sha256": "2" * 64,
        "accepted_coupling_result_sha256": "2" * 64,
    }
    generation = "adjoint-sensitivity-221"
    summary[
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    ] = {
        "sensitivity_generation": generation,
        **{
            key: generation
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
        },
        "objective_tag": "torque_ripple_rms",
        "result_objective_tag": "torque_ripple_rms",
        "design_variable": "magnet_arc_rad",
        "result_design_variable": "magnet_arc_rad",
        "design_scale": 0.1,
        "result_design_scale": 0.1,
        "active_constraint": "magnet_volume_constant",
        "result_active_constraint": "magnet_volume_constant",
        "adjoint_gradient": 2.5,
        "chainrule_gradient": 2.5,
        "finite_difference_gradient": 2.500001,
        "gradient_tolerance": 1.0e-4,
        "result_gradient_tolerance": 1.0e-4,
        "fd_perturbation": 1.0e-5,
        "result_fd_perturbation": 1.0e-5,
        "sensitivity_mesh_sha256": "3" * 64,
        "result_sensitivity_mesh_sha256": "3" * 64,
        "primal_solution_sha256": "4" * 64,
        "result_primal_solution_sha256": "4" * 64,
        "gradient_result_sha256": "5" * 64,
        "accepted_gradient_result_sha256": "5" * 64,
    }
    return summary


def _modal_response(
    response_frequencies: list[float],
    mode_frequencies: list[float],
    damping: list[float],
    participation: list[float],
    probe_factors: list[float],
) -> list[list[float]]:
    values = []
    for frequency in response_frequencies:
        omega = 2.0 * math.pi * frequency
        response = 0.0j
        for mode_hz, zeta, factor, probe in zip(
            mode_frequencies, damping, participation, probe_factors, strict=True
        ):
            omega_mode = 2.0 * math.pi * mode_hz
            response += probe * factor / complex(
                omega_mode**2 - omega**2,
                2.0 * zeta * omega_mode * omega,
            )
        values.append([response.real, response.imag])
    return values


def _with_v36_force_and_modal_identity(summary: dict) -> dict:
    summary = _with_v35_multirate_and_adjoint_identity(summary)
    generation = "virtual-work-coenergy-231"
    displacement = [-1.0e-4, 0.0, 1.0e-4]
    coenergy = [0.4997, 0.5, 0.5003]
    force = (coenergy[2] - coenergy[0]) / (displacement[2] - displacement[0])
    summary[
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    ] = {
        "force_generation": generation,
        **{
            key: generation
            for key in (
                "displacement_generation",
                "coenergy_generation",
                "current_generation",
                "mesh_generation",
                "frame_generation",
                "solution_generation",
                "result_generation",
            )
        },
        "displacement_m": displacement,
        "result_displacement_m": displacement,
        "coenergy_j": coenergy,
        "result_coenergy_j": coenergy,
        "held_source_convention": "constant_current",
        "result_held_source_convention": "constant_current",
        "force_sign_convention": "positive_dcoenergy_dx",
        "result_force_sign_convention": "positive_dcoenergy_dx",
        "central_coenergy_force_n": force,
        "result_force_n": force,
        "force_tolerance_n": 1.0e-9,
        "result_force_tolerance_n": 1.0e-9,
        "coordinate_frame": "stationary_cartesian_x",
        "result_coordinate_frame": "stationary_cartesian_x",
        "displaced_mesh_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "result_displaced_mesh_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "force_solution_owner": "std1/sol1:parametric_displacement",
        "result_force_solution_owner": "std1/sol1:parametric_displacement",
        "force_result_sha256": "4" * 64,
        "accepted_force_result_sha256": "4" * 64,
    }

    generation = "acoustic-modal-participation-231"
    frequencies = [100.0, 160.0]
    masses = [2.0, 1.5]
    participation = [0.5, 0.4]
    damping = [0.01, 0.02]
    probes = [1.0, 0.8]
    response_frequencies = [90.0, 120.0, 180.0]
    response = _modal_response(
        response_frequencies, frequencies, damping, participation, probes
    )
    effective_masses = [
        factor**2 * mass
        for factor, mass in zip(participation, masses, strict=True)
    ]
    summary[
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    ] = {
        "modal_generation": generation,
        **{
            key: generation
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
        },
        "normalization": "peak_displacement",
        "result_normalization": "peak_displacement",
        "mode_frequency_hz": frequencies,
        "result_mode_frequency_hz": frequencies,
        "modal_mass_kg": masses,
        "result_modal_mass_kg": masses,
        "participation_factor": participation,
        "result_participation_factor": participation,
        "effective_modal_mass_kg": effective_masses,
        "result_effective_modal_mass_kg": effective_masses,
        "damping_ratio": damping,
        "result_damping_ratio": damping,
        "probe_mode_factor": probes,
        "result_probe_mode_factor": probes,
        "response_frequency_hz": response_frequencies,
        "result_response_frequency_hz": response_frequencies,
        "probe_response_complex": response,
        "result_probe_response_complex": response,
        "response_tolerance": 1.0e-12,
        "result_response_tolerance": 1.0e-12,
        "modal_mesh_sha256": "5" * 64,
        "result_modal_mesh_sha256": "5" * 64,
        "modal_result_sha256": "6" * 64,
        "accepted_modal_result_sha256": "6" * 64,
    }
    return summary


def _with_v37_capacitance_and_thermoelastic_identity(summary: dict) -> dict:
    summary = _with_v36_force_and_modal_identity(summary)
    generation = "capacitance-closure-241"
    matrix = [[2.0e-12, -2.0e-12], [-2.0e-12, 2.0e-12]]
    summary[
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    ] = {
        "capacitance_generation": generation,
        **{
            key: generation
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
        },
        "terminal_names": ["terminal_1", "reference_2"],
        "result_terminal_names": ["terminal_1", "reference_2"],
        "reference_terminal": "reference_2",
        "result_reference_terminal": "reference_2",
        "capacitance_matrix_f": matrix,
        "result_capacitance_matrix_f": matrix,
        "terminal_potential_v": [1.0, 0.0],
        "result_terminal_potential_v": [1.0, 0.0],
        "terminal_charge_c": [2.0e-12, -2.0e-12],
        "result_terminal_charge_c": [2.0e-12, -2.0e-12],
        "stored_energy_j": 1.0e-12,
        "result_stored_energy_j": 1.0e-12,
        "reciprocity_tolerance": 1.0e-12,
        "result_reciprocity_tolerance": 1.0e-12,
        "terminal_owner": "comp1/es/terminals-241",
        "accepted_terminal_owner": "comp1/es/terminals-241",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "thermoelastic-harmonic-241"
    heat = [10.0, 2.0]
    summary[
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    ] = {
        "thermoelastic_generation": generation,
        **{
            key: generation
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
        },
        "frequency_hz": 1000.0,
        "result_frequency_hz": 1000.0,
        "heat_source_complex_w": heat,
        "result_heat_source_complex_w": heat,
        "heat_source_phase_rad": math.atan2(heat[1], heat[0]),
        "result_heat_source_phase_rad": math.atan2(heat[1], heat[0]),
        "temperature_complex_k": [5.0, 1.0],
        "result_temperature_complex_k": [5.0, 1.0],
        "displacement_complex_m": [1.0e-6, -2.0e-7],
        "result_displacement_complex_m": [1.0e-6, -2.0e-7],
        "thermal_expansion_work_j": 2.0e-3,
        "result_thermal_expansion_work_j": 2.0e-3,
        "mechanical_loss_j": 1.0e-4,
        "result_mechanical_loss_j": 1.0e-4,
        "loss_convention": "positive_dissipated_per_cycle",
        "result_loss_convention": "positive_dissipated_per_cycle",
        "mesh_owner": "comp1/mesh1:thermoelastic-241",
        "accepted_mesh_owner": "comp1/mesh1:thermoelastic-241",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _with_v38_thermoviscous_and_piezoelectric_identity(summary: dict) -> dict:
    summary = _with_v37_capacitance_and_thermoelastic_identity(summary)
    generation = "thermoviscous-interface-258"
    summary[
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    ] = {
        "interface_generation": generation,
        **{
            key: generation
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
        },
        "frequency_hz": 1.0e4,
        "result_frequency_hz": 1.0e4,
        "interface_area_m2": 1.0e-2,
        "result_interface_area_m2": 1.0e-2,
        "normal_velocity_complex_m_per_s": [1.0e-2, 2.0e-3],
        "result_normal_velocity_complex_m_per_s": [1.0e-2, 2.0e-3],
        "pressure_complex_pa": [200.0, 40.0],
        "result_pressure_complex_pa": [200.0, 40.0],
        "traction_sign": "minus_pressure_times_outward_normal",
        "result_traction_sign": "minus_pressure_times_outward_normal",
        "normal_orientation": "thermoviscous_to_pressure_acoustics",
        "result_normal_orientation": "thermoviscous_to_pressure_acoustics",
        "interface_power_w": 1.04e-2,
        "result_interface_power_w": 1.04e-2,
        "viscous_loss_w": 3.0e-3,
        "result_viscous_loss_w": 3.0e-3,
        "thermal_loss_w": 2.0e-3,
        "result_thermal_loss_w": 2.0e-3,
        "outgoing_acoustic_power_w": 5.4e-3,
        "result_outgoing_acoustic_power_w": 5.4e-3,
        "mesh_owner": "comp1/mesh1:thermoviscous-258",
        "accepted_mesh_owner": "comp1/mesh1:thermoviscous-258",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "piezoelectric-reciprocity-258"
    summary[
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    ] = {
        "piezo_generation": generation,
        **{
            key: generation
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
        },
        "direct_coefficient_c_per_n": 2.0e-10,
        "result_direct_coefficient_c_per_n": 2.0e-10,
        "converse_coefficient_m_per_v": 2.0e-10,
        "result_converse_coefficient_m_per_v": 2.0e-10,
        "electric_field_v_per_m": 1.0e5,
        "result_electric_field_v_per_m": 1.0e5,
        "mechanical_stress_pa": 1.0e6,
        "result_mechanical_stress_pa": 1.0e6,
        "induced_strain": 2.0e-5,
        "result_induced_strain": 2.0e-5,
        "induced_charge_density_c_per_m2": 2.0e-4,
        "result_induced_charge_density_c_per_m2": 2.0e-4,
        "terminal_charge_c": 4.0e-6,
        "result_terminal_charge_c": 4.0e-6,
        "electrical_work_j": 3.0e-2,
        "result_electrical_work_j": 3.0e-2,
        "elastic_energy_j": 2.0e-2,
        "result_elastic_energy_j": 2.0e-2,
        "coupling_energy_j": 1.0e-2,
        "result_coupling_energy_j": 1.0e-2,
        "total_stored_energy_j": 4.0e-2,
        "result_total_stored_energy_j": 4.0e-2,
        "polarization_frame": "material_axis_3",
        "result_polarization_frame": "material_axis_3",
        "mesh_owner": "comp1/mesh1:piezo-258",
        "accepted_mesh_owner": "comp1/mesh1:piezo-258",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _with_v39_poroelastic_and_induction_identity(summary: dict) -> dict:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(summary)
    generation = "poroelastic-biot-377"
    alpha, pressure, strain, modulus = 0.8, 1.0e5, 1.0e-3, 1.0e9
    permeability, viscosity, gradient = 1.0e-12, 1.0e-3, 1.0e6
    volume, timestep = 1.0e-2, 0.1
    flux = -permeability * gradient / viscosity
    mirrored = {
        "biot_coefficient": alpha, "pore_pressure_pa": pressure,
        "volumetric_strain": strain, "biot_modulus_pa": modulus,
        "fluid_content_increment": alpha * strain + pressure / modulus,
        "permeability_m2": permeability, "dynamic_viscosity_pa_s": viscosity,
        "pressure_gradient_pa_per_m": gradient, "darcy_flux_m_per_s": flux,
        "domain_volume_m3": volume, "time_step_s": timestep,
        "interface_traction_pa": -alpha * pressure,
        "storage_energy_j": 0.5 * pressure * pressure / modulus * volume,
        "skeleton_coupling_work_j": alpha * pressure * strain * volume,
        "fluid_dissipation_j": viscosity / permeability * flux * flux * volume * timestep,
    }
    summary["poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"] = {
        "poroelastic_generation": generation,
        **{key: generation for key in (
            "biot_generation", "pressure_generation", "displacement_generation",
            "flux_generation", "storage_generation", "dissipation_generation",
            "interface_generation", "mesh_generation", "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "interface_normal": "porous_skeleton_to_free_fluid",
        "result_interface_normal": "porous_skeleton_to_free_fluid",
        "mesh_owner": "comp1/mesh1:poroelastic-377",
        "accepted_mesh_owner": "comp1/mesh1:poroelastic-377",
        "mesh_sha256": "1" * 64, "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "rotating-induction-377"
    frequency, pole_pairs, rotor_speed, torque = 50.0, 2, 150.0, 20.0
    synchronous = 2.0 * math.pi * frequency / pole_pairs
    slip = (synchronous - rotor_speed) / synchronous
    airgap_power, resistance = torque * synchronous, 0.2
    current = math.sqrt((slip * airgap_power) / (3.0 * resistance))
    mirrored = {
        "supply_frequency_hz": frequency, "pole_pairs": pole_pairs,
        "synchronous_speed_rad_per_s": synchronous,
        "rotor_speed_rad_per_s": rotor_speed, "slip": slip,
        "rotor_electrical_frequency_hz": slip * frequency,
        "rotor_phase_current_a_rms": current,
        "rotor_phase_resistance_ohm": resistance,
        "rotor_copper_loss_w": 3.0 * current * current * resistance,
        "airgap_torque_nm": torque, "airgap_power_w": airgap_power,
        "mechanical_power_w": torque * rotor_speed,
    }
    summary["rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"] = {
        "induction_generation": generation,
        **{key: generation for key in (
            "slip_generation", "frequency_generation", "current_generation",
            "loss_generation", "torque_generation", "power_generation",
            "frame_generation", "mesh_generation", "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "rotating_frame": "rotor_mechanical_frame",
        "result_rotating_frame": "rotor_mechanical_frame",
        "mesh_owner": "comp1/mesh1:induction-377",
        "accepted_mesh_owner": "comp1/mesh1:induction-377",
        "mesh_sha256": "3" * 64, "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


_THERMOACOUSTIC_KEY = (
    "thermoacoustic_meanflow_convected_wavenumber_flux_impedance_power_mesh_"
    "result_generation_identity"
)


_BATTERY_KEY = (
    "battery_electrothermal_soc_current_heat_temperature_energy_safety_mesh_"
    "result_generation_identity"
)


def _with_v40_thermoacoustic_and_battery_identity(summary: dict) -> dict:
    summary = _with_v39_poroelastic_and_induction_identity(summary)
    generation = "thermoacoustic-719"
    frequency, sound_speed, mach = 1000.0, 343.0, 0.1
    density, pressure, area = 1.2, 1.0, 0.1
    flow_speed = mach * sound_speed
    particle_velocity = pressure / (density * sound_speed)
    intensity = pressure * particle_velocity
    power = intensity * area
    mirrored = {
        "frequency_hz": frequency, "sound_speed_m_per_s": sound_speed,
        "mean_flow_mach": mach, "mean_flow_speed_m_per_s": flow_speed,
        "downstream_wavenumber_rad_per_m": 2.0 * math.pi * frequency / (sound_speed + flow_speed),
        "upstream_wavenumber_rad_per_m": 2.0 * math.pi * frequency / (sound_speed - flow_speed),
        "density_kg_per_m3": density, "pressure_rms_pa": pressure,
        "particle_velocity_rms_m_per_s": particle_velocity,
        "acoustic_intensity_w_per_m2": intensity, "boundary_area_m2": area,
        "boundary_impedance_pa_s_per_m": density * sound_speed,
        "boundary_flux_power_w": power, "impedance_work_w": power,
        "dissipated_power_w": power, "power_balance_residual_w": 0.0,
    }
    summary[_THERMOACOUSTIC_KEY] = {
        "thermoacoustic_generation": generation,
        **{key: generation for key in (
            "meanflow_generation", "wavenumber_generation", "flux_generation",
            "impedance_generation", "power_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:thermoacoustic-719",
        "accepted_mesh_owner": "comp1/mesh1:thermoacoustic-719",
        "mesh_sha256": "1" * 64, "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }

    generation = "battery-electrothermal-719"
    capacity, initial_soc, current, voltage = 3600.0, 0.8, 2.0, 3.7
    timestep, resistance = 10.0, 0.05
    irreversible_heat, reversible_heat = current * current * resistance * timestep, 0.5
    thermal_energy = irreversible_heat + reversible_heat
    mass, heat_capacity, initial_temperature = 0.05, 1000.0, 298.15
    mirrored = {
        "capacity_c": capacity, "initial_state_of_charge": initial_soc,
        "terminal_current_a": current, "terminal_voltage_v": voltage,
        "time_step_s": timestep,
        "final_state_of_charge": initial_soc - current * timestep / capacity,
        "internal_resistance_ohm": resistance,
        "irreversible_heat_j": irreversible_heat, "reversible_heat_j": reversible_heat,
        "thermal_energy_j": thermal_energy,
        "electrical_energy_j": voltage * current * timestep,
        "cell_mass_kg": mass, "specific_heat_j_per_kg_k": heat_capacity,
        "initial_temperature_k": initial_temperature,
        "final_temperature_k": initial_temperature + thermal_energy / (mass * heat_capacity),
        "maximum_safe_temperature_k": 333.15, "thermal_balance_residual_j": 0.0,
    }
    summary[_BATTERY_KEY] = {
        "battery_generation": generation,
        **{key: generation for key in (
            "soc_generation", "current_generation", "heat_generation",
            "temperature_generation", "energy_generation", "safety_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored, **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:battery-719",
        "accepted_mesh_owner": "comp1/mesh1:battery-719",
        "mesh_sha256": "3" * 64, "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary


_PIEZO_KEY = (
    "piezoelectric_admittance_resonance_antiresonance_coupling_energy_phase_"
    "mesh_result_generation_identity"
)


_BEARING_KEY = (
    "fluidfilm_bearing_reynolds_pressure_load_friction_temperature_power_mesh_"
    "result_generation_identity"
)


def _with_v41_piezoelectric_and_bearing_identity(summary: dict) -> dict:
    summary = _with_v40_thermoacoustic_and_battery_identity(summary)
    generation = "piezoelectric-724"
    resonance, antiresonance = 100_000.0, 105_000.0
    voltage, current, phase_deg = 10.0, 0.02, -30.0
    apparent_power = voltage * current
    real_power = apparent_power * math.cos(math.radians(phase_deg))
    mechanical_power = 0.15
    mirrored = {
        "resonance_frequency_hz": resonance,
        "antiresonance_frequency_hz": antiresonance,
        "electromechanical_coupling_squared": 1.0 - (resonance / antiresonance) ** 2,
        "voltage_rms_v": voltage,
        "current_rms_a": current,
        "admittance_magnitude_s": current / voltage,
        "admittance_phase_deg": phase_deg,
        "real_electrical_power_w": real_power,
        "reactive_electrical_power_var": apparent_power * math.sin(math.radians(phase_deg)),
        "mechanical_output_power_w": mechanical_power,
        "dielectric_loss_w": real_power - mechanical_power,
        "mechanical_stored_energy_j": 2.0e-6,
        "electric_stored_energy_j": 3.0e-6,
        "power_balance_residual_w": 0.0,
    }
    summary[_PIEZO_KEY] = {
        "piezoelectric_generation": generation,
        **{key: generation for key in (
            "admittance_generation", "resonance_generation", "coupling_generation",
            "phase_generation", "energy_generation", "power_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:piezoelectric-724",
        "accepted_mesh_owner": "comp1/mesh1:piezoelectric-724",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "fluidfilm-bearing-724"
    clearance, eccentricity = 50.0e-6, 0.6
    angular_speed, friction_torque = 100.0, 2.0
    shaft_power, viscous_power = angular_speed * friction_torque, 150.0
    mirrored = {
        "journal_radius_m": 0.025,
        "bearing_length_m": 0.05,
        "radial_clearance_m": clearance,
        "eccentricity_ratio": eccentricity,
        "minimum_film_thickness_m": clearance * (1.0 - eccentricity),
        "maximum_pressure_pa": 4.0e6,
        "integrated_load_n": 1000.0,
        "attitude_angle_deg": 55.0,
        "angular_speed_rad_per_s": angular_speed,
        "friction_torque_nm": friction_torque,
        "shaft_power_w": shaft_power,
        "viscous_dissipation_w": viscous_power,
        "removed_heat_w": shaft_power - viscous_power,
        "inlet_temperature_k": 313.15,
        "maximum_temperature_k": 333.15,
        "power_balance_residual_w": 0.0,
    }
    summary[_BEARING_KEY] = {
        "fluidfilm_generation": generation,
        **{key: generation for key in (
            "film_generation", "pressure_generation", "load_generation",
            "friction_generation", "temperature_generation", "power_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:fluidfilm-724",
        "accepted_mesh_owner": "comp1/mesh1:fluidfilm-724",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


_INDUCTION_KEY = (
    "inductionheating_skin_proximity_joule_thermal_flux_temperature_energy_"
    "mesh_result_generation_identity"
)


_SPECIES_KEY = (
    "species_transport_reaction_diffusion_flux_massbalance_rate_temperature_"
    "mesh_result_generation_identity"
)


def _with_v42_induction_and_species_identity(summary: dict) -> dict:
    summary = _with_v41_piezoelectric_and_bearing_identity(summary)
    generation = "induction-heating-725"
    frequency = 10_000.0
    conductivity = 5.8e7
    permeability = 4.0e-7 * math.pi
    joule_loss, magnetic_loss = 780.0, 20.0
    input_power = joule_loss + magnetic_loss
    ambient, maximum = 293.15, 373.15
    mirrored = {
        "frequency_hz": frequency,
        "conductivity_s_per_m": conductivity,
        "relative_permeability": 1.0,
        "skin_depth_m": math.sqrt(2.0 / (2.0 * math.pi * frequency * permeability * conductivity)),
        "surface_current_density_a_per_m": [1200.0, 1500.0, 1100.0],
        "proximity_current_density_a_per_m2": [4.0e7, 6.0e7, 3.5e7],
        "joule_loss_w": joule_loss,
        "magnetic_loss_w": magnetic_loss,
        "electromagnetic_input_power_w": input_power,
        "outward_thermal_flux_w": input_power,
        "ambient_temperature_k": ambient,
        "maximum_temperature_k": maximum,
        "temperature_rise_k": maximum - ambient,
        "electromagnetic_power_balance_residual_w": 0.0,
        "thermal_power_balance_residual_w": 0.0,
    }
    summary[_INDUCTION_KEY] = {
        "induction_generation": generation,
        **{key: generation for key in (
            "skin_generation", "proximity_generation", "joule_generation",
            "thermal_generation", "temperature_generation", "energy_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:induction-725",
        "accepted_mesh_owner": "component/mesh:induction-725",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "reacting-species-725"
    temperature, gas_constant = 350.0, 8.314462618
    preexponential, activation_energy = 1.0e6, 40_000.0
    rate = preexponential * math.exp(-activation_energy / (gas_constant * temperature))
    concentration, volume = 2.0, 0.01
    consumption = rate * concentration * volume
    mirrored = {
        "diffusivity_m2_per_s": 2.0e-9,
        "temperature_k": temperature,
        "gas_constant_j_per_mol_k": gas_constant,
        "preexponential_factor_per_s": preexponential,
        "activation_energy_j_per_mol": activation_energy,
        "reaction_rate_constant_per_s": rate,
        "mean_concentration_mol_per_m3": concentration,
        "domain_volume_m3": volume,
        "integrated_species_mol": concentration * volume,
        "integrated_consumption_mol_per_s": consumption,
        "inward_boundary_flux_mol_per_s": consumption,
        "mass_balance_residual_mol_per_s": 0.0,
    }
    summary[_SPECIES_KEY] = {
        "species_generation": generation,
        **{key: generation for key in (
            "diffusion_generation", "reaction_generation", "flux_generation",
            "mass_generation", "temperature_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:species-725",
        "accepted_mesh_owner": "component/mesh:species-725",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


_MICROWAVE_KEY = (
    "microwaveheating_sparameter_absorbedpower_jouleheat_temperature_energy_"
    "mesh_result_generation_identity"
)


_POROELASTIC_KEY = (
    "poroelastic_wave_pressure_displacement_flux_dissipation_mass_energy_"
    "mesh_result_generation_identity"
)


def _with_v43_microwave_and_poroelastic_identity(summary: dict) -> dict:
    summary = _with_v42_induction_and_species_identity(summary)
    generation = "microwave-heating-726"
    incident, reflected, transmitted = 100.0, 10.0, 5.0
    absorbed = incident - reflected - transmitted
    mirrored = {
        "frequency_hz": 2.45e9,
        "reference_impedance_ohm": 50.0,
        "s11_magnitude": math.sqrt(reflected / incident),
        "s21_magnitude": math.sqrt(transmitted / incident),
        "incident_power_w": incident,
        "reflected_power_w": reflected,
        "transmitted_power_w": transmitted,
        "absorbed_power_w": absorbed,
        "joule_heat_w": 80.0,
        "dielectric_heat_w": 5.0,
        "electromagnetic_power_residual_w": 0.0,
        "outward_thermal_flux_w": absorbed,
        "ambient_temperature_k": 293.15,
        "maximum_temperature_k": 335.65,
        "temperature_rise_k": 42.5,
        "thermal_power_residual_w": 0.0,
    }
    summary[_MICROWAVE_KEY] = {
        "microwave_generation": generation,
        **{key: generation for key in (
            "sparameter_generation", "power_generation", "heat_generation",
            "temperature_generation", "energy_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:microwave-726",
        "accepted_mesh_owner": "component/mesh:microwave-726",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "poroelastic-wave-726"
    mirrored = {
        "frequency_hz": 125.0,
        "porosity": 0.32,
        "solid_displacement_amplitude_m": 2.5e-6,
        "pore_pressure_amplitude_pa": 1250.0,
        "darcy_flux_amplitude_m_per_s": 3.0e-5,
        "pressure_displacement_phase_deg": -35.0,
        "fluid_mass_kg": 0.032,
        "fluid_mass_rate_kg_per_s": 0.004,
        "net_inward_mass_flux_kg_per_s": 0.004,
        "solid_energy_j": 1.2,
        "fluid_energy_j": 0.8,
        "dissipated_power_w": 0.4,
        "input_power_w": 0.4,
        "mass_balance_residual_kg_per_s": 0.0,
        "energy_balance_residual_w": 0.0,
    }
    summary[_POROELASTIC_KEY] = {
        "poroelastic_generation": generation,
        **{key: generation for key in (
            "pressure_generation", "displacement_generation", "flux_generation",
            "mass_generation", "energy_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:poroelastic-726",
        "accepted_mesh_owner": "component/mesh:poroelastic-726",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def _with_v44(summary: dict) -> dict:
    summary = _with_v43_microwave_and_poroelastic_identity(summary)
    generation = "microwave-port-744"
    values = {
        "frequency_hz": 2.45e9, "port_normalization_ohm": 50.0,
        "incident_power_w": 100.0, "reflected_power_w": 10.0,
        "transmitted_power_w": 5.0, "absorbed_power_w": 85.0,
        "s11_power_fraction": 0.10, "s21_power_fraction": 0.05,
        "thermal_coupling_power_w": 85.0, "temperature_rise_k": 42.5,
    }
    key = "microwave_boundaryport_sparameter_power_normalization_temperature_coupling_restart_owner_result_identity"
    summary[key] = {
        "microwave_port_generation": generation, "port_generation": generation,
        "power_generation": generation, "thermal_generation": generation,
        "restart_generation": generation, "owner_generation": generation,
        "result_generation": generation, **values,
        **{f"result_{name}": value for name, value in values.items()},
        "restart_checkpoint_id": "chk-microwave-744", "result_restart_checkpoint_id": "chk-microwave-744",
        "mesh_owner": "component/mesh:microwave-744", "accepted_mesh_owner": "component/mesh:microwave-744",
        "mesh_sha256": "7" * 64, "accepted_mesh_sha256": "7" * 64,
        "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
    }
    generation = "acoustic-poroelastic-744"
    values = {
        "frequency_hz": 125.0, "impedance_magnitude_pa_s_per_m": 2.5e5,
        "impedance_phase_deg": -35.0, "pressure_flux_w": 0.4,
        "mechanical_energy_j": 1.2, "fluid_energy_j": 0.8,
        "dissipated_power_w": 0.4, "time_window_s": 0.008,
        "energy_balance_residual_w": 0.0,
    }
    key = "acoustics_poroelastic_impedance_phase_flux_energy_timewindow_dataset_owner_result_identity"
    summary[key] = {
        "acoustic_generation": generation, "impedance_generation": generation,
        "phase_generation": generation, "flux_generation": generation,
        "energy_generation": generation, "dataset_generation": generation,
        "owner_generation": generation, "result_generation": generation,
        **values, **{f"result_{name}": value for name, value in values.items()},
        "dataset_tag": "dset-acoustic-744", "result_dataset_tag": "dset-acoustic-744",
        "boundary_owner": "component/boundary:acoustic-744", "accepted_boundary_owner": "component/boundary:acoustic-744",
        "boundary_sha256": "9" * 64, "accepted_boundary_sha256": "9" * 64,
        "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
    }
    return summary


def _with_v45(summary: dict) -> dict:
    summary["microwave_sparameter_port_reference_plane_deembed_complex_power_mesh_result_identity"] = {
        "generation": "microwave-port-v45-812", "reference_plane_m_generation": "microwave-port-v45-812",
        "deembed_length_m_generation": "microwave-port-v45-812", "complex_power_w_generation": "microwave-port-v45-812",
        "mesh_generation_generation": "mesh-microwave-v45-812", "port_mode_generation": "microwave-port-v45-812",
        "s11_complex_generation": "microwave-port-v45-812", "s21_complex_generation": "microwave-port-v45-812",
        "frequency_hz": 2.45e9, "reference_plane_m": 0.012, "deembed_length_m": 0.008, "complex_power_w": 85.0,
        "port_mode": "TE10",
        "s11_complex": {"real": -0.2, "imag": 0.1}, "s21_complex": {"real": 0.7, "imag": -0.1},
        "result_reference_plane_m": 0.012, "result_deembed_length_m": 0.008, "result_complex_power_w": 85.0,
        "result_port_mode": "TE10", "result_s11_complex": {"real": -0.2, "imag": 0.1}, "result_s21_complex": {"real": 0.7, "imag": -0.1},
        "owner": "model/microwave-v45-812", "accepted_owner": "model/microwave-v45-812", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
    }
    summary["acoustic_impedance_absorption_phase_energy_flux_farfield_window_dataset_result_identity"] = {
        "generation": "acoustic-farfield-v45-812", "frequency_hz_generation": "acoustic-farfield-v45-812",
        "impedance_magnitude_pa_s_per_m_generation": "acoustic-farfield-v45-812", "impedance_phase_deg_generation": "acoustic-farfield-v45-812",
        "absorption_coefficient_generation": "acoustic-farfield-v45-812", "normal_energy_flux_w_generation": "acoustic-farfield-v45-812",
        "farfield_radius_m_generation": "acoustic-farfield-v45-812", "time_window_s_generation": "acoustic-farfield-v45-812",
        "dataset_tag_generation": "acoustic-farfield-v45-812", "frequency_hz": 125.0, "impedance_magnitude_pa_s_per_m": 2.5e5,
        "impedance_phase_deg": -35.0, "absorption_coefficient": 0.72, "normal_energy_flux_w": 0.4, "farfield_radius_m": 2.0,
        "time_window_s": 0.008, "dataset_tag": "dset-acoustic-v45-812", "result_frequency_hz": 125.0,
        "result_impedance_magnitude_pa_s_per_m": 2.5e5, "result_impedance_phase_deg": -35.0, "result_absorption_coefficient": 0.72,
        "result_normal_energy_flux_w": 0.4, "result_farfield_radius_m": 2.0, "result_time_window_s": 0.008, "result_dataset_tag": "dset-acoustic-v45-812",
        "owner": "model/acoustic-v45-812", "accepted_owner": "model/acoustic-v45-812", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    return summary


def _records_v47() -> dict[str, object]:
    rows_generation = "multi-output-rows-v47-901"
    cache_generation = "cache-chain-v47-901"
    keys = ["speed=1000|current=5", "speed=2000|current=5", "speed=3000|current=5"]
    chain = ["model:m1", "mesh:mesh1", "study:std1", "solution:sol1", "result:r1"]
    return {
        "force_torque_energy_parameter_row_key_identity": {
            "generation": rows_generation,
            "force_generation": rows_generation,
            "torque_generation": rows_generation,
            "energy_generation": rows_generation,
            "result_generation": rows_generation,
            "parameter_row_keys": keys,
            "force_parameter_row_keys": keys,
            "torque_parameter_row_keys": keys,
            "energy_parameter_row_keys": keys,
            "parameter_row_order_sha256": "1" * 64,
            "result_parameter_row_order_sha256": "1" * 64,
            "owner": "result/multi-output-v47-901",
            "accepted_owner": "result/multi-output-v47-901",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        "model_mesh_study_result_cache_owner_chain_identity": {
            "generation": cache_generation,
            "model_generation": cache_generation,
            "mesh_generation": cache_generation,
            "study_generation": cache_generation,
            "solution_generation": cache_generation,
            "result_generation": cache_generation,
            "cache_generation": cache_generation,
            "owner_chain": chain,
            "cached_result_owner_chain": chain,
            "model_mesh_study_result_sha256": "3" * 64,
            "cached_owner_chain_sha256": "3" * 64,
            "owner": "cache/result-v47-901",
            "accepted_owner": "cache/result-v47-901",
            "result_sha256": "4" * 64,
            "accepted_result_sha256": "4" * 64,
        },
    }


def _records_v48() -> dict[str, object]:
    ale_generation = "ale-force-v48"
    segregated_generation = "segregated-v48"
    groups = ["magnetic_vector_potential", "temperature", "displacement"]
    scaling = {"magnetic_vector_potential": 1.0, "temperature": 300.0, "displacement": 1.0e-3}
    iterations = ["iter=1|group=magnetic_vector_potential", "iter=1|group=temperature", "iter=1|group=displacement"]
    return {
        "ale_reference_current_force_quadrature_owner_identity": {
            "generation": ale_generation,
            "reference_mesh_generation": ale_generation,
            "current_mesh_generation": ale_generation,
            "quadrature_generation": ale_generation,
            "normal_generation": ale_generation,
            "result_generation": ale_generation,
            "reference_configuration_id": "ale/reference-v48",
            "result_reference_configuration_id": "ale/reference-v48",
            "current_configuration_id": "ale/current-v48",
            "result_current_configuration_id": "ale/current-v48",
            "quadrature_rule": "gauss-surface-order-4",
            "result_quadrature_rule": "gauss-surface-order-4",
            "normal_orientation_sha256": "1" * 64,
            "result_normal_orientation_sha256": "1" * 64,
            "body_owner": "body:moving-domain-v48",
            "result_body_owner": "body:moving-domain-v48",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        "segregated_variable_scaling_residual_iteration_solution_identity": {
            "generation": segregated_generation,
            "variable_group_generation": segregated_generation,
            "scaling_generation": segregated_generation,
            "residual_generation": segregated_generation,
            "iteration_generation": segregated_generation,
            "solution_generation": segregated_generation,
            "result_generation": segregated_generation,
            "variable_groups": groups,
            "result_variable_groups": groups,
            "variable_scaling": scaling,
            "result_variable_scaling": scaling,
            "residual_norm": "scaled_l2",
            "result_residual_norm": "scaled_l2",
            "iteration_rows": iterations,
            "result_iteration_rows": iterations,
            "solution_owner": "solution:segregated-v48",
            "result_solution_owner": "solution:segregated-v48",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def _records_v49() -> dict[str, object]:
    material_generation = "nonlinear-material-v49"
    sliding_generation = "sliding-interface-v49"
    return {
        "nonlinear_material_interpolation_branch_unit_temperature_extrapolation_owner_identity": {
            "generation": material_generation,
            "material_generation": material_generation,
            "branch_generation": material_generation,
            "temperature_generation": material_generation,
            "interpolation_generation": material_generation,
            "result_generation": material_generation,
            "interpolation_branch": "ascending-major-loop",
            "result_interpolation_branch": "ascending-major-loop",
            "input_units": {"magnetic_flux_density": "T", "magnetic_field": "A/m", "temperature": "K"},
            "result_input_units": {"magnetic_flux_density": "T", "magnetic_field": "A/m", "temperature": "K"},
            "temperature_value": 353.15,
            "result_temperature_value": 353.15,
            "extrapolation_policy": "reject",
            "result_extrapolation_policy": "reject",
            "material_owner": "material:nonlinear-steel-v49",
            "result_material_owner": "material:nonlinear-steel-v49",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        "moving_mesh_sliding_interface_frame_time_remesh_solution_owner_identity": {
            "generation": sliding_generation,
            "mesh_generation": sliding_generation,
            "interface_generation": sliding_generation,
            "frame_generation": sliding_generation,
            "time_generation": sliding_generation,
            "remesh_generation": sliding_generation,
            "solution_generation": sliding_generation,
            "result_generation": sliding_generation,
            "sliding_interface_map_sha256": "2" * 64,
            "result_sliding_interface_map_sha256": "2" * 64,
            "coordinate_frame": "spatial",
            "result_coordinate_frame": "spatial",
            "time_value_s": 0.0125,
            "result_time_value_s": 0.0125,
            "remesh_revision": "remesh-v49-r3",
            "result_remesh_revision": "remesh-v49-r3",
            "solution_owner": "solution:moving-mesh-v49",
            "result_solution_owner": "solution:moving-mesh-v49",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def _records_v50() -> dict[str, object]:
    frequency_generation = "frequency-sweep-v50-1001"
    contact_generation = "contact-v50-1001"
    frequencies = [100.0, 1000.0, 10000.0]
    gaps = [0.0, 1e-6, 2e-6]
    pressure = [2e6, 1e6, 0.0]
    return {
        "frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity": {
            "generation": frequency_generation,
            **{
                name: frequency_generation
                for name in (
                    "frequency_generation",
                    "branch_generation",
                    "phase_generation",
                    "unit_generation",
                    "dataset_generation",
                    "interpolation_generation",
                    "solution_generation",
                    "result_generation",
                )
            },
            "frequency_hz": frequencies,
            "result_frequency_hz": frequencies,
            "complex_branch": "positive_frequency",
            "result_complex_branch": "positive_frequency",
            "phase_convention": "exp(+jomega_t)",
            "result_phase_convention": "exp(+jomega_t)",
            "field_units": {"electric_field": "V/m", "magnetic_field": "A/m"},
            "result_field_units": {"electric_field": "V/m", "magnetic_field": "A/m"},
            "dataset_tag": "dataset:v50-frequency",
            "result_dataset_tag": "dataset:v50-frequency",
            "dataset_interpolation": "linear_complex",
            "result_dataset_interpolation": "linear_complex",
            "solution_owner": "solution:frequency-v50-1001",
            "result_solution_owner": "solution:frequency-v50-1001",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        "contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity": {
            "generation": contact_generation,
            **{
                name: contact_generation
                for name in (
                    "pair_generation",
                    "method_generation",
                    "penalty_generation",
                    "gap_generation",
                    "pressure_generation",
                    "frame_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "contact_pair_id": "pair:source-destination-v50",
            "result_contact_pair_id": "pair:source-destination-v50",
            "contact_method": "augmented_lagrange",
            "result_contact_method": "augmented_lagrange",
            "penalty_factor": 1e9,
            "result_penalty_factor": 1e9,
            "gap_m": gaps,
            "result_gap_m": gaps,
            "pressure_pa": pressure,
            "result_pressure_pa": pressure,
            "coordinate_frame": "spatial",
            "result_coordinate_frame": "spatial",
            "contact_owner": "contact:v50-1001",
            "result_contact_owner": "contact:v50-1001",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def _records_v51() -> dict[str, object]:
    eigen = "eigenmode-v51-1001"
    continuation = "continuation-v51-1001"
    frequencies = [1240.0, 1240.0]
    phase = {"dof": 17, "component": "real", "sign": "positive"}
    states = ["predictor:42", "corrector:42"]
    load_path = [0.0, 0.35, 0.7, 0.93, 0.88]
    return {
        "eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity": {
            "generation": eigen,
            **{name: eigen for name in ("frequency_generation", "normalization_generation", "phase_generation", "subspace_generation", "mesh_generation", "owner_generation", "result_generation")},
            "frequency_hz": frequencies,
            "result_frequency_hz": frequencies,
            "normalization": "unit_generalized_mass",
            "result_normalization": "unit_generalized_mass",
            "phase_anchor": phase,
            "result_phase_anchor": phase,
            "degenerate_subspace_basis_sha256": "1" * 64,
            "result_degenerate_subspace_basis_sha256": "1" * 64,
            "mesh_revision": "mesh:v51-r3",
            "result_mesh_revision": "mesh:v51-r3",
            "mode_owner": "mode-set:v51-1001",
            "result_mode_owner": "mode-set:v51-1001",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        "continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity": {
            "generation": continuation,
            **{name: continuation for name in ("branch_generation", "state_generation", "loadpath_generation", "turningpoint_generation", "owner_generation", "result_generation")},
            "branch_id": "branch:v51-primary",
            "result_branch_id": "branch:v51-primary",
            "predictor_corrector_states": states,
            "result_predictor_corrector_states": states,
            "load_path": load_path,
            "result_load_path": load_path,
            "turning_point_index": 3,
            "result_turning_point_index": 3,
            "solution_owner": "solution:continuation-v51-1001",
            "result_solution_owner": "solution:continuation-v51-1001",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def _generation(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": prefix, **{name: prefix for name in names}}


def _records_v52():
    adjoint = {
        **_generation("adj-v52-test", ("objective_generation", "scaling_generation", "conjugation_generation", "design_generation", "gradient_generation", "owner_generation", "result_generation")),
        "objective_tag": "obj_loss", "result_objective_tag": "obj_loss",
        "objective_scale": 0.01, "result_objective_scale": 0.01,
        "complex_adjoint_convention": "hermitian_conjugate", "result_complex_adjoint_convention": "hermitian_conjugate",
        "design_variable_order": ["x", "y"], "result_design_variable_order": ["x", "y"],
        "scaled_gradient": [1.0, -2.0], "result_scaled_gradient": [1.0, -2.0],
        "solution_owner": "solution:adj-v52-test", "result_solution_owner": "solution:adj-v52-test",
        "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
    }
    weak = {
        **_generation("weak-v52-test", ("testfunction_generation", "sign_generation", "orientation_generation", "measure_generation", "term_generation", "owner_generation", "result_generation")),
        "test_function": "test(u)", "result_test_function": "test(u)",
        "residual_sign": "lhs_minus_rhs", "result_residual_sign": "lhs_minus_rhs",
        "boundary_orientation": "outward_normal", "result_boundary_orientation": "outward_normal",
        "integration_measure": "surface_jacobian", "result_integration_measure": "surface_jacobian",
        "weak_terms": ["test(u)*u", "dot(grad(test(u)),grad(u))"],
        "result_weak_terms": ["test(u)*u", "dot(grad(test(u)),grad(u))"],
        "form_owner": "weak-form:v52-test", "result_form_owner": "weak-form:v52-test",
        "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
    }
    return {ADJOINT: adjoint, WEAK_FORM: weak}


def _records_v53() -> dict[str, object]:
    generation = "multiphysics-public-v53"
    generations = lambda names: {name: generation for name in names}
    pair = {"source": "left", "destination": "right"}
    return {
        ELECTROTHERMAL: {
            "generation": generation, **generations(("contact_generation", "electric_generation", "thermal_generation", "time_generation", "owner_generation", "result_generation")),
            "contact_resistance_ohm": 0.02, "result_contact_resistance_ohm": 0.02,
            "contact_current_a": 3.0, "result_contact_current_a": 3.0,
            "electric_power_w": 0.18, "result_electric_power_w": 0.18,
            "deposited_heat_w": 0.18, "result_deposited_heat_w": 0.18,
            "time_s": 0.5, "result_time_s": 0.5,
            "solution_owner": "solution:electrothermal-v53", "result_solution_owner": "solution:electrothermal-v53",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        FLOQUET: {
            "generation": generation, **generations(("phase_generation", "wavevector_generation", "pair_generation", "orientation_generation", "owner_generation", "result_generation")),
            "phase_rad": 0.6, "result_phase_rad": 0.6,
            "wave_vector_per_m": [20.0, 0.0, 0.0], "result_wave_vector_per_m": [20.0, 0.0, 0.0],
            "translation_m": [0.03, 0.0, 0.0], "result_translation_m": [0.03, 0.0, 0.0],
            "boundary_pair": pair, "result_boundary_pair": pair,
            "orientation": "source_to_destination", "result_orientation": "source_to_destination",
            "field_owner": "field:floquet-v53", "result_field_owner": "field:floquet-v53",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def _records_v54() -> dict[str, object]:
    piezo_generation = "multiphysics-public-piezo-v54"
    species_generation = "multiphysics-public-species-v54"
    fractions = {"H2": 0.10, "O2": 0.20, "H2O": 0.15, "N2": 0.55}
    rates = {"H2": -2.0, "O2": -1.0, "H2O": 2.0, "N2": 0.0}
    return {
        PIEZO: {
            "generation": piezo_generation,
            **{name: piezo_generation for name in ("electric_generation", "mechanical_generation", "phase_generation", "owner_generation", "result_generation")},
            "voltage_v": 10.0, "result_voltage_v": 10.0,
            "charge_c": 2.0e-6, "result_charge_c": 2.0e-6,
            "electric_work_j": 1.0e-5, "result_electric_work_j": 1.0e-5,
            "mechanical_work_j": 1.0e-5, "result_mechanical_work_j": 1.0e-5,
            "harmonic_phase_rad": 0.25, "result_harmonic_phase_rad": 0.25,
            "solution_owner": "solution:piezo-v54", "result_solution_owner": "solution:piezo-v54",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        SPECIES: {
            "generation": species_generation,
            **{name: species_generation for name in ("fraction_generation", "rate_generation", "flux_generation", "time_generation", "owner_generation", "result_generation")},
            "species_mass_fraction": fractions, "result_species_mass_fraction": fractions,
            "stoichiometric_rate_mol_m3_s": rates, "result_stoichiometric_rate_mol_m3_s": rates,
            "total_species_mass_rate_kg_s": 2.0e-5, "result_total_species_mass_rate_kg_s": 2.0e-5,
            "boundary_mass_flux_kg_s": -2.0e-5, "result_boundary_mass_flux_kg_s": -2.0e-5,
            "time_s": 0.02, "result_time_s": 0.02,
            "solution_owner": "solution:reacting-v54", "result_solution_owner": "solution:reacting-v54",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def _records_v55() -> dict[str, object]:
    thermo_generation = "multiphysics-public-thermo-v55"
    electrochem_generation = "multiphysics-public-electrochem-v55"
    frequency_hz = 12500.0
    decay_hz = -2.5
    quality_factor = frequency_hz / (-2.0 * decay_hz)
    stored_energy_j = 4.0e-6
    cycle_dissipation_j = 2.0 * 3.141592653589793 * stored_energy_j / quality_factor
    rates = {"Li": -1.0e-6, "Li_plus": 1.0e-6, "electron": 1.0e-6}
    charges = {"Li": 0, "Li_plus": 1, "electron": -1}
    return {
        THERMO: {
            "generation": thermo_generation,
            **{
                name: thermo_generation
                for name in (
                    "eigenfrequency_generation",
                    "energy_generation",
                    "dissipation_generation",
                    "normalization_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "complex_eigenfrequency_hz": [frequency_hz, decay_hz],
            "result_complex_eigenfrequency_hz": [frequency_hz, decay_hz],
            "stored_energy_j": stored_energy_j,
            "result_stored_energy_j": stored_energy_j,
            "cycle_dissipation_j": cycle_dissipation_j,
            "result_cycle_dissipation_j": cycle_dissipation_j,
            "quality_factor": quality_factor,
            "result_quality_factor": quality_factor,
            "modal_normalization": "unit_total_stored_energy",
            "result_modal_normalization": "unit_total_stored_energy",
            "solution_owner": "solution:thermoelastic-v55",
            "result_solution_owner": "solution:thermoelastic-v55",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        ELECTROCHEM: {
            "generation": electrochem_generation,
            **{
                name: electrochem_generation
                for name in (
                    "current_generation",
                    "species_generation",
                    "stoichiometry_generation",
                    "flux_generation",
                    "time_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "terminal_current_a": 0.09648533212,
            "result_terminal_current_a": 0.09648533212,
            "species_rate_mol_s": rates,
            "result_species_rate_mol_s": rates,
            "species_charge_number": charges,
            "result_species_charge_number": charges,
            "boundary_species_flux_mol_s": {"Li_plus": -1.0e-6},
            "result_boundary_species_flux_mol_s": {"Li_plus": -1.0e-6},
            "time_s": 0.5,
            "result_time_s": 0.5,
            "solution_owner": "solution:electrochem-v55",
            "result_solution_owner": "solution:electrochem-v55",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def _records_v56() -> dict[str, object]:
    modal_generation = "modal-public-v56"
    induction_generation = "induction-public-v56"
    modes = [{"mode": "mode:1", "frequency_hz": 500.0, "participation": 0.8, "normalized_energy": 0.8}, {"mode": "mode:2", "frequency_hz": 1000.0, "participation": 0.2, "normalized_energy": 0.2}]
    response = [{"frequency_hz": 500.0, "pressure_pa": 2.0}, {"frequency_hz": 1000.0, "pressure_pa": 0.5}]
    times = [0.0, 0.5, 1.0]
    temperatures = [293.15, 303.15, 313.15]
    return {
        MODAL: {
            "generation": modal_generation, **{field: modal_generation for field in ("participation_generation", "response_generation", "normalization_generation", "energy_generation", "owner_generation", "result_generation")},
            "modal_rows": modes, "result_modal_rows": modes,
            "frequency_response": response, "result_frequency_response": response,
            "normalization": "unit_total_modal_energy", "result_normalization": "unit_total_modal_energy",
            "total_normalized_energy": 1.0, "result_total_normalized_energy": 1.0,
            "solution_owner": "solution:modal-v56", "result_solution_owner": "solution:modal-v56",
            "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        INDUCTION: {
            "generation": induction_generation, **{field: induction_generation for field in ("input_generation", "joule_generation", "thermal_generation", "time_generation", "owner_generation", "result_generation")},
            "coil_input_energy_j": 100.0, "result_coil_input_energy_j": 100.0,
            "joule_heat_energy_j": 80.0, "result_joule_heat_energy_j": 80.0,
            "stored_thermal_energy_j": 15.0, "result_stored_thermal_energy_j": 15.0,
            "boundary_heat_loss_j": 5.0, "result_boundary_heat_loss_j": 5.0,
            "time_s": times, "result_time_s": times,
            "average_temperature_k": temperatures, "result_average_temperature_k": temperatures,
            "solution_owner": "solution:induction-v56", "result_solution_owner": "solution:induction-v56",
            "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }
