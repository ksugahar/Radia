from __future__ import annotations

import math
from copy import deepcopy

from radia_mcp.radia_ngsolve.adjoint_weakform_identity_v52 import (
    ADJOINT,
    WEAK_FORM,
    validate_public_v52_identity,
)
from radia_mcp.radia_ngsolve.conservation_identity_v54 import (
    PIEZO,
    SPECIES,
    validate_public_v54_identity,
)
from radia_mcp.radia_ngsolve.coupled_periodic_identity_v53 import (
    ELECTROTHERMAL,
    FLOQUET,
    validate_public_v53_identity,
)
from radia_mcp.radia_ngsolve.cross_artifact_lineage_v47 import (
    validate_public_identity,
)
from radia_mcp.radia_ngsolve.dissipation_reaction_identity_v55 import (
    ELECTROCHEM,
    THERMO,
    validate_public_v55_identity,
)
from radia_mcp.radia_ngsolve.frequency_contact_identity_v50 import (
    validate_public_v50_identity,
)
from radia_mcp.radia_ngsolve.modal_continuation_identity_v51 import (
    validate_public_v51_identity,
)
from radia_mcp.radia_ngsolve.multiphysics_identity_v56 import (
    INDUCTION,
    MODAL,
    validate_public_v56_identity,
)
from radia_mcp.radia_ngsolve.rotational_eddy_brake_energy_gate import (
    rotational_eddy_brake_energy_gate as gate,
)
from radia_mcp.radia_ngsolve.solver_state_identity_v49 import (
    validate_public_v49_identity,
)
from radia_mcp.radia_ngsolve.transform_normalization_v48 import (
    validate_public_v48_identity,
)
from test_rotational_eddy_brake_energy_gate import (
    _summary,
)

from _multiphysics_generalization_payloads import (
    _BATTERY_KEY,
    _BEARING_KEY,
    _INDUCTION_KEY,
    _MICROWAVE_KEY,
    _PIEZO_KEY,
    _POROELASTIC_KEY,
    _SPECIES_KEY,
    _THERMOACOUSTIC_KEY,
    _records_v47,
    _records_v48,
    _records_v49,
    _records_v50,
    _records_v51,
    _records_v52,
    _records_v53,
    _records_v54,
    _records_v55,
    _records_v56,
    _with_v29_sliding_and_radiation_identity,
    _with_v30_joule_and_eigenmode_identity,
    _with_v31_transform_and_force_identity,
    _with_v32_nonlinear_and_eigenmode_identity,
    _with_v33_contact_and_dae_identity,
    _with_v34_arclength_and_electrochemical_identity,
    _with_v35_multirate_and_adjoint_identity,
    _with_v36_force_and_modal_identity,
    _with_v37_capacitance_and_thermoelastic_identity,
    _with_v38_thermoviscous_and_piezoelectric_identity,
    _with_v39_poroelastic_and_induction_identity,
    _with_v40_thermoacoustic_and_battery_identity,
    _with_v41_piezoelectric_and_bearing_identity,
    _with_v42_induction_and_species_identity,
    _with_v43_microwave_and_poroelastic_identity,
    _with_v44,
    _with_v45,
)


def test_v29_public_rotating_sliding_interface_sector_pitch_azimuth_interpolation_torque_periodicity_mismatch() -> None:
    summary = _with_v29_sliding_and_radiation_identity(_summary())
    summary[
        "rotating_sliding_interface_sector_pitch_azimuth_interpolation_frame_periodicity_mesh_torque_generation_identity"
    ].update(
        {
            "sector_sliding_generation": "sliding-interface-160",
            "interpolation_sliding_generation": "sliding-interface-159",
            "result_sliding_generation": "sliding-interface-158",
            "result_sector_pitch_deg": 60.0,
            "result_sector_count": 10,
            "result_azimuth_origin_deg": 15.0,
            "result_source_interface_tag": "stator_if",
            "result_target_interface_tag": "rotor_if",
            "result_interpolation": "nearest_neighbor",
            "result_rotor_frame": "global_cartesian",
            "result_periodic_phase_deg": 180.0,
            "result_azimuth_samples_deg": [0.0, 10.0, 20.0],
            "result_torque_nm": [10.0, 8.0, 12.0],
            "result_sliding_mesh_sha256": "8" * 64,
            "accepted_torque_result_sha256": "9" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "rotating_sliding_interface_uses_current_sector_azimuth_interpolation_frame_periodicity_mesh_and_torque"
    ]


def test_v29_public_acoustic_radiation_impedance_modal_projection_reference_area_power_flux_mismatch() -> None:
    summary = _with_v29_sliding_and_radiation_identity(_summary())
    summary[
        "acoustic_radiation_impedance_modal_trace_reference_area_pressure_velocity_power_frequency_result_generation_identity"
    ].update(
        {
            "modal_radiation_generation": "acoustic-radiation-160",
            "convention_radiation_generation": "acoustic-radiation-159",
            "result_radiation_generation": "acoustic-radiation-158",
            "result_modal_basis_id": "p0_unnormalized",
            "result_mode_indices": [3, 2, 1],
            "result_trace_projection": "point_sample",
            "result_reference_area_m2": 1.0,
            "result_pressure_velocity_convention": "inward_positive_velocity",
            "result_frequency_grid_hz": [100.0, 300.0],
            "result_radiation_impedance_ri": [[-20.0, 5.0]],
            "result_outward_power_flux_w": [-1.0],
            "result_radiation_mesh_sha256": "a" * 64,
            "accepted_radiation_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "acoustic_radiation_uses_current_modes_trace_area_convention_frequency_power_and_result"
    ]


def test_v30_public_joule_heat_source_mapping_resistivity_temperature_time_average_energy_balance_mismatch() -> None:
    summary = _with_v30_joule_and_eigenmode_identity(_summary())
    summary[
        "joule_heat_source_current_density_resistivity_temperature_frame_time_average_energy_mesh_result_generation_identity"
    ].update(
        {
            "mapping_joule_generation": "joule-heat-closure-170",
            "averaging_joule_generation": "joule-heat-closure-169",
            "result_joule_generation": "joule-heat-closure-168",
            "result_current_density_field_id": "ec.J-old",
            "result_temperature_field_id": "ht.T-old",
            "result_resistivity_model_id": "rho-constant",
            "result_source_frame": "global_spatial",
            "result_averaging_window_s": [0.0, 0.01],
            "result_time_average_method": "final_sample",
            "result_electric_loss_w": 10.0,
            "result_heat_source_integral_w": 20.0,
            "result_coupled_mesh_sha256": "9" * 64,
            "accepted_joule_heat_result_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "joule_heat_uses_current_mapping_resistivity_temperature_frame_average_energy_mesh_and_result"
    ]


def test_v30_public_nonlinear_eigenmode_mac_branch_normalization_parameter_continuation_mismatch() -> None:
    summary = _with_v30_joule_and_eigenmode_identity(_summary())
    summary[
        "nonlinear_eigenmode_continuation_parameter_normalization_phase_mac_branch_eigenvalue_mesh_result_generation_identity"
    ].update(
        {
            "continuation_eigenmode_generation": "nonlinear-eigenmode-170",
            "mac_eigenmode_generation": "nonlinear-eigenmode-169",
            "result_eigenmode_generation": "nonlinear-eigenmode-168",
            "result_continuation_parameter_name": "temperature",
            "result_continuation_parameter_values": [1.0, 0.5, 0.0],
            "result_mode_normalization": "unit_max",
            "result_phase_anchor_dof": "base-x",
            "result_mac_reference_branch_ids": [2, 1],
            "result_mode_branch_ids": [[2, 1], [1, 2]],
            "result_eigenvalues_ri": [[[100.0, 0.0], [151.0, 0.0]]],
            "result_mac_assignment_sha256": "b" * 64,
            "result_eigenmode_mesh_sha256": "c" * 64,
            "accepted_eigenmode_result_sha256": "d" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_eigenmodes_use_current_continuation_normalization_phase_mac_branch_eigenvalues_mesh_and_result"
    ]


def test_v31_public_frequency_time_reconstruction_hermitian_window_group_delay_parseval_mismatch() -> None:
    summary = _with_v31_transform_and_force_identity(_summary())
    summary[
        "frequency_time_hermitian_spacing_window_group_delay_parseval_mesh_result_generation_identity"
    ].update(
        {
            "spectrum_transform_generation": "frequency-time-closure-180",
            "delay_transform_generation": "frequency-time-closure-179",
            "result_transform_generation": "frequency-time-closure-178",
            "result_frequencies_hz": [0.0, 90.0, 200.0],
            "result_frequency_spacing_hz": 90.0,
            "result_spectrum_ri": [[1.0, 0.0], [0.5, 0.25]],
            "result_hermitian_completion": "copy_without_conjugation",
            "result_window_name": "rectangular",
            "result_window_coherent_gain": 1.0,
            "result_group_delay_s": 0.0,
            "result_time_origin_s": -0.001,
            "time_energy": 3.0,
            "result_transform_mesh_sha256": "9" * 64,
            "accepted_time_trace_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "frequency_time_reconstruction_uses_current_hermitian_spacing_window_delay_parseval_mesh_and_result"
    ]


def test_v31_public_rotating_force_virtual_work_stress_tensor_phase_frame_torque_balance_mismatch() -> None:
    summary = _with_v31_transform_and_force_identity(_summary())
    summary[
        "rotating_force_virtual_work_stress_phase_frame_lever_angle_power_mesh_result_generation_identity"
    ].update(
        {
            "virtual_work_force_generation": "rotating-force-balance-180",
            "phase_force_generation": "rotating-force-balance-179",
            "result_force_generation": "rotating-force-balance-178",
            "result_phasor_convention": "exp_negative_jwt_peak",
            "result_coordinate_frame": "global_spatial",
            "result_lever_arm_m": [0.0, 0.05, 0.0],
            "result_mechanical_angles_rad": [0.0, 0.2, 0.1],
            "stress_tensor_torque_nm": [-1.0, -1.2, -1.1],
            "airgap_power_w": 95.0,
            "result_force_mesh_sha256": "b" * 64,
            "accepted_force_result_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "rotating_force_uses_current_virtual_work_stress_phase_frame_angle_power_mesh_and_result"
    ]


def test_v32_public_nonlinear_segregated_iteration_relaxation_residual_jacobian_continuation_mismatch() -> None:
    summary = _with_v32_nonlinear_and_eigenmode_identity(_summary())
    summary[
        "nonlinear_segregated_group_relaxation_residual_jacobian_continuation_mesh_result_generation_identity"
    ].update(
        {
            "segregated_group_generation": "nonlinear-segregated-190",
            "jacobian_generation": "nonlinear-segregated-189",
            "result_generation": "nonlinear-segregated-188",
            "result_segregated_group_order": ["thermal", "magnetic", "structural"],
            "result_relaxation_schedule": [1.0, 1.0, 1.0],
            "accepted_residual_norms": [1.0e-2, 4.0e-4, 9.0e-5],
            "accepted_continuation_parameter": 0.8,
            "accepted_continuation_unit": "percent",
            "accepted_jacobian_sha256": "c" * 64,
            "accepted_nonlinear_mesh_sha256": "d" * 64,
            "accepted_nonlinear_solution_sha256": "e" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_segregated_solutions_use_current_groups_relaxation_residual_jacobian_continuation_mesh_and_result"
    ]


def test_v32_public_degenerate_eigenmode_subspace_phase_normalization_participation_mass_mismatch() -> None:
    summary = _with_v32_nonlinear_and_eigenmode_identity(_summary())
    summary[
        "degenerate_eigenmode_subspace_phase_normalization_participation_mass_mesh_owner_result_generation_identity"
    ].update(
        {
            "subspace_generation": "degenerate-eigenmode-190",
            "normalization_generation": "degenerate-eigenmode-189",
            "result_generation": "degenerate-eigenmode-188",
            "result_eigenvalues": [100.0, 101.0],
            "result_degenerate_subspace_sha256": "f" * 64,
            "result_phase_anchor_dofs": [37, 12],
            "result_normalization": "max_component",
            "result_participation_factors": [0.6, -0.8],
            "result_effective_masses_kg": [0.36, 0.80],
            "result_mode_owner": "component2/acpr/eig1",
            "result_eigenmode_mesh_sha256": "0" * 64,
            "accepted_eigenmode_result_sha256": "1" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "degenerate_eigenmodes_use_current_subspace_phase_normalization_participation_mass_mesh_owner_and_result"
    ]


def test_v33_public_contact_complementarity_gap_pressure_active_set_friction_dissipation_mismatch() -> None:
    summary = _with_v33_contact_and_dae_identity(_summary())
    summary[
        "contact_gap_pressure_active_set_friction_dissipation_normal_mesh_result_generation_identity"
    ].update(
        {
            "gap_generation": "contact-complementarity-200",
            "active_set_generation": "contact-complementarity-199",
            "result_generation": "contact-complementarity-198",
            "result_contact_pair": "contact_old",
            "result_active_contact_ids": [7, 9],
            "result_normal_gap_m": [1.0e-3, -2.0e-4],
            "result_normal_pressure_pa": [-1.0e6, 2.0e6],
            "result_tangential_slip_m": [1.0e-3, 0.0],
            "result_friction_traction_pa": [4.0e5, 0.0],
            "result_contact_area_m2": [2.0e-4, 1.0e-4],
            "result_friction_coefficient": 0.1,
            "result_friction_dissipation_j": -0.01,
            "result_normal_orientation": "outward_master_to_slave",
            "result_contact_mesh_sha256": "9" * 64,
            "accepted_contact_result_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "contact_results_satisfy_current_complementarity_active_set_friction_dissipation_normal_mesh_and_result"
    ]


def test_v33_public_field_circuit_dae_charge_current_event_energy_constraint_mismatch() -> None:
    summary = _with_v33_contact_and_dae_identity(_summary())
    summary[
        "field_circuit_dae_charge_current_event_energy_time_dataset_result_generation_identity"
    ].update(
        {
            "charge_generation": "field-circuit-dae-200",
            "event_generation": "field-circuit-dae-199",
            "result_generation": "field-circuit-dae-198",
            "result_time_s": [0.0, 0.6e-3, 1.0e-3],
            "result_switch_event_time_s": 0.7e-3,
            "result_event_side": "left_limit_before_event",
            "result_charge_c": [0.0, -1.0e-6, 0.5e-6],
            "result_integrated_current_c": [0.0, 0.4e-6, 1.5e-6],
            "accepted_algebraic_residual_c": [0.0, 1.0e-5, 0.0],
            "result_current_sign_convention": "positive_out_of_field_device",
            "result_stored_energy_after_j": 0.0022,
            "result_switch_dissipation_j": -0.0002,
            "result_dataset_owner": "dset_old/sol1",
            "result_dae_dataset_sha256": "b" * 64,
            "accepted_dae_result_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "field_circuit_dae_results_use_current_charge_current_event_energy_time_dataset_and_result"
    ]


def test_v34_public_nonlinear_arclength_tangent_branch_turning_point_residual_owner_mismatch() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    summary[
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    ].update(
        {
            "arclength_generation": "arclength-210",
            "branch_generation": "arclength-209",
            "result_generation": "arclength-208",
            "result_previous_augmented_state": [0.0, 1.0, 1.0],
            "result_predictor_tangent": [-0.8, 0.0, 0.6],
            "result_arclength_step": -0.05,
            "result_predictor_augmented_state": [1.2, 0.0, 0.9],
            "result_corrected_augmented_state": [1.5, 0.0, 0.5],
            "result_branch_id": "lower_branch",
            "result_turning_point_side": "post_turn_negative_parameter_tangent",
            "result_corrected_residual_norm": 1.0e-2,
            "result_residual_tolerance": 1.0e-12,
            "result_continuation_mesh_sha256": "9" * 64,
            "accepted_continuation_result_sha256": "a" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_continuation_uses_current_arclength_tangent_branch_turning_residual_mesh_and_result"
    ]


def test_v34_public_electrochemical_species_flux_charge_mass_reaction_boundary_energy_mismatch() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    summary[
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    ].update(
        {
            "species_generation": "electrochemical-210",
            "reaction_generation": "electrochemical-209",
            "result_generation": "electrochemical-208",
            "result_species_order": ["C_neutral", "B_minus", "A_plus"],
            "result_charge_numbers": [0, -1, 1],
            "result_molar_mass_basis": [2.0, 1.0, 1.0],
            "result_reaction_stoichiometry": [1.0, -1.0, 0.0],
            "result_reaction_extent_mol": 0.8,
            "result_initial_inventory_mol": [0.0, 1.0, 1.0],
            "result_final_inventory_mol": [2.0, 0.1, 0.1],
            "result_integrated_boundary_flux_mol": [1.0, 0.0, 0.0],
            "result_integrated_electric_current_c": 96485.0,
            "result_initial_free_energy_j": 1.0,
            "result_final_free_energy_j": 2.0,
            "result_dissipated_free_energy_j": -1.0,
            "result_time_s": [1.0, 0.0],
            "result_electrochemical_mesh_sha256": "b" * 64,
            "accepted_electrochemical_result_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrochemical_results_use_current_species_flux_charge_mass_reaction_energy_time_mesh_and_result"
    ]


def test_v34_public_rejects_self_consistent_nonunit_arclength_tangent() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    identity = summary[
        "nonlinear_arclength_tangent_branch_turning_residual_mesh_result_generation_identity"
    ]
    identity["predictor_tangent"] = [0.6, 0.0, 0.6]
    identity["result_predictor_tangent"] = [0.6, 0.0, 0.6]
    assert gate(summary)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_nonstoichiometric_inventory() -> None:
    summary = _with_v34_arclength_and_electrochemical_identity(_summary())
    identity = summary[
        "electrochemical_species_flux_charge_mass_reaction_energy_time_mesh_result_generation_identity"
    ]
    identity["final_inventory_mol"] = [0.4, 0.5, 1.0]
    identity["result_final_inventory_mol"] = [0.4, 0.5, 1.0]
    assert gate(summary)["status"] == "needs_attention"


def test_v35_public_multirate_electromechanical_event_interpolation_work_power_timegrid_mismatch() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    summary[
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    ].update(
        {
            "event_generation": "multirate-coupling-220",
            "timegrid_generation": "multirate-coupling-219",
            "result_generation": "multirate-coupling-218",
            "result_electrical_time_s": [0.0, 0.0004, 0.001],
            "result_mechanical_time_s": [0.0, 0.0005, 0.001],
            "result_event_time_s": 0.0006,
            "result_event_interpolation_side": "left_continuous_before_event",
            "result_substep_owner": "mechanical_only",
            "result_coordinate_frame": "rotor_cylindrical",
            "result_electrical_input_energy_j": 0.010,
            "result_mechanical_output_work_j": 0.011,
            "result_dissipated_energy_j": -0.001,
            "result_energy_balance_tolerance_j": 1.0e-15,
            "result_coupling_mesh_sha256": "a" * 64,
            "accepted_coupling_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "multirate_electromechanical_results_use_current_event_timegrids_work_power_frame_mesh_and_result"
    ]


def test_v35_public_adjoint_sensitivity_objective_chainrule_constraint_fd_mesh_owner_mismatch() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    summary[
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    ].update(
        {
            "objective_generation": "adjoint-sensitivity-220",
            "fd_generation": "adjoint-sensitivity-219",
            "result_generation": "adjoint-sensitivity-218",
            "result_objective_tag": "mean_torque",
            "result_design_variable": "magnet_arc_deg",
            "result_design_scale": 10.0,
            "result_active_constraint": "none",
            "chainrule_gradient": -2.5,
            "finite_difference_gradient": 0.25,
            "result_gradient_tolerance": 1.0e-12,
            "result_fd_perturbation": 0.1,
            "result_sensitivity_mesh_sha256": "c" * 64,
            "result_primal_solution_sha256": "d" * 64,
            "accepted_gradient_result_sha256": "e" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adjoint_sensitivities_use_current_objective_design_chainrule_constraint_fd_mesh_solution_and_result"
    ]


def test_v35_public_rejects_self_consistent_multirate_energy_creation() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    identity = summary[
        "multirate_electromechanical_event_interpolation_work_power_timegrid_frame_mesh_result_generation_identity"
    ]
    identity["mechanical_output_work_j"] = 0.013
    identity["result_mechanical_output_work_j"] = 0.013
    assert gate(summary)["status"] == "needs_attention"


def test_v35_public_rejects_self_consistent_adjoint_fd_disagreement() -> None:
    summary = _with_v35_multirate_and_adjoint_identity(_summary())
    identity = summary[
        "adjoint_objective_design_chainrule_constraint_fd_mesh_solution_gradient_generation_identity"
    ]
    identity["finite_difference_gradient"] = 3.0
    assert gate(summary)["status"] == "needs_attention"


def test_v36_public_magnetostatic_virtual_work_coenergy_force_displacement_mesh_owner_mismatch() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    summary[
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    ].update(
        {
            "coenergy_generation": "virtual-work-coenergy-230",
            "mesh_generation": "virtual-work-coenergy-229",
            "result_generation": "virtual-work-coenergy-228",
            "result_displacement_m": [-2.0e-4, 0.0, 1.0e-4],
            "result_coenergy_j": [0.4999, 0.5, 0.5001],
            "result_held_source_convention": "constant_voltage",
            "result_force_sign_convention": "negative_dcoenergy_dx",
            "result_force_n": -3.0,
            "result_coordinate_frame": "rotor_cylindrical",
            "result_displaced_mesh_sha256": ["a" * 64, "b" * 64, "c" * 64],
            "result_force_solution_owner": "std_old/sol0",
            "accepted_force_result_sha256": "d" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetostatic_force_uses_current_displacement_coenergy_source_sign_mesh_frame_owner_and_result"
    ]


def test_v36_public_acoustic_modal_participation_effective_mass_damping_reconstruction_mismatch() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    summary[
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    ].update(
        {
            "normalization_generation": "acoustic-modal-participation-230",
            "result_generation": "acoustic-modal-participation-228",
            "result_normalization": "mass_normalized",
            "result_mode_frequency_hz": [160.0, 100.0],
            "result_modal_mass_kg": [1.0, 1.0],
            "result_participation_factor": [0.4, -0.5],
            "result_effective_modal_mass_kg": [0.16, 0.25],
            "result_damping_ratio": [-0.02, 0.01],
            "result_response_frequency_hz": [90.0, 125.0, 180.0],
            "result_probe_response_complex": [[0.0, 0.0]] * 3,
            "result_modal_mesh_sha256": "e" * 64,
            "accepted_modal_result_sha256": "f" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "acoustic_modes_use_current_normalization_mass_participation_damping_reconstruction_mesh_and_result"
    ]


def test_v36_public_rejects_self_consistent_wrong_coenergy_derivative() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    identity = summary[
        "magnetostatic_virtual_work_coenergy_force_displacement_current_mesh_frame_solution_result_generation_identity"
    ]
    identity["coenergy_j"] = [0.4999, 0.5, 0.5001]
    identity["result_coenergy_j"] = [0.4999, 0.5, 0.5001]
    assert gate(summary)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_wrong_effective_mass() -> None:
    summary = _with_v36_force_and_modal_identity(_summary())
    identity = summary[
        "acoustic_modal_normalization_effective_mass_participation_damping_frequency_reconstruction_mesh_result_generation_identity"
    ]
    identity["effective_modal_mass_kg"] = [0.6, 0.24]
    identity["result_effective_modal_mass_kg"] = [0.6, 0.24]
    assert gate(summary)["status"] == "needs_attention"


def test_v37_public_capacitance_matrix_charge_energy_gauge_reciprocity_terminal_owner_mismatch() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "matrix_generation": "capacitance-closure-240",
            "result_reference_terminal": "terminal_1",
            "result_capacitance_matrix_f": [[2.0e-12, 1.0e-12], [-2.0e-12, 1.0e-12]],
            "result_terminal_charge_c": [1.0e-12, 1.0e-12],
            "result_stored_energy_j": -1.0e-12,
            "accepted_terminal_owner": "comp1/es/old-terminals",
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "capacitance_results_use_current_matrix_charge_energy_gauge_reciprocity_terminals_mesh_and_result"
    ]


def test_v37_public_thermoelastic_harmonic_phase_loss_work_temperature_displacement_mesh_mismatch() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "phase_generation": "thermoelastic-harmonic-240",
            "result_frequency_hz": -1000.0,
            "result_heat_source_phase_rad": -2.0,
            "result_temperature_complex_k": [-5.0, -1.0],
            "result_displacement_complex_m": [-1.0e-6, 2.0e-7],
            "result_thermal_expansion_work_j": -2.0e-3,
            "result_mechanical_loss_j": -1.0e-4,
            "accepted_mesh_sha256": "c" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thermoelastic_harmonics_use_current_heat_phase_temperature_displacement_work_loss_frequency_mesh_and_result"
    ]


def test_v37_public_rejects_self_consistent_nonzero_capacitance_row_sum() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "capacitance_matrix_charge_energy_gauge_reciprocity_terminal_mesh_result_generation_identity"
    ]
    wrong = [[2.0e-12, -1.0e-12], [-1.0e-12, 2.0e-12]]
    identity["capacitance_matrix_f"] = wrong
    identity["result_capacitance_matrix_f"] = wrong
    assert gate(summary)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_negative_thermoelastic_loss() -> None:
    summary = _with_v37_capacitance_and_thermoelastic_identity(_summary())
    identity = summary[
        "thermoelastic_harmonic_heat_phase_temperature_displacement_work_loss_frequency_mesh_result_generation_identity"
    ]
    identity["mechanical_loss_j"] = -1.0e-4
    identity["result_mechanical_loss_j"] = -1.0e-4
    assert gate(summary)["status"] == "needs_attention"


def test_v38_public_thermoviscous_pressure_acoustics_interface_velocity_traction_dissipation_power_mismatch() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "velocity_generation": "thermoviscous-interface-257",
            "power_generation": "thermoviscous-interface-256",
            "result_generation": "thermoviscous-interface-255",
            "result_frequency_hz": -1.0e4,
            "result_normal_velocity_complex_m_per_s": [-1.0e-2, 2.0e-3],
            "result_pressure_complex_pa": [200.0, -40.0],
            "result_traction_sign": "plus_pressure_times_normal",
            "result_normal_orientation": "pressure_to_thermoviscous",
            "result_interface_power_w": -1.04e-2,
            "result_viscous_loss_w": -3.0e-3,
            "result_outgoing_acoustic_power_w": 1.0e-1,
            "accepted_mesh_owner": "comp1/mesh0:old",
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thermoviscous_interfaces_use_current_velocity_traction_dissipation_power_normal_mesh_and_result"
    ]


def test_v38_public_piezoelectric_charge_strain_reciprocity_electromechanical_energy_terminal_mismatch() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "reciprocity_generation": "piezoelectric-reciprocity-257",
            "coupling_energy_generation": "piezoelectric-reciprocity-256",
            "result_generation": "piezoelectric-reciprocity-255",
            "result_direct_coefficient_c_per_n": -2.0e-10,
            "result_converse_coefficient_m_per_v": 4.0e-10,
            "result_induced_strain": -2.0e-5,
            "result_induced_charge_density_c_per_m2": 1.0e-3,
            "result_terminal_charge_c": -4.0e-6,
            "result_coupling_energy_j": -1.0e-2,
            "result_total_stored_energy_j": 6.0e-2,
            "result_polarization_frame": "global_z",
            "accepted_mesh_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "piezoelectric_results_use_current_charge_strain_reciprocity_energy_polarization_mesh_and_result"
    ]


def test_v38_public_rejects_self_consistent_negative_thermal_loss() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    ]
    identity["thermal_loss_w"] = identity["result_thermal_loss_w"] = -2.0e-3
    identity["outgoing_acoustic_power_w"] = identity[
        "result_outgoing_acoustic_power_w"
    ] = 9.4e-3
    assert gate(summary)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_nonreciprocal_piezo_coefficients() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    ]
    identity["converse_coefficient_m_per_v"] = identity[
        "result_converse_coefficient_m_per_v"
    ] = 4.0e-10
    identity["induced_strain"] = identity["result_induced_strain"] = 4.0e-5
    assert gate(summary)["status"] == "needs_attention"


def test_v39_public_poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mismatch() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"]
    identity.update({
        "flux_generation": "poroelastic-biot-376", "result_generation": "poroelastic-biot-375",
        "result_fluid_content_increment": -1.0, "result_darcy_flux_m_per_s": 1.0e-3,
        "result_interface_traction_pa": 8.0e4, "result_storage_energy_j": -1.0,
        "result_skeleton_coupling_work_j": -1.0, "result_fluid_dissipation_j": -1.0,
        "result_interface_normal": "free_fluid_to_skeleton",
        "accepted_mesh_owner": "comp1/mesh0:old", "accepted_result_sha256": "9" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["poroelastic_results_use_current_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_and_result"]


def test_v39_public_rotating_induction_slip_frequency_current_loss_torque_power_frame_mismatch() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"]
    identity.update({
        "slip_generation": "rotating-induction-376", "result_generation": "rotating-induction-375",
        "result_slip": -1.0, "result_rotor_electrical_frequency_hz": -50.0,
        "result_rotor_phase_current_a_rms": -1.0, "result_rotor_copper_loss_w": -1.0,
        "result_airgap_torque_nm": -20.0, "result_mechanical_power_w": 100.0,
        "result_rotating_frame": "stator_frame", "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["rotating_induction_results_use_current_slip_frequency_current_loss_torque_power_frame_mesh_and_result"]


def test_v39_public_rejects_self_consistent_wrong_darcy_flux() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["poroelastic_biot_pressure_displacement_flux_storage_dissipation_interface_mesh_result_generation_identity"]
    identity["darcy_flux_m_per_s"] = identity["result_darcy_flux_m_per_s"] = 1.0e-3
    assert gate(summary)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_induction_power_creation() -> None:
    summary = _with_v39_poroelastic_and_induction_identity(_summary())
    identity = summary["rotating_induction_slip_frequency_current_loss_torque_power_frame_mesh_result_generation_identity"]
    identity["airgap_power_w"] = identity["result_airgap_power_w"] = 1.0
    assert gate(summary)["status"] == "needs_attention"


def test_v40_public_thermoacoustic_meanflow_convected_wavenumber_flux_impedance_power_mismatch() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_THERMOACOUSTIC_KEY]
    identity.update({
        "wavenumber_generation": "thermoacoustic-718",
        "result_generation": "thermoacoustic-717",
        "result_mean_flow_mach": -2.0,
        "result_downstream_wavenumber_rad_per_m": -1.0,
        "result_acoustic_intensity_w_per_m2": -1.0,
        "result_boundary_impedance_pa_s_per_m": -1.0,
        "result_boundary_flux_power_w": -1.0,
        "result_power_balance_residual_w": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["thermoacoustic_results_use_current_meanflow_wavenumber_flux_impedance_power_mesh_and_result"]


def test_v40_public_battery_electrothermal_soc_current_heat_temperature_energy_safety_mismatch() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_BATTERY_KEY]
    identity.update({
        "heat_generation": "battery-electrothermal-718",
        "result_generation": "battery-electrothermal-717",
        "result_final_state_of_charge": 2.0,
        "result_terminal_current_a": -2.0,
        "result_irreversible_heat_j": -1.0,
        "result_thermal_energy_j": -1.0,
        "result_final_temperature_k": 500.0,
        "result_thermal_balance_residual_j": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["battery_results_use_current_soc_current_heat_temperature_energy_safety_mesh_and_result"]


def test_v40_public_rejects_self_consistent_nonconvected_wavenumber() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_THERMOACOUSTIC_KEY]
    value = 2.0 * math.pi * identity["frequency_hz"] / identity["sound_speed_m_per_s"]
    identity["downstream_wavenumber_rad_per_m"] = value
    identity["result_downstream_wavenumber_rad_per_m"] = value
    assert gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_battery_heat_creation() -> None:
    summary = _with_v40_thermoacoustic_and_battery_identity(_summary())
    identity = summary[_BATTERY_KEY]
    identity["thermal_energy_j"] = identity["result_thermal_energy_j"] = 1.0
    identity["thermal_balance_residual_j"] = identity["result_thermal_balance_residual_j"] = 0.0
    assert gate(summary)["status"] == "needs_attention"


def test_v41_public_piezoelectric_admittance_resonance_antiresonance_coupling_energy_phase_mismatch() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    summary[_PIEZO_KEY].update({
        "coupling_generation": "piezoelectric-723",
        "result_generation": "piezoelectric-722",
        "result_antiresonance_frequency_hz": 90_000.0,
        "result_electromechanical_coupling_squared": -1.0,
        "result_admittance_phase_deg": 150.0,
        "result_real_electrical_power_w": -2.0,
        "result_power_balance_residual_w": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["piezoelectric_admittance_results_use_current_resonance_coupling_phase_energy_power_mesh_and_result"]


def test_v41_public_fluidfilm_bearing_reynolds_pressure_load_friction_temperature_power_mismatch() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    summary[_BEARING_KEY].update({
        "pressure_generation": "fluidfilm-bearing-723",
        "result_generation": "fluidfilm-bearing-722",
        "result_minimum_film_thickness_m": -1.0,
        "result_maximum_pressure_pa": -1.0,
        "result_integrated_load_n": -1.0,
        "result_friction_torque_nm": -2.0,
        "result_maximum_temperature_k": 250.0,
        "result_power_balance_residual_w": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fluidfilm_bearing_results_use_current_film_pressure_load_friction_temperature_power_mesh_and_result"]


def test_v41_public_rejects_self_consistent_antiresonance_below_resonance() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    identity = summary[_PIEZO_KEY]
    antiresonance = 90_000.0
    coupling = 1.0 - (identity["resonance_frequency_hz"] / antiresonance) ** 2
    identity["antiresonance_frequency_hz"] = identity["result_antiresonance_frequency_hz"] = antiresonance
    identity["electromechanical_coupling_squared"] = identity["result_electromechanical_coupling_squared"] = coupling
    assert gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_overclosed_bearing_film() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    identity = summary[_BEARING_KEY]
    identity["eccentricity_ratio"] = identity["result_eccentricity_ratio"] = 1.2
    thickness = identity["radial_clearance_m"] * (1.0 - 1.2)
    identity["minimum_film_thickness_m"] = identity["result_minimum_film_thickness_m"] = thickness
    assert gate(summary)["status"] == "needs_attention"


def test_v42_public_inductionheating_skin_proximity_joule_thermal_flux_temperature_energy_mismatch() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    summary[_INDUCTION_KEY].update({
        "skin_generation": "induction-heating-724",
        "result_generation": "induction-heating-723",
        "result_skin_depth_m": -1.0,
        "result_proximity_current_density_a_per_m2": [-1.0],
        "result_joule_loss_w": -10.0,
        "result_outward_thermal_flux_w": 10.0,
        "result_maximum_temperature_k": 250.0,
        "result_electromagnetic_power_balance_residual_w": 99.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["induction_heating_results_use_current_skin_proximity_joule_thermal_temperature_energy_mesh_and_result"]


def test_v42_public_species_transport_reaction_diffusion_flux_massbalance_rate_temperature_mismatch() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    summary[_SPECIES_KEY].update({
        "reaction_generation": "reacting-species-724",
        "result_generation": "reacting-species-723",
        "result_diffusivity_m2_per_s": -1.0,
        "result_temperature_k": 0.0,
        "result_reaction_rate_constant_per_s": -1.0,
        "result_integrated_species_mol": -1.0,
        "result_inward_boundary_flux_mol_per_s": -1.0,
        "result_mass_balance_residual_mol_per_s": 1.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["reacting_species_results_use_current_diffusion_rate_flux_mass_temperature_mesh_and_result"]


def test_v42_public_rejects_self_consistent_wrong_skin_depth() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    identity = summary[_INDUCTION_KEY]
    wrong = 2.0 * identity["skin_depth_m"]
    identity["skin_depth_m"] = identity["result_skin_depth_m"] = wrong
    assert gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_species_mass_leak() -> None:
    summary = _with_v42_induction_and_species_identity(_summary())
    identity = summary[_SPECIES_KEY]
    consumption = identity["integrated_consumption_mol_per_s"]
    identity["inward_boundary_flux_mol_per_s"] = 2.0 * consumption
    identity["result_inward_boundary_flux_mol_per_s"] = 2.0 * consumption
    identity["mass_balance_residual_mol_per_s"] = consumption
    identity["result_mass_balance_residual_mol_per_s"] = consumption
    assert gate(summary)["status"] == "needs_attention"


def test_v43_public_microwaveheating_sparameter_absorbedpower_jouleheat_temperature_energy_mismatch() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    summary[_MICROWAVE_KEY].update({
        "power_generation": "microwave-heating-725",
        "result_generation": "microwave-heating-724",
        "result_s11_magnitude": 1.5,
        "result_absorbed_power_w": -10.0,
        "result_electromagnetic_power_residual_w": 95.0,
        "result_thermal_power_residual_w": 85.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["microwave_heating_results_use_current_sparameters_power_heat_temperature_energy_mesh_and_result"]


def test_v43_public_poroelastic_wave_pressure_displacement_flux_dissipation_mass_energy_mismatch() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    summary[_POROELASTIC_KEY].update({
        "flux_generation": "poroelastic-wave-725",
        "result_generation": "poroelastic-wave-724",
        "result_porosity": 1.4,
        "result_mass_balance_residual_kg_per_s": 1.0,
        "result_energy_balance_residual_w": 2.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["poroelastic_wave_results_use_current_pressure_displacement_flux_mass_energy_mesh_and_result"]


def test_v43_public_rejects_self_consistent_nonpassive_sparameters() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    identity = summary[_MICROWAVE_KEY]
    identity["s11_magnitude"] = identity["result_s11_magnitude"] = 1.1
    assert gate(summary)["status"] == "needs_attention"


def test_v43_public_rejects_self_consistent_invalid_poroelastic_phase() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    identity = summary[_POROELASTIC_KEY]
    identity["pressure_displacement_phase_deg"] = 220.0
    identity["result_pressure_displacement_phase_deg"] = 220.0
    assert gate(summary)["status"] == "needs_attention"


def test_v44_public_positive_port_and_acoustic_closure() -> None:
    assert gate(_with_v44(_summary()))["status"] == "ok"


def test_v44_public_rejects_port_power_mismatch() -> None:
    summary = _with_v44(_summary())
    summary["microwave_boundaryport_sparameter_power_normalization_temperature_coupling_restart_owner_result_identity"]["result_absorbed_power_w"] = -1.0
    assert gate(summary)["status"] == "needs_attention"


def test_v44_public_rejects_acoustic_phase_mismatch() -> None:
    summary = _with_v44(_summary())
    summary["acoustics_poroelastic_impedance_phase_flux_energy_timewindow_dataset_owner_result_identity"]["result_impedance_phase_deg"] = 45.0
    assert gate(summary)["status"] == "needs_attention"


def test_v45_public_positive_identity() -> None:
    assert gate(_with_v45(_summary()))["status"] == "ok"


def test_v45_public_rejects_port_deembed_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["microwave_sparameter_port_reference_plane_deembed_complex_power_mesh_result_identity"]["result_deembed_length_m"] = -1.0
    assert gate(summary)["status"] == "needs_attention"


def test_v45_public_rejects_acoustic_window_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["acoustic_impedance_absorption_phase_energy_flux_farfield_window_dataset_result_identity"]["result_time_window_s"] = -1.0
    assert gate(summary)["status"] == "needs_attention"


def test_v47_positive_replays_are_accepted() -> None:
    result = validate_public_identity(_records_v47())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v47_row_permutation_is_rejected() -> None:
    records = _records_v47()
    row = records["force_torque_energy_parameter_row_key_identity"]
    row["torque_parameter_row_keys"] = list(reversed(row["parameter_row_keys"]))
    assert validate_public_identity(records)["status"] == "needs_attention"


def test_v47_stale_cache_chain_is_rejected() -> None:
    records = _records_v47()
    row = records["model_mesh_study_result_cache_owner_chain_identity"]
    row["study_generation"] = "old"
    row["cached_result_owner_chain"] = ["model:m1", "mesh:mesh0", "study:std0", "solution:sol0", "result:r1"]
    assert validate_public_identity(records)["status"] == "needs_attention"


def test_v48_positive_replays_are_accepted() -> None:
    result = validate_public_v48_identity(_records_v48())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v48_mixed_ale_configuration_is_rejected() -> None:
    records = deepcopy(_records_v48())
    row = records["ale_reference_current_force_quadrature_owner_identity"]
    row["result_current_configuration_id"] = "ale/current-old"
    row["result_quadrature_rule"] = "gauss-surface-order-2"
    row["result_normal_orientation_sha256"] = "a" * 64
    row["result_body_owner"] = "body:fixed-old"
    assert validate_public_v48_identity(records)["status"] == "needs_attention"


def test_v48_mixed_segregated_solver_rows_are_rejected() -> None:
    records = deepcopy(_records_v48())
    row = records["segregated_variable_scaling_residual_iteration_solution_identity"]
    row["result_variable_groups"] = list(reversed(row["variable_groups"]))
    row["result_variable_scaling"]["temperature"] = 1.0
    row["result_residual_norm"] = "unscaled_linf"
    row["result_iteration_rows"] = row["iteration_rows"][:-1]
    row["result_solution_owner"] = "solution:segregated-old"
    assert validate_public_v48_identity(records)["status"] == "needs_attention"


def test_v49_public_positive_replay_is_accepted() -> None:
    result = validate_public_v49_identity(_records_v49())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v49_public_mixed_material_state_is_rejected() -> None:
    records = deepcopy(_records_v49())
    row = records["nonlinear_material_interpolation_branch_unit_temperature_extrapolation_owner_identity"]
    row["result_interpolation_branch"] = "descending-recoil-loop"
    row["result_input_units"] = {"temperature": "degC"}
    row["result_temperature_value"] = 80.0
    row["result_extrapolation_policy"] = "linear"
    row["result_material_owner"] = "material:old"
    assert validate_public_v49_identity(records)["status"] == "needs_attention"


def test_v49_public_mixed_sliding_mesh_state_is_rejected() -> None:
    records = deepcopy(_records_v49())
    row = records["moving_mesh_sliding_interface_frame_time_remesh_solution_owner_identity"]
    row["result_sliding_interface_map_sha256"] = "a" * 64
    row["result_coordinate_frame"] = "material"
    row["result_time_value_s"] = 0.01
    row["result_remesh_revision"] = "remesh-v49-r2"
    row["result_solution_owner"] = "solution:old"
    assert validate_public_v49_identity(records)["status"] == "needs_attention"


def test_v50_public_positive_replay_is_accepted() -> None:
    result = validate_public_v50_identity(_records_v50())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v50_public_mixed_frequency_sweep_is_rejected() -> None:
    records = deepcopy(_records_v50())
    row = records["frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity"]
    row.update(
        {
            "result_frequency_hz": [10000.0, 1000.0, 100.0],
            "result_complex_branch": "negative_frequency",
            "result_phase_convention": "exp(-jomega_t)",
            "result_field_units": {"electric_field": "mV/m", "magnetic_field": "A/m"},
            "result_dataset_tag": "dataset:old",
            "result_dataset_interpolation": "nearest_magnitude",
            "result_solution_owner": "solution:old",
        }
    )
    assert validate_public_v50_identity(records)["status"] == "needs_attention"


def test_v50_public_mixed_contact_state_is_rejected() -> None:
    records = deepcopy(_records_v50())
    row = records["contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity"]
    row.update(
        {
            "result_contact_pair_id": "pair:reversed",
            "result_contact_method": "penalty",
            "result_penalty_factor": 1e6,
            "result_gap_m": [2e-6, 1e-6, 0.0],
            "result_pressure_pa": [0.0, 1e6, 2e6],
            "result_coordinate_frame": "material",
            "result_contact_owner": "contact:old",
        }
    )
    assert validate_public_v50_identity(records)["status"] == "needs_attention"


def test_v50_public_invalid_canonical_sequences_are_rejected() -> None:
    records = deepcopy(_records_v50())
    records["frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity"]["frequency_hz"] = [1000.0, 100.0]
    records["contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity"]["pressure_pa"] = [0.0, 1e6, 2e6]
    assert validate_public_v50_identity(records)["status"] == "needs_attention"


def test_v51_public_positive_replay_is_accepted() -> None:
    result = validate_public_v51_identity(_records_v51())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v51_public_mixed_eigenmode_identity_is_rejected() -> None:
    records = deepcopy(_records_v51())
    row = records["eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity"]
    row.update({"result_frequency_hz": [1235.0, 1245.0], "result_normalization": "unit_peak", "result_phase_anchor": {"dof": 21, "component": "imag", "sign": "negative"}, "result_degenerate_subspace_basis_sha256": "a" * 64, "result_mesh_revision": "mesh:v51-r2", "result_mode_owner": "mode-set:stale"})
    assert validate_public_v51_identity(records)["status"] == "needs_attention"


def test_v51_public_mixed_continuation_identity_is_rejected() -> None:
    records = deepcopy(_records_v51())
    row = records["continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity"]
    row.update({"result_branch_id": "branch:v51-secondary", "result_predictor_corrector_states": ["predictor:41", "corrector:41"], "result_load_path": [0.0, 0.35, 0.7, 0.88], "result_turning_point_index": 2, "result_solution_owner": "solution:continuation-stale"})
    assert validate_public_v51_identity(records)["status"] == "needs_attention"


def test_v51_public_invalid_canonical_records_are_rejected() -> None:
    records = deepcopy(_records_v51())
    records["eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity"]["phase_anchor"] = {"dof": -1, "component": "real", "sign": "positive"}
    records["continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity"]["turning_point_index"] = 2
    assert validate_public_v51_identity(records)["status"] == "needs_attention"


def test_v52_public_positive_replay_is_accepted():
    assert validate_public_v52_identity(_records_v52())["status"] == "ok"


def test_v52_public_mixed_adjoint_identity_is_rejected():
    value = deepcopy(_records_v52())
    value[ADJOINT]["result_complex_adjoint_convention"] = "transpose_without_conjugation"
    assert validate_public_v52_identity(value)["status"] == "needs_attention"


def test_v52_public_mixed_weakform_identity_is_rejected():
    value = deepcopy(_records_v52())
    value[WEAK_FORM]["result_boundary_orientation"] = "inward_normal"
    assert validate_public_v52_identity(value)["status"] == "needs_attention"


def test_v53_public_positive_replay_is_accepted() -> None:
    assert validate_public_v53_identity(_records_v53())["status"] == "ok"


def test_v53_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records_v53())
    value[ELECTROTHERMAL].update({"result_electric_power_w": 0.36, "result_deposited_heat_w": 0.12, "result_solution_owner": "solution:stale"})
    value[FLOQUET].update({"result_phase_rad": -0.6, "result_boundary_pair": {"source": "right", "destination": "left"}, "result_field_owner": "field:stale"})
    assert validate_public_v53_identity(value)["status"] == "needs_attention"


def test_v53_self_consistent_energy_and_phase_errors_are_rejected() -> None:
    value = deepcopy(_records_v53())
    value[ELECTROTHERMAL]["electric_power_w"] = value[ELECTROTHERMAL]["result_electric_power_w"] = 0.12
    value[ELECTROTHERMAL]["deposited_heat_w"] = value[ELECTROTHERMAL]["result_deposited_heat_w"] = 0.12
    value[FLOQUET]["phase_rad"] = value[FLOQUET]["result_phase_rad"] = -0.6
    assert validate_public_v53_identity(value)["status"] == "needs_attention"


def test_v54_public_positive_replay_is_accepted() -> None:
    assert validate_public_v54_identity(_records_v54())["status"] == "ok"


def test_v54_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records_v54())
    value[PIEZO].update({"result_voltage_v": 5.0, "result_electric_work_j": 2.0e-5, "result_mechanical_work_j": 8.0e-6, "result_harmonic_phase_rad": -0.25, "result_solution_owner": "solution:stale"})
    value[SPECIES].update({"result_species_mass_fraction": {"H2": 0.8, "O2": 0.8}, "result_boundary_mass_flux_kg_s": 2.0e-5, "result_time_s": 0.01, "result_solution_owner": "solution:stale"})
    assert validate_public_v54_identity(value)["status"] == "needs_attention"


def test_v54_self_consistent_conservation_errors_are_rejected() -> None:
    value = deepcopy(_records_v54())
    value[PIEZO]["electric_work_j"] = value[PIEZO]["result_electric_work_j"] = 2.0e-5
    value[SPECIES]["species_mass_fraction"] = value[SPECIES]["result_species_mass_fraction"] = {"H2": 0.8, "O2": 0.8}
    value[SPECIES]["boundary_mass_flux_kg_s"] = value[SPECIES]["result_boundary_mass_flux_kg_s"] = 2.0e-5
    assert validate_public_v54_identity(value)["status"] == "needs_attention"


def test_v55_public_positive_replay_is_accepted() -> None:
    assert validate_public_v55_identity(_records_v55())["status"] == "ok"


def test_v55_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records_v55())
    value[THERMO].update(
        {
            "result_complex_eigenfrequency_hz": [12500.0, 2.5],
            "result_stored_energy_j": 8.0e-6,
            "result_cycle_dissipation_j": 0.0,
            "result_modal_normalization": "peak_displacement",
            "result_solution_owner": "solution:stale",
        }
    )
    value[ELECTROCHEM].update(
        {
            "result_terminal_current_a": 1.0,
            "result_species_rate_mol_s": {"Li_plus": -1.0},
            "result_species_charge_number": {"Li_plus": -1},
            "result_boundary_species_flux_mol_s": {"Li_plus": 1.0},
            "result_time_s": 0.25,
            "result_solution_owner": "solution:stale",
        }
    )
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v55_self_consistent_nonphysical_damping_is_rejected() -> None:
    value = deepcopy(_records_v55())
    value[THERMO]["complex_eigenfrequency_hz"] = value[THERMO][
        "result_complex_eigenfrequency_hz"
    ] = [12500.0, 2.5]
    value[THERMO]["cycle_dissipation_j"] = value[THERMO][
        "result_cycle_dissipation_j"
    ] = 0.0
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v55_self_consistent_charge_or_flux_imbalance_is_rejected() -> None:
    value = deepcopy(_records_v55())
    bad_rates = {"Li": -1.0e-6, "Li_plus": 1.0e-6, "electron": 0.5e-6}
    value[ELECTROCHEM]["species_rate_mol_s"] = value[ELECTROCHEM][
        "result_species_rate_mol_s"
    ] = bad_rates
    value[ELECTROCHEM]["boundary_species_flux_mol_s"] = value[ELECTROCHEM][
        "result_boundary_species_flux_mol_s"
    ] = {"Li_plus": 1.0e-6}
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v55_numeric_sha256_values_are_rejected() -> None:
    value = _records_v55()
    numeric_digest = int("9" * 64)
    for row in value.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert validate_public_v55_identity(value)["status"] == "needs_attention"


def test_v56_public_positive_replay_is_accepted() -> None:
    assert validate_public_v56_identity(_records_v56())["status"] == "ok"


def test_v56_frozen_public_mutations_are_rejected() -> None:
    value = deepcopy(_records_v56())
    value[MODAL].update({"result_modal_rows": [], "result_normalization": "peak_pressure", "result_solution_owner": "solution:stale"})
    value[INDUCTION].update({"result_coil_input_energy_j": 50.0, "result_time_s": [0.0, 2.0], "result_solution_owner": "solution:stale"})
    assert validate_public_v56_identity(value)["status"] == "needs_attention"


def test_v56_self_consistent_energy_and_history_contradictions_are_rejected() -> None:
    value = deepcopy(_records_v56())
    value[MODAL]["total_normalized_energy"] = value[MODAL]["result_total_normalized_energy"] = 2.0
    value[INDUCTION]["average_temperature_k"] = value[INDUCTION]["result_average_temperature_k"] = [313.15, 303.15, 293.15]
    assert validate_public_v56_identity(value)["status"] == "needs_attention"


def test_v56_malformed_values_reject_without_raising() -> None:
    value = deepcopy(_records_v56())
    value[MODAL]["modal_rows"] = [{"mode": ["mode:1"]}]
    value[INDUCTION]["time_s"] = [[0.0]]
    assert validate_public_v56_identity(value)["status"] == "needs_attention"


def test_v56_numeric_digests_are_rejected() -> None:
    value = deepcopy(_records_v56())
    numeric_digest = int("1" * 64)
    for identity in (MODAL, INDUCTION):
        value[identity]["result_sha256"] = numeric_digest
        value[identity]["accepted_result_sha256"] = numeric_digest
    assert validate_public_v56_identity(value)["status"] == "needs_attention"
