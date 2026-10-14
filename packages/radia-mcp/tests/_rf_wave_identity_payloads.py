"""Payload builders for the RF/wave-port identity generalization gate tests.

Data only; not collected by pytest. Chain roots: nonlinear_inductance_sweep_gate
summaries (_summary_v23 .. _summary_v43) and independent _payload_vNN builders for
the standalone v44/v51-v55 identity validators."""
from __future__ import annotations

import cmath
import math
from copy import deepcopy
from radia_mcp.radia_ngsolve.network_artifact_identity_v51 import GROUP_DELAY, S_PARAMETER
from radia_mcp.radia_ngsolve.wave_energy_identity_v52 import ENERGY_BALANCE, EIGENMODE_Q
from radia_mcp.radia_ngsolve.wave_energy_identity_v55 import ANTENNA, RESONATOR
from radia_mcp.radia_ngsolve.wave_port_identity_v53 import FARFIELD, WAVEGUIDE
from radia_mcp.radia_ngsolve.wave_sar_identity_v54 import CUTOFF, SAR
from test_nonlinear_inductance_sweep_gate import _summary_v22


def _summary_v23():
    summary = _summary_v22()
    for index, row in enumerate(summary["runs"]):
        generation = f"broadband-sparam-{51 + index}"
        row["broadband_adaptive_mesh_sparam_renormalization_port_generation_identity"] = {
            "sweep_generation": generation,
            "adaptive_mesh_sweep_generation": generation,
            "frequency_interpolation_sweep_generation": generation,
            "port_mode_sweep_generation": generation,
            "renormalization_sweep_generation": generation,
            "sparameter_result_sweep_generation": generation,
            "adaptive_mesh_sha256": "1" * 64,
            "result_adaptive_mesh_sha256": "1" * 64,
            "frequency_samples_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_samples_hz": [1.0e9, 1.5e9, 2.0e9],
            "frequency_interpolation": "vector_fitting",
            "result_frequency_interpolation": "vector_fitting",
            "port_mode_ids": ["P1:M1", "P2:M1"],
            "result_port_mode_ids": ["P1:M1", "P2:M1"],
            "renormalization_impedance_ohm": [[50.0, 0.0], [50.0, 0.0]],
            "result_renormalization_impedance_ohm": [[50.0, 0.0], [50.0, 0.0]],
            "sparameter_table_sha256": "2" * 64,
            "result_sparameter_table_sha256": "2" * 64,
        }
        transient_generation = f"transient-monitor-{51 + index}"
        row["transient_monitor_time_origin_excitation_waveform_mesh_generation_identity"] = {
            "transient_generation": transient_generation,
            "time_origin_transient_generation": transient_generation,
            "excitation_waveform_transient_generation": transient_generation,
            "monitor_frame_transient_generation": transient_generation,
            "mesh_transient_generation": transient_generation,
            "field_result_transient_generation": transient_generation,
            "time_origin_s": 0.0,
            "result_time_origin_s": 0.0,
            "excitation_waveform_sha256": "3" * 64,
            "result_excitation_waveform_sha256": "3" * 64,
            "monitor_coordinate_frame": "global_xyz",
            "result_monitor_coordinate_frame": "global_xyz",
            "monitor_ids": [101, 102],
            "result_monitor_ids": [101, 102],
            "mesh_sha256": "4" * 64,
            "result_mesh_sha256": "4" * 64,
            "time_samples_s": [0.0, 1.0e-10, 2.0e-10],
            "result_time_samples_s": [0.0, 1.0e-10, 2.0e-10],
            "monitor_field_table_sha256": "5" * 64,
            "result_monitor_field_table_sha256": "5" * 64,
        }
    return summary


def _summary_v24():
    summary = _summary_v23()
    for index, row in enumerate(summary["runs"]):
        deembed_generation = f"deembed-{101 + index}"
        row[
            "deembedding_reference_plane_phase_causality_passivity_grid_generation_identity"
        ] = {
            "deembedding_generation": deembed_generation,
            "reference_plane_deembedding_generation": deembed_generation,
            "phase_deembedding_generation": deembed_generation,
            "causality_deembedding_generation": deembed_generation,
            "passivity_deembedding_generation": deembed_generation,
            "frequency_grid_deembedding_generation": deembed_generation,
            "result_deembedding_generation": deembed_generation,
            "port_mode_ids": ["P1:M1", "P2:M1"],
            "result_port_mode_ids": ["P1:M1", "P2:M1"],
            "reference_plane_offsets_m": [0.001, 0.001],
            "result_reference_plane_offsets_m": [0.001, 0.001],
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "unwrapped_phase_rad": [-0.2, -0.35, -0.5],
            "result_unwrapped_phase_rad": [-0.2, -0.35, -0.5],
            "causality_check_passed": True,
            "result_causality_check_passed": True,
            "passivity_max_singular_values": [0.82, 0.84, 0.86],
            "result_passivity_max_singular_values": [0.82, 0.84, 0.86],
            "deembedded_network_sha256": "1" * 64,
            "result_deembedded_network_sha256": "1" * 64,
        }
        cosim_generation = f"field-circuit-{101 + index}"
        row[
            "field_circuit_cosim_port_sign_impedance_power_balance_generation_identity"
        ] = {
            "cosim_generation": cosim_generation,
            "field_port_cosim_generation": cosim_generation,
            "circuit_port_cosim_generation": cosim_generation,
            "sign_cosim_generation": cosim_generation,
            "impedance_cosim_generation": cosim_generation,
            "power_balance_cosim_generation": cosim_generation,
            "result_cosim_generation": cosim_generation,
            "port_id": "P1:M1",
            "result_port_id": "P1:M1",
            "current_sign_convention": "positive_into_field_port",
            "result_current_sign_convention": "positive_into_field_port",
            "voltage_reference": "positive_to_negative_terminal",
            "result_voltage_reference": "positive_to_negative_terminal",
            "port_voltage_ri_v": [10.0, 2.0],
            "result_port_voltage_ri_v": [10.0, 2.0],
            "port_current_ri_a": [0.2, -0.04],
            "result_port_current_ri_a": [0.2, -0.04],
            "port_impedance_ri_ohm": [46.15384615384615, 19.23076923076923],
            "result_port_impedance_ri_ohm": [46.15384615384615, 19.23076923076923],
            "phasor_amplitude_convention": "rms",
            "result_phasor_amplitude_convention": "rms",
            "field_absorbed_power_w": 1.92,
            "circuit_delivered_power_w": 1.92,
            "power_balance_residual_w": 0.0,
            "result_power_balance_residual_w": 0.0,
            "cosim_result_sha256": "2" * 64,
            "reported_cosim_result_sha256": "2" * 64,
        }
    return summary


def _summary_v25():
    summary = _summary_v24()
    for index, row in enumerate(summary["runs"]):
        generation = f"adaptive-pass-{201 + index}"
        row[
            "adaptive_mesh_pass_sparameter_energy_convergence_grid_generation_identity"
        ] = {
            "adaptive_generation": generation,
            "mesh_pass_adaptive_generation": generation,
            "sparameter_adaptive_generation": generation,
            "energy_adaptive_generation": generation,
            "frequency_grid_adaptive_generation": generation,
            "stopping_rule_adaptive_generation": generation,
            "result_adaptive_generation": generation,
            "mesh_pass_ids": [0, 1, 2],
            "result_mesh_pass_ids": [0, 1, 2],
            "mesh_cell_counts": [10000, 18000, 29000],
            "result_mesh_cell_counts": [10000, 18000, 29000],
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "maximum_sparameter_delta": [0.1, 0.03, 0.005],
            "result_maximum_sparameter_delta": [0.1, 0.03, 0.005],
            "stored_energy_closure_residual": [0.05, 0.01, 0.001],
            "result_stored_energy_closure_residual": [0.05, 0.01, 0.001],
            "sparameter_delta_tolerance": 0.01,
            "energy_closure_tolerance": 0.005,
            "converged_pass_id": 2,
            "result_converged_pass_id": 2,
            "adaptive_result_sha256": "1" * 64,
            "reported_adaptive_result_sha256": "1" * 64,
        }
        generation = f"eigenmode-track-{201 + index}"
        row[
            "eigenmode_tracking_phase_normalization_port_coupling_mesh_generation_identity"
        ] = {
            "tracking_generation": generation,
            "modal_subspace_tracking_generation": generation,
            "phase_tracking_generation": generation,
            "normalization_tracking_generation": generation,
            "port_coupling_tracking_generation": generation,
            "mesh_tracking_generation": generation,
            "result_tracking_generation": generation,
            "sweep_parameters": [0.0, 0.5, 1.0],
            "result_sweep_parameters": [0.0, 0.5, 1.0],
            "tracked_mode_ids": ["mode-1", "mode-2"],
            "result_tracked_mode_ids": ["mode-1", "mode-2"],
            "modal_subspace_sha256": ["2" * 64, "3" * 64, "4" * 64],
            "result_modal_subspace_sha256": ["2" * 64, "3" * 64, "4" * 64],
            "phase_anchor_ids": ["probe-ez", "probe-hy"],
            "result_phase_anchor_ids": ["probe-ez", "probe-hy"],
            "normalization": "stored_energy_1j",
            "result_normalization": "stored_energy_1j",
            "port_coupling_magnitudes": [[0.8, 0.1], [0.78, 0.12], [0.75, 0.15]],
            "result_port_coupling_magnitudes": [[0.8, 0.1], [0.78, 0.12], [0.75, 0.15]],
            "mesh_sha256": ["5" * 64, "6" * 64, "7" * 64],
            "result_mesh_sha256": ["5" * 64, "6" * 64, "7" * 64],
            "eigenmode_track_sha256": "8" * 64,
            "reported_eigenmode_track_sha256": "8" * 64,
        }
    return summary


def _summary_v26():
    summary = _summary_v25()
    for index, row in enumerate(summary["runs"]):
        generation = f"port-network-{301 + index}"
        row[
            "port_deembedding_reference_plane_impedance_mode_normalization_smatrix_generation_identity"
        ] = {
            "network_generation": generation,
            "port_mode_network_generation": generation,
            "deembedding_network_generation": generation,
            "reference_impedance_network_generation": generation,
            "normalization_network_generation": generation,
            "frequency_grid_network_generation": generation,
            "result_network_generation": generation,
            "port_mode_ids": ["P1:M1", "P2:M1"],
            "result_port_mode_ids": ["P1:M1", "P2:M1"],
            "reference_plane_offsets_m": [0.001, 0.0015],
            "result_reference_plane_offsets_m": [0.001, 0.0015],
            "reference_impedance_ri_ohm": [[50.0, 0.0], [50.0, 0.0]],
            "result_reference_impedance_ri_ohm": [[50.0, 0.0], [50.0, 0.0]],
            "wave_normalization": "power_wave",
            "result_wave_normalization": "power_wave",
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "smatrix_ri": [[[0.1, -0.01], [0.8, -0.1]], [[0.79, -0.12], [0.11, -0.02]]],
            "result_smatrix_ri": [[[0.1, -0.01], [0.8, -0.1]], [[0.79, -0.12], [0.11, -0.02]]],
            "smatrix_sha256": "1" * 64,
            "reported_smatrix_sha256": "1" * 64,
        }
        generation = f"farfield-{301 + index}"
        row[
            "farfield_angular_grid_polarization_coordinate_power_normalization_mesh_generation_identity"
        ] = {
            "farfield_generation": generation,
            "angular_grid_farfield_generation": generation,
            "polarization_farfield_generation": generation,
            "coordinate_farfield_generation": generation,
            "power_farfield_generation": generation,
            "mesh_farfield_generation": generation,
            "result_farfield_generation": generation,
            "theta_deg": [0.0, 45.0, 90.0],
            "result_theta_deg": [0.0, 45.0, 90.0],
            "phi_deg": [0.0, 90.0],
            "result_phi_deg": [0.0, 90.0],
            "polarization_basis": "ludwig3_co_cross",
            "result_polarization_basis": "ludwig3_co_cross",
            "coordinate_frame": "global_xyz_z_up",
            "result_coordinate_frame": "global_xyz_z_up",
            "radiated_power_w": 0.8,
            "result_radiated_power_w": 0.8,
            "field_normalization": "sqrt_radiated_power",
            "result_field_normalization": "sqrt_radiated_power",
            "mesh_sha256": "2" * 64,
            "result_mesh_sha256": "2" * 64,
            "farfield_sha256": "3" * 64,
            "reported_farfield_sha256": "3" * 64,
        }
    return summary


def _summary_v27():
    summary = _summary_v26()
    for index, row in enumerate(summary["runs"]):
        generation = f"td-port-{311 + index}"
        row[
            "time_domain_port_waveform_normalization_fft_window_grid_deembedding_smatrix_generation_identity"
        ] = {
            "time_domain_generation": generation,
            "waveform_time_domain_generation": generation,
            "normalization_time_domain_generation": generation,
            "fft_time_domain_generation": generation,
            "grid_time_domain_generation": generation,
            "deembedding_time_domain_generation": generation,
            "smatrix_time_domain_generation": generation,
            "result_time_domain_generation": generation,
            "port_mode_ids": ["P1:M1", "P2:M1"],
            "result_port_mode_ids": ["P1:M1", "P2:M1"],
            "time_grid_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12, 4.0e-12, 5.0e-12],
            "result_time_grid_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12, 4.0e-12, 5.0e-12],
            "incident_waveform": [0.0, 0.2, 1.0, 0.2, 0.0, 0.0],
            "result_incident_waveform": [0.0, 0.2, 1.0, 0.2, 0.0, 0.0],
            "wave_normalization": "power-wave",
            "result_wave_normalization": "power-wave",
            "reference_impedance_ohm": [50.0, 50.0],
            "result_reference_impedance_ohm": [50.0, 50.0],
            "fft_window": "tukey-0.2",
            "result_fft_window": "tukey-0.2",
            "frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "result_frequency_grid_hz": [1.0e9, 1.5e9, 2.0e9],
            "deembedding_offsets_m": [0.001, 0.0015],
            "result_deembedding_offsets_m": [0.001, 0.0015],
            "smatrix_ri": [
                [[0.1, -0.01], [0.8, -0.1]],
                [[0.79, -0.12], [0.11, -0.02]],
            ],
            "result_smatrix_ri": [
                [[0.1, -0.01], [0.8, -0.1]],
                [[0.79, -0.12], [0.11, -0.02]],
            ],
            "time_result_sha256": "1" * 64,
            "accepted_time_result_sha256": "1" * 64,
            "smatrix_sha256": "2" * 64,
            "accepted_smatrix_sha256": "2" * 64,
        }
        generation = f"huygens-{311 + index}"
        row[
            "huygens_box_orientation_phase_center_frequency_mesh_near_far_transform_generation_identity"
        ] = {
            "huygens_generation": generation,
            "orientation_huygens_generation": generation,
            "phase_center_huygens_generation": generation,
            "frequency_huygens_generation": generation,
            "mesh_huygens_generation": generation,
            "transform_huygens_generation": generation,
            "result_huygens_generation": generation,
            "box_face_ids": ["-x", "+x", "-y", "+y", "-z", "+z"],
            "result_box_face_ids": ["-x", "+x", "-y", "+y", "-z", "+z"],
            "outward_orientation_sign": [-1, 1, -1, 1, -1, 1],
            "result_outward_orientation_sign": [-1, 1, -1, 1, -1, 1],
            "phase_center_m": [0.0, 0.0, 0.0],
            "result_phase_center_m": [0.0, 0.0, 0.0],
            "frequency_hz": 10.0e9,
            "result_frequency_hz": 10.0e9,
            "near_far_transform": "equivalent-current-near-to-far",
            "result_near_far_transform": "equivalent-current-near-to-far",
            "encloses_all_sources": True,
            "result_encloses_all_sources": True,
            "enclosing_mesh_sha256": "3" * 64,
            "result_enclosing_mesh_sha256": "3" * 64,
            "near_field_sha256": "4" * 64,
            "accepted_near_field_sha256": "4" * 64,
            "far_field_sha256": "5" * 64,
            "accepted_far_field_sha256": "5" * 64,
        }
    return summary


def _summary_v28():
    summary = _summary_v27()
    for index, row in enumerate(summary["runs"]):
        generation = f"waveguide-port-{321 + index}"
        row[
            "waveguide_port_mode_cutoff_normalization_reference_plane_mesh_field_result_generation_identity"
        ] = {
            "port_generation": generation,
            "mode_port_generation": generation,
            "cutoff_port_generation": generation,
            "normalization_port_generation": generation,
            "reference_plane_port_generation": generation,
            "mesh_port_generation": generation,
            "field_port_generation": generation,
            "result_port_generation": generation,
            "port_id": "P1",
            "result_port_id": "P1",
            "mode_id": "TE10",
            "result_mode_id": "TE10",
            "cutoff_frequency_hz": 6.557e9,
            "result_cutoff_frequency_hz": 6.557e9,
            "evaluation_frequency_hz": 10.0e9,
            "result_evaluation_frequency_hz": 10.0e9,
            "normalization": "unit-power-wave",
            "result_normalization": "unit-power-wave",
            "reference_plane_m": 0.005,
            "result_reference_plane_m": 0.005,
            "port_mesh_sha256": "1" * 64,
            "result_port_mesh_sha256": "1" * 64,
            "field_eigenvector_sha256": "2" * 64,
            "result_field_eigenvector_sha256": "2" * 64,
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        }
        generation = f"wake-impedance-{321 + index}"
        row[
            "wake_impedance_bunch_profile_time_grid_frequency_transform_normalization_mesh_result_generation_identity"
        ] = {
            "wake_generation": generation,
            "bunch_wake_generation": generation,
            "time_wake_generation": generation,
            "transform_wake_generation": generation,
            "frequency_wake_generation": generation,
            "normalization_wake_generation": generation,
            "mesh_wake_generation": generation,
            "result_wake_generation": generation,
            "bunch_profile": "gaussian",
            "result_bunch_profile": "gaussian",
            "bunch_sigma_s": 1.0e-12,
            "result_bunch_sigma_s": 1.0e-12,
            "bunch_charge_c": 1.0e-9,
            "result_bunch_charge_c": 1.0e-9,
            "time_grid_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12],
            "result_time_grid_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12],
            "wake_potential_v_c": [0.0, 2.0e12, 1.0e12, 0.0],
            "result_wake_potential_v_c": [0.0, 2.0e12, 1.0e12, 0.0],
            "fft_convention": "exp-minus-i-omega-t",
            "result_fft_convention": "exp-minus-i-omega-t",
            "frequency_grid_hz": [0.0, 1.0e9, 2.0e9],
            "result_frequency_grid_hz": [0.0, 1.0e9, 2.0e9],
            "impedance_normalization": "longitudinal-v-per-coulomb",
            "result_impedance_normalization": "longitudinal-v-per-coulomb",
            "mesh_sha256": "4" * 64,
            "result_mesh_sha256": "4" * 64,
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
        }
    return summary


def _summary_v29():
    summary = _summary_v28()
    for index, row in enumerate(summary["runs"]):
        generation = f"dispersive-fit-{341 + index}"
        row["dispersive_vector_fit_passivity_causality_temperature_generation_identity"] = {
            "fit_generation": generation, "pole_fit_generation": generation,
            "residue_fit_generation": generation, "passivity_fit_generation": generation,
            "causality_fit_generation": generation, "temperature_fit_generation": generation,
            "frequency_fit_generation": generation, "material_fit_generation": generation,
            "result_fit_generation": generation,
            "temperature_c": 25.0, "result_temperature_c": 25.0,
            "frequency_grid_hz": [1.0e9, 2.0e9, 4.0e9, 8.0e9],
            "result_frequency_grid_hz": [1.0e9, 2.0e9, 4.0e9, 8.0e9],
            "poles_rad_s": [[-1.0e9, 1.0e10], [-2.0e9, 2.0e10]],
            "result_poles_rad_s": [[-1.0e9, 1.0e10], [-2.0e9, 2.0e10]],
            "residues": [[1.0e9, 2.0e8], [5.0e8, 1.0e8]],
            "result_residues": [[1.0e9, 2.0e8], [5.0e8, 1.0e8]],
            "passivity_enforced": True, "result_passivity_enforced": True,
            "minimum_dissipation": 0.01, "result_minimum_dissipation": 0.01,
            "causality_residual": 1.0e-8, "result_causality_residual": 1.0e-8,
            "causality_residual_limit": 1.0e-6,
            "material_table_sha256": "1" * 64, "result_material_table_sha256": "1" * 64,
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        }
        generation = f"array-scan-{341 + index}"
        row["array_embedded_pattern_feed_phase_active_reflection_scan_generation_identity"] = {
            "array_generation": generation, "pattern_array_generation": generation,
            "element_array_generation": generation, "phase_array_generation": generation,
            "reflection_array_generation": generation, "scan_array_generation": generation,
            "power_array_generation": generation, "mesh_array_generation": generation,
            "result_array_generation": generation,
            "element_order": [1, 2, 3, 4], "result_element_order": [1, 2, 3, 4],
            "embedded_pattern_sha256": ["3" * 64, "4" * 64, "5" * 64, "6" * 64],
            "result_embedded_pattern_sha256": ["3" * 64, "4" * 64, "5" * 64, "6" * 64],
            "scan_angles_deg": [-30.0, 0.0, 30.0], "result_scan_angles_deg": [-30.0, 0.0, 30.0],
            "feed_phase_deg": [[0.0, -45.0, -90.0, -135.0], [0.0, 0.0, 0.0, 0.0], [0.0, 45.0, 90.0, 135.0]],
            "result_feed_phase_deg": [[0.0, -45.0, -90.0, -135.0], [0.0, 0.0, 0.0, 0.0], [0.0, 45.0, 90.0, 135.0]],
            "active_reflection_magnitude": [0.2, 0.1, 0.2],
            "result_active_reflection_magnitude": [0.2, 0.1, 0.2],
            "accepted_power_fraction": [0.96, 0.99, 0.96],
            "result_accepted_power_fraction": [0.96, 0.99, 0.96],
            "array_mesh_sha256": "7" * 64, "result_array_mesh_sha256": "7" * 64,
            "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
        }
    return summary


def _summary_v30():
    summary = _summary_v29()
    for index, row in enumerate(summary["runs"]):
        generation = f"waveguide-port-{351 + index}"
        row["waveguide_port_mode_power_deembed_impedance_frequency_port_smatrix_result_identity"] = {
            "port_generation": generation, "mode_port_generation": generation, "power_port_generation": generation,
            "deembed_port_generation": generation, "impedance_port_generation": generation,
            "frequency_port_generation": generation, "order_port_generation": generation, "result_port_generation": generation,
            "mode_ids": ["port1:TE10", "port2:TE10"], "result_mode_ids": ["port1:TE10", "port2:TE10"],
            "power_normalization_w": [1.0, 1.0], "result_power_normalization_w": [1.0, 1.0],
            "deembed_plane_m": [0.0, 0.1], "result_deembed_plane_m": [0.0, 0.1],
            "reference_impedance_ohm": [50.0, 50.0], "result_reference_impedance_ohm": [50.0, 50.0],
            "frequency_hz": [8.0e9, 9.0e9, 10.0e9], "result_frequency_hz": [8.0e9, 9.0e9, 10.0e9],
            "port_order": [1, 2], "result_port_order": [1, 2],
            "smatrix_ri": [[[0.1, 0.0], [0.9, 0.0]], [[0.9, 0.0], [0.1, 0.0]]],
            "result_smatrix_ri": [[[0.1, 0.0], [0.9, 0.0]], [[0.9, 0.0], [0.1, 0.0]]],
            "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        }
        generation = f"sar-mass-{351 + index}"
        row["sar_mass_density_voxel_frequency_field_mesh_result_identity"] = {
            "sar_generation": generation, "mass_sar_generation": generation, "density_sar_generation": generation,
            "voxel_sar_generation": generation, "frequency_sar_generation": generation,
            "field_sar_generation": generation, "mesh_sar_generation": generation, "result_sar_generation": generation,
            "averaging_mass_kg": 0.01, "result_averaging_mass_kg": 0.01,
            "tissue_density_kg_m3": 1000.0, "result_tissue_density_kg_m3": 1000.0,
            "voxel_ids": [101, 102, 103], "result_voxel_ids": [101, 102, 103],
            "voxel_mass_kg": [0.003, 0.004, 0.003], "result_voxel_mass_kg": [0.003, 0.004, 0.003],
            "frequency_hz": 2.45e9, "result_frequency_hz": 2.45e9,
            "field_normalization": "accepted_power_1w", "result_field_normalization": "accepted_power_1w",
            "sar_w_kg": 1.2, "result_sar_w_kg": 1.2,
            "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        }
    return summary


def _summary_v31():
    summary = _summary_v30()
    for index, row in enumerate(summary["runs"]):
        generation = f"wave-port-reference-{361 + index}"
        row["wave_port_modal_power_impedance_deembed_phase_balance_result_identity"] = {
            "port_generation": generation,
            "mode_port_generation": generation,
            "power_port_generation": generation,
            "impedance_port_generation": generation,
            "deembed_port_generation": generation,
            "phase_port_generation": generation,
            "owner_port_generation": generation,
            "balance_port_generation": generation,
            "result_port_generation": generation,
            "mode_ids": ["port1:TE10", "port2:TE10"],
            "result_mode_ids": ["port1:TE10", "port2:TE10"],
            "modal_power_normalization_w": [1.0, 1.0],
            "result_modal_power_normalization_w": [1.0, 1.0],
            "reference_impedance_ohm": [50.0, 50.0],
            "result_reference_impedance_ohm": [50.0, 50.0],
            "deembed_plane_m": [0.0, 0.1],
            "result_deembed_plane_m": [0.0, 0.1],
            "phase_reference_rad": [0.0, 0.0],
            "result_phase_reference_rad": [0.0, 0.0],
            "port_mode_owner_ids": ["project-361:port-1:TE10", "project-361:port-2:TE10"],
            "result_port_mode_owner_ids": ["project-361:port-1:TE10", "project-361:port-2:TE10"],
            "incident_power_w": 1.0,
            "reflected_power_w": 0.04,
            "transmitted_power_w": 0.94,
            "dissipated_power_w": 0.02,
            "result_power_balance_w": 1.0,
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        }
        generation = f"farfield-basis-{361 + index}"
        row["farfield_spherical_basis_handedness_polarization_phase_power_result_identity"] = {
            "farfield_generation": generation,
            "basis_farfield_generation": generation,
            "handedness_farfield_generation": generation,
            "order_farfield_generation": generation,
            "polarization_farfield_generation": generation,
            "phase_farfield_generation": generation,
            "weights_farfield_generation": generation,
            "power_farfield_generation": generation,
            "owner_farfield_generation": generation,
            "result_farfield_generation": generation,
            "spherical_basis": "e_theta_e_phi",
            "result_spherical_basis": "e_theta_e_phi",
            "coordinate_handedness": "right_handed",
            "result_coordinate_handedness": "right_handed",
            "angular_order": "theta_major_phi_minor",
            "result_angular_order": "theta_major_phi_minor",
            "polarization_phase_convention": "exp_plus_j_phase",
            "result_polarization_phase_convention": "exp_plus_j_phase",
            "theta_weights": [0.25, 0.5, 0.25],
            "result_theta_weights": [0.25, 0.5, 0.25],
            "phi_weights": [0.5, 0.5],
            "result_phi_weights": [0.5, 0.5],
            "radiated_power_w": 0.94,
            "integrated_radiated_power_w": 0.94,
            "farfield_owner_id": "project-361:monitor-ff-1",
            "accepted_farfield_owner_id": "project-361:monitor-ff-1",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        }
    return summary


def _summary_v32():
    summary = _summary_v31()
    for index, row in enumerate(summary["runs"]):
        generation = f"dispersive-port-{371 + index}"
        row[
            "dispersive_port_mode_branch_cutoff_normalization_beta_phase_group_delay_mesh_result_identity"
        ] = {
            "port_generation": generation,
            **{
                key: generation
                for key in (
                    "mode_port_generation",
                    "branch_port_generation",
                    "cutoff_port_generation",
                    "normalization_port_generation",
                    "beta_port_generation",
                    "phase_port_generation",
                    "delay_port_generation",
                    "mesh_port_generation",
                    "result_port_generation",
                )
            },
            "mode_id": "port1:TE10",
            "result_mode_id": "port1:TE10",
            "tracked_branch_id": "branch:TE10:forward",
            "result_tracked_branch_id": "branch:TE10:forward",
            "cutoff_frequency_hz": 6.56e9,
            "result_cutoff_frequency_hz": 6.56e9,
            "frequency_hz": [8.0e9, 9.0e9, 10.0e9],
            "result_frequency_hz": [8.0e9, 9.0e9, 10.0e9],
            "modal_normalization": "unit_forward_power",
            "result_modal_normalization": "unit_forward_power",
            "propagation_constant_sign": "positive_forward",
            "result_propagation_constant_sign": "positive_forward",
            "propagation_constant_rad_per_m": [101.0, 132.0, 158.0],
            "result_propagation_constant_rad_per_m": [101.0, 132.0, 158.0],
            "deembedded_phase_rad": [0.4, 0.2, 0.0],
            "result_deembedded_phase_rad": [0.4, 0.2, 0.0],
            "group_delay_s": 3.183098861837907e-11,
            "result_group_delay_s": 3.183098861837907e-11,
            "mesh_sha256": "1" * 64,
            "result_mesh_sha256": "1" * 64,
            "result_owner": "waveguide/case-371/port1-te10",
            "accepted_result_owner": "waveguide/case-371/port1-te10",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        }
        generation = f"transient-farfield-{371 + index}"
        row[
            "transient_farfield_time_gate_fft_window_phase_center_angular_energy_monitor_result_identity"
        ] = {
            "farfield_generation": generation,
            **{
                key: generation
                for key in (
                    "gate_farfield_generation",
                    "fft_farfield_generation",
                    "phase_center_farfield_generation",
                    "angular_farfield_generation",
                    "energy_farfield_generation",
                    "monitor_farfield_generation",
                    "owner_farfield_generation",
                    "result_farfield_generation",
                )
            },
            "time_gate_s": [2.0e-9, 8.0e-9],
            "result_time_gate_s": [2.0e-9, 8.0e-9],
            "fft_window": "hann",
            "result_fft_window": "hann",
            "fft_normalization": "one_sided_energy_preserving",
            "result_fft_normalization": "one_sided_energy_preserving",
            "phase_center_m": [0.0, 0.0, 0.0],
            "result_phase_center_m": [0.0, 0.0, 0.0],
            "theta_deg": [0.0, 45.0, 90.0],
            "result_theta_deg": [0.0, 45.0, 90.0],
            "phi_deg": [0.0, 90.0, 180.0, 270.0],
            "result_phi_deg": [0.0, 90.0, 180.0, 270.0],
            "accepted_energy_j": 1.0,
            "result_accepted_energy_j": 1.0,
            "radiated_energy_j": 0.82,
            "result_radiated_energy_j": 0.82,
            "monitor_owner": "project-371:monitor-transient-ff",
            "result_monitor_owner": "project-371:monitor-transient-ff",
            "result_owner": "farfield/case-371/time-domain",
            "accepted_result_owner": "farfield/case-371/time-domain",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        }
    return summary


def _summary_v33():
    summary = _summary_v32()
    for index, row in enumerate(summary["runs"]):
        generation = f"eigenmode-q-{381 + index}"
        q_conductor = 10000.0
        q_dielectric = 20000.0
        q_radiation = 50000.0
        q_total = 1.0 / (
            1.0 / q_conductor + 1.0 / q_dielectric + 1.0 / q_radiation
        )
        row[
            "eigenmode_frequency_branch_energy_conductor_dielectric_radiation_q_mesh_owner_result_identity"
        ] = {
            "eigenmode_generation": generation,
            **{
                key: generation
                for key in (
                    "frequency_generation",
                    "branch_generation",
                    "energy_generation",
                    "conductor_q_generation",
                    "dielectric_q_generation",
                    "radiation_q_generation",
                    "inverse_sum_generation",
                    "mesh_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "mode_id": "mode-1",
            "result_mode_id": "mode-1",
            "mode_branch": "fundamental",
            "result_mode_branch": "fundamental",
            "frequency_hz": 10.0e9,
            "result_frequency_hz": 10.0e9,
            "electric_energy_j": 0.5,
            "result_electric_energy_j": 0.5,
            "magnetic_energy_j": 0.5,
            "result_magnetic_energy_j": 0.5,
            "stored_energy_j": 1.0,
            "result_stored_energy_j": 1.0,
            "q_conductor": q_conductor,
            "result_q_conductor": q_conductor,
            "q_dielectric": q_dielectric,
            "result_q_dielectric": q_dielectric,
            "q_radiation": q_radiation,
            "result_q_radiation": q_radiation,
            "q_total": q_total,
            "result_q_total": q_total,
            "mesh_sha256": "1" * 64,
            "result_mesh_sha256": "1" * 64,
            "mode_owner": "cavity/case-381/mode-1",
            "accepted_mode_owner": "cavity/case-381/mode-1",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        }
        generation = f"tdr-reference-{381 + index}"
        times = [0.0, 1.0e-9, 2.0e-9, 3.0e-9, 4.0e-9]
        waveform = [0.0, 0.0, 0.0, 0.1, 0.05]
        row[
            "tdr_reference_plane_velocity_time_zero_impedance_arrival_window_causality_energy_owner_result_identity"
        ] = {
            "tdr_generation": generation,
            **{
                key: generation
                for key in (
                    "reference_generation",
                    "velocity_generation",
                    "time_zero_generation",
                    "impedance_generation",
                    "arrival_generation",
                    "window_generation",
                    "causality_generation",
                    "energy_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "reference_plane_m": 0.0,
            "result_reference_plane_m": 0.0,
            "propagation_velocity_m_per_s": 2.0e8,
            "result_propagation_velocity_m_per_s": 2.0e8,
            "time_zero_s": 1.0e-9,
            "result_time_zero_s": 1.0e-9,
            "characteristic_impedance_ohm": 50.0,
            "result_characteristic_impedance_ohm": 50.0,
            "reflection_distance_m": 0.2,
            "result_reflection_distance_m": 0.2,
            "reflection_arrival_s": 3.0e-9,
            "result_reflection_arrival_s": 3.0e-9,
            "time_window_s": [0.0, 4.0e-9],
            "result_time_window_s": [0.0, 4.0e-9],
            "time_samples_s": times,
            "result_time_samples_s": list(times),
            "reflection_waveform": waveform,
            "result_reflection_waveform": list(waveform),
            "pre_arrival_max_abs": 0.0,
            "result_pre_arrival_max_abs": 0.0,
            "incident_energy_j": 1.0,
            "result_incident_energy_j": 1.0,
            "reflected_energy_j": 0.1,
            "result_reflected_energy_j": 0.1,
            "accepted_energy_j": 0.9,
            "result_accepted_energy_j": 0.9,
            "waveform_owner": "tdr/case-381/port-1",
            "accepted_waveform_owner": "tdr/case-381/port-1",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        }
    return summary


def _summary_v34():
    summary = _summary_v33()
    for index, row in enumerate(summary["runs"]):
        generation = f"sparameter-gated-{391 + index}"
        row[
            "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
        ] = {
            "sparameter_generation": generation,
            **{
                key: generation
                for key in (
                    "reference_generation",
                    "gate_generation",
                    "causality_generation",
                    "passivity_generation",
                    "energy_generation",
                    "port_generation",
                    "frequency_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "reference_plane_shift_m": 0.01,
            "result_reference_plane_shift_m": 0.01,
            "time_gate_window_s": [0.0, 4.0e-9],
            "result_time_gate_window_s": [0.0, 4.0e-9],
            "impulse_time_s": [-1.0e-9, 0.0, 1.0e-9, 2.0e-9, 3.0e-9],
            "result_impulse_time_s": [-1.0e-9, 0.0, 1.0e-9, 2.0e-9, 3.0e-9],
            "impulse_response": [0.0, 0.2, 0.1, 0.05, 0.0],
            "result_impulse_response": [0.0, 0.2, 0.1, 0.05, 0.0],
            "pre_zero_max_abs": 0.0,
            "result_pre_zero_max_abs": 0.0,
            "maximum_singular_values": [0.8, 0.9, 0.85],
            "result_maximum_singular_values": [0.8, 0.9, 0.85],
            "incident_energy_j": 1.0,
            "result_incident_energy_j": 1.0,
            "reflected_energy_j": 0.1,
            "result_reflected_energy_j": 0.1,
            "transmitted_energy_j": 0.8,
            "result_transmitted_energy_j": 0.8,
            "absorbed_energy_j": 0.1,
            "result_absorbed_energy_j": 0.1,
            "port_impedance_ohm": [50.0, 50.0],
            "result_port_impedance_ohm": [50.0, 50.0],
            "frequency_grid_hz": [1.0e9, 2.0e9, 3.0e9],
            "result_frequency_grid_hz": [1.0e9, 2.0e9, 3.0e9],
            "sparameter_owner": "network/case-391/gated",
            "accepted_sparameter_owner": "network/case-391/gated",
            "sparameter_sha256": "1" * 64,
            "accepted_sparameter_sha256": "1" * 64,
        }

        generation = f"degenerate-eigenmode-{391 + index}"
        row[
            "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
        ] = {
            "degenerate_mode_generation": generation,
            **{
                key: generation
                for key in (
                    "subspace_generation",
                    "principal_angle_generation",
                    "mass_generation",
                    "phase_generation",
                    "tracking_generation",
                    "residual_generation",
                    "mesh_generation",
                    "owner_generation",
                    "result_generation",
                )
            },
            "mode_frequencies_hz": [10.0e9, 10.0e9 + 1000.0],
            "result_mode_frequencies_hz": [10.0e9, 10.0e9 + 1000.0],
            "principal_angles_rad": [0.0, 0.001],
            "result_principal_angles_rad": [0.0, 0.001],
            "mass_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "result_mass_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "mass_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "result_mass_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "phase_anchor_complex": [[1.0, 0.0], [1.0, 0.0]],
            "result_phase_anchor_complex": [[1.0, 0.0], [1.0, 0.0]],
            "tracking_subspace_ids": ["subspace-391", "subspace-391"],
            "result_tracking_subspace_ids": ["subspace-391", "subspace-391"],
            "residual_norms": [1.0e-9, 2.0e-9],
            "result_residual_norms": [1.0e-9, 2.0e-9],
            "mesh_cell_counts": [1000, 8000, 64000],
            "result_mesh_cell_counts": [1000, 8000, 64000],
            "mesh_converged_frequency_hz": [9.8e9, 9.98e9, 10.0e9],
            "result_mesh_converged_frequency_hz": [9.8e9, 9.98e9, 10.0e9],
            "eigenmode_mesh_sha256": "2" * 64,
            "result_eigenmode_mesh_sha256": "2" * 64,
            "field_owner": "eigenmode/case-391/subspace",
            "accepted_field_owner": "eigenmode/case-391/subspace",
            "field_sha256": "3" * 64,
            "accepted_field_sha256": "3" * 64,
        }
    return summary


def _summary_v35():
    summary = _summary_v34()
    c0 = 299_792_458.0
    cutoff = 6.5e9
    frequencies = [8.0e9, 10.0e9, 12.0e9]
    impedance = [377.0 / math.sqrt(1.0 - (cutoff / frequency) ** 2) for frequency in frequencies]
    beta = [2.0 * math.pi / c0 * math.sqrt(frequency**2 - cutoff**2) for frequency in frequencies]
    for index, row in enumerate(summary["runs"]):
        generation = f"waveguide-port-{411 + index}"
        row[
            "waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_frequency_owner_result_identity"
        ] = {
            "waveguide_port_generation": generation,
            **{
                key: generation
                for key in (
                    "mode_generation", "power_generation", "orthogonality_generation",
                    "impedance_generation", "deembed_generation", "cutoff_generation",
                    "frequency_generation", "mesh_generation", "owner_generation",
                    "result_generation",
                )
            },
            "mode_name": "TE10", "result_mode_name": "TE10",
            "normalization": "accepted_power_1w", "result_normalization": "accepted_power_1w",
            "modal_power_w": [1.0, 1.0, 1.0], "result_modal_power_w": [1.0, 1.0, 1.0],
            "mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "result_mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "result_mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "impedance_definition": "te_wave_impedance",
            "result_impedance_definition": "te_wave_impedance",
            "modal_impedance_ohm": impedance, "result_modal_impedance_ohm": list(impedance),
            "frequency_grid_hz": frequencies, "result_frequency_grid_hz": list(frequencies),
            "cutoff_frequency_hz": cutoff, "result_cutoff_frequency_hz": cutoff,
            "propagation_constant_rad_m": beta, "result_propagation_constant_rad_m": list(beta),
            "reference_plane_m": 0.0, "result_reference_plane_m": 0.0,
            "deembedded_reference_plane_m": 0.01, "result_deembedded_reference_plane_m": 0.01,
            "deembed_phase_rad": [-item * 0.01 for item in beta],
            "result_deembed_phase_rad": [-item * 0.01 for item in beta],
            "port_mesh_sha256": "1" * 64, "result_port_mesh_sha256": "1" * 64,
            "port_owner": "waveguide-port/case-411",
            "accepted_port_owner": "waveguide-port/case-411",
            "port_result_sha256": "2" * 64, "accepted_port_result_sha256": "2" * 64,
        }
        generation = f"nearfar-{411 + index}"
        weights = [math.pi / 2.0] * 8
        intensity = [8.0 / math.pi, 8.0 / math.pi] + [0.0] * 6
        row[
            "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
        ] = {
            "nearfar_generation": generation,
            **{
                key: generation
                for key in (
                    "sphere_generation", "power_generation", "directivity_generation",
                    "gain_generation", "efficiency_generation", "polarization_generation",
                    "quadrature_generation", "mesh_generation", "owner_generation",
                    "result_generation",
                )
            },
            "frequency_hz": 10.0e9, "result_frequency_hz": 10.0e9,
            "accepted_power_w": 10.0, "result_accepted_power_w": 10.0,
            "enclosing_sphere_power_w": 8.0, "result_enclosing_sphere_power_w": 8.0,
            "radiated_power_w": 8.0, "result_radiated_power_w": 8.0,
            "radiation_efficiency": 0.8, "result_radiation_efficiency": 0.8,
            "maximum_directivity_linear": 4.0, "result_maximum_directivity_linear": 4.0,
            "realized_gain_linear": 3.2, "result_realized_gain_linear": 3.2,
            "polarization_basis": "theta_phi_right_handed",
            "result_polarization_basis": "theta_phi_right_handed",
            "copolar_definition": "ludwig3", "result_copolar_definition": "ludwig3",
            "angular_quadrature_weights_sr": weights,
            "result_angular_quadrature_weights_sr": list(weights),
            "radiation_intensity_w_sr": intensity,
            "result_radiation_intensity_w_sr": list(intensity),
            "farfield_mesh_sha256": "3" * 64, "result_farfield_mesh_sha256": "3" * 64,
            "nearfar_owner": "nearfar/case-411", "accepted_nearfar_owner": "nearfar/case-411",
            "nearfar_result_sha256": "4" * 64, "accepted_nearfar_result_sha256": "4" * 64,
        }
    return summary


def _summary_v36():
    summary = _summary_v35()
    for index, row in enumerate(summary["runs"]):
        generation = f"cavity-mode-{412 + index}"
        q_external, q_dielectric, q_conductor = 20000.0, 50000.0, 40000.0
        q_total = 1.0 / (1.0 / q_external + 1.0 / q_dielectric + 1.0 / q_conductor)
        row["cavity_eigenmode_frequency_q_energy_orthogonality_degeneracy_mesh_owner_result_identity"] = {
            "cavity_generation": generation,
            **{key: generation for key in (
                "frequency_generation", "q_generation", "energy_generation",
                "orthogonality_generation", "degeneracy_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )},
            "mode_ids": [1, 2], "result_mode_ids": [1, 2],
            "mode_frequency_hz": [10.0e9, 10.5e9], "result_mode_frequency_hz": [10.0e9, 10.5e9],
            "q_external": q_external, "result_q_external": q_external,
            "q_dielectric": q_dielectric, "result_q_dielectric": q_dielectric,
            "q_conductor": q_conductor, "result_q_conductor": q_conductor,
            "q_total": q_total, "result_q_total": q_total,
            "electric_energy_j": [0.5, 0.5], "result_electric_energy_j": [0.5, 0.5],
            "magnetic_energy_j": [0.5, 0.5], "result_magnetic_energy_j": [0.5, 0.5],
            "normalization": "unit_total_energy_1j", "result_normalization": "unit_total_energy_1j",
            "mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "result_mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "result_mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "degeneracy_order": [1, 2], "result_degeneracy_order": [1, 2],
            "mesh_dof": [10000, 20000, 40000], "result_mesh_dof": [10000, 20000, 40000],
            "mesh_frequency_hz": [9.9e9, 9.98e9, 10.0e9],
            "result_mesh_frequency_hz": [9.9e9, 9.98e9, 10.0e9],
            "cavity_mesh_sha256": "1" * 64, "result_cavity_mesh_sha256": "1" * 64,
            "cavity_owner": "cavity/case-412", "accepted_cavity_owner": "cavity/case-412",
            "cavity_result_sha256": "2" * 64, "accepted_cavity_result_sha256": "2" * 64,
        }
        generation = f"active-array-{412 + index}"
        excitation = [1.0 + 0.0j, cmath.exp(1j * math.pi / 2.0)]
        matrix = [[0.1 + 0.0j, 0.2 + 0.0j], [0.2 + 0.0j, 0.1 + 0.0j]]
        active = [sum(matrix[r][c] * excitation[c] for c in range(2)) / excitation[r] for r in range(2)]
        pair = lambda value: [value.real, value.imag]
        row["active_sparameter_embedded_pattern_scan_impedance_power_frequency_mesh_owner_result_identity"] = {
            "array_generation": generation,
            **{key: generation for key in (
                "sparameter_generation", "pattern_generation", "scan_generation",
                "impedance_generation", "power_generation", "frequency_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )},
            "frequency_hz": 10.0e9, "result_frequency_hz": 10.0e9,
            "port_impedance_ohm": [50.0, 50.0], "result_port_impedance_ohm": [50.0, 50.0],
            "scan_phase_rad": [0.0, math.pi / 2.0], "result_scan_phase_rad": [0.0, math.pi / 2.0],
            "excitation_complex": [pair(value) for value in excitation],
            "result_excitation_complex": [pair(value) for value in excitation],
            "s_matrix_complex": [[pair(value) for value in item] for item in matrix],
            "result_s_matrix_complex": [[pair(value) for value in item] for item in matrix],
            "active_s_complex": [pair(value) for value in active],
            "result_active_s_complex": [pair(value) for value in active],
            "embedded_pattern_complex": [[[1.0, 0.0], [0.0, 1.0]], [[0.5, 0.0], [0.0, 0.5]]],
            "result_embedded_pattern_complex": [[[1.0, 0.0], [0.0, 1.0]], [[0.5, 0.0], [0.0, 0.5]]],
            "incident_power_w": 2.0, "result_incident_power_w": 2.0,
            "reflected_power_w": 0.2, "result_reflected_power_w": 0.2,
            "accepted_power_w": 1.8, "result_accepted_power_w": 1.8,
            "radiated_power_w": 1.5, "result_radiated_power_w": 1.5,
            "dissipated_power_w": 0.3, "result_dissipated_power_w": 0.3,
            "array_mesh_sha256": "3" * 64, "result_array_mesh_sha256": "3" * 64,
            "array_owner": "array/case-412", "accepted_array_owner": "array/case-412",
            "array_result_sha256": "4" * 64, "accepted_array_result_sha256": "4" * 64,
        }
    return summary


C0 = 299_792_458.0


ETA0 = 376.730313668


def _waveguide_v37(index: int):
    generation = f"waveguide-broadband-{513 + index}"
    width, length = 0.02, 0.12
    cutoff = C0 / (2.0 * width)
    frequencies = [8.0e9, 10.0e9, 12.0e9]
    factors = [math.sqrt(1.0 - (cutoff / frequency) ** 2) for frequency in frequencies]
    mirrored = {
        "waveguide_width_m": width,
        "waveguide_length_m": length,
        "mode": "TE10",
        "cutoff_frequency_hz": cutoff,
        "frequency_hz": frequencies,
        "modal_impedance_ohm": [ETA0 / factor for factor in factors],
        "propagation_constant_rad_m": [2.0 * math.pi * frequency * factor / C0 for frequency, factor in zip(frequencies, factors)],
        "group_delay_s": [length / (C0 * factor) for factor in factors],
        "power_normalization_w": 1.0,
        "mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
        "mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
        "mesh_dof": [12000, 24000, 48000],
        "mesh_cutoff_hz": [cutoff * 1.01, cutoff * 1.002, cutoff],
        "waveguide_mesh_sha256": "1" * 64,
    }
    return {
        "waveguide_generation": generation,
        **{key: generation for key in (
            "cutoff_generation", "impedance_generation", "propagation_generation",
            "group_delay_generation", "power_generation", "orthogonality_generation",
            "mesh_generation", "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "waveguide_owner": f"waveguide/case-{513 + index}",
        "accepted_waveguide_owner": f"waveguide/case-{513 + index}",
        "waveguide_result_sha256": "2" * 64,
        "accepted_waveguide_result_sha256": "2" * 64,
    }


def _probe(index: int):
    generation = f"emc-probe-{513 + index}"
    mirrored = {
        "coordinate_frame": "global_cartesian_xyz_m",
        "probe_xyz_m": [0.01, 0.0, 0.02],
        "interpolation_scheme": "trilinear_terminal_mesh",
        "interpolation_node_ids": [101, 102, 103, 104],
        "interpolation_weights": [0.25, 0.25, 0.25, 0.25],
        "time_s": [0.0, 1.0e-9, 2.0e-9, 3.0e-9],
        "field_trace_v_m": [1.0, 0.0, -1.0, 0.0],
        "fft_window": "rectangular",
        "fft_scaling": "unscaled_forward_inverse_over_n",
        "frequency_hz": [0.0, 250.0e6, 500.0e6, 750.0e6],
        "fft_complex_v_m": [[0.0, 0.0], [2.0, 0.0], [0.0, 0.0], [2.0, 0.0]],
        "selected_frequency_hz": 250.0e6,
        "selected_fft_complex_v_m": [2.0, 0.0],
        "time_energy": 2.0,
        "frequency_energy_over_n": 2.0,
        "probe_mesh_sha256": "3" * 64,
    }
    return {
        "probe_generation": generation,
        **{key: generation for key in (
            "coordinate_generation", "interpolation_generation", "time_generation",
            "window_generation", "fft_generation", "frequency_generation",
            "parseval_generation", "mesh_generation", "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "probe_owner": f"emc/probe-{513 + index}",
        "accepted_probe_owner": f"emc/probe-{513 + index}",
        "probe_result_sha256": "4" * 64,
        "accepted_probe_result_sha256": "4" * 64,
    }


def _summary_v37():
    summary = _summary_v36()
    for index, row in enumerate(summary["runs"]):
        row["waveguide_cutoff_mode_impedance_group_delay_power_orthogonality_mesh_owner_result_identity"] = _waveguide_v37(index)
        row["emc_probe_coordinate_interpolation_time_fft_window_parseval_mesh_owner_result_identity"] = _probe(index)
    return summary


def _microstrip(index: int):
    generation = f"microstrip-quasitem-{614 + index}"
    width, height, copper, relative_permittivity, length = (
        2.0e-3, 1.0e-3, 35.0e-6, 4.0, 0.1
    )
    ratio = width / height
    effective_permittivity = (
        (relative_permittivity + 1.0) / 2.0
        + (relative_permittivity - 1.0)
        / (2.0 * math.sqrt(1.0 + 12.0 / ratio))
    )
    impedance = 120.0 * math.pi / (
        math.sqrt(effective_permittivity)
        * (ratio + 1.393 + 0.667 * math.log(ratio + 1.444))
    )
    conductor_loss, dielectric_loss, s11 = 0.02, 0.01, 0.1
    s21 = math.sqrt(1.0 - s11**2 - conductor_loss - dielectric_loss)
    mirrored = {
        "trace_width_m": width,
        "substrate_height_m": height,
        "copper_thickness_m": copper,
        "relative_permittivity": relative_permittivity,
        "line_length_m": length,
        "quasitem_model": "hammerstad_zero_thickness_core",
        "effective_permittivity": effective_permittivity,
        "characteristic_impedance_ohm": impedance,
        "propagation_delay_s": length * math.sqrt(effective_permittivity) / C0,
        "frequency_hz": [1.0e9, 2.0e9, 3.0e9],
        "conductor_loss_fraction": conductor_loss,
        "dielectric_loss_fraction": dielectric_loss,
        "s11_magnitude": s11,
        "s21_magnitude": s21,
        "incident_power_w": 1.0,
        "mesh_dof": [20_000, 40_000, 80_000],
        "mesh_impedance_ohm": [impedance * 1.02, impedance * 1.004, impedance],
        "maximum_final_relative_change": 0.01,
        "microstrip_mesh_sha256": "1" * 64,
    }
    return {
        "microstrip_generation": generation,
        **{key: generation for key in (
            "geometry_generation", "material_generation", "impedance_generation",
            "delay_generation", "loss_generation", "sparameter_generation",
            "frequency_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "microstrip_owner": f"pcb/microstrip-{614 + index}",
        "accepted_microstrip_owner": f"pcb/microstrip-{614 + index}",
        "microstrip_result_sha256": "2" * 64,
        "accepted_microstrip_result_sha256": "2" * 64,
    }


def _shielding(index: int):
    generation = f"shield-aperture-{614 + index}"
    incident = [1.0, 1.0, 1.0]
    transmitted = [0.1, 0.05, 0.02]
    incident_power = [item**2 for item in incident]
    transmitted_power = [item**2 for item in transmitted]
    field_se = [20.0 * math.log10(a / b) for a, b in zip(incident, transmitted)]
    power_se = [
        10.0 * math.log10(a / b)
        for a, b in zip(incident_power, transmitted_power)
    ]
    mirrored = {
        "aperture_plane": "xy_normal_positive_z",
        "incident_polarization": "x_linear",
        "frequency_hz": [1.0e9, 2.0e9, 3.0e9],
        "incident_field_v_m": incident,
        "transmitted_field_v_m": transmitted,
        "incident_power_density_normalized": incident_power,
        "transmitted_power_density_normalized": transmitted_power,
        "shielding_effectiveness_field_db": field_se,
        "shielding_effectiveness_power_db": power_se,
        "probe_frame": "global_cartesian_xyz_m",
        "mesh_dof": [30_000, 60_000, 120_000],
        "mesh_selected_se_db": [power_se[1] - 1.0, power_se[1] - 0.2, power_se[1]],
        "maximum_final_se_change_db": 0.5,
        "shield_mesh_sha256": "3" * 64,
    }
    return {
        "shield_generation": generation,
        **{key: generation for key in (
            "aperture_generation", "polarization_generation", "field_generation",
            "power_generation", "se_generation", "frequency_generation",
            "probe_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "shield_owner": f"emc/shield-{614 + index}",
        "accepted_shield_owner": f"emc/shield-{614 + index}",
        "shield_result_sha256": "4" * 64,
        "accepted_shield_result_sha256": "4" * 64,
    }


def _summary_v38():
    summary = _summary_v37()
    for index, row in enumerate(summary["runs"]):
        row[
            "microstrip_quasitem_geometry_permittivity_impedance_delay_loss_sparameter_frequency_mesh_owner_result_identity"
        ] = _microstrip(index)
        row[
            "shielding_aperture_orientation_polarization_field_power_se_frequency_probe_mesh_owner_result_identity"
        ] = _shielding(index)
    return summary


_WAVEGUIDE_KEY = (
    "waveguide_mode_cutoff_impedance_power_orthogonality_propagation_deembed_"
    "mesh_owner_result_identity"
)


_ANTENNA_KEY = (
    "antenna_nearfar_directivity_gain_efficiency_polarization_power_mesh_owner_"
    "result_identity"
)


def _waveguide_v39(index: int) -> dict:
    generation = f"waveguide-mode-{715 + index}"
    width, height, frequency = 0.02, 0.01, 12.0e9
    cutoff = C0 / (2.0 * width)
    factor = math.sqrt(1.0 - (cutoff / frequency) ** 2)
    impedance = ETA0 / factor
    beta = 2.0 * math.pi * frequency * factor / C0
    distance, raw_phase = 0.012, -1.2
    mirrored = {
        "waveguide_width_m": width,
        "waveguide_height_m": height,
        "mode_name": "TE10",
        "mode_index": [1, 0],
        "cutoff_frequency_hz": cutoff,
        "frequency_hz": frequency,
        "modal_impedance_ohm": impedance,
        "propagation_constant_rad_m": beta,
        "normalized_forward_power_w": 1.0,
        "mode_overlap_real": [[1.0, 0.0], [0.0, 1.0]],
        "mode_overlap_imag": [[0.0, 0.0], [0.0, 0.0]],
        "orthogonality_tolerance": 1.0e-9,
        "port_plane_m": 0.0,
        "reference_plane_m": distance,
        "deembed_distance_m": distance,
        "deembed_convention": "port_to_reference_add_beta_l",
        "raw_s21_phase_rad": raw_phase,
        "deembedded_s21_phase_rad": raw_phase + beta * distance,
        "waveguide_mesh_sha256": "1" * 64,
    }
    return {
        "waveguide_mode_generation": generation,
        **{
            key: generation
            for key in (
                "cutoff_generation",
                "impedance_generation",
                "power_generation",
                "orthogonality_generation",
                "propagation_generation",
                "deembed_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "waveguide_mode_owner": f"waveguide/mode-{715 + index}",
        "accepted_waveguide_mode_owner": f"waveguide/mode-{715 + index}",
        "waveguide_mode_result_sha256": "2" * 64,
        "accepted_waveguide_mode_result_sha256": "2" * 64,
    }


def _antenna(index: int) -> dict:
    generation = f"antenna-nearfar-{715 + index}"
    directivity, efficiency, mismatch = 4.0, 0.8, 0.9
    gain = directivity * efficiency
    realized = gain * mismatch
    mirrored = {
        "frequency_hz": 2.4e9,
        "nearfield_surface_closed": True,
        "near_to_far_transform": "equivalence_surface_stratton_chu",
        "polarization_basis": "ieee_theta_phi",
        "co_polar_component": "theta",
        "cross_polar_component": "phi",
        "directivity_linear": directivity,
        "directivity_dbi": 10.0 * math.log10(directivity),
        "radiation_efficiency": efficiency,
        "gain_linear": gain,
        "mismatch_efficiency": mismatch,
        "realized_gain_linear": realized,
        "realized_gain_dbi": 10.0 * math.log10(realized),
        "accepted_power_w": 1.0,
        "radiated_power_w": efficiency,
        "loss_power_w": 1.0 - efficiency,
        "power_balance_residual_w": 0.0,
        "farfield_sphere_samples": 2592,
        "antenna_mesh_sha256": "3" * 64,
    }
    return {
        "antenna_generation": generation,
        **{
            key: generation
            for key in (
                "nearfield_generation",
                "farfield_generation",
                "directivity_generation",
                "gain_generation",
                "efficiency_generation",
                "polarization_generation",
                "power_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "antenna_owner": f"antenna/nearfar-{715 + index}",
        "accepted_antenna_owner": f"antenna/nearfar-{715 + index}",
        "antenna_result_sha256": "4" * 64,
        "accepted_antenna_result_sha256": "4" * 64,
    }


def _summary_v39() -> dict:
    summary = _summary_v38()
    for index, row in enumerate(summary["runs"]):
        row[_WAVEGUIDE_KEY] = _waveguide_v39(index)
        row[_ANTENNA_KEY] = _antenna(index)
    return summary


_VIA = "pcb_via_stub_geometry_impedance_resonance_sparameter_current_loss_power_mesh_owner_result_identity"


_CAVITY = "cavity_degenerate_mode_frequency_orthogonality_energy_symmetry_quality_mesh_owner_result_identity"


def _via(index: int) -> dict:
    generation = f"pcb-via-{724 + index}"
    board, diameter, stub, epsilon_r = 1.6e-3, 0.30e-3, 0.80e-3, 4.0
    impedance = 60.0 / math.sqrt(epsilon_r) * math.log(4.0 * board / diameter)
    resonance = C0 / (4.0 * stub * math.sqrt(epsilon_r))
    s11, s21 = [0.2, -0.1], [0.85, -0.05]
    incident = 1.0
    reflected = incident * sum(item * item for item in s11)
    accepted = incident - reflected
    transmitted = incident * sum(item * item for item in s21)
    conductor_loss = 0.12
    values = {
        "board_thickness_m": board, "via_diameter_m": diameter,
        "stub_length_m": stub, "relative_permittivity": epsilon_r,
        "characteristic_impedance_ohm": impedance,
        "quarterwave_resonance_hz": resonance, "s11_complex": s11,
        "s21_complex": s21, "incident_power_w": incident,
        "reflected_power_w": reflected, "accepted_power_w": accepted,
        "transmitted_power_w": transmitted,
        "barrel_current_rms_a": math.sqrt(accepted / impedance),
        "conductor_loss_w": conductor_loss,
        "dielectric_loss_w": accepted - transmitted - conductor_loss,
        "power_balance_residual_w": 0.0, "pcb_mesh_sha256": "1" * 64,
    }
    return {
        "via_generation": generation,
        **{key: generation for key in ("geometry_generation", "impedance_generation", "resonance_generation", "sparameter_generation", "current_generation", "loss_generation", "power_generation", "mesh_generation", "owner_generation", "result_generation")},
        **values, **{f"result_{key}": value for key, value in values.items()},
        "via_owner": f"pcb/via-{724 + index}",
        "accepted_via_owner": f"pcb/via-{724 + index}",
        "via_result_sha256": "2" * 64, "accepted_via_result_sha256": "2" * 64,
    }


def _cavity(index: int) -> dict:
    generation = f"cavity-degenerate-{724 + index}"
    side = 0.03
    frequency = C0 / 2.0 * math.sqrt(2.0 / side**2)
    values = {
        "cavity_dimensions_m": [side, side, side],
        "mode_names": ["TE101", "TE011"],
        "mode_frequencies_hz": [frequency, frequency],
        "degeneracy_tolerance_relative": 1.0e-9,
        "modal_overlap_matrix": [[1.0, 0.0], [0.0, 1.0]],
        "electric_energy_j": [0.5, 0.5], "magnetic_energy_j": [0.5, 0.5],
        "symmetry_classes": ["even_x_odd_y", "odd_x_even_y"],
        "quality_factors": [5000.0, 5000.0], "cavity_mesh_sha256": "3" * 64,
    }
    return {
        "cavity_generation": generation,
        **{key: generation for key in ("frequency_generation", "orthogonality_generation", "energy_generation", "symmetry_generation", "quality_generation", "mesh_generation", "owner_generation", "result_generation")},
        **values, **{f"result_{key}": value for key, value in values.items()},
        "cavity_owner": f"cavity/degenerate-{724 + index}",
        "accepted_cavity_owner": f"cavity/degenerate-{724 + index}",
        "cavity_result_sha256": "4" * 64, "accepted_cavity_result_sha256": "4" * 64,
    }


def _summary_v40() -> dict:
    summary = _summary_v39()
    for index, run in enumerate(summary["runs"]):
        run[_VIA] = _via(index)
        run[_CAVITY] = _cavity(index)
    return summary


_WAVEGUIDE_v41 = "waveguide_cutoff_propagation_impedance_mode_orthogonality_sparameter_power_owner_result_identity"


_ANTENNA = "antenna_farfield_directivity_gain_efficiency_power_polarization_mesh_owner_result_identity"


MU0 = 4.0e-7 * math.pi


def _summary_v41() -> dict:
    summary = _summary_v40()
    for index, run in enumerate(summary["runs"]):
        generation = f"waveguide-te10-{724 + index}"
        width, height, frequency = 22.86e-3, 10.16e-3, 10.0e9
        omega = 2.0 * math.pi * frequency
        cutoff = C0 / (2.0 * width)
        propagation = math.sqrt((omega / C0) ** 2 - (math.pi / width) ** 2)
        s11, s21 = [0.1, 0.0], [math.sqrt(0.97), 0.0]
        incident = 1.0
        reflected = incident * sum(item * item for item in s11)
        transmitted = incident * sum(item * item for item in s21)
        values = {
            "waveguide_width_m": width, "waveguide_height_m": height,
            "mode_name": "TE10", "frequency_hz": frequency,
            "cutoff_frequency_hz": cutoff,
            "propagation_constant_rad_per_m": propagation,
            "guide_wavelength_m": 2.0 * math.pi / propagation,
            "guide_impedance_ohm": omega * MU0 / propagation,
            "modal_overlap_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "s11_complex": s11, "s21_complex": s21,
            "incident_power_w": incident, "reflected_power_w": reflected,
            "transmitted_power_w": transmitted,
            "wall_loss_w": incident - reflected - transmitted,
            "power_balance_residual_w": 0.0,
        }
        run[_WAVEGUIDE_v41] = {
            "waveguide_generation": generation,
            **{key: generation for key in ("geometry_generation", "mode_generation", "cutoff_generation", "propagation_generation", "impedance_generation", "orthogonality_generation", "sparameter_generation", "power_generation", "owner_generation", "result_generation")},
            **values, **{f"result_{key}": value for key, value in values.items()},
            "waveguide_owner": f"waveguide/te10-{724 + index}",
            "accepted_waveguide_owner": f"waveguide/te10-{724 + index}",
            "waveguide_result_sha256": "5" * 64,
            "accepted_waveguide_result_sha256": "5" * 64,
        }

        generation = f"antenna-farfield-{724 + index}"
        incident, reflected = 10.0, 1.0
        accepted, conductor_loss, dielectric_loss = 9.0, 0.5, 0.5
        radiated = accepted - conductor_loss - dielectric_loss
        efficiency, directivity = radiated / accepted, 6.0
        values = {
            "incident_power_w": incident, "reflected_power_w": reflected,
            "accepted_power_w": accepted, "radiated_power_w": radiated,
            "conductor_loss_w": conductor_loss, "dielectric_loss_w": dielectric_loss,
            "radiation_efficiency": efficiency, "directivity_linear": directivity,
            "gain_linear": directivity * efficiency,
            "maximum_radiation_intensity_w_per_sr": directivity * radiated / (4.0 * math.pi),
            "polarization_basis": "linear_xy", "co_polar_fraction": 1.0,
            "cross_polar_fraction": 0.0, "power_balance_residual_w": 0.0,
        }
        run[_ANTENNA] = {
            "antenna_generation": generation,
            **{key: generation for key in ("excitation_generation", "farfield_generation", "directivity_generation", "gain_generation", "efficiency_generation", "power_generation", "polarization_generation", "mesh_generation", "owner_generation", "result_generation")},
            **values, **{f"result_{key}": value for key, value in values.items()},
            "antenna_owner": f"antenna/farfield-{724 + index}",
            "accepted_antenna_owner": f"antenna/farfield-{724 + index}",
            "antenna_result_sha256": "6" * 64,
            "accepted_antenna_result_sha256": "6" * 64,
        }
    return summary


_EMC_v43 = "emc_shielding_incident_transmitted_reflected_poynting_se_power_energy_generation_identity"


_CONNECTOR = "differential_connector_sdd_scc_modeconversion_passivity_loss_power_generation_identity"


def _summary_v43() -> dict:
    summary = _summary_v41()
    for index, run in enumerate(summary["runs"]):
        generation = f"emc-shield-843-{index}"
        values = {
            "incident_power_w": 100.0,
            "transmitted_power_w": 1.0,
            "reflected_power_w": 9.0,
            "absorbed_power_w": 90.0,
            "shielding_effectiveness_db": 20.0,
            "poynting_flux_orientation": "outward",
            "power_closure_residual_w": 0.0,
            "mesh_owner": f"mesh:{generation}",
        }
        run[_EMC_v43] = {
            "emc_generation": generation,
            "generations": {name: generation for name in ("incident", "transmitted", "reflected", "poynting", "shielding", "power", "mesh", "result")},
            "values": values,
            "result_values": dict(values),
            "emc_result_sha256": "7" * 64,
            "accepted_emc_result_sha256": "7" * 64,
        }
        generation = f"connector-mixed-843-{index}"
        values = {
            "port_order": ["P1+", "P1-", "P2+", "P2-"],
            "mixed_mode_order": ["Sdd11", "Sdd21", "Scc11", "Scc21"],
            "sdd_ri": [[0.1, 0.0], [0.8, 0.0]],
            "scc_ri": [[0.2, 0.0], [0.7, 0.0]],
            "return_loss_db": [20.0, 1.938200260161128],
            "insertion_loss_db": [1.938200260161128, 3.098039199714863],
            "passivity_margin": 0.1,
            "dissipated_power_w": 0.2,
            "reference_impedance_ohm": 100.0,
            "mesh_owner": f"mesh:{generation}",
        }
        run[_CONNECTOR] = {
            "connector_generation": generation,
            "generations": {name: generation for name in ("port", "mixed_mode", "sdd", "scc", "loss", "passivity", "power", "impedance", "mesh", "result")},
            "values": values,
            "result_values": dict(values),
            "connector_result_sha256": "9" * 64,
            "accepted_connector_result_sha256": "9" * 64,
        }
    return summary


_WAVEGUIDE_v44 = "waveguide_modal_cutoff_impedance_groupdelay_power_orthogonality_mesh_result_identity"


_EMC_v44 = "emc_probe_interpolation_timewindow_fft_parseval_coordinate_monitor_owner_identity"


def _payload_v44() -> dict:
    wave = "waveguide-mode-844-0"
    emc = "emc-probe-844-0"
    return {"runs": [{
        _WAVEGUIDE_v44: {
            "waveguide_generation": wave,
            **{key: wave for key in ("frequency_generation", "cutoff_generation", "impedance_generation", "groupdelay_generation", "power_generation", "orthogonality_generation", "mesh_generation", "result_generation")},
            "frequency_hz": [8.0e9, 10.0e9, 12.0e9], "result_frequency_hz": [8.0e9, 10.0e9, 12.0e9], "cutoff_frequency_hz": 6.0e9, "result_cutoff_frequency_hz": 6.0e9,
            "modal_impedance_ohm": [50.0, 55.0, 60.0], "result_modal_impedance_ohm": [50.0, 55.0, 60.0], "group_delay_s": [1.0e-9, 1.1e-9, 1.2e-9], "result_group_delay_s": [1.0e-9, 1.1e-9, 1.2e-9],
            "power_normalization_w": [1.0, 1.0, 1.0], "result_power_normalization_w": [1.0, 1.0, 1.0], "orthogonality_matrix": [[1.0, 0.0], [0.0, 1.0]], "result_orthogonality_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "mesh_owner": "mesh:" + wave, "result_mesh_owner": "mesh:" + wave, "waveguide_result_sha256": "d" * 64, "accepted_waveguide_result_sha256": "d" * 64,
        },
        _EMC_v44: {
            "emc_probe_generation": emc,
            **{key: emc for key in ("coordinate_generation", "interpolation_generation", "timewindow_generation", "fft_generation", "parseval_generation", "monitor_generation", "result_generation")},
            "coordinate_system": "global_cartesian", "result_coordinate_system": "global_cartesian", "interpolation_order": 2, "result_interpolation_order": 2, "time_window_s": [0.0, 1.0e-9], "result_time_window_s": [0.0, 1.0e-9],
            "fft_normalization": "parseval_unitary", "result_fft_normalization": "parseval_unitary", "parseval_time_energy_j": 1.0, "result_parseval_time_energy_j": 1.0, "parseval_frequency_energy_j": 1.0, "result_parseval_frequency_energy_j": 1.0,
            "monitor_owner": "monitor:" + emc, "result_monitor_owner": "monitor:" + emc, "emc_probe_result_sha256": "e" * 64, "accepted_emc_probe_result_sha256": "e" * 64,
        },
    }]}


def _payload_v51() -> dict[str, object]:
    generation = "network-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    frequencies = [1.0e9, 1.5e9, 2.0e9, 2.5e9, 3.0e9]
    phase = [0.0, -18.0, -36.0, -54.0, -72.0]
    delay = [1.0e-10] * 5
    run = {
        S_PARAMETER: {
            "generation": generation, "renormalization_generation": generation, "zref_generation": generation,
            "modal_generation": generation, "wave_generation": generation, "owner_generation": generation,
            "result_generation": generation, "renormalized": True, "result_renormalized": True,
            "complex_reference_impedance_ohm": [50.0, 5.0], "result_complex_reference_impedance_ohm": [50.0, 5.0],
            "modal_impedances_ohm": {"port1_mode1": [48.0, 2.0], "port2_mode1": [52.0, -1.0]},
            "result_modal_impedances_ohm": {"port1_mode1": [48.0, 2.0], "port2_mode1": [52.0, -1.0]},
            "wave_basis": "power_waves", "result_wave_basis": "power_waves", "port_owner": "port:network-v51",
            "result_port_owner": "port:network-v51", **result,
        },
        GROUP_DELAY: {
            "generation": generation, "unwrap_generation": generation, "frequency_generation": generation,
            "derivative_generation": generation, "smoothing_generation": generation, "owner_generation": generation,
            "result_generation": generation, "frequency_hz": frequencies, "result_frequency_hz": frequencies,
            "unwrapped_phase_deg": phase, "result_unwrapped_phase_deg": phase, "group_delay_s": delay,
            "result_group_delay_s": delay, "derivative_definition": "minus_dphi_rad_domega",
            "result_derivative_definition": "minus_dphi_rad_domega", "smoothing_window_points": 5,
            "result_smoothing_window_points": 5, "trace_owner": "trace:s21-v51", "result_trace_owner": "trace:s21-v51", **result,
        },
    }
    return {"runs": [deepcopy(run), deepcopy(run)]}


def _generations_v52(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload_v52() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    dissipated = {"ohmic": 40.0, "dielectric": 20.0}
    q_factor = 2.0 * math.pi * 1.0e9 * 0.01 / 3.0
    run = {
        ENERGY_BALANCE: {
            **_generations_v52("energy-balance-v52", ("incident_generation", "scattered_generation", "absorbed_generation", "dissipated_generation", "owner_generation", "result_generation")),
            "incident_power_w": 100.0, "result_incident_power_w": 100.0,
            "reflected_power_w": 10.0, "result_reflected_power_w": 10.0,
            "transmitted_power_w": 30.0, "result_transmitted_power_w": 30.0,
            "absorbed_power_w": 60.0, "result_absorbed_power_w": 60.0,
            "dissipated_power_w": dissipated, "result_dissipated_power_w": dissipated,
            "run_owner": "run:energy-balance-v52", "result_run_owner": "run:energy-balance-v52", **result,
        },
        EIGENMODE_Q: {
            **_generations_v52("eigenmode-q-v52", ("frequency_generation", "energy_generation", "loss_generation", "normalization_generation", "owner_generation", "result_generation")),
            "frequency_hz": 1.0e9, "result_frequency_hz": 1.0e9,
            "stored_energy_j": 0.01, "result_stored_energy_j": 0.01,
            "boundary_loss_w": 2.0, "result_boundary_loss_w": 2.0,
            "volume_loss_w": 1.0, "result_volume_loss_w": 1.0,
            "q_factor": q_factor, "result_q_factor": q_factor,
            "normalization": "physical_stored_energy", "result_normalization": "physical_stored_energy",
            "mode_owner": "mode:eigenmode-v52", "result_mode_owner": "mode:eigenmode-v52", **result,
        },
    }
    return {"runs": [deepcopy(run), deepcopy(run)]}


def _generations_v53(generation: str, fields: tuple[str, ...]) -> dict[str, str]: return {"generation": generation, **{field: generation for field in fields}}


def _payload_v53():
    width = 0.02286; cutoff = 299792458.0 / (2.0 * width)
    waveguide = {**_generations_v53("wave-v53", ("mode_generation", "cutoff_generation", "normalization_generation", "plane_generation", "owner_generation", "result_generation")), "mode_id": "TE10", "result_mode_id": "TE10", "waveguide_width_m": width, "result_waveguide_width_m": width, "cutoff_frequency_hz": cutoff, "result_cutoff_frequency_hz": cutoff, "sample_frequency_hz": 10.0e9, "result_sample_frequency_hz": 10.0e9, "normalization": "unit_incident_power_w", "result_normalization": "unit_incident_power_w", "reference_plane_offset_m": -0.01, "result_reference_plane_offset_m": -0.01, "port_owner": "port:v53", "result_port_owner": "port:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    theta = [0.0, 90.0, 180.0]; phi = [0.0, 90.0, 180.0, 270.0]; gain = [[1.0] * 4, [8.0] * 4, [1.0] * 4]
    farfield = {**_generations_v53("far-v53", ("gain_generation", "polarization_generation", "grid_generation", "owner_generation", "result_generation")), "theta_deg": theta, "result_theta_deg": theta, "phi_deg": phi, "result_phi_deg": phi, "realized_gain_dbi": gain, "result_realized_gain_dbi": gain, "polarization_basis": "ludwig3", "result_polarization_basis": "ludwig3", "monitor_owner": "monitor:v53", "result_monitor_owner": "monitor:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {"runs": [{WAVEGUIDE: waveguide, FARFIELD: farfield}]}


_C0 = 299792458.0


_ETA0 = 376.730313668


def _generations_v54(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload_v54():
    width = 0.02286
    cutoff_hz = _C0 / (2.0 * width)
    frequency_hz = 10.0e9
    impedance = _ETA0 / (1.0 - (cutoff_hz / frequency_hz) ** 2) ** 0.5
    cutoff = {
        **_generations_v54("cutoff-v54", ("cutoff_generation", "mode_generation", "normalization_generation", "power_generation", "impedance_generation", "owner_generation", "result_generation")),
        "mode_id": "TE10", "result_mode_id": "TE10",
        "waveguide_width_m": width, "result_waveguide_width_m": width,
        "cutoff_frequency_hz": cutoff_hz, "result_cutoff_frequency_hz": cutoff_hz,
        "sample_frequency_hz": frequency_hz, "result_sample_frequency_hz": frequency_hz,
        "field_normalization": "unit_modal_power_w", "result_field_normalization": "unit_modal_power_w",
        "modal_power_w": 1.0, "result_modal_power_w": 1.0,
        "mode_impedance_ohm": impedance, "result_mode_impedance_ohm": impedance,
        "port_owner": "port:v54", "result_port_owner": "port:v54",
        "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
    }
    voxels = [1.0e-6] * 10
    fields = [12.0 + value for value in range(10)]
    sar = {
        **_generations_v54("sar-v54", ("mass_generation", "density_generation", "voxel_generation", "frequency_generation", "field_generation", "owner_generation", "result_generation")),
        "averaging_mass_kg": 0.01, "result_averaging_mass_kg": 0.01,
        "tissue_density_kg_m3": 1000.0, "result_tissue_density_kg_m3": 1000.0,
        "voxel_support_m3": voxels, "result_voxel_support_m3": voxels,
        "electric_field_rms_v_m": fields, "result_electric_field_rms_v_m": fields,
        "frequency_hz": 2.45e9, "result_frequency_hz": 2.45e9,
        "field_solution_sha256": "2" * 64, "result_field_solution_sha256": "2" * 64,
        "monitor_owner": "monitor:v54", "result_monitor_owner": "monitor:v54",
        "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
    }
    return {"runs": [{CUTOFF: cutoff, SAR: sar}]}


CASE_IDS = {"v55_public_resonator_loaded_unloaded_q_coupling_linewidth_energy_owner_mismatch", "v55_public_antenna_efficiency_accepted_radiated_loss_gain_directivity_owner_mismatch"}


def _payload_v55():
    gen=lambda name,fields:{"generation":name,**{field:name for field in fields}}
    f0=2.4e9; ql=1200.0; q0=3000.0; qe=2000.0; power=1.0
    resonator={**gen("resonator-v55",("q_generation","coupling_generation","linewidth_generation","energy_generation","owner_generation","result_generation")),"resonance_frequency_hz":f0,"result_resonance_frequency_hz":f0,"loaded_q":ql,"result_loaded_q":ql,"unloaded_q":q0,"result_unloaded_q":q0,"external_q":qe,"result_external_q":qe,"coupling_beta":q0/qe,"result_coupling_beta":q0/qe,"linewidth_hz":f0/ql,"result_linewidth_hz":f0/ql,"dissipated_power_w":power,"result_dissipated_power_w":power,"stored_energy_j":ql*power/(2*math.pi*f0),"result_stored_energy_j":ql*power/(2*math.pi*f0),"monitor_owner":"monitor:v55","result_monitor_owner":"monitor:v55","result_sha256":"d"*64,"accepted_result_sha256":"d"*64}
    accepted=1.0;radiated=.75;loss=.25;eta=.75;directivity=6.0;gain=directivity+10*math.log10(eta)
    antenna={**gen("antenna-v55",("power_generation","efficiency_generation","gain_generation","directivity_generation","owner_generation","result_generation")),"accepted_power_w":accepted,"result_accepted_power_w":accepted,"radiated_power_w":radiated,"result_radiated_power_w":radiated,"loss_power_w":loss,"result_loss_power_w":loss,"radiation_efficiency":eta,"result_radiation_efficiency":eta,"directivity_dbi":directivity,"result_directivity_dbi":directivity,"gain_dbi":gain,"result_gain_dbi":gain,"farfield_owner":"farfield:v55","result_farfield_owner":"farfield:v55","result_sha256":"e"*64,"accepted_result_sha256":"e"*64}
    return {"runs":[{RESONATOR:resonator,ANTENNA:antenna}]}
