"""RF/wave-port, S-parameter, far-field and EMC identity generalization gates.

One cumulative positive closure per gate; every negative keeps its own failure signal.
"""
from __future__ import annotations

import math
from copy import deepcopy
from radia_mcp.radia_ngsolve.network_artifact_identity_v51 import (
    GROUP_DELAY,
    S_PARAMETER,
    validate_public_v51_identity,
)
from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import nonlinear_inductance_sweep_gate
from radia_mcp.radia_ngsolve.wave_energy_identity_v52 import (
    ENERGY_BALANCE,
    EIGENMODE_Q,
    validate_public_v52_identity,
)
from radia_mcp.radia_ngsolve.wave_energy_identity_v55 import (
    ANTENNA,
    RESONATOR,
    validate_public_v55_identity,
)
from radia_mcp.radia_ngsolve.wave_port_identity_v53 import (
    FARFIELD,
    WAVEGUIDE,
    validate_public_v53_identity,
)
from radia_mcp.radia_ngsolve.wave_sar_identity_v54 import CUTOFF, SAR, validate_public_v54_identity
from radia_mcp.radia_ngsolve.waveguide_emc_v44_identity import validate_public_identity

from _rf_wave_identity_payloads import (
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
    _WAVEGUIDE_KEY,
    _ANTENNA_KEY,
    _summary_v39,
    _VIA,
    _CAVITY,
    _summary_v40,
    _WAVEGUIDE_v41,
    _ANTENNA,
    _summary_v41,
    _EMC_v43,
    _CONNECTOR,
    _summary_v43,
    _WAVEGUIDE_v44,
    _EMC_v44,
    _payload_v44,
    _payload_v51,
    _payload_v52,
    _payload_v53,
    _payload_v54,
    _payload_v55,
)


def test_v23_public_broadband_adaptive_mesh_sparam_renormalization_port_generation_mismatch() -> None:
    summary = _summary_v23()
    summary["runs"][0][
        "broadband_adaptive_mesh_sparam_renormalization_port_generation_identity"
    ].update(
        {
            "adaptive_mesh_sweep_generation": "broadband-sparam-50",
            "frequency_interpolation_sweep_generation": "broadband-sparam-49",
            "port_mode_sweep_generation": "broadband-sparam-48",
            "renormalization_sweep_generation": "broadband-sparam-47",
            "sparameter_result_sweep_generation": "broadband-sparam-46",
            "result_adaptive_mesh_sha256": "a" * 64,
            "result_frequency_samples_hz": [1.0e9, 1.6e9, 2.0e9],
            "result_frequency_interpolation": "linear",
            "result_port_mode_ids": ["P2:M1", "P1:M2"],
            "result_renormalization_impedance_ohm": [[75.0, 0.0], [50.0, 5.0]],
            "result_sparameter_table_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "broadband_sparameters_use_current_mesh_interpolation_modes_and_impedance"
    ]


def test_v23_public_transient_monitor_time_origin_excitation_waveform_mesh_generation_mismatch() -> None:
    summary = _summary_v23()
    summary["runs"][0][
        "transient_monitor_time_origin_excitation_waveform_mesh_generation_identity"
    ].update(
        {
            "time_origin_transient_generation": "transient-monitor-50",
            "excitation_waveform_transient_generation": "transient-monitor-49",
            "monitor_frame_transient_generation": "transient-monitor-48",
            "mesh_transient_generation": "transient-monitor-47",
            "field_result_transient_generation": "transient-monitor-46",
            "result_time_origin_s": 1.0e-9,
            "result_excitation_waveform_sha256": "c" * 64,
            "result_monitor_coordinate_frame": "port_local",
            "result_monitor_ids": [102, 103],
            "result_mesh_sha256": "d" * 64,
            "result_time_samples_s": [1.0e-9, 1.1e-9, 1.2e-9],
            "result_monitor_field_table_sha256": "e" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "transient_monitors_use_current_time_waveform_frame_and_mesh"
    ]


def test_v24_public_deembedding_phase_causality_passivity_grid_mismatch() -> None:
    summary = _summary_v24()
    summary["runs"][0][
        "deembedding_reference_plane_phase_causality_passivity_grid_generation_identity"
    ].update(
        {
            "reference_plane_deembedding_generation": "deembed-100",
            "phase_deembedding_generation": "deembed-99",
            "causality_deembedding_generation": "deembed-98",
            "passivity_deembedding_generation": "deembed-97",
            "frequency_grid_deembedding_generation": "deembed-96",
            "result_reference_plane_offsets_m": [0.002, -0.001],
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "result_unwrapped_phase_rad": [-0.2, 5.9, -0.5],
            "result_causality_check_passed": False,
            "result_passivity_max_singular_values": [0.82, 1.12, 0.86],
            "result_deembedded_network_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "deembedded_network_uses_current_planes_phase_causality_passivity_and_grid"
    ]


def test_v24_public_field_circuit_cosim_sign_impedance_power_mismatch() -> None:
    summary = _summary_v24()
    summary["runs"][0][
        "field_circuit_cosim_port_sign_impedance_power_balance_generation_identity"
    ].update(
        {
            "field_port_cosim_generation": "field-circuit-100",
            "circuit_port_cosim_generation": "field-circuit-99",
            "sign_cosim_generation": "field-circuit-98",
            "impedance_cosim_generation": "field-circuit-97",
            "power_balance_cosim_generation": "field-circuit-96",
            "result_current_sign_convention": "positive_out_of_field_port",
            "result_voltage_reference": "negative_to_positive_terminal",
            "result_port_impedance_ri_ohm": [-46.15384615384615, -19.23076923076923],
            "circuit_delivered_power_w": -1.2,
            "result_power_balance_residual_w": 3.16,
            "reported_cosim_result_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "field_circuit_cosim_uses_current_sign_impedance_and_power_balance"
    ]


def test_v25_public_adaptive_mesh_convergence_generation_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["runs"][0][
        "adaptive_mesh_pass_sparameter_energy_convergence_grid_generation_identity"
    ]
    identity.update(
        {
            "mesh_pass_adaptive_generation": "adaptive-pass-200",
            "sparameter_adaptive_generation": "adaptive-pass-199",
            "energy_adaptive_generation": "adaptive-pass-198",
            "frequency_grid_adaptive_generation": "adaptive-pass-197",
            "stopping_rule_adaptive_generation": "adaptive-pass-196",
            "result_mesh_pass_ids": [0, 2, 1],
            "result_mesh_cell_counts": [10000, 29000, 18000],
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "result_maximum_sparameter_delta": [0.1, 0.03, 0.02],
            "result_stored_energy_closure_residual": [0.05, 0.01, 0.02],
            "result_converged_pass_id": 1,
            "reported_adaptive_result_sha256": "d" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "adaptive_results_use_current_mesh_pass_sparameter_energy_grid_and_stop_rule"
    ]


def test_v25_public_eigenmode_tracking_generation_mismatch() -> None:
    summary = _summary_v25()
    identity = summary["runs"][0][
        "eigenmode_tracking_phase_normalization_port_coupling_mesh_generation_identity"
    ]
    identity.update(
        {
            "modal_subspace_tracking_generation": "eigenmode-track-200",
            "phase_tracking_generation": "eigenmode-track-199",
            "normalization_tracking_generation": "eigenmode-track-198",
            "port_coupling_tracking_generation": "eigenmode-track-197",
            "mesh_tracking_generation": "eigenmode-track-196",
            "result_sweep_parameters": [1.0, 0.5, 0.0],
            "result_tracked_mode_ids": ["mode-2", "mode-1"],
            "result_modal_subspace_sha256": ["4" * 64, "3" * 64, "2" * 64],
            "result_phase_anchor_ids": ["probe-hy", "probe-ez"],
            "result_normalization": "peak_field_1",
            "result_port_coupling_magnitudes": [[0.1, 0.8]],
            "result_mesh_sha256": ["7" * 64, "6" * 64, "5" * 64],
            "reported_eigenmode_track_sha256": "e" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "eigenmodes_use_current_subspace_phase_normalization_ports_and_mesh"
    ]


def test_v26_public_rejects_port_network_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary["runs"][0][
        "port_deembedding_reference_plane_impedance_mode_normalization_smatrix_generation_identity"
    ]
    identity.update(
        {
            "port_mode_network_generation": "port-network-300",
            "result_port_mode_ids": ["P2:M1", "P1:M2"],
            "result_reference_plane_offsets_m": [0.0, 0.002],
            "result_reference_impedance_ri_ohm": [[75.0, 0.0], [50.0, 5.0]],
            "result_wave_normalization": "voltage_wave",
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "reported_smatrix_sha256": "8" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "sparameters_use_current_port_modes_planes_impedances_normalization_grid_and_result"
    ]


def test_v26_public_rejects_farfield_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary["runs"][0][
        "farfield_angular_grid_polarization_coordinate_power_normalization_mesh_generation_identity"
    ]
    identity.update(
        {
            "angular_grid_farfield_generation": "farfield-300",
            "result_theta_deg": [90.0, 45.0, 0.0],
            "result_polarization_basis": "spherical_theta_phi",
            "result_coordinate_frame": "local_wcs_y_up",
            "result_radiated_power_w": 1.2,
            "result_field_normalization": "peak_field",
            "result_mesh_sha256": "9" * 64,
            "reported_farfield_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "farfields_use_current_angular_grid_polarization_coordinates_power_mesh_and_result"
    ]


def test_v27_public_rejects_time_domain_port_identity_mismatch() -> None:
    summary = _summary_v27()
    identity = summary["runs"][0][
        "time_domain_port_waveform_normalization_fft_window_grid_deembedding_smatrix_generation_identity"
    ]
    identity.update(
        {
            "waveform_time_domain_generation": "td-port-310",
            "fft_time_domain_generation": "td-port-309",
            "smatrix_time_domain_generation": "td-port-308",
            "result_port_mode_ids": ["P2:M1", "P1:M2"],
            "result_time_grid_s": [0.0, 2.0e-12, 4.0e-12],
            "result_incident_waveform": [0.0, 1.0, 0.0],
            "result_wave_normalization": "voltage-wave",
            "result_reference_impedance_ohm": [75.0, 50.0],
            "result_fft_window": "rectangular",
            "result_frequency_grid_hz": [1.0e9, 1.6e9, 2.0e9],
            "result_deembedding_offsets_m": [0.0, 0.002],
            "result_smatrix_ri": [[[0.8, -0.1]]],
            "accepted_time_result_sha256": "b" * 64,
            "accepted_smatrix_sha256": "c" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "time_domain_sparameters_use_current_waveform_normalization_fft_grid_deembedding_and_result"
    ]


def test_v27_public_rejects_huygens_identity_mismatch() -> None:
    summary = _summary_v27()
    identity = summary["runs"][0][
        "huygens_box_orientation_phase_center_frequency_mesh_near_far_transform_generation_identity"
    ]
    identity.update(
        {
            "orientation_huygens_generation": "huygens-310",
            "phase_center_huygens_generation": "huygens-309",
            "mesh_huygens_generation": "huygens-308",
            "result_box_face_ids": ["+x", "-x", "+y"],
            "result_outward_orientation_sign": [1, 1, 1],
            "result_phase_center_m": [0.1, 0.0, 0.0],
            "result_frequency_hz": 9.0e9,
            "result_near_far_transform": "direct-field-copy",
            "result_encloses_all_sources": False,
            "result_enclosing_mesh_sha256": "d" * 64,
            "accepted_near_field_sha256": "e" * 64,
            "accepted_far_field_sha256": "f" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "near_to_far_results_use_current_huygens_orientation_phase_center_frequency_mesh_and_transform"
    ]


def test_v28_public_rejects_waveguide_port_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary["runs"][0][
        "waveguide_port_mode_cutoff_normalization_reference_plane_mesh_field_result_generation_identity"
    ]
    identity.update(
        {
            "mode_port_generation": "waveguide-port-320",
            "mesh_port_generation": "waveguide-port-319",
            "result_port_id": "P2",
            "result_mode_id": "TM11",
            "result_cutoff_frequency_hz": 12.0e9,
            "result_evaluation_frequency_hz": 9.0e9,
            "result_normalization": "unit-voltage-wave",
            "result_reference_plane_m": 0.0,
            "result_port_mesh_sha256": "a" * 64,
            "result_field_eigenvector_sha256": "b" * 64,
            "accepted_result_sha256": "c" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "waveguide_ports_use_current_mode_cutoff_normalization_plane_mesh_field_and_result"
    ]


def test_v28_public_rejects_wake_impedance_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary["runs"][0][
        "wake_impedance_bunch_profile_time_grid_frequency_transform_normalization_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "bunch_wake_generation": "wake-impedance-320",
            "frequency_wake_generation": "wake-impedance-319",
            "result_bunch_profile": "rectangular",
            "result_bunch_sigma_s": 2.0e-12,
            "result_bunch_charge_c": 2.0e-9,
            "result_time_grid_s": [0.0, 2.0e-12, 5.0e-12],
            "result_wake_potential_v_c": [0.0, -1.0e12, 0.0],
            "result_fft_convention": "exp-plus-i-omega-t",
            "result_frequency_grid_hz": [0.0, 1.4e9, 3.0e9],
            "result_impedance_normalization": "transverse-ohm-per-metre",
            "result_mesh_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "wake_impedance_uses_current_bunch_time_grid_transform_frequency_normalization_mesh_and_result"
    ]


def test_v29_public_rejects_dispersive_fit_identity_mismatch():
    summary = _summary_v29()
    identity = summary["runs"][0]["dispersive_vector_fit_passivity_causality_temperature_generation_identity"]
    identity.update({
        "pole_fit_generation": "dispersive-fit-340", "temperature_fit_generation": "dispersive-fit-339",
        "result_temperature_c": 100.0, "result_frequency_grid_hz": [1.0e9, 3.0e9, 9.0e9],
        "result_poles_rad_s": [[1.0e9, 1.0e10]], "result_residues": [[-1.0e9, 0.0]],
        "result_passivity_enforced": False, "result_minimum_dissipation": -0.1,
        "result_causality_residual": 1.0e-2, "result_material_table_sha256": "e" * 64,
        "accepted_result_sha256": "f" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "dispersive_vector_fit_uses_current_stable_poles_residues_passivity_causality_temperature_and_result"
    ]


def test_v29_public_rejects_array_scan_identity_mismatch():
    summary = _summary_v29()
    identity = summary["runs"][0]["array_embedded_pattern_feed_phase_active_reflection_scan_generation_identity"]
    identity.update({
        "pattern_array_generation": "array-scan-340", "phase_array_generation": "array-scan-339",
        "result_element_order": [4, 3, 2, 1], "result_embedded_pattern_sha256": ["e" * 64],
        "result_scan_angles_deg": [30.0, 0.0, -30.0], "result_feed_phase_deg": [[0.0, 90.0]],
        "result_active_reflection_magnitude": [1.2, 0.1], "result_accepted_power_fraction": [0.2, 1.1],
        "result_array_mesh_sha256": "f" * 64, "accepted_result_sha256": "a" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "array_scan_uses_current_embedded_patterns_element_order_phases_reflection_power_mesh_and_result"
    ]


def test_v30_public_waveguide_port_mode_power_normalization_deembed_plane_smatrix_reference_mismatch():
    summary = _summary_v30(); identity = summary["runs"][0]["waveguide_port_mode_power_deembed_impedance_frequency_port_smatrix_result_identity"]
    identity.update({
        "mode_port_generation": "waveguide-port-350", "frequency_port_generation": "waveguide-port-349",
        "result_mode_ids": ["port1:TM01"], "result_power_normalization_w": [0.5, 2.0],
        "result_deembed_plane_m": [0.02, 0.08], "result_reference_impedance_ohm": [75.0, 50.0],
        "result_frequency_hz": [8.5e9, 9.5e9], "result_port_order": [2, 1],
        "result_smatrix_ri": [[[1.2, 0.0]]], "result_mesh_sha256": "7" * 64,
        "accepted_result_sha256": "8" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary); assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["waveguide_ports_use_current_modes_power_deembed_impedance_frequency_order_smatrix_mesh_and_result"]


def test_v30_public_sar_mass_average_density_voxel_frequency_field_tissue_result_mismatch():
    summary = _summary_v30(); identity = summary["runs"][0]["sar_mass_density_voxel_frequency_field_mesh_result_identity"]
    identity.update({
        "mass_sar_generation": "sar-mass-350", "voxel_sar_generation": "sar-mass-349",
        "result_averaging_mass_kg": 0.001, "result_tissue_density_kg_m3": 800.0,
        "result_voxel_ids": [103, 101], "result_voxel_mass_kg": [0.01, 0.01],
        "result_frequency_hz": 5.8e9, "result_field_normalization": "peak_field_1v_m",
        "result_sar_w_kg": -1.0, "result_mesh_sha256": "9" * 64,
        "accepted_result_sha256": "a" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary); assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["sar_uses_current_average_mass_density_voxels_frequency_field_mesh_and_result"]


def test_v31_public_wave_port_modal_power_impedance_deembed_phase_reference_balance_mismatch():
    summary = _summary_v31()
    identity = summary["runs"][0][
        "wave_port_modal_power_impedance_deembed_phase_balance_result_identity"
    ]
    identity.update(
        {
            "phase_port_generation": "wave-port-reference-360",
            "balance_port_generation": "wave-port-reference-359",
            "result_mode_ids": ["port1:TM01", "port2:TE10"],
            "result_modal_power_normalization_w": [0.5, 2.0],
            "result_reference_impedance_ohm": [75.0, 50.0],
            "result_deembed_plane_m": [0.02, 0.08],
            "result_phase_reference_rad": [3.141592653589793, 0.0],
            "result_port_mode_owner_ids": ["old:port-1", "old:port-2"],
            "result_power_balance_w": 1.2,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "wave_ports_use_current_modal_power_impedance_deembed_phase_owner_balance_and_result"
    ]


def test_v31_public_farfield_spherical_basis_handedness_polarization_phase_radiated_power_mismatch():
    summary = _summary_v31()
    identity = summary["runs"][0][
        "farfield_spherical_basis_handedness_polarization_phase_power_result_identity"
    ]
    identity.update(
        {
            "basis_farfield_generation": "farfield-basis-360",
            "power_farfield_generation": "farfield-basis-359",
            "result_spherical_basis": "e_phi_e_theta",
            "result_coordinate_handedness": "left_handed",
            "result_angular_order": "phi_major_theta_minor",
            "result_polarization_phase_convention": "exp_minus_j_phase",
            "result_theta_weights": [1.0, -1.0, 1.0],
            "result_phi_weights": [1.0],
            "integrated_radiated_power_w": 1.25,
            "accepted_farfield_owner_id": "old-project:old-monitor",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "farfields_use_current_spherical_basis_handedness_order_polarization_weights_power_owner_and_result"
    ]


def test_v32_public_dispersive_port_mode_branch_cutoff_normalization_group_delay_mismatch():
    summary = _summary_v32()
    identity = summary["runs"][0][
        "dispersive_port_mode_branch_cutoff_normalization_beta_phase_group_delay_mesh_result_identity"
    ]
    identity.update(
        {
            "branch_port_generation": "dispersive-port-370",
            "mesh_port_generation": "dispersive-port-369",
            "result_port_generation": "dispersive-port-368",
            "result_mode_id": "port1:TM01",
            "result_tracked_branch_id": "branch:TM01:backward",
            "result_cutoff_frequency_hz": 7.2e9,
            "result_modal_normalization": "unit_voltage",
            "result_propagation_constant_sign": "negative_forward",
            "result_propagation_constant_rad_per_m": [-101.0, -132.0, -158.0],
            "result_deembedded_phase_rad": [0.0, 0.2, 0.4],
            "result_group_delay_s": -3.183098861837907e-11,
            "result_mesh_sha256": "8" * 64,
            "accepted_result_owner": "waveguide/old-port",
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "dispersive_ports_use_current_mode_branch_cutoff_power_normalization_beta_phase_group_delay_mesh_and_result"
    ]


def test_v32_public_transient_farfield_time_gate_fft_window_phase_center_energy_mismatch():
    summary = _summary_v32()
    identity = summary["runs"][0][
        "transient_farfield_time_gate_fft_window_phase_center_angular_energy_monitor_result_identity"
    ]
    identity.update(
        {
            "gate_farfield_generation": "transient-farfield-370",
            "monitor_farfield_generation": "transient-farfield-369",
            "result_farfield_generation": "transient-farfield-368",
            "result_time_gate_s": [0.0, 4.0e-9],
            "result_fft_window": "rectangular",
            "result_fft_normalization": "raw_fft",
            "result_phase_center_m": [0.01, 0.0, 0.0],
            "result_theta_deg": [90.0, 45.0, 0.0],
            "result_phi_deg": [0.0, 180.0],
            "result_accepted_energy_j": 0.5,
            "result_radiated_energy_j": 1.2,
            "result_monitor_owner": "project-old:monitor-old",
            "accepted_result_owner": "farfield/old-result",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "transient_farfields_use_current_time_gate_fft_phase_center_angles_energy_monitor_owner_and_result"
    ]


def test_v33_public_eigenmode_q_energy_conductor_dielectric_radiation_inverse_sum_branch_mismatch():
    summary = _summary_v33()
    identity = summary["runs"][0][
        "eigenmode_frequency_branch_energy_conductor_dielectric_radiation_q_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "branch_generation": "eigenmode-q-380",
            "energy_generation": "eigenmode-q-379",
            "result_generation": "eigenmode-q-378",
            "result_mode_id": "mode-2",
            "result_mode_branch": "stale-crossing-branch",
            "result_frequency_hz": 9.0e9,
            "result_electric_energy_j": 0.1,
            "result_magnetic_energy_j": 0.2,
            "result_stored_energy_j": 2.0,
            "result_q_conductor": 1000.0,
            "result_q_dielectric": 2000.0,
            "result_q_radiation": -5000.0,
            "result_q_total": 10000.0,
            "result_mesh_sha256": "8" * 64,
            "accepted_mode_owner": "cavity/old-mode",
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "eigenmodes_use_current_frequency_branch_energy_q_inverse_sum_mesh_owner_and_result"
    ]


def test_v33_public_tdr_reference_plane_velocity_time_zero_impedance_reflection_causality_energy_mismatch():
    summary = _summary_v33()
    identity = summary["runs"][0][
        "tdr_reference_plane_velocity_time_zero_impedance_arrival_window_causality_energy_owner_result_identity"
    ]
    identity.update(
        {
            "reference_generation": "tdr-reference-380",
            "arrival_generation": "tdr-reference-379",
            "result_generation": "tdr-reference-378",
            "result_reference_plane_m": 0.1,
            "result_propagation_velocity_m_per_s": 3.0e8,
            "result_time_zero_s": -1.0e-9,
            "result_characteristic_impedance_ohm": 75.0,
            "result_reflection_distance_m": 0.5,
            "result_reflection_arrival_s": 1.0e-9,
            "result_time_window_s": [2.0e-9, 1.0e-9],
            "result_time_samples_s": [4.0e-9, 3.0e-9, 2.0e-9],
            "result_reflection_waveform": [0.2, 0.1, 0.0],
            "result_pre_arrival_max_abs": 0.2,
            "result_incident_energy_j": 0.5,
            "result_reflected_energy_j": 0.8,
            "result_accepted_energy_j": -0.3,
            "accepted_waveform_owner": "tdr/old-port",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "tdr_uses_current_reference_velocity_time_zero_impedance_arrival_causality_energy_owner_and_result"
    ]


def test_v34_public_sparameter_reference_plane_time_gate_passivity_causality_energy_mismatch():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
    ]
    identity.update(
        {
            "reference_generation": "sparameter-gated-390",
            "energy_generation": "sparameter-gated-389",
            "result_generation": "sparameter-gated-388",
            "result_reference_plane_shift_m": -0.02,
            "result_time_gate_window_s": [5.0e-9, 1.0e-9],
            "result_impulse_time_s": [2.0e-9, 1.0e-9, -1.0e-9],
            "result_impulse_response": [0.1, 0.2, 0.5],
            "result_pre_zero_max_abs": 0.5,
            "result_maximum_singular_values": [1.2, 1.5],
            "result_incident_energy_j": 0.5,
            "result_reflected_energy_j": 0.8,
            "result_transmitted_energy_j": 0.6,
            "result_absorbed_energy_j": -0.9,
            "result_port_impedance_ohm": [75.0],
            "result_frequency_grid_hz": [3.0e9, 2.0e9, 1.0e9],
            "accepted_sparameter_owner": "network/old",
            "accepted_sparameter_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "sparameters_use_current_reference_gate_causality_passivity_energy_ports_frequency_owner_and_result"
    ]


def test_v34_public_eigenmode_degeneracy_subspace_orthogonality_tracking_mesh_convergence_mismatch():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "subspace_generation": "degenerate-eigenmode-390",
            "tracking_generation": "degenerate-eigenmode-389",
            "result_generation": "degenerate-eigenmode-388",
            "result_mode_frequencies_hz": [9.0e9, 11.0e9],
            "result_principal_angles_rad": [0.5, 1.0],
            "result_mass_gram_real": [[1.0, 0.5], [0.2, 0.1]],
            "result_mass_gram_imag": [[0.0, 1.0], [-1.0, 0.0]],
            "result_phase_anchor_complex": [[-1.0, 1.0], [0.0, 0.0]],
            "result_tracking_subspace_ids": ["mode-a", "mode-b"],
            "result_residual_norms": [1.0, 2.0],
            "result_mesh_cell_counts": [64000, 8000, 1000],
            "result_mesh_converged_frequency_hz": [10.0e9, 9.0e9, 11.0e9],
            "result_eigenmode_mesh_sha256": "b" * 64,
            "accepted_field_owner": "eigenmode/old",
            "accepted_field_sha256": "c" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "degenerate_eigenmodes_use_current_subspace_angles_mass_orthogonality_phase_tracking_residual_mesh_owner_and_result"
    ]


def test_v34_public_self_consistent_nonpassive_sparameters_are_rejected():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "sparameter_reference_plane_time_gate_causality_passivity_energy_port_frequency_owner_result_identity"
    ]
    identity["maximum_singular_values"] = [0.8, 1.01, 0.85]
    identity["result_maximum_singular_values"] = [0.8, 1.01, 0.85]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v34_public_self_consistent_vector_labels_are_not_subspace_tracking():
    summary = _summary_v34()
    identity = summary["runs"][0][
        "eigenmode_degenerate_subspace_principal_angle_mass_orthogonality_phase_tracking_residual_mesh_owner_result_identity"
    ]
    identity["tracking_subspace_ids"] = ["mode-a", "mode-b"]
    identity["result_tracking_subspace_ids"] = ["mode-a", "mode-b"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v35_public_waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_mismatch():
    summary = _summary_v35()
    record = summary["runs"][0][
        "waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_frequency_owner_result_identity"
    ]
    record.update(
        {
            "mode_generation": "waveguide-port-410", "deembed_generation": "waveguide-port-409",
            "result_generation": "waveguide-port-408", "result_mode_name": "TM01",
            "result_normalization": "peak_field", "result_modal_power_w": [-1.0, 2.0],
            "result_mode_gram_real": [[1.0, 0.5], [0.5, 0.1]],
            "result_mode_gram_imag": [[0.0, 1.0], [-1.0, 0.0]],
            "result_impedance_definition": "lumped_50_ohm",
            "result_modal_impedance_ohm": [50.0, -50.0],
            "result_frequency_grid_hz": [12.0e9, 6.0e9],
            "result_cutoff_frequency_hz": 15.0e9,
            "result_propagation_constant_rad_m": [-1.0, -2.0],
            "result_deembedded_reference_plane_m": -0.01,
            "result_deembed_phase_rad": [1.0, 2.0],
            "result_port_mesh_sha256": "a" * 64,
            "accepted_port_owner": "waveguide-port/old",
            "accepted_port_result_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "waveguide_port_modes_use_current_power_orthogonality_impedance_deembed_cutoff_frequency_mesh_owner_and_result"
    ]


def test_v35_public_nearfar_directivity_gain_efficiency_polarization_sphere_power_mismatch():
    summary = _summary_v35()
    record = summary["runs"][0][
        "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
    ]
    record.update(
        {
            "power_generation": "nearfar-410", "polarization_generation": "nearfar-409",
            "result_generation": "nearfar-408", "result_frequency_hz": 9.0e9,
            "result_accepted_power_w": -10.0, "result_enclosing_sphere_power_w": 20.0,
            "result_radiated_power_w": 4.0, "result_radiation_efficiency": 1.5,
            "result_maximum_directivity_linear": -4.0, "result_realized_gain_linear": 9.0,
            "result_polarization_basis": "left_handed_xy", "result_copolar_definition": "unknown",
            "result_angular_quadrature_weights_sr": [1.0, -1.0],
            "result_radiation_intensity_w_sr": [100.0],
            "result_farfield_mesh_sha256": "c" * 64,
            "accepted_nearfar_owner": "nearfar/old",
            "accepted_nearfar_result_sha256": "d" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "nearfar_results_use_current_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_and_result"
    ]


def test_v35_public_rejects_self_consistent_nearfar_gain_over_efficiency():
    summary = _summary_v35()
    for row in summary["runs"]:
        record = row[
            "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
        ]
        record["realized_gain_linear"] = 5.0
        record["result_realized_gain_linear"] = 5.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v36_public_cavity_eigenmode_frequency_q_energy_normalization_orthogonality_mesh_mismatch():
    summary = _summary_v36()
    record = summary["runs"][0]["cavity_eigenmode_frequency_q_energy_orthogonality_degeneracy_mesh_owner_result_identity"]
    record.update({
        "frequency_generation": "cavity-mode-411", "q_generation": "cavity-mode-410",
        "result_generation": "cavity-mode-409", "result_mode_frequency_hz": [1.0, 1.0],
        "result_q_total": -1.0, "result_electric_energy_j": [2.0],
        "result_magnetic_energy_j": [-1.0], "result_normalization": "peak_field",
        "result_mode_gram_real": [[1.0, 1.0], [1.0, 1.0]],
        "result_degeneracy_order": [2, 1], "result_mesh_dof": [40000, 10000],
        "result_mesh_frequency_hz": [8.0e9, 7.0e9], "result_cavity_mesh_sha256": "a" * 64,
        "accepted_cavity_owner": "cavity/old", "accepted_cavity_result_sha256": "b" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["cavity_eigenmodes_use_current_frequency_q_energy_orthogonality_degeneracy_mesh_owner_and_result"]


def test_v36_public_active_sparameter_embedded_pattern_scan_impedance_power_closure_mismatch():
    summary = _summary_v36()
    record = summary["runs"][0]["active_sparameter_embedded_pattern_scan_impedance_power_frequency_mesh_owner_result_identity"]
    record.update({
        "sparameter_generation": "active-array-411", "power_generation": "active-array-410",
        "result_generation": "active-array-409", "result_frequency_hz": 9.0e9,
        "result_port_impedance_ohm": [-50.0], "result_scan_phase_rad": [math.pi],
        "result_excitation_complex": [[0.0, 0.0]], "result_s_matrix_complex": [[[2.0, 0.0]]],
        "result_active_s_complex": [[9.0, 9.0]], "result_embedded_pattern_complex": [],
        "result_incident_power_w": -2.0, "result_reflected_power_w": 3.0,
        "result_accepted_power_w": -1.0, "result_radiated_power_w": 2.0,
        "result_dissipated_power_w": -3.0, "result_array_mesh_sha256": "c" * 64,
        "accepted_array_owner": "array/old", "accepted_array_result_sha256": "d" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["active_arrays_use_current_sparameters_patterns_scan_impedance_power_frequency_mesh_owner_and_result"]


def test_v36_public_rejects_self_consistent_wrong_cavity_q_sum():
    summary = _summary_v36()
    for row in summary["runs"]:
        record = row["cavity_eigenmode_frequency_q_energy_orthogonality_degeneracy_mesh_owner_result_identity"]
        record["q_total"] = 100.0
        record["result_q_total"] = 100.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_active_array_power_imbalance():
    summary = _summary_v36()
    for row in summary["runs"]:
        record = row["active_sparameter_embedded_pattern_scan_impedance_power_frequency_mesh_owner_result_identity"]
        record["radiated_power_w"] = 1.0
        record["result_radiated_power_w"] = 1.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v37_public_waveguide_cutoff_mode_impedance_group_delay_power_orthogonality_mesh_mismatch():
    summary = _summary_v37()
    record = summary["runs"][0]["waveguide_cutoff_mode_impedance_group_delay_power_orthogonality_mesh_owner_result_identity"]
    record.update({
        "cutoff_generation": "waveguide-broadband-512", "group_delay_generation": "waveguide-broadband-511",
        "result_generation": "waveguide-broadband-510", "result_cutoff_frequency_hz": 1.0,
        "result_modal_impedance_ohm": [-1.0], "result_propagation_constant_rad_m": [0.0],
        "result_group_delay_s": [-1.0],
        "result_power_normalization_w": -1.0, "result_mode_gram_real": [[1.0, 1.0], [1.0, 1.0]],
        "result_mesh_cutoff_hz": [9.0e9, 8.0e9],
        "accepted_waveguide_owner": "waveguide/old", "accepted_waveguide_result_sha256": "8" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["broadband_waveguides_use_current_cutoff_impedance_group_delay_power_orthogonality_mesh_owner_and_result"]


def test_v37_public_emc_probe_coordinate_interpolation_time_frequency_window_parseval_owner_mismatch():
    summary = _summary_v37()
    record = summary["runs"][0]["emc_probe_coordinate_interpolation_time_fft_window_parseval_mesh_owner_result_identity"]
    record.update({
        "coordinate_generation": "emc-probe-512", "fft_generation": "emc-probe-511",
        "result_generation": "emc-probe-510", "result_coordinate_frame": "local_spherical",
        "result_interpolation_weights": [2.0, -1.0], "result_time_s": [0.0, 2.0e-9],
        "result_field_trace_v_m": [99.0], "result_fft_window": "hann",
        "result_fft_scaling": "unknown", "result_selected_frequency_hz": 1.0,
        "result_selected_fft_complex_v_m": [99.0, 99.0], "result_frequency_energy_over_n": 99.0,
        "accepted_probe_owner": "emc/old",
        "accepted_probe_result_sha256": "9" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"]["emc_probes_use_current_coordinates_interpolation_time_fft_window_parseval_mesh_owner_and_result"]


def test_v37_public_rejects_self_consistent_wrong_cutoff():
    summary = _summary_v37()
    for row in summary["runs"]:
        record = row["waveguide_cutoff_mode_impedance_group_delay_power_orthogonality_mesh_owner_result_identity"]
        record["cutoff_frequency_hz"] *= 0.5
        record["result_cutoff_frequency_hz"] = record["cutoff_frequency_hz"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_parseval_error():
    summary = _summary_v37()
    for row in summary["runs"]:
        record = row["emc_probe_coordinate_interpolation_time_fft_window_parseval_mesh_owner_result_identity"]
        record["time_energy"] = 3.0
        record["result_time_energy"] = 3.0
        record["frequency_energy_over_n"] = 3.0
        record["result_frequency_energy_over_n"] = 3.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v38_public_microstrip_quasitem_impedance_effective_permittivity_delay_loss_sparameter_mismatch():
    summary = _summary_v38()
    row = summary["runs"][0][
        "microstrip_quasitem_geometry_permittivity_impedance_delay_loss_sparameter_frequency_mesh_owner_result_identity"
    ]
    row.update({
        "impedance_generation": "microstrip-quasitem-613",
        "loss_generation": "microstrip-quasitem-612",
        "result_generation": "microstrip-quasitem-611",
        "result_effective_permittivity": -1.0,
        "result_characteristic_impedance_ohm": -50.0,
        "result_propagation_delay_s": -1.0,
        "result_conductor_loss_fraction": -0.2,
        "result_dielectric_loss_fraction": 2.0,
        "result_s11_magnitude": 2.0,
        "result_s21_magnitude": 2.0,
        "result_mesh_impedance_ohm": [10.0, 100.0],
        "accepted_microstrip_owner": "pcb/old",
        "accepted_microstrip_result_sha256": "a" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "microstrip_results_use_current_quasitem_impedance_permittivity_delay_loss_passivity_mesh_owner_and_result"
    ]


def test_v38_public_shielding_aperture_incident_transmitted_field_power_se_frequency_mesh_mismatch():
    summary = _summary_v38()
    row = summary["runs"][0][
        "shielding_aperture_orientation_polarization_field_power_se_frequency_probe_mesh_owner_result_identity"
    ]
    row.update({
        "aperture_generation": "shield-aperture-613",
        "se_generation": "shield-aperture-612",
        "result_generation": "shield-aperture-611",
        "result_aperture_plane": "unknown",
        "result_incident_polarization": "left_circular",
        "result_frequency_hz": [3.0e9, 1.0e9],
        "result_incident_field_v_m": [-1.0],
        "result_transmitted_field_v_m": [2.0],
        "result_transmitted_power_density_normalized": [4.0],
        "result_shielding_effectiveness_field_db": [-6.0],
        "result_shielding_effectiveness_power_db": [6.0],
        "result_probe_frame": "local_spherical",
        "result_mesh_selected_se_db": [40.0, 20.0],
        "accepted_shield_owner": "emc/old",
        "accepted_shield_result_sha256": "b" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "shielding_results_use_current_aperture_polarization_field_power_se_frequency_probe_mesh_owner_and_result"
    ]


def test_v38_public_rejects_self_consistent_wrong_microstrip_impedance():
    summary = _summary_v38()
    for row in summary["runs"]:
        record = row[
            "microstrip_quasitem_geometry_permittivity_impedance_delay_loss_sparameter_frequency_mesh_owner_result_identity"
        ]
        record["characteristic_impedance_ohm"] *= 2.0
        record["result_characteristic_impedance_ohm"] = record["characteristic_impedance_ohm"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_field_power_se_sign_error():
    summary = _summary_v38()
    for row in summary["runs"]:
        record = row[
            "shielding_aperture_orientation_polarization_field_power_se_frequency_probe_mesh_owner_result_identity"
        ]
        wrong = [-item for item in record["shielding_effectiveness_power_db"]]
        record["shielding_effectiveness_power_db"] = wrong
        record["result_shielding_effectiveness_power_db"] = wrong
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v39_public_waveguide_mode_cutoff_impedance_power_orthogonality_deembed_mismatch() -> None:
    summary = _summary_v39()
    row = summary["runs"][0][_WAVEGUIDE_KEY]
    row.update(
        {
            "cutoff_generation": "waveguide-mode-714",
            "deembed_generation": "waveguide-mode-713",
            "result_generation": "waveguide-mode-712",
            "result_cutoff_frequency_hz": -1.0,
            "result_modal_impedance_ohm": -50.0,
            "result_normalized_forward_power_w": -1.0,
            "result_mode_overlap_real": [[1.0, 1.0], [1.0, 1.0]],
            "result_propagation_constant_rad_m": -1.0,
            "result_deembed_distance_m": -0.01,
            "result_deembedded_s21_phase_rad": 9.0,
            "accepted_waveguide_mode_owner": "waveguide/old",
            "accepted_waveguide_mode_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "waveguide_modes_use_current_cutoff_impedance_power_orthogonality_propagation_deembed_mesh_owner_and_result"
    ]


def test_v39_public_antenna_nearfar_directivity_gain_efficiency_polarization_power_mismatch() -> None:
    summary = _summary_v39()
    row = summary["runs"][0][_ANTENNA_KEY]
    row.update(
        {
            "farfield_generation": "antenna-nearfar-714",
            "power_generation": "antenna-nearfar-713",
            "result_generation": "antenna-nearfar-712",
            "result_nearfield_surface_closed": False,
            "result_near_to_far_transform": "unknown",
            "result_directivity_linear": -1.0,
            "result_realized_gain_linear": 9.0,
            "result_radiation_efficiency": 2.0,
            "result_polarization_basis": "left_handed_local",
            "result_accepted_power_w": -1.0,
            "result_radiated_power_w": 2.0,
            "result_power_balance_residual_w": 3.0,
            "accepted_antenna_owner": "antenna/old",
            "accepted_antenna_result_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "antenna_nearfar_uses_current_transform_directivity_gain_efficiency_polarization_power_mesh_owner_and_result"
    ]


def test_v39_public_rejects_self_consistent_wrong_waveguide_impedance() -> None:
    summary = _summary_v39()
    for run in summary["runs"]:
        row = run[_WAVEGUIDE_KEY]
        row["modal_impedance_ohm"] *= 2.0
        row["result_modal_impedance_ohm"] = row["modal_impedance_ohm"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_unknown_nearfar_transform() -> None:
    summary = _summary_v39()
    for run in summary["runs"]:
        row = run[_ANTENNA_KEY]
        row["near_to_far_transform"] = "unknown"
        row["result_near_to_far_transform"] = "unknown"
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_pcb_via_stub_impedance_resonance_sparameter_current_loss_power_mismatch() -> None:
    summary = _summary_v40()
    summary["runs"][0][_VIA].update({"geometry_generation": "pcb-via-723", "result_characteristic_impedance_ohm": -1.0, "result_s11_complex": [2.0, 0.0], "result_accepted_power_w": -1.0, "accepted_via_owner": "pcb/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_cavity_degenerate_modes_frequency_orthogonality_energy_symmetry_q_mismatch() -> None:
    summary = _summary_v40()
    summary["runs"][0][_CAVITY].update({"frequency_generation": "cavity-degenerate-723", "result_mode_frequencies_hz": [1.0, 2.0], "result_modal_overlap_matrix": [[1.0, 1.0], [1.0, 1.0]], "result_quality_factors": [-1.0, -1.0], "accepted_cavity_owner": "cavity/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_via_power_gap() -> None:
    summary = _summary_v40()
    for run in summary["runs"]:
        row = run[_VIA]
        row["dielectric_loss_w"] = row["result_dielectric_loss_w"] = 0.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_cavity_energy_imbalance() -> None:
    summary = _summary_v40()
    for run in summary["runs"]:
        row = run[_CAVITY]
        row["magnetic_energy_j"] = row["result_magnetic_energy_j"] = [0.4, 0.4]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_waveguide_cutoff_propagation_impedance_mode_orthogonality_sparameter_power_mismatch() -> None:
    summary = _summary_v41()
    summary["runs"][0][_WAVEGUIDE_v41].update({"cutoff_generation": "waveguide-te10-723", "result_cutoff_frequency_hz": -1.0, "result_modal_overlap_matrix": [[1.0, 1.0], [1.0, 1.0]], "accepted_waveguide_owner": "waveguide/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_antenna_farfield_directivity_gain_efficiency_power_polarization_mismatch() -> None:
    summary = _summary_v41()
    summary["runs"][0][_ANTENNA].update({"farfield_generation": "antenna-farfield-723", "result_radiation_efficiency": 2.0, "result_gain_linear": 20.0, "result_polarization_basis": "unknown", "accepted_antenna_owner": "antenna/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_waveguide_impedance() -> None:
    summary = _summary_v41()
    for run in summary["runs"]:
        row = run[_WAVEGUIDE_v41]
        row["guide_impedance_ohm"] = row["result_guide_impedance_ohm"] = 50.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_antenna_gain() -> None:
    summary = _summary_v41()
    for run in summary["runs"]:
        row = run[_ANTENNA]
        row["gain_linear"] = row["result_gain_linear"] = row["directivity_linear"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v43_public_positive_emc_and_differential_connector() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v43())["status"] == "ok"


def test_v43_public_emc_rejects_power_closure_and_digest_mismatch() -> None:
    summary = _summary_v43()
    row = summary["runs"][0][_EMC_v43]
    row["result_values"]["transmitted_power_w"] = 20.0
    row["result_values"]["power_closure_residual_w"] = -19.0
    row["accepted_emc_result_sha256"] = "d" * 64
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v43_public_connector_rejects_mode_order_and_passivity_mismatch() -> None:
    summary = _summary_v43()
    row = summary["runs"][0][_CONNECTOR]
    row["result_values"]["port_order"] = ["P1-", "P1+", "P2+", "P2-"]
    row["result_values"]["passivity_margin"] = -0.2
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v43_public_rejects_self_consistent_wrong_shielding_effectiveness() -> None:
    summary = _summary_v43()
    for run in summary["runs"]:
        row = run[_EMC_v43]
        row["values"]["shielding_effectiveness_db"] = 10.0
        row["result_values"]["shielding_effectiveness_db"] = 10.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v44_public_waveguide_and_emc_identity_positive() -> None:
    assert validate_public_identity(_payload_v44()) == {"waveguide_v44_modal_identity": True, "waveguide_v44_emc_probe_identity": True}


def test_v44_public_identity_rejects_cutoff_and_parseval_mutations() -> None:
    payload = _payload_v44()
    payload["runs"][0][_WAVEGUIDE_v44]["result_cutoff_frequency_hz"] = 12.0e9
    payload["runs"][0][_EMC_v44]["result_parseval_frequency_energy_j"] = 2.0
    result = validate_public_identity(payload)
    assert result["waveguide_v44_modal_identity"] is False
    assert result["waveguide_v44_emc_probe_identity"] is False


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_v51_identity(_payload_v51()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    payload = _payload_v51()
    payload["runs"][0][S_PARAMETER].update({"result_wave_basis": "pseudo_waves", "result_port_owner": "port:stale"})
    payload["runs"][0][GROUP_DELAY].update({"result_derivative_definition": "dphi_deg_df", "result_trace_owner": "trace:stale"})
    assert not all(validate_public_v51_identity(payload).values())


def test_v51_self_consistent_wrong_network_semantics_are_rejected() -> None:
    payload = _payload_v51()
    for run in payload["runs"]:
        run[S_PARAMETER]["wave_basis"] = run[S_PARAMETER]["result_wave_basis"] = "pseudo_waves"
        run[GROUP_DELAY]["smoothing_window_points"] = run[GROUP_DELAY]["result_smoothing_window_points"] = 4
    assert not all(validate_public_v51_identity(payload).values())


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_v52_identity(_payload_v52()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    payload = _payload_v52()
    payload["runs"][0][ENERGY_BALANCE]["result_absorbed_power_w"] = 80.0
    payload["runs"][0][EIGENMODE_Q]["result_q_factor"] = 1.0
    assert not all(validate_public_v52_identity(payload).values())


def test_v52_self_consistent_wrong_physics_are_rejected() -> None:
    payload = _payload_v52()
    for run in payload["runs"]:
        run[ENERGY_BALANCE]["absorbed_power_w"] = run[ENERGY_BALANCE]["result_absorbed_power_w"] = 80.0
        run[EIGENMODE_Q]["q_factor"] = run[EIGENMODE_Q]["result_q_factor"] = 1.0
    assert not all(validate_public_v52_identity(payload).values())


def test_v53_positive_public_artifacts_are_accepted(): assert all(validate_public_v53_identity(_payload_v53()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    payload = deepcopy(_payload_v53()); payload["runs"][0][WAVEGUIDE]["result_mode_id"] = "TE20"; payload["runs"][0][FARFIELD]["result_polarization_basis"] = "spherical"; assert not all(validate_public_v53_identity(payload).values())


def test_v53_self_consistent_invalid_semantics_are_rejected():
    payload = deepcopy(_payload_v53()); payload["runs"][0][WAVEGUIDE]["cutoff_frequency_hz"] = payload["runs"][0][WAVEGUIDE]["result_cutoff_frequency_hz"] = 1.0; payload["runs"][0][FARFIELD]["realized_gain_dbi"] = payload["runs"][0][FARFIELD]["result_realized_gain_dbi"] = [[1.0]]; assert not all(validate_public_v53_identity(payload).values())


def test_v54_positive_public_artifacts_are_accepted():
    assert all(validate_public_v54_identity(_payload_v54()).values())


def test_v54_frozen_counterfactuals_are_rejected():
    payload = deepcopy(_payload_v54())
    payload["runs"][0][CUTOFF]["result_mode_impedance_ohm"] = 50.0
    payload["runs"][0][SAR]["result_field_solution_sha256"] = "9" * 64
    assert not all(validate_public_v54_identity(payload).values())


def test_v54_self_consistent_nonphysical_artifacts_are_rejected():
    payload = deepcopy(_payload_v54())
    payload["runs"][0][CUTOFF]["modal_power_w"] = payload["runs"][0][CUTOFF]["result_modal_power_w"] = 2.0
    payload["runs"][0][SAR]["voxel_support_m3"] = payload["runs"][0][SAR]["result_voxel_support_m3"] = [1.0e-6]
    assert not all(validate_public_v54_identity(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload_v54())
    payload["runs"][0][CUTOFF]["sample_frequency_hz"] = [10.0e9]
    payload["runs"][0][SAR]["electric_field_rms_v_m"] = [[12.0]]
    assert not all(validate_public_v54_identity(payload).values())


def test_v55_positive_public_artifacts_are_accepted(): assert all(validate_public_v55_identity(_payload_v55()).values())


def test_v55_frozen_counterfactuals_are_rejected():
    p=deepcopy(_payload_v55());p["runs"][0][RESONATOR]["result_monitor_owner"]="monitor:stale";p["runs"][0][ANTENNA]["result_farfield_owner"]="farfield:stale";assert not all(validate_public_v55_identity(p).values())


def test_v55_self_consistent_nonphysical_artifacts_are_rejected():
    p=deepcopy(_payload_v55());p["runs"][0][RESONATOR]["stored_energy_j"]=p["runs"][0][RESONATOR]["result_stored_energy_j"]=-1.0;p["runs"][0][ANTENNA]["radiation_efficiency"]=p["runs"][0][ANTENNA]["result_radiation_efficiency"]=1.5;assert not all(validate_public_v55_identity(p).values())


def test_v55_malformed_values_reject_without_raising():
    p=deepcopy(_payload_v55());p["runs"][0][RESONATOR]["loaded_q"]=[1200.0];p["runs"][0][ANTENNA]["accepted_power_w"]=[1.0];assert not all(validate_public_v55_identity(p).values())


def test_v55_numeric_digests_are_rejected():
    p = deepcopy(_payload_v55())
    numeric_digest = int("1" * 64)
    for identity in (RESONATOR, ANTENNA):
        p["runs"][0][identity]["result_sha256"] = numeric_digest
        p["runs"][0][identity]["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_v55_identity(p).values())
