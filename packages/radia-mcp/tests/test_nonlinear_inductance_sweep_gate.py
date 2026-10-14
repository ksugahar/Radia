import copy
import json
import math

import pytest

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from radia_mcp.radia_ngsolve.server import nonlinear_inductance_sweep_gate as mcp_gate


def _summary():
    runs = []
    levels = (
        (0.5, 2.0, 2.5),
        (2.0, 3.0, 3.6),
        (10.0, 2.0, 0.8),
        (50.0, 0.8, 0.2),
    )
    for current, apparent, incremental in levels:
        for replay in (1, 2):
            matrix = [[apparent, -0.5 * apparent], [-0.5 * apparent, apparent]]
            tangent = [
                [incremental, -0.5 * incremental],
                [-0.5 * incremental, incremental],
            ]
            flux = [apparent * current, -0.5 * apparent * current]
            energy = 0.4 * current * flux[0]
            runs.append(
                {
                    "current_A_requested": current,
                    "replay": replay,
                    "apparent_inductance_H": matrix,
                    "incremental_inductance_H": tangent,
                    "current_A": [current, 0.0],
                    "flux_linkage_Vs": flux,
                    "energy_J": energy,
                    "coenergy_J": current * flux[0] - energy,
                    "final_nonlinear_residual_log10": -7.0,
                    "result_metadata": {
                        "energy": {"run_id": 0},
                        "coenergy": {"run_id": 0},
                        "residual": {"run_id": 0},
                    },
                }
            )
    return {"runs": runs}


def _with_v8_generations(summary):
    for row in summary["runs"]:
        row.update(
            {
                "solve_sweep_generation": "nonlinear-sweep-42",
                "apparent_matrix_sweep_generation": "nonlinear-sweep-42",
                "incremental_matrix_sweep_generation": "nonlinear-sweep-42",
            }
        )
    summary["energy_history_segments"] = [
        {
            "segment_generation": "segment-1",
            "start_run_index": 0,
            "end_run_index": 3,
            "coenergy_offset_in_J": 0.0,
            "coenergy_offset_out_J": 1.25,
        },
        {
            "segment_generation": "segment-2",
            "start_run_index": 4,
            "end_run_index": 7,
            "coenergy_offset_in_J": 1.25,
            "coenergy_offset_out_J": 4.0,
        },
    ]
    return summary


def _with_v9_bindings(summary):
    summary = _with_v8_generations(summary)
    for row in summary["runs"]:
        row["matrix_port_order"] = {
            "run_current": ["primary", "secondary"],
            "flux_linkage": ["primary", "secondary"],
            "apparent_rows": ["primary", "secondary"],
            "apparent_columns": ["primary", "secondary"],
            "incremental_rows": ["primary", "secondary"],
            "incremental_columns": ["primary", "secondary"],
        }
        row["energy_loss_basis"] = {
            "stored_energy_unit": "J",
            "coenergy_unit": "J",
            "loss_series_unit": "J",
            "loss_series_scale_to_J": 1.0,
            "shared_accumulation_basis": "J",
        }
    return summary


def _with_v10_bindings(summary):
    summary = _with_v9_bindings(summary)
    for row in summary["runs"]:
        row["sparameter_reference_impedance"] = {
            "port_order": ["primary", "secondary"],
            "solver_reference_impedance_ohm_complex": [
                [50.0, 0.0],
                [50.0, 0.0],
            ],
            "comparison_reference_impedance_ohm_complex": [
                [50.0, 0.0],
                [50.0, 0.0],
            ],
            "renormalization_applied": False,
            "reference_impedance_generation": "zref-42",
        }
        row["frequency_axis_identity"] = {
            "numeric_axis_unit": "GHz",
            "metadata_axis_unit": "GHz",
            "scale_to_hz": 1.0e9,
            "normalized_axis_unit": "Hz",
            "normalization_applied_once": True,
            "frequency_axis_generation": "frequency-grid-42",
        }
    return summary


def _with_v11_bindings(summary):
    summary = _with_v10_bindings(summary)
    for row in summary["runs"]:
        row["sparameter_reference_plane_identity"] = {
            "port_order": ["primary", "secondary"],
            "original_reference_plane_ids": [
                "rp-primary-42",
                "rp-secondary-42",
            ],
            "target_reference_plane_ids": ["rp-primary-43", "rp-secondary-43"],
            "compared_port_mode_reference_plane_ids": [
                "rp-primary-43",
                "rp-secondary-43",
            ],
            "deembedding_applied": True,
            "deembedding_generation": "deembed-43",
            "sparameter_generation": "deembed-43",
        }
        row["energy_q_frequency_identity"] = {
            "q_frequency_hz": 1.0e9,
            "stored_energy_frequency_hz": 1.0e9,
            "loss_frequency_hz": 1.0e9,
            "adaptive_sample_id": "adaptive-f-43",
            "stored_energy_sample_id": "adaptive-f-43",
            "loss_sample_id": "adaptive-f-43",
        }
    return summary


def _with_v12_bindings(summary):
    summary = _with_v11_bindings(summary)
    for row in summary["runs"]:
        row["mixed_mode_sparameter_basis_identity"] = {
            "single_ended_port_order": ["P1+", "P1-", "P2+", "P2-"],
            "sparameter_port_order": ["P1+", "P1-", "P2+", "P2-"],
            "basis_matrix_port_order": ["P1+", "P1-", "P2+", "P2-"],
            "port_order_generation": "port-order-14",
            "basis_matrix_port_order_generation": "port-order-14",
            "basis_matrix_sha256": "6" * 64,
        }
        row["farfield_realized_gain_power_frequency_identity"] = {
            "farfield_frequency_hz": 2.45e9,
            "accepted_power_frequency_hz": 2.45e9,
            "farfield_adaptive_sample_id": "adaptive-frequency-14",
            "accepted_power_adaptive_sample_id": "adaptive-frequency-14",
            "farfield_result_generation": "farfield-result-14",
            "accepted_power_result_generation": "farfield-result-14",
        }
    return summary


def _with_v13_bindings(summary):
    summary = _with_v12_bindings(summary)
    for row in summary["runs"]:
        row["field_monitor_interpolation_mesh_pass_identity"] = {
            "active_adaptive_pass_id": "adaptive-pass-15",
            "field_monitor_adaptive_pass_id": "adaptive-pass-15",
            "interpolation_weight_adaptive_pass_id": "adaptive-pass-15",
            "integral_adaptive_pass_id": "adaptive-pass-15",
            "active_mesh_generation": "adaptive-mesh-15",
            "field_monitor_mesh_generation": "adaptive-mesh-15",
            "interpolation_weight_mesh_generation": "adaptive-mesh-15",
            "interpolation_weight_sha256": "9" * 64,
            "integral_weight_sha256": "9" * 64,
        }
        row["port_deembed_reference_plane_unit_identity"] = {
            "model_length_unit": "mm",
            "model_length_scale_to_m": 0.001,
            "reference_plane_offset_numeric": 2.5,
            "reference_plane_offset_unit": "mm",
            "reference_plane_offset_scale_to_m": 0.001,
            "result_reference_plane_offset_numeric": 2.5,
            "result_reference_plane_offset_unit": "mm",
            "result_reference_plane_offset_scale_to_m": 0.001,
            "port_setup_generation": "port-setup-15",
            "sparameter_result_generation": "port-setup-15",
        }
    return summary


def _with_v14_bindings(summary):
    summary = _with_v13_bindings(summary)
    for row in summary["runs"]:
        row["sparameter_reference_impedance_renormalization_identity"] = {
            "port_setup_generation": "port-setup-16",
            "sparameter_result_generation": "port-setup-16",
            "reference_impedance_basis": "complex_ohm",
            "result_reference_impedance_basis": "complex_ohm",
            "reference_impedance_real_ohm": 50.0,
            "reference_impedance_imag_ohm": 0.0,
            "result_reference_impedance_real_ohm": 50.0,
            "result_reference_impedance_imag_ohm": 0.0,
            "renormalization_applied": True,
            "renormalization_generation": "renormalization-16",
            "result_renormalization_generation": "renormalization-16",
            "sparameter_array_sha256": "c" * 64,
            "renormalized_array_sha256": "c" * 64,
        }
        row["farfield_ludwig_polarization_basis_identity"] = {
            "farfield_result_generation": "farfield-result-16",
            "co_cross_result_generation": "farfield-result-16",
            "coordinate_frame_id": "global-spherical",
            "co_cross_coordinate_frame_id": "global-spherical",
            "ludwig_basis_definition": "ludwig_3",
            "co_cross_ludwig_basis_definition": "ludwig_3",
            "polarization_reference_axis": "global-x",
            "co_cross_polarization_reference_axis": "global-x",
            "polarization_basis_sha256": "e" * 64,
            "co_cross_polarization_basis_sha256": "e" * 64,
        }
        row["sparameter_power_wave_normalization_identity"] = {
            "port_setup_generation": "port-setup-17",
            "modal_result_port_setup_generation": "port-setup-17",
            "sparameter_result_port_setup_generation": "port-setup-17",
            "incident_modal_amplitude_normalization": "power_wave",
            "reflected_modal_amplitude_normalization": "power_wave",
            "sparameter_normalization": "power_wave",
            "reference_impedance_basis": "complex_ohm",
            "reference_impedance_real_ohm": 50.0,
            "reference_impedance_imag_ohm": 5.0,
            "modal_normalization_generation": "modal-normalization-17",
            "sparameter_normalization_generation": "modal-normalization-17",
            "modal_normalization_sha256": "1" * 64,
            "sparameter_normalization_sha256": "1" * 64,
        }
        row["time_domain_fft_window_coherent_gain_identity"] = {
            "time_trace_generation": "time-trace-17",
            "fft_input_trace_generation": "time-trace-17",
            "window_generation": "fft-window-17",
            "coherent_gain_window_generation": "fft-window-17",
            "fft_result_window_generation": "fft-window-17",
            "window_definition": "periodic_hann",
            "sample_count": 1024,
            "coherent_gain": 0.5,
            "fft_coherent_gain_correction": 2.0,
            "coherent_gain_application_count": 1,
            "window_coefficients_sha256": "2" * 64,
            "fft_window_coefficients_sha256": "2" * 64,
        }
        row["sparameter_complex_impedance_renormalization_identity"] = {
            "port_calibration_generation": "port-calibration-18",
            "sparameter_result_port_calibration_generation": "port-calibration-18",
            "renormalization_port_calibration_generation": "port-calibration-18",
            "source_reference_impedance_basis": "complex_ohm",
            "renormalization_reference_impedance_basis": "complex_ohm",
            "source_reference_impedance_real_ohm": 50.0,
            "source_reference_impedance_imag_ohm": 5.0,
            "target_reference_impedance_real_ohm": 75.0,
            "target_reference_impedance_imag_ohm": -2.0,
            "renormalization_transform_generation": "renormalization-18",
            "result_renormalization_transform_generation": "renormalization-18",
            "renormalization_transform_sha256": "5" * 64,
            "result_renormalization_transform_sha256": "5" * 64,
        }
        row["farfield_polarization_basis_transform_identity"] = {
            "farfield_solve_generation": "farfield-solve-18",
            "farfield_result_generation": "farfield-solve-18",
            "source_polarization_basis": "spherical_theta_phi",
            "comparison_polarization_basis": "ludwig_3",
            "basis_transform_applied": True,
            "basis_transform_generation": "polarization-transform-18",
            "comparison_transform_generation": "polarization-transform-18",
            "angular_coordinate_frame": "global_spherical",
            "comparison_angular_coordinate_frame": "global_spherical",
            "basis_transform_sha256": "6" * 64,
            "comparison_basis_transform_sha256": "6" * 64,
        }
        row["mixed_mode_sparameter_port_pair_order_identity"] = {
            "port_calibration_generation": "port-calibration-19",
            "single_ended_result_port_calibration_generation": "port-calibration-19",
            "mixed_mode_pairing_port_calibration_generation": "port-calibration-19",
            "single_ended_port_order": ["P1+", "P1-", "P2+", "P2-"],
            "mixed_mode_pair_order": [["P1+", "P1-"], ["P2+", "P2-"]],
            "mixed_mode_pair_polarity": [1, 1],
            "transform_pair_order": [["P1+", "P1-"], ["P2+", "P2-"]],
            "transform_pair_polarity": [1, 1],
            "pairing_sha256": "1" * 64,
            "transform_pairing_sha256": "1" * 64,
            "basis_matrix_sha256": "6" * 64,
            "transform_basis_matrix_sha256": "6" * 64,
        }
        row["nearfield_farfield_phase_center_coordinate_frame_identity"] = {
            "nearfield_result_generation": "nearfield-19",
            "farfield_transform_nearfield_generation": "nearfield-19",
            "phase_center_coordinate_frame_generation": "phase-frame-19",
            "farfield_phase_center_frame_generation": "phase-frame-19",
            "farfield_result_generation": "farfield-solve-18",
            "phase_center_farfield_result_generation": "farfield-solve-18",
            "phase_center_coordinate_frame": "global_cartesian",
            "farfield_phase_center_coordinate_frame": "global_cartesian",
            "phase_center_coordinate_unit": "m",
            "farfield_phase_center_coordinate_unit": "m",
            "phase_center_coordinates_m": [0.01, -0.02, 0.03],
            "farfield_phase_center_coordinates_m": [0.01, -0.02, 0.03],
            "phase_center_sha256": "2" * 64,
            "farfield_phase_center_sha256": "2" * 64,
        }
        row[
            "sparameter_deembed_reference_plane_per_port_generation_identity"
        ] = {
            "sparameter_generation": "sparameter-20",
            "deembedded_result_sparameter_generation": "sparameter-20",
            "port_generation": "port-20",
            "reference_plane_port_generation": "port-20",
            "deembedded_result_port_generation": "port-20",
            "port_ids": ["P1", "P2"],
            "reference_plane_port_ids": ["P1", "P2"],
            "reference_plane_offsets_m": [0.001, 0.002],
            "applied_reference_plane_offsets_m": [0.001, 0.002],
            "reference_plane_map_sha256": "5" * 64,
            "deembedded_reference_plane_map_sha256": "5" * 64,
        }
        row["time_domain_port_signal_gate_window_generation_identity"] = {
            "signal_generation": "port-signal-20",
            "gate_window_signal_generation": "port-signal-20",
            "transform_signal_generation": "port-signal-20",
            "gate_generation": "gate-20",
            "transform_gate_generation": "gate-20",
            "signal_sample_start": 0,
            "signal_sample_end": 1023,
            "gate_window_start_sample": 120,
            "gate_window_end_sample": 880,
            "transform_gate_window": [120, 880],
            "normalization_basis": "incident_wave_peak",
            "transform_normalization_basis": "incident_wave_peak",
            "gate_window_sha256": "6" * 64,
            "transform_gate_window_sha256": "6" * 64,
        }
        row[
            "sparameter_port_renormalization_reference_impedance_generation_identity"
        ] = {
            "sparameter_generation": "sparameter-21",
            "renormalized_result_sparameter_generation": "sparameter-21",
            "port_calibration_generation": "port-calibration-21",
            "reference_impedance_port_calibration_generation": "port-calibration-21",
            "renormalization_port_calibration_generation": "port-calibration-21",
            "renormalized_result_port_calibration_generation": "port-calibration-21",
            "port_ids": ["P1", "P2"],
            "reference_impedance_port_ids": ["P1", "P2"],
            "reference_impedances_ohm": [50.0, 75.0],
            "applied_reference_impedances_ohm": [50.0, 75.0],
            "frequency_grid_sha256": "a" * 64,
            "renormalized_frequency_grid_sha256": "a" * 64,
            "reference_impedance_map_sha256": "b" * 64,
            "renormalized_reference_impedance_map_sha256": "b" * 64,
        }
        row[
            "realized_gain_accepted_power_port_excitation_generation_identity"
        ] = {
            "realized_gain_generation": "realized-gain-21",
            "result_realized_gain_generation": "realized-gain-21",
            "excitation_generation": "excitation-21",
            "accepted_power_excitation_generation": "excitation-21",
            "port_coefficient_excitation_generation": "excitation-21",
            "realized_gain_excitation_generation": "excitation-21",
            "port_ids": ["P1", "P2"],
            "accepted_power_port_ids": ["P1", "P2"],
            "excitation_coefficient_port_ids": ["P1", "P2"],
            "accepted_power_w": [0.6, 0.4],
            "realized_gain_accepted_power_w": [0.6, 0.4],
            "excitation_coefficients_re_im": [[1.0, 0.0], [0.0, 0.0]],
            "realized_gain_excitation_coefficients_re_im": [
                [1.0, 0.0],
                [0.0, 0.0],
            ],
            "accepted_power_unit": "W",
            "gain_unit": "dBi",
            "excitation_table_sha256": "c" * 64,
            "realized_gain_excitation_table_sha256": "c" * 64,
        }
        row[
            "farfield_polarization_basis_phase_center_coordinate_generation_identity"
        ] = {
            "farfield_generation": "farfield-22",
            "result_farfield_generation": "farfield-22",
            "monitor_coordinate_generation": "farfield-monitor-22",
            "theta_basis_monitor_coordinate_generation": "farfield-monitor-22",
            "phi_basis_monitor_coordinate_generation": "farfield-monitor-22",
            "phase_center_monitor_coordinate_generation": "farfield-monitor-22",
            "result_monitor_coordinate_generation": "farfield-monitor-22",
            "sample_ids": [1, 2, 3],
            "result_sample_ids": [1, 2, 3],
            "theta_basis_sha256": ["1" * 64, "2" * 64, "3" * 64],
            "result_theta_basis_sha256": ["1" * 64, "2" * 64, "3" * 64],
            "phi_basis_sha256": ["4" * 64, "5" * 64, "6" * 64],
            "result_phi_basis_sha256": ["4" * 64, "5" * 64, "6" * 64],
            "phase_center_m": [0.0, 0.0, 0.0],
            "result_phase_center_m": [0.0, 0.0, 0.0],
            "polarization_table_sha256": "7" * 64,
            "result_polarization_table_sha256": "7" * 64,
        }
        row["broadband_energy_q_port_loss_normalization_generation_identity"] = {
            "analysis_generation": "broadband-analysis-22",
            "energy_analysis_generation": "broadband-analysis-22",
            "port_power_analysis_generation": "broadband-analysis-22",
            "loss_analysis_generation": "broadband-analysis-22",
            "q_result_analysis_generation": "broadband-analysis-22",
            "frequency_grid_generation": "frequency-grid-22",
            "energy_frequency_grid_generation": "frequency-grid-22",
            "port_power_frequency_grid_generation": "frequency-grid-22",
            "loss_frequency_grid_generation": "frequency-grid-22",
            "q_frequency_grid_generation": "frequency-grid-22",
            "excitation_generation": "excitation-22",
            "energy_excitation_generation": "excitation-22",
            "port_power_excitation_generation": "excitation-22",
            "loss_excitation_generation": "excitation-22",
            "frequencies_hz": [1.0e9, 1.5e9, 2.0e9],
            "q_frequencies_hz": [1.0e9, 1.5e9, 2.0e9],
            "stored_energy_j": [1.0e-6, 1.1e-6, 1.2e-6],
            "q_stored_energy_j": [1.0e-6, 1.1e-6, 1.2e-6],
            "port_power_w": [1.0, 1.0, 1.0],
            "q_port_power_w": [1.0, 1.0, 1.0],
            "loss_power_w": [0.1, 0.11, 0.12],
            "q_loss_power_w": [0.1, 0.11, 0.12],
            "energy_q_input_sha256": "8" * 64,
            "result_energy_q_input_sha256": "8" * 64,
        }
    return summary


def test_nonlinear_inductance_sweep_accepts_crossover_duality_and_replay():
    result = nonlinear_inductance_sweep_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["differential_to_apparent_primary_ratios"][0] > 1.0
    assert result["differential_to_apparent_primary_ratios"][-1] < 1.0


def test_nonlinear_inductance_sweep_rejects_wrong_global_order_and_duality():
    bad = copy.deepcopy(_summary())
    for row in bad["runs"]:
        row["incremental_inductance_H"] = row["apparent_inductance_H"]
    bad["runs"][0]["coenergy_J"] *= 0.5
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["initial_magnetization_rise_is_observed"] is False
    assert result["checks"]["differential_to_apparent_crossover_is_observed"] is False
    assert result["checks"]["all_run_identities_and_matrices_close"] is False


def test_nonlinear_inductance_sweep_rejects_asymmetric_incremental_matrix():
    bad = copy.deepcopy(_summary())
    bad["runs"][3]["incremental_inductance_H"][0][1] *= 0.5
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False


def test_nonlinear_inductance_sweep_mcp_dispatches_and_rejects_bad_shape():
    result = json.loads(mcp_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
    invalid = json.loads(mcp_gate('{"runs": []}'))
    assert invalid["status"] == "invalid_input"


def test_nonlinear_inductance_sweep_rejects_indefinite_tangent_matrix():
    bad = copy.deepcopy(_summary())
    bad["runs"][0]["incremental_inductance_H"][0][0] *= -1.0
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False


@pytest.mark.parametrize(
    "case_id",
    ["matrix_symmetry", "requested_current", "legendre_duality", "replay_matrix", "nonlinear_residual"],
)
def test_counterfactual_curriculum90_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "matrix_symmetry":
        bad["runs"][0]["incremental_inductance_H"][0][1] *= 0.5
    elif case_id == "requested_current":
        bad["runs"][0]["current_A"][0] *= 0.5
    elif case_id == "legendre_duality":
        bad["runs"][0]["coenergy_J"] *= 0.5
    elif case_id == "replay_matrix":
        bad["runs"][1]["apparent_inductance_H"][0][0] *= 1.1
    else:
        bad["runs"][0]["final_nonlinear_residual_log10"] = -2.0
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


def test_generalization_v3s_rejects_nonzero_open_secondary_current():
    bad = copy.deepcopy(_summary())
    bad["runs"][0]["current_A"][1] = 1.0
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v4_apparent_symmetry", "v4_flux_identity", "v4_negative_energy", "v4_duplicate_current_level", "v4_saturation_reversal"],
)
def test_counterfactual_curriculum90_v4_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v4_apparent_symmetry":
        bad["runs"][0]["apparent_inductance_H"][0][1] *= 0.25
    elif case_id == "v4_flux_identity":
        bad["runs"][0]["flux_linkage_Vs"][0] *= 1.2
    elif case_id == "v4_negative_energy":
        bad["runs"][0]["energy_J"] = -1.0
    elif case_id == "v4_duplicate_current_level":
        bad["runs"][2]["current_A_requested"] = bad["runs"][0]["current_A_requested"]
    else:
        bad["runs"][6]["apparent_inductance_H"][0][0] *= 10.0
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


def test_generalization_v5_rejects_noncanonical_replay_index():
    bad = copy.deepcopy(_summary())
    bad["runs"][0]["replay"] = 3
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_incremental_matrix_asymmetry", "v6_public_result_metadata_run_mismatch"],
)
def test_generalization_v6_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v6_public_incremental_matrix_asymmetry":
        bad["runs"][0]["incremental_inductance_H"][0][1] *= 0.50
    else:
        bad["runs"][0]["result_metadata"]["energy"]["run_id"] = "wrong-run"
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_operating_point_matrix_mix",
        "v7_public_coenergy_unit_shadowing",
    ],
)
def test_generalization_v7_public(case_id):
    bad = copy.deepcopy(_summary())
    for index, row in enumerate(bad["runs"]):
        operating_point_id = f"op-{index // 2}"
        row.update(
            {
                "operating_point_id": operating_point_id,
                "apparent_matrix_operating_point_id": operating_point_id,
                "incremental_matrix_operating_point_id": operating_point_id,
                "apparent_matrix_current_A": list(row["current_A"]),
                "incremental_matrix_current_A": list(row["current_A"]),
                "reported_units": {
                    "current": "A",
                    "flux_linkage": "Vs",
                    "inductance": "H",
                    "energy": "J",
                    "coenergy": "J",
                },
                "artifact_units": {
                    "current": "A",
                    "flux_linkage": "Vs",
                    "inductance": "H",
                    "energy": "J",
                    "coenergy": "J",
                },
            }
        )
    if case_id == "v7_public_operating_point_matrix_mix":
        for row in bad["runs"][:2]:
            row["incremental_matrix_operating_point_id"] = "op-1"
            row["incremental_matrix_current_A"] = [2.0, 0.0]
    else:
        bad["runs"][0]["artifact_units"]["coenergy"] = "mJ"
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


def test_accepts_v8_sweep_generation_and_restart_offset_chain():
    result = nonlinear_inductance_sweep_gate(_with_v8_generations(_summary()))
    assert result["status"] == "ok"


def test_v8_public_inductance_matrix_prior_sweep_generation():
    bad = _with_v8_generations(_summary())
    bad["runs"][0]["incremental_matrix_sweep_generation"] = "nonlinear-sweep-41"
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False
    assert (
        result["runs"][0]["checks"][
            "inductance_matrices_share_solve_sweep_generation"
        ]
        is False
    )


def test_v8_public_energy_history_restart_offset_lost():
    bad = _with_v8_generations(_summary())
    bad["energy_history_segments"][1]["coenergy_offset_in_J"] = 0.0
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["restart_energy_history_offsets_are_continuous"] is False


def test_v9_public_inductance_matrix_port_order_permuted():
    bad = _with_v9_bindings(_summary())
    bad["runs"][0]["matrix_port_order"].update(
        {
            "incremental_rows": ["secondary", "primary"],
            "incremental_columns": ["secondary", "primary"],
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False
    assert (
        result["runs"][0]["checks"][
            "matrix_rows_columns_and_vectors_share_port_order"
        ]
        is False
    )


def test_v9_public_energy_loss_unit_scale_mismatch():
    bad = _with_v9_bindings(_summary())
    bad["runs"][0]["energy_loss_basis"].update(
        {"loss_series_unit": "nJ", "loss_series_scale_to_J": 1.0}
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False
    assert (
        result["runs"][0]["checks"][
            "stored_energy_and_loss_series_share_si_basis"
        ]
        is False
    )


def test_v10_public_sparameter_reference_impedance_mismatch():
    bad = _with_v10_bindings(_summary())
    bad["runs"][0]["sparameter_reference_impedance"].update(
        {
            "comparison_reference_impedance_ohm_complex": [
                [50.0, 0.0],
                [75.0, 10.0],
            ],
            "renormalization_applied": False,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameters_share_complex_reference_impedance_or_renormalization"
        ]
        is False
    )


def test_v10_public_frequency_axis_unit_metadata_mismatch():
    bad = _with_v10_bindings(_summary())
    bad["runs"][0]["frequency_axis_identity"]["metadata_axis_unit"] = "Hz"
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "frequency_axis_unit_and_hz_scale_share_identity"
        ]
        is False
    )


def test_v11_public_sparameter_port_mode_reference_plane_mismatch():
    bad = _with_v11_bindings(_summary())
    bad["runs"][0]["sparameter_reference_plane_identity"][
        "compared_port_mode_reference_plane_ids"
    ] = ["rp-primary-43", "rp-secondary-42"]
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameter_port_modes_share_deembedded_reference_planes"
        ]
        is False
    )


def test_v11_public_energy_q_factor_frequency_sample_mismatch():
    bad = _with_v11_bindings(_summary())
    bad["runs"][0]["energy_q_frequency_identity"].update(
        {
            "loss_frequency_hz": 1.001e9,
            "loss_sample_id": "adaptive-f-44",
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"]["energy_and_loss_share_q_frequency_sample"]
        is False
    )


def test_v12_public_mixed_mode_sparameter_basis_generation_mismatch():
    bad = _with_v12_bindings(_summary())
    bad["runs"][0]["mixed_mode_sparameter_basis_identity"].update(
        {
            "basis_matrix_port_order": ["P1+", "P1-", "P2-", "P2+"],
            "basis_matrix_port_order_generation": "port-order-13",
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "mixed_mode_basis_matches_current_single_ended_port_order"
        ]
        is False
    )


def test_v12_public_farfield_realized_gain_power_frequency_sample_mismatch():
    bad = _with_v12_bindings(_summary())
    bad["runs"][0]["farfield_realized_gain_power_frequency_identity"].update(
        {
            "accepted_power_frequency_hz": 2.451e9,
            "accepted_power_adaptive_sample_id": "adaptive-frequency-15",
            "accepted_power_result_generation": "farfield-result-15",
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "realized_gain_and_accepted_power_share_frequency_sample"
        ]
        is False
    )


def test_v13_public_field_monitor_interpolation_mesh_pass_mismatch():
    bad = _with_v13_bindings(_summary())
    bad["runs"][0]["field_monitor_interpolation_mesh_pass_identity"].update(
        {
            "interpolation_weight_adaptive_pass_id": "adaptive-pass-14",
            "interpolation_weight_mesh_generation": "adaptive-mesh-14",
            "interpolation_weight_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "field_monitor_interpolation_matches_current_mesh_pass"
        ]
        is False
    )


def test_v13_public_port_deembed_reference_plane_length_unit_mismatch():
    bad = _with_v13_bindings(_summary())
    bad["runs"][0]["port_deembed_reference_plane_unit_identity"].update(
        {
            "result_reference_plane_offset_unit": "m",
            "result_reference_plane_offset_scale_to_m": 1.0,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "port_deembed_reference_plane_uses_explicit_length_unit"
        ]
        is False
    )


def test_accepts_v14_sparameter_and_farfield_basis_lineage():
    result = nonlinear_inductance_sweep_gate(_with_v14_bindings(_summary()))
    assert result["status"] == "ok"


def test_v14_public_sparameter_reference_impedance_renormalization_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0][
        "sparameter_reference_impedance_renormalization_identity"
    ].update(
        {
            "result_reference_impedance_real_ohm": 75.0,
            "result_reference_impedance_imag_ohm": 5.0,
            "result_renormalization_generation": "renormalization-15",
            "renormalized_array_sha256": "d" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameter_renormalization_matches_complex_reference_impedance"
        ]
        is False
    )


def test_v14_public_farfield_ludwig_polarization_basis_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0]["farfield_ludwig_polarization_basis_identity"].update(
        {
            "co_cross_result_generation": "farfield-result-15",
            "co_cross_ludwig_basis_definition": "ludwig_2",
            "co_cross_polarization_reference_axis": "global-y",
            "co_cross_polarization_basis_sha256": "0" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "farfield_co_cross_uses_current_ludwig_polarization_basis"
        ]
        is False
    )


def test_v15_public_sparameter_power_wave_pseudo_wave_normalization_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0]["sparameter_power_wave_normalization_identity"].update(
        {
            "reflected_modal_amplitude_normalization": "pseudo_wave",
            "sparameter_normalization": "pseudo_wave",
            "sparameter_normalization_generation": "modal-normalization-16",
            "sparameter_normalization_sha256": "5" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameters_use_one_power_wave_normalization_generation"
        ]
        is False
    )


def test_v15_public_time_domain_fft_window_coherent_gain_normalization_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0]["time_domain_fft_window_coherent_gain_identity"].update(
        {
            "coherent_gain_window_generation": "fft-window-16",
            "coherent_gain": 1.0,
            "fft_coherent_gain_correction": 1.0,
            "fft_window_coefficients_sha256": "5" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "fft_window_uses_current_coherent_gain_correction"
        ]
        is False
    )


def test_v16_public_sparameter_complex_reference_impedance_renormalization_generation_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0][
        "sparameter_complex_impedance_renormalization_identity"
    ].update(
        {
            "renormalization_port_calibration_generation": "port-calibration-17",
            "result_renormalization_transform_generation": "renormalization-17",
            "result_renormalization_transform_sha256": "4" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameter_renormalization_uses_current_complex_impedance_calibration"
        ]
        is False
    )


def test_v16_public_farfield_polarization_theta_phi_ludwig_basis_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0]["farfield_polarization_basis_transform_identity"].update(
        {
            "basis_transform_applied": False,
            "comparison_transform_generation": "polarization-transform-17",
            "comparison_basis_transform_sha256": "4" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "farfield_comparison_uses_explicit_current_polarization_transform"
        ]
        is False
    )


def test_v17_public_mixed_mode_sparameter_port_pair_order_generation_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0]["mixed_mode_sparameter_port_pair_order_identity"].update(
        {
            "mixed_mode_pairing_port_calibration_generation": "port-calibration-18",
            "transform_pair_order": [["P1+", "P2-"], ["P2+", "P1-"]],
            "transform_pair_polarity": [1, -1],
            "transform_pairing_sha256": "5" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "mixed_mode_sparameters_use_current_port_pair_order_and_polarity"
        ]
        is False
    )


def test_v17_public_nearfield_farfield_phase_center_coordinate_frame_mismatch():
    bad = _with_v14_bindings(_summary())
    bad["runs"][0][
        "nearfield_farfield_phase_center_coordinate_frame_identity"
    ].update(
        {
            "farfield_phase_center_coordinate_frame": "antenna_local",
            "farfield_phase_center_coordinates_m": [0.0, 0.0, 0.0],
            "farfield_phase_center_sha256": "5" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "nearfield_farfield_phase_center_uses_one_global_coordinate_frame"
        ]
        is False
    )


def test_v18_public_sparameter_deembed_reference_plane_per_port_generation_mismatch():
    bad = _with_v14_bindings(_summary())
    deembed = bad["runs"][0][
        "sparameter_deembed_reference_plane_per_port_generation_identity"
    ]
    deembed.update(
        {
            "reference_plane_port_generation": "port-19",
            "deembedded_result_port_generation": "port-19",
            "reference_plane_port_ids": ["P2", "P1"],
            "applied_reference_plane_offsets_m": [0.002, 0.001],
            "deembedded_reference_plane_map_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameter_deembed_uses_current_per_port_reference_planes"
        ]
        is False
    )


def test_v18_public_time_domain_port_signal_gate_window_generation_mismatch():
    bad = _with_v14_bindings(_summary())
    gate_window = bad["runs"][0][
        "time_domain_port_signal_gate_window_generation_identity"
    ]
    gate_window.update(
        {
            "gate_window_signal_generation": "port-signal-19",
            "transform_gate_generation": "gate-19",
            "transform_gate_window": [80, 760],
            "transform_normalization_basis": "total_signal_peak",
            "transform_gate_window_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "time_domain_port_transform_uses_current_gate_window"
        ]
        is False
    )


def test_v19_public_sparameter_port_renormalization_reference_impedance_generation_mismatch():
    summary = _with_v14_bindings(_summary())
    identity = summary["runs"][0][
        "sparameter_port_renormalization_reference_impedance_generation_identity"
    ]
    identity.update(
        {
            "reference_impedance_port_calibration_generation": "port-calibration-20",
            "renormalized_result_port_calibration_generation": "port-calibration-20",
            "reference_impedance_port_ids": ["P2", "P1"],
            "applied_reference_impedances_ohm": [75.0, 50.0],
            "renormalized_reference_impedance_map_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "sparameter_renormalization_uses_current_port_reference_impedances"
        ]
        is False
    )


def test_v19_public_realized_gain_accepted_power_port_excitation_generation_mismatch():
    summary = _with_v14_bindings(_summary())
    identity = summary["runs"][0][
        "realized_gain_accepted_power_port_excitation_generation_identity"
    ]
    identity.update(
        {
            "accepted_power_excitation_generation": "excitation-20",
            "realized_gain_excitation_generation": "excitation-20",
            "accepted_power_port_ids": ["P2", "P1"],
            "realized_gain_accepted_power_w": [0.4, 0.6],
            "realized_gain_excitation_coefficients_re_im": [
                [0.0, 0.0],
                [1.0, 0.0],
            ],
            "realized_gain_excitation_table_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "realized_gain_uses_current_excitation_and_accepted_power"
        ]
        is False
    )


def test_v20_public_farfield_polarization_basis_phase_center_coordinate_generation_mismatch():
    summary = _with_v14_bindings(_summary())
    identity = summary["runs"][0][
        "farfield_polarization_basis_phase_center_coordinate_generation_identity"
    ]
    identity.update(
        {
            "theta_basis_monitor_coordinate_generation": "farfield-monitor-21",
            "phase_center_monitor_coordinate_generation": "farfield-monitor-21",
            "result_sample_ids": [3, 2, 1],
            "result_theta_basis_sha256": ["3" * 64, "2" * 64, "1" * 64],
            "result_phi_basis_sha256": ["6" * 64, "5" * 64, "4" * 64],
            "result_phase_center_m": [0.01, 0.0, 0.0],
            "result_polarization_table_sha256": "f" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "farfield_polarization_basis_and_phase_center_use_current_coordinates"
        ]
        is False
    )


def test_v20_public_broadband_energy_q_port_loss_normalization_generation_mismatch():
    summary = _with_v14_bindings(_summary())
    identity = summary["runs"][0][
        "broadband_energy_q_port_loss_normalization_generation_identity"
    ]
    identity.update(
        {
            "port_power_analysis_generation": "broadband-analysis-21",
            "loss_frequency_grid_generation": "frequency-grid-21",
            "energy_excitation_generation": "excitation-21",
            "q_frequencies_hz": [1.0e9, 2.0e9, 1.5e9],
            "q_stored_energy_j": [1.2e-6, 1.0e-6, 1.1e-6],
            "q_port_power_w": [0.9, 1.0, 1.0],
            "q_loss_power_w": [0.12, 0.1, 0.11],
            "result_energy_q_input_sha256": "f" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "broadband_energy_q_uses_current_port_and_loss_normalization"
        ]
        is False
    )


def _summary_v21():
    summary = _with_v14_bindings(_summary())
    for index, row in enumerate(summary["runs"]):
        port_generation = f"mixed-mode-port-{31 + index}"
        row["mixed_mode_pair_impedance_reference_plane_generation_identity"] = {
            "result_generation": f"mixed-mode-result-{31 + index}",
            "decoded_result_generation": f"mixed-mode-result-{31 + index}",
            "port_generation": port_generation,
            "pair_map_port_generation": port_generation,
            "modal_impedance_port_generation": port_generation,
            "polarity_port_generation": port_generation,
            "reference_plane_port_generation": port_generation,
            "pair_ids": [["P1", "P2"], ["P3", "P4"]],
            "result_pair_ids": [["P1", "P2"], ["P3", "P4"]],
            "pair_polarities": [[1, -1], [1, -1]],
            "result_pair_polarities": [[1, -1], [1, -1]],
            "differential_impedances_ohm": [100.0, 100.0],
            "result_differential_impedances_ohm": [100.0, 100.0],
            "common_impedances_ohm": [25.0, 25.0],
            "result_common_impedances_ohm": [25.0, 25.0],
            "reference_planes_m": [0.0, 0.01],
            "result_reference_planes_m": [0.0, 0.01],
            "mixed_mode_port_table_sha256": "1" * 64,
            "result_mixed_mode_port_table_sha256": "1" * 64,
        }
        monitor_generation = f"time-farfield-monitor-{31 + index}"
        row["time_farfield_fft_window_phase_center_generation_identity"] = {
            "farfield_generation": f"time-farfield-{31 + index}",
            "result_farfield_generation": f"time-farfield-{31 + index}",
            "monitor_generation": monitor_generation,
            "time_grid_monitor_generation": monitor_generation,
            "window_monitor_generation": monitor_generation,
            "fft_scaling_monitor_generation": monitor_generation,
            "phase_center_monitor_generation": monitor_generation,
            "time_samples_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12],
            "fft_time_samples_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12],
            "window_samples": [0.0, 0.75, 0.75, 0.0],
            "fft_window_samples": [0.0, 0.75, 0.75, 0.0],
            "fft_scaling": "one_sided_amplitude",
            "result_fft_scaling": "one_sided_amplitude",
            "phase_center_m": [0.0, 0.0, 0.01],
            "result_phase_center_m": [0.0, 0.0, 0.01],
            "time_farfield_input_sha256": "2" * 64,
            "result_time_farfield_input_sha256": "2" * 64,
        }
    return summary


def test_v21_public_positive_mixed_mode_and_time_farfield_identity():
    result = nonlinear_inductance_sweep_gate(_summary_v21())
    assert result["status"] == "ok"
    assert all(
        row["checks"][
            "mixed_mode_uses_current_pairs_impedances_polarities_and_planes"
        ]
        for row in result["runs"]
    )
    assert all(
        row["checks"][
            "time_farfield_fft_uses_current_grid_window_scaling_and_phase_center"
        ]
        for row in result["runs"]
    )


def test_v21_public_mixed_mode_pair_impedance_reference_plane_generation_mismatch():
    summary = _summary_v21()
    summary["runs"][0][
        "mixed_mode_pair_impedance_reference_plane_generation_identity"
    ].update(
        {
            "pair_map_port_generation": "mixed-mode-port-30",
            "modal_impedance_port_generation": "mixed-mode-port-29",
            "polarity_port_generation": "mixed-mode-port-28",
            "reference_plane_port_generation": "mixed-mode-port-27",
            "result_pair_ids": [["P2", "P1"], ["P4", "P3"]],
            "result_pair_polarities": [[-1, 1], [-1, 1]],
            "result_differential_impedances_ohm": [90.0, 110.0],
            "result_common_impedances_ohm": [30.0, 20.0],
            "result_reference_planes_m": [0.01, 0.0],
            "result_mixed_mode_port_table_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["runs"][0]["checks"][
        "mixed_mode_uses_current_pairs_impedances_polarities_and_planes"
    ] is False


def test_v21_public_time_farfield_fft_window_phase_center_generation_mismatch():
    summary = _summary_v21()
    summary["runs"][0][
        "time_farfield_fft_window_phase_center_generation_identity"
    ].update(
        {
            "time_grid_monitor_generation": "time-farfield-monitor-30",
            "window_monitor_generation": "time-farfield-monitor-29",
            "fft_scaling_monitor_generation": "time-farfield-monitor-28",
            "phase_center_monitor_generation": "time-farfield-monitor-27",
            "fft_time_samples_s": [0.0, 2.0e-12, 1.0e-12, 3.0e-12],
            "fft_window_samples": [1.0, 1.0, 1.0, 1.0],
            "result_fft_scaling": "two_sided_power",
            "result_phase_center_m": [0.01, 0.0, 0.0],
            "result_time_farfield_input_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert result["runs"][0]["checks"][
        "time_farfield_fft_uses_current_grid_window_scaling_and_phase_center"
    ] is False


def _summary_v22():
    summary = _summary_v21()
    for index, row in enumerate(summary["runs"]):
        sweep_generation = f"mode-sweep-{41 + index}"
        row["waveguide_degenerate_mode_phase_order_overlap_generation_identity"] = {
            "sweep_generation": sweep_generation,
            "mesh_sweep_generation": sweep_generation,
            "phase_sweep_generation": sweep_generation,
            "modal_order_sweep_generation": sweep_generation,
            "overlap_sweep_generation": sweep_generation,
            "result_sweep_generation": sweep_generation,
            "mode_ids": ["degenerate-a", "degenerate-b"],
            "result_mode_ids": ["degenerate-a", "degenerate-b"],
            "modal_order": [1, 2],
            "result_modal_order": [1, 2],
            "phase_reference_deg": [0.0, 90.0],
            "result_phase_reference_deg": [0.0, 90.0],
            "overlap_vectors": [[1.0, 0.0], [0.0, 1.0]],
            "result_overlap_vectors": [[1.0, 0.0], [0.0, 1.0]],
            "mode_tracking_table_sha256": "1" * 64,
            "result_mode_tracking_table_sha256": "1" * 64,
        }
        fit_generation = f"dispersive-fit-{41 + index}"
        row["dispersive_causal_pole_fit_temperature_unit_generation_identity"] = {
            "fit_generation": fit_generation,
            "pole_fit_generation": fit_generation,
            "causal_convention_fit_generation": fit_generation,
            "temperature_fit_generation": fit_generation,
            "frequency_unit_fit_generation": fit_generation,
            "field_result_fit_generation": fit_generation,
            "causal_convention": "exp(-iwt)",
            "result_causal_convention": "exp(-iwt)",
            "temperature_c": 85.0,
            "result_temperature_c": 85.0,
            "frequency_unit": "Hz",
            "result_frequency_unit": "Hz",
            "pole_pairs_rad_per_s": [[-1.0e9, 2.0e10], [-1.0e9, -2.0e10]],
            "result_pole_pairs_rad_per_s": [
                [-1.0e9, 2.0e10],
                [-1.0e9, -2.0e10],
            ],
            "residues": [[2.0e9, 1.0e8], [2.0e9, -1.0e8]],
            "result_residues": [[2.0e9, 1.0e8], [2.0e9, -1.0e8]],
            "pole_fit_sha256": "2" * 64,
            "result_pole_fit_sha256": "2" * 64,
        }
    return summary


def test_v22_public_positive_mode_tracking_and_dispersive_fit_identity():
    assert nonlinear_inductance_sweep_gate(_summary_v22())["status"] == "ok"


def test_v22_public_waveguide_degenerate_mode_tracking_phase_order_generation_mismatch():
    summary = _summary_v22()
    identity = summary["runs"][0][
        "waveguide_degenerate_mode_phase_order_overlap_generation_identity"
    ]
    identity.update(
        {
            "mesh_sweep_generation": "mode-sweep-40",
            "phase_sweep_generation": "mode-sweep-39",
            "modal_order_sweep_generation": "mode-sweep-38",
            "overlap_sweep_generation": "mode-sweep-37",
            "result_mode_ids": ["degenerate-b", "degenerate-a"],
            "result_modal_order": [2, 1],
            "result_phase_reference_deg": [90.0, 0.0],
            "result_overlap_vectors": [[0.0, 1.0], [1.0, 0.0]],
            "result_mode_tracking_table_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "degenerate_modes_use_current_mesh_phase_order_and_overlap"
    ]


def test_v22_public_dispersive_material_causal_pole_fit_temperature_unit_generation_mismatch():
    summary = _summary_v22()
    identity = summary["runs"][0][
        "dispersive_causal_pole_fit_temperature_unit_generation_identity"
    ]
    identity.update(
        {
            "pole_fit_generation": "dispersive-fit-40",
            "causal_convention_fit_generation": "dispersive-fit-39",
            "temperature_fit_generation": "dispersive-fit-38",
            "frequency_unit_fit_generation": "dispersive-fit-37",
            "result_causal_convention": "exp(+iwt)",
            "result_temperature_c": 25.0,
            "result_frequency_unit": "GHz",
            "result_pole_pairs_rad_per_s": [
                [-1.0e8, 2.0e9],
                [-1.0e8, -2.0e9],
            ],
            "result_residues": [[2.0e8, -1.0e7], [2.0e8, 1.0e7]],
            "result_pole_fit_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "dispersive_fields_use_current_causal_poles_temperature_and_units"
    ]


def test_v17_public_mixed_mode_identity_rejects_malformed_types() -> None:
    bad = _with_v14_bindings(_summary())
    identity = bad["runs"][0]["mixed_mode_sparameter_port_pair_order_identity"]
    identity["single_ended_port_order"] = [["P1+"], "P1-"]
    identity["mixed_mode_pair_order"] = [["P1+", "P1-"]]
    identity["mixed_mode_pair_polarity"] = [True]
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["runs"][0]["checks"][
            "mixed_mode_sparameters_use_current_port_pair_order_and_polarity"
        ]
        is False
    )


def test_v17_public_mixed_mode_pairs_may_use_nonadjacent_port_indices() -> None:
    summary = _with_v14_bindings(_summary())
    row = summary["runs"][0]
    interleaved = ["P1+", "P2+", "P1-", "P2-"]
    basis = row["mixed_mode_sparameter_basis_identity"]
    basis["single_ended_port_order"] = interleaved
    basis["sparameter_port_order"] = interleaved
    basis["basis_matrix_port_order"] = interleaved
    row["mixed_mode_sparameter_port_pair_order_identity"][
        "single_ended_port_order"
    ] = interleaved
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "ok"


def test_v17_public_phase_center_accepts_roundoff_in_global_meters() -> None:
    summary = _with_v14_bindings(_summary())
    identity = summary["runs"][0][
        "nearfield_farfield_phase_center_coordinate_frame_identity"
    ]
    identity["farfield_phase_center_coordinates_m"] = [
        math.nextafter(value, math.inf)
        for value in identity["phase_center_coordinates_m"]
    ]
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "ok"
