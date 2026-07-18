from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v36 import _summary_v36


_PROMOTED_CASE_IDS = (
    "v37_public_waveguide_cutoff_mode_impedance_group_delay_power_orthogonality_mesh_mismatch",
    "v37_public_emc_probe_coordinate_interpolation_time_frequency_window_parseval_owner_mismatch",
)
C0 = 299_792_458.0
ETA0 = 376.730313668


def _waveguide(index: int):
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
        row["waveguide_cutoff_mode_impedance_group_delay_power_orthogonality_mesh_owner_result_identity"] = _waveguide(index)
        row["emc_probe_coordinate_interpolation_time_fft_window_parseval_mesh_owner_result_identity"] = _probe(index)
    return summary


def test_v37_public_positive_waveguide_and_emc_probe_closure():
    assert nonlinear_inductance_sweep_gate(_summary_v37())["status"] == "ok"


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
