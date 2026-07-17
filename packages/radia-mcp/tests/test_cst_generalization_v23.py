from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
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


def test_v23_public_positive_broadband_and_transient_identity() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v23())["status"] == "ok"


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
