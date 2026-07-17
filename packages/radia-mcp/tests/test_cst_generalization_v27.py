from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v26 import _summary_v26


_PROMOTED_CASE_IDS = (
    "v27_public_time_domain_port_waveform_normalization_fft_window_grid_smatrix_generation_mismatch",
    "v27_public_huygens_box_orientation_phase_center_frequency_mesh_near_far_transform_mismatch",
)


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


def test_v27_public_positive_time_domain_and_huygens_identities() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v27())["status"] == "ok"


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
