from __future__ import annotations

from test_comsol_generalization_v36 import _summary, gate


_PROMOTED_CASE_IDS = (
    "v45_public_microwave_sparameter_port_reference_plane_deembed_complex_power_mesh_result_mismatch",
    "v45_public_acoustic_impedance_absorption_phase_energy_flux_farfield_window_dataset_mismatch",
)


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


def test_v45_public_positive_identity() -> None:
    assert gate(_with_v45(_summary()))["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v45_public_rejects_port_deembed_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["microwave_sparameter_port_reference_plane_deembed_complex_power_mesh_result_identity"]["result_deembed_length_m"] = -1.0
    assert gate(summary)["status"] == "needs_attention"


def test_v45_public_rejects_acoustic_window_mismatch() -> None:
    summary = _with_v45(_summary())
    summary["acoustic_impedance_absorption_phase_energy_flux_farfield_window_dataset_result_identity"]["result_time_window_s"] = -1.0
    assert gate(summary)["status"] == "needs_attention"
