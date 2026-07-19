from __future__ import annotations

from radia_mcp.radia_ngsolve.waveguide_emc_v44_identity import validate_public_identity


_WAVEGUIDE = "waveguide_modal_cutoff_impedance_groupdelay_power_orthogonality_mesh_result_identity"
_EMC = "emc_probe_interpolation_timewindow_fft_parseval_coordinate_monitor_owner_identity"
_PROMOTED_CASE_IDS = (
    "v44_public_waveguide_modal_cutoff_impedance_groupdelay_power_orthogonality_mesh_mismatch",
    "v44_public_emc_probe_interpolation_timewindow_fft_parseval_coordinate_monitor_owner_mismatch",
)


def _payload() -> dict:
    wave = "waveguide-mode-844-0"
    emc = "emc-probe-844-0"
    return {"runs": [{
        _WAVEGUIDE: {
            "waveguide_generation": wave,
            **{key: wave for key in ("frequency_generation", "cutoff_generation", "impedance_generation", "groupdelay_generation", "power_generation", "orthogonality_generation", "mesh_generation", "result_generation")},
            "frequency_hz": [8.0e9, 10.0e9, 12.0e9], "result_frequency_hz": [8.0e9, 10.0e9, 12.0e9], "cutoff_frequency_hz": 6.0e9, "result_cutoff_frequency_hz": 6.0e9,
            "modal_impedance_ohm": [50.0, 55.0, 60.0], "result_modal_impedance_ohm": [50.0, 55.0, 60.0], "group_delay_s": [1.0e-9, 1.1e-9, 1.2e-9], "result_group_delay_s": [1.0e-9, 1.1e-9, 1.2e-9],
            "power_normalization_w": [1.0, 1.0, 1.0], "result_power_normalization_w": [1.0, 1.0, 1.0], "orthogonality_matrix": [[1.0, 0.0], [0.0, 1.0]], "result_orthogonality_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "mesh_owner": "mesh:" + wave, "result_mesh_owner": "mesh:" + wave, "waveguide_result_sha256": "d" * 64, "accepted_waveguide_result_sha256": "d" * 64,
        },
        _EMC: {
            "emc_probe_generation": emc,
            **{key: emc for key in ("coordinate_generation", "interpolation_generation", "timewindow_generation", "fft_generation", "parseval_generation", "monitor_generation", "result_generation")},
            "coordinate_system": "global_cartesian", "result_coordinate_system": "global_cartesian", "interpolation_order": 2, "result_interpolation_order": 2, "time_window_s": [0.0, 1.0e-9], "result_time_window_s": [0.0, 1.0e-9],
            "fft_normalization": "parseval_unitary", "result_fft_normalization": "parseval_unitary", "parseval_time_energy_j": 1.0, "result_parseval_time_energy_j": 1.0, "parseval_frequency_energy_j": 1.0, "result_parseval_frequency_energy_j": 1.0,
            "monitor_owner": "monitor:" + emc, "result_monitor_owner": "monitor:" + emc, "emc_probe_result_sha256": "e" * 64, "accepted_emc_probe_result_sha256": "e" * 64,
        },
    }]}


def test_v44_public_waveguide_and_emc_identity_positive() -> None:
    assert validate_public_identity(_payload()) == {"waveguide_v44_modal_identity": True, "waveguide_v44_emc_probe_identity": True}


def test_v44_public_identity_rejects_cutoff_and_parseval_mutations() -> None:
    payload = _payload()
    payload["runs"][0][_WAVEGUIDE]["result_cutoff_frequency_hz"] = 12.0e9
    payload["runs"][0][_EMC]["result_parseval_frequency_energy_j"] = 2.0
    result = validate_public_identity(payload)
    assert result["waveguide_v44_modal_identity"] is False
    assert result["waveguide_v44_emc_probe_identity"] is False
