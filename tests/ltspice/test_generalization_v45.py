from radia.ltspice.ltspice_v45_gates import validate_ltspice_v45_identity

PROMOTED_CASE_IDS = (
    "v45_public_buck_startup_softstart_inductor_current_ripple_efficiency_energy_waveform_mismatch",
    "v45_public_noise_ac_transfer_psd_sidedness_bandwidth_correlation_measure_result_mismatch",
)


def _positive():
    return {
        "buck_v45_startup_softstart_inductor_current_ripple_efficiency_energy_waveform_identity": {
            "generation_id": "buck-v45",
            "release_id": "LTSpice-2026.1",
            "result_release_id": "LTSpice-2026.1",
            "waveform_owner": "waveform/current",
            "accepted_waveform_owner": "waveform/current",
            "result_waveform_owner": "waveform/current",
            "result_sha256": "9" * 64,
            "accepted_result_sha256": "9" * 64,
        },
        "noise_v45_ac_transfer_psd_sidedness_bandwidth_correlation_measure_identity": {
            "generation_id": "noise-v45",
            "release_id": "LTSpice-2026.1",
            "result_release_id": "LTSpice-2026.1",
            "measure_owner": "measure/current",
            "accepted_measure_owner": "measure/current",
            "result_measure_owner": "measure/current",
            "result_sha256": "a" * 64,
            "accepted_result_sha256": "a" * 64,
        },
    }


def test_v45_public_release_owner_digest_identity():
    assert validate_ltspice_v45_identity(_positive()) is True


def test_v45_public_rejects_release_and_owner_mutations():
    identity = _positive()
    identity["buck_v45_startup_softstart_inductor_current_ripple_efficiency_energy_waveform_identity"]["result_release_id"] = "LTSpice-2025.4"
    assert validate_ltspice_v45_identity(identity) is False
