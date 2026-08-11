import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v31 import _v31


_PROMOTED_CASE_IDS = (
    "v32_public_smps_startup_softstart_uvlo_switching_cycle_energy_timestep_result_mismatch",
    "v32_public_noise_source_correlation_spectral_density_bandwidth_integration_rms_mismatch",
)


def _v32():
    payload = _v31()
    positive = payload["metrics"]["positive"]
    generation = "smps-startup-191"
    positive[
        "smps_startup_softstart_uvlo_switch_cycle_timestep_energy_waveform_result_identity"
    ] = {
        "startup_generation_id": generation,
        **{
            key: generation
            for key in (
                "softstart_startup_generation_id",
                "uvlo_startup_generation_id",
                "switch_startup_generation_id",
                "timestep_startup_generation_id",
                "energy_startup_generation_id",
                "waveform_startup_generation_id",
                "result_startup_generation_id",
            )
        },
        "softstart_time_s": [0.0, 1.0e-3, 2.0e-3],
        "result_softstart_time_s": [0.0, 1.0e-3, 2.0e-3],
        "softstart_command": [0.0, 0.5, 1.0],
        "result_softstart_command": [0.0, 0.5, 1.0],
        "uvlo_on_v": 8.0,
        "result_uvlo_on_v": 8.0,
        "uvlo_off_v": 7.0,
        "result_uvlo_off_v": 7.0,
        "first_switching_cycle_s": [2.1e-3, 2.11e-3],
        "result_first_switching_cycle_s": [2.1e-3, 2.11e-3],
        "aligned_timestep_grid_s": [2.1e-3, 2.105e-3, 2.11e-3],
        "result_aligned_timestep_grid_s": [2.1e-3, 2.105e-3, 2.11e-3],
        "input_energy_j": 0.01,
        "result_input_energy_j": 0.01,
        "output_energy_j": 0.007,
        "result_output_energy_j": 0.007,
        "stored_energy_j": 0.002,
        "result_stored_energy_j": 0.002,
        "loss_energy_j": 0.001,
        "result_loss_energy_j": 0.001,
        "waveform_sha256": "1" * 64,
        "result_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "noise-band-191"
    integrated_noise = 990.0e-12
    positive[
        "noise_source_correlation_psd_grid_bandwidth_transfer_integration_rms_model_result_identity"
    ] = {
        "noise_generation_id": generation,
        **{
            key: generation
            for key in (
                "correlation_noise_generation_id",
                "psd_noise_generation_id",
                "grid_noise_generation_id",
                "bandwidth_noise_generation_id",
                "transfer_noise_generation_id",
                "integration_noise_generation_id",
                "model_noise_generation_id",
                "result_noise_generation_id",
            )
        },
        "source_order": ["R1", "R2"],
        "result_source_order": ["R1", "R2"],
        "source_correlation": [[1.0, 0.0], [0.0, 1.0]],
        "result_source_correlation": [[1.0, 0.0], [0.0, 1.0]],
        "psd_convention": "one_sided_v2_per_hz",
        "result_psd_convention": "one_sided_v2_per_hz",
        "frequency_hz": [10.0, 100.0, 1000.0],
        "result_frequency_hz": [10.0, 100.0, 1000.0],
        "integration_bandwidth_hz": [10.0, 1000.0],
        "result_integration_bandwidth_hz": [10.0, 1000.0],
        "transfer_magnitude": [1.0, 1.0, 1.0],
        "result_transfer_magnitude": [1.0, 1.0, 1.0],
        "output_psd_v2_per_hz": [1.0e-12, 1.0e-12, 1.0e-12],
        "result_output_psd_v2_per_hz": [1.0e-12, 1.0e-12, 1.0e-12],
        "integrated_noise_v2": integrated_noise,
        "result_integrated_noise_v2": integrated_noise,
        "rms_noise_v": math.sqrt(integrated_noise),
        "result_rms_noise_v": math.sqrt(integrated_noise),
        "noise_model_sha256": "3" * 64,
        "result_noise_model_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v32_public_positive_smps_startup_and_noise_band_identities():
    assert ideal_transformer_identity_gate(_v32())["status"] == "ok"


def test_v32_public_smps_startup_softstart_uvlo_switching_cycle_energy_timestep_result_mismatch():
    payload = _v32()
    identity = payload["metrics"]["positive"][
        "smps_startup_softstart_uvlo_switch_cycle_timestep_energy_waveform_result_identity"
    ]
    identity.update(
        {
            "softstart_startup_generation_id": "smps-startup-190",
            "energy_startup_generation_id": "smps-startup-189",
            "result_startup_generation_id": "smps-startup-188",
            "result_softstart_time_s": [0.0, 0.5e-3],
            "result_softstart_command": [1.0, 0.0],
            "result_uvlo_on_v": 6.0,
            "result_uvlo_off_v": 9.0,
            "result_first_switching_cycle_s": [1.0e-3, 1.01e-3],
            "result_aligned_timestep_grid_s": [0.0, 1.0e-3],
            "result_input_energy_j": 0.005,
            "result_output_energy_j": 0.009,
            "result_stored_energy_j": -0.001,
            "result_loss_energy_j": 0.0,
            "result_waveform_sha256": "9" * 64,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "smps_startup_uses_current_softstart_uvlo_switch_cycle_timestep_energy_waveform_and_result"
    ]


def test_v32_public_noise_source_correlation_spectral_density_bandwidth_integration_rms_mismatch():
    payload = _v32()
    identity = payload["metrics"]["positive"][
        "noise_source_correlation_psd_grid_bandwidth_transfer_integration_rms_model_result_identity"
    ]
    identity.update(
        {
            "correlation_noise_generation_id": "noise-band-190",
            "model_noise_generation_id": "noise-band-189",
            "result_noise_generation_id": "noise-band-188",
            "result_source_order": ["R2", "R1"],
            "result_source_correlation": [[1.0, 1.2], [0.0, 1.0]],
            "result_psd_convention": "two_sided_v2_per_hz",
            "result_frequency_hz": [1000.0, 100.0, 10.0],
            "result_integration_bandwidth_hz": [100.0, 500.0],
            "result_transfer_magnitude": [0.0, 0.0],
            "result_output_psd_v2_per_hz": [-1.0e-12],
            "result_integrated_noise_v2": -1.0,
            "result_rms_noise_v": 1.0,
            "result_noise_model_sha256": "b" * 64,
            "accepted_result_sha256": "c" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "noise_bands_use_current_sources_correlation_psd_grid_bandwidth_transfer_integration_model_and_result"
    ]


def test_v32_public_noise_rejects_symmetric_non_psd_correlation_matrix():
    payload = _v32()
    identity = payload["metrics"]["positive"][
        "noise_source_correlation_psd_grid_bandwidth_transfer_integration_rms_model_result_identity"
    ]
    correlation = [
        [1.0, -0.9, -0.9],
        [-0.9, 1.0, -0.9],
        [-0.9, -0.9, 1.0],
    ]
    identity.update(
        {
            "source_order": ["R1", "R2", "R3"],
            "result_source_order": ["R1", "R2", "R3"],
            "source_correlation": correlation,
            "result_source_correlation": correlation,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "noise_bands_use_current_sources_correlation_psd_grid_bandwidth_transfer_integration_model_and_result"
    ]
