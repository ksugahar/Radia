from copy import deepcopy
import math

from radia.ltspice.ltspice_v53_gates import NOISE, STEADY, validate_ltspice_v53_identity


PROMOTED_CASE_IDS = {
    "v53_public_noise_spectraldensity_integrated_bandwidth_source_owner_mismatch",
    "v53_public_switchmode_sampled_steadystate_cycle_phase_ripple_waveform_owner_mismatch",
}


def _positive() -> dict[str, object]:
    noise_generation = "noise-public-v53"
    steady_generation = "steady-public-v53"
    frequency = [10.0, 100.0, 1000.0, 10000.0]
    density = [4.0e-18, 2.0e-18, 1.2e-18, 1.0e-18]
    variance = sum(0.5 * (right_density + left_density) * (right_frequency - left_frequency) for left_frequency, right_frequency, left_density, right_density in zip(frequency, frequency[1:], density, density[1:]))
    waveforms = [[4.94, 4.99, 5.03, 5.01, 4.94], [4.945, 4.992, 5.031, 5.012, 4.945], [4.947, 4.993, 5.032, 5.013, 4.947]]
    sampled = [row[2] for row in waveforms]
    ripple = [max(row) - min(row) for row in waveforms]
    return {
        NOISE: {
            "generation_id": noise_generation,
            **{name: noise_generation for name in ("spectral_generation_id", "bandwidth_generation_id", "source_generation_id", "owner_generation_id", "result_generation_id")},
            "frequency_hz": frequency, "result_frequency_hz": frequency,
            "spectral_density_v2_per_hz": density, "result_spectral_density_v2_per_hz": density,
            "integration_bandwidth_hz": [10.0, 10000.0], "result_integration_bandwidth_hz": [10.0, 10000.0],
            "integrated_noise_v_rms": math.sqrt(variance), "result_integrated_noise_v_rms": math.sqrt(variance),
            "contributing_source_fraction": {"source:R1": 0.65, "source:M1": 0.35}, "result_contributing_source_fraction": {"source:R1": 0.65, "source:M1": 0.35},
            "trace_owner": "trace:onoise-v53", "result_trace_owner": "trace:onoise-v53",
            "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
        },
        STEADY: {
            "generation_id": steady_generation,
            **{name: steady_generation for name in ("cycle_generation_id", "phase_generation_id", "ripple_generation_id", "waveform_generation_id", "owner_generation_id", "result_generation_id")},
            "cycle_index": [98, 99, 100], "result_cycle_index": [98, 99, 100],
            "sample_phase_fraction": 0.5, "result_sample_phase_fraction": 0.5,
            "cycle_waveform_v": waveforms, "result_cycle_waveform_v": waveforms,
            "sampled_output_v": sampled, "result_sampled_output_v": sampled,
            "ripple_pp_v": ripple, "result_ripple_pp_v": ripple,
            "waveform_owner": "waveform:buck-steady-v53", "result_waveform_owner": "waveform:buck-steady-v53",
            "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
        },
    }


def test_v53_positive_identity_is_accepted() -> None:
    assert validate_ltspice_v53_identity(_positive()) is True


def test_v53_frozen_public_counterfactuals_are_rejected() -> None:
    value = deepcopy(_positive())
    value[NOISE].update({"result_spectral_density_v2_per_hz": [1.0e-12] * 4, "result_integration_bandwidth_hz": [100.0, 1000.0], "result_integrated_noise_v_rms": 1.0, "result_contributing_source_fraction": {"source:foreign": 1.0}, "result_trace_owner": "trace:foreign"})
    value[STEADY].update({"result_cycle_index": [1, 2, 3], "result_sample_phase_fraction": 0.0, "result_sampled_output_v": [0.0] * 3, "result_ripple_pp_v": [9.0] * 3, "result_waveform_owner": "waveform:foreign"})
    assert validate_ltspice_v53_identity(value) is False


def test_v53_self_consistent_wrong_noise_integral_is_rejected() -> None:
    value = deepcopy(_positive())
    value[NOISE]["integrated_noise_v_rms"] = value[NOISE]["result_integrated_noise_v_rms"] = 1.0
    value[NOISE]["contributing_source_fraction"] = value[NOISE]["result_contributing_source_fraction"] = {"source:R1": 0.8, "source:M1": 0.8}
    assert validate_ltspice_v53_identity(value) is False


def test_v53_self_consistent_wrong_sampling_is_rejected() -> None:
    value = deepcopy(_positive())
    value[STEADY]["sampled_output_v"] = value[STEADY]["result_sampled_output_v"] = [0.0, 0.0, 0.0]
    value[STEADY]["ripple_pp_v"] = value[STEADY]["result_ripple_pp_v"] = [1.0, 1.0, 1.0]
    assert validate_ltspice_v53_identity(value) is False
