from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v32 import _v32


_PROMOTED_CASE_IDS = (
    "v33_public_transmission_line_z0_delay_reflection_arrival_polarity_energy_causality_mismatch",
    "v33_public_sampled_loop_sideband_injection_period_fft_bin_nyquist_crossover_mismatch",
)


def _v33():
    payload = _v32()
    positive = payload["metrics"]["positive"]
    generation = "tline-transient-201"
    positive[
        "transmission_line_z0_delay_reflection_arrival_polarity_causality_energy_waveform_result_identity"
    ] = {
        "transmission_line_generation_id": generation,
        **{
            key: generation
            for key in (
                "impedance_transmission_line_generation_id",
                "delay_transmission_line_generation_id",
                "reflection_transmission_line_generation_id",
                "arrival_transmission_line_generation_id",
                "polarity_transmission_line_generation_id",
                "causality_transmission_line_generation_id",
                "energy_transmission_line_generation_id",
                "waveform_transmission_line_generation_id",
                "result_transmission_line_generation_id",
            )
        },
        "characteristic_impedance_ohm": 50.0,
        "result_characteristic_impedance_ohm": 50.0,
        "source_impedance_ohm": 50.0,
        "result_source_impedance_ohm": 50.0,
        "load_impedance_ohm": 100.0,
        "result_load_impedance_ohm": 100.0,
        "source_reflection_coefficient": 0.0,
        "result_source_reflection_coefficient": 0.0,
        "load_reflection_coefficient": 1.0 / 3.0,
        "result_load_reflection_coefficient": 1.0 / 3.0,
        "one_way_delay_s": 2.0e-9,
        "result_one_way_delay_s": 2.0e-9,
        "incident_arrival_s": 2.0e-9,
        "result_incident_arrival_s": 2.0e-9,
        "reflected_source_arrival_s": 4.0e-9,
        "result_reflected_source_arrival_s": 4.0e-9,
        "incident_pulse_polarity": "positive",
        "result_incident_pulse_polarity": "positive",
        "reflected_pulse_polarity": "positive",
        "result_reflected_pulse_polarity": "positive",
        "time_s": [0.0, 1.0e-9, 2.0e-9, 3.0e-9, 4.0e-9, 5.0e-9],
        "result_time_s": [0.0, 1.0e-9, 2.0e-9, 3.0e-9, 4.0e-9, 5.0e-9],
        "source_observation_v": [0.0, 0.0, 1.0, 1.0, 4.0 / 3.0, 4.0 / 3.0],
        "result_source_observation_v": [
            0.0,
            0.0,
            1.0,
            1.0,
            4.0 / 3.0,
            4.0 / 3.0,
        ],
        "pre_arrival_max_abs_v": 0.0,
        "result_pre_arrival_max_abs_v": 0.0,
        "incident_energy_j": 9.0,
        "result_incident_energy_j": 9.0,
        "reflected_energy_j": 1.0,
        "result_reflected_energy_j": 1.0,
        "accepted_energy_j": 8.0,
        "result_accepted_energy_j": 8.0,
        "waveform_sha256": "1" * 64,
        "result_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "sampled-loop-201"
    positive[
        "sampled_loop_sideband_injection_period_fft_bin_phase_nyquist_crossover_waveform_result_identity"
    ] = {
        "sampled_loop_generation_id": generation,
        **{
            key: generation
            for key in (
                "sideband_sampled_loop_generation_id",
                "injection_sampled_loop_generation_id",
                "period_sampled_loop_generation_id",
                "fft_sampled_loop_generation_id",
                "phase_sampled_loop_generation_id",
                "nyquist_sampled_loop_generation_id",
                "crossover_sampled_loop_generation_id",
                "waveform_sampled_loop_generation_id",
                "result_sampled_loop_generation_id",
            )
        },
        "sideband_order": 0,
        "result_sideband_order": 0,
        "sideband_selection": "fundamental",
        "result_sideband_selection": "fundamental",
        "injection_point": "control_to_duty_break",
        "result_injection_point": "control_to_duty_break",
        "switching_period_s": 10.0e-6,
        "result_switching_period_s": 10.0e-6,
        "switching_frequency_hz": 100.0e3,
        "result_switching_frequency_hz": 100.0e3,
        "record_duration_s": 200.0e-6,
        "result_record_duration_s": 200.0e-6,
        "injection_frequency_hz": 5.0e3,
        "result_injection_frequency_hz": 5.0e3,
        "coherent_fft_bin": 1,
        "result_coherent_fft_bin": 1,
        "loop_frequency_hz": [2.0e3, 5.0e3, 10.0e3],
        "result_loop_frequency_hz": [2.0e3, 5.0e3, 10.0e3],
        "loop_magnitude": [2.0, 1.0, 0.5],
        "result_loop_magnitude": [2.0, 1.0, 0.5],
        "loop_phase_deg": [-90.0, -135.0, -180.0],
        "result_loop_phase_deg": [-90.0, -135.0, -180.0],
        "crossover_frequency_hz": 5.0e3,
        "result_crossover_frequency_hz": 5.0e3,
        "nyquist_clockwise_encirclements_minus_one": 0,
        "result_nyquist_clockwise_encirclements_minus_one": 0,
        "open_loop_right_half_plane_poles": 0,
        "result_open_loop_right_half_plane_poles": 0,
        "closed_loop_right_half_plane_poles": 0,
        "result_closed_loop_right_half_plane_poles": 0,
        "waveform_sha256": "3" * 64,
        "result_waveform_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v33_public_positive_transmission_line_and_sampled_loop_identities():
    assert ideal_transformer_identity_gate(_v33())["status"] == "ok"


def test_v33_public_transmission_line_z0_delay_reflection_arrival_polarity_energy_causality_mismatch():
    payload = _v33()
    identity = payload["metrics"]["positive"][
        "transmission_line_z0_delay_reflection_arrival_polarity_causality_energy_waveform_result_identity"
    ]
    identity.update(
        {
            "impedance_transmission_line_generation_id": "tline-transient-200",
            "energy_transmission_line_generation_id": "tline-transient-199",
            "result_transmission_line_generation_id": "tline-transient-198",
            "result_characteristic_impedance_ohm": 75.0,
            "result_source_impedance_ohm": 25.0,
            "result_load_impedance_ohm": 25.0,
            "result_source_reflection_coefficient": -0.5,
            "result_load_reflection_coefficient": -0.5,
            "result_one_way_delay_s": 3.0e-9,
            "result_incident_arrival_s": 1.0e-9,
            "result_reflected_source_arrival_s": 3.0e-9,
            "result_incident_pulse_polarity": "negative",
            "result_reflected_pulse_polarity": "negative",
            "result_time_s": [0.0, 1.0e-9, 0.5e-9],
            "result_source_observation_v": [0.5, 1.0, -1.0],
            "result_pre_arrival_max_abs_v": 0.5,
            "result_incident_energy_j": 4.0,
            "result_reflected_energy_j": 5.0,
            "result_accepted_energy_j": -1.0,
            "result_waveform_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "transmission_lines_use_current_z0_delay_reflections_arrivals_polarity_causality_energy_waveform_and_result"
    ]


def test_v33_public_sampled_loop_sideband_injection_period_fft_bin_nyquist_crossover_mismatch():
    payload = _v33()
    identity = payload["metrics"]["positive"][
        "sampled_loop_sideband_injection_period_fft_bin_phase_nyquist_crossover_waveform_result_identity"
    ]
    identity.update(
        {
            "sideband_sampled_loop_generation_id": "sampled-loop-200",
            "waveform_sampled_loop_generation_id": "sampled-loop-199",
            "result_sampled_loop_generation_id": "sampled-loop-198",
            "result_sideband_order": 1,
            "result_sideband_selection": "upper_first",
            "result_injection_point": "output_shunt",
            "result_switching_period_s": 20.0e-6,
            "result_switching_frequency_hz": 40.0e3,
            "result_record_duration_s": 190.0e-6,
            "result_injection_frequency_hz": 6.0e3,
            "result_coherent_fft_bin": 2,
            "result_loop_frequency_hz": [10.0e3, 5.0e3, 2.0e3],
            "result_loop_magnitude": [0.5, 1.2, 2.0],
            "result_loop_phase_deg": [180.0, 135.0, 90.0],
            "result_crossover_frequency_hz": 8.0e3,
            "result_nyquist_clockwise_encirclements_minus_one": 1,
            "result_open_loop_right_half_plane_poles": 1,
            "result_closed_loop_right_half_plane_poles": 2,
            "result_waveform_sha256": "c" * 64,
            "accepted_result_sha256": "d" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sampled_loops_use_current_sideband_injection_period_fft_phase_nyquist_crossover_waveform_and_result"
    ]


def test_v33_public_rejects_self_consistent_but_nonphysical_reflection_energy():
    payload = _v33()
    identity = payload["metrics"]["positive"][
        "transmission_line_z0_delay_reflection_arrival_polarity_causality_energy_waveform_result_identity"
    ]
    identity["load_reflection_coefficient"] = 0.5
    identity["result_load_reflection_coefficient"] = 0.5
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"


def test_v33_public_rejects_self_consistent_noncoherent_fft_bin():
    payload = _v33()
    identity = payload["metrics"]["positive"][
        "sampled_loop_sideband_injection_period_fft_bin_phase_nyquist_crossover_waveform_result_identity"
    ]
    identity["coherent_fft_bin"] = 2
    identity["result_coherent_fft_bin"] = 2
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
