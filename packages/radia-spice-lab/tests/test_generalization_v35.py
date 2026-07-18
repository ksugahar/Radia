import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v34 import _v34

_PROMOTED_CASE_IDS = (
    "v35_public_switched_converter_poincare_limitcycle_floquet_event_energy_mismatch",
    "v35_public_noise_input_output_referred_density_correlation_integration_gain_mismatch",
)


def _v35():
    payload = _v34()
    positive = payload["metrics"]["positive"]
    generation = "switched-limit-cycle-221"
    positive[
        "switched_converter_poincare_limitcycle_floquet_event_flux_charge_cycle_energy_step_waveform_owner_result_identity"
    ] = {
        "limit_cycle_generation_id": generation,
        **{
            key: generation
            for key in (
                "poincare_limit_cycle_generation_id",
                "floquet_limit_cycle_generation_id",
                "event_limit_cycle_generation_id",
                "flux_limit_cycle_generation_id",
                "charge_limit_cycle_generation_id",
                "energy_limit_cycle_generation_id",
                "step_limit_cycle_generation_id",
                "waveform_limit_cycle_generation_id",
                "owner_limit_cycle_generation_id",
                "result_limit_cycle_generation_id",
            )
        },
        "poincare_state_start": [2.0, 12.0],
        "poincare_state_end": [2.0, 12.0],
        "result_poincare_state_start": [2.0, 12.0],
        "result_poincare_state_end": [2.0, 12.0],
        "poincare_tolerance": 1e-8,
        "result_poincare_tolerance": 1e-8,
        "event_sequence": ["switch_on", "switch_off", "switch_on"],
        "result_event_sequence": ["switch_on", "switch_off", "switch_on"],
        "event_times_s": [0.0, 5e-6, 1e-5],
        "result_event_times_s": [0.0, 5e-6, 1e-5],
        "floquet_multipliers_ri": [[0.8, 0.0], [0.6, 0.1]],
        "result_floquet_multipliers_ri": [[0.8, 0.0], [0.6, 0.1]],
        "inductor_flux_start_wb": [2e-4],
        "inductor_flux_end_wb": [2e-4],
        "result_inductor_flux_start_wb": [2e-4],
        "result_inductor_flux_end_wb": [2e-4],
        "capacitor_charge_start_c": [1.2e-5],
        "capacitor_charge_end_c": [1.2e-5],
        "result_capacitor_charge_start_c": [1.2e-5],
        "result_capacitor_charge_end_c": [1.2e-5],
        "input_cycle_energy_j": 1e-3,
        "output_cycle_energy_j": 8e-4,
        "loss_cycle_energy_j": 2e-4,
        "result_input_cycle_energy_j": 1e-3,
        "result_output_cycle_energy_j": 8e-4,
        "result_loss_cycle_energy_j": 2e-4,
        "step_owner": "converter/step-221",
        "accepted_step_owner": "converter/step-221",
        "step_sha256": "1" * 64,
        "accepted_step_sha256": "1" * 64,
        "waveform_sha256": "2" * 64,
        "accepted_waveform_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    generation = "noise-referral-221"
    positive[
        "noise_input_output_referred_density_correlation_integration_gain_circuit_result_identity"
    ] = {
        "noise_generation_id": generation,
        **{
            key: generation
            for key in (
                "input_noise_generation_id",
                "output_noise_generation_id",
                "correlation_noise_generation_id",
                "integration_noise_generation_id",
                "gain_noise_generation_id",
                "circuit_noise_generation_id",
                "result_noise_generation_id",
            )
        },
        "frequency_hz": [1e3, 2e3],
        "result_frequency_hz": [1e3, 2e3],
        "input_referred_density_v_per_sqrt_hz": [1e-9, 1e-9],
        "result_input_referred_density_v_per_sqrt_hz": [1e-9, 1e-9],
        "output_referred_density_v_per_sqrt_hz": [2e-9, 3e-9],
        "result_output_referred_density_v_per_sqrt_hz": [2e-9, 3e-9],
        "transfer_gain_magnitude": [2.0, 3.0],
        "result_transfer_gain_magnitude": [2.0, 3.0],
        "source_correlation_matrix": [[1.0, 0.2], [0.2, 1.0]],
        "result_source_correlation_matrix": [[1.0, 0.2], [0.2, 1.0]],
        "spectral_density_unit": "V/sqrt(Hz)",
        "result_spectral_density_unit": "V/sqrt(Hz)",
        "integration_bandwidth_hz": [1e3, 2e3],
        "result_integration_bandwidth_hz": [1e3, 2e3],
        "integration_bin_width_hz": [500.0, 500.0],
        "result_integration_bin_width_hz": [500.0, 500.0],
        "total_output_rms_noise_v": math.sqrt(6.5e-15),
        "result_total_output_rms_noise_v": math.sqrt(6.5e-15),
        "circuit_owner": "noise/circuit-221",
        "accepted_circuit_owner": "noise/circuit-221",
        "circuit_sha256": "4" * 64,
        "accepted_circuit_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return payload


def test_v35_public_positive_limit_cycle_and_noise_referral_closure():
    assert ideal_transformer_identity_gate(_v35())["status"] == "ok"


def test_v35_public_switched_converter_poincare_limitcycle_floquet_event_energy_mismatch():
    payload = _v35()
    identity = payload["metrics"]["positive"][
        "switched_converter_poincare_limitcycle_floquet_event_flux_charge_cycle_energy_step_waveform_owner_result_identity"
    ]
    identity.update(
        {
            "poincare_limit_cycle_generation_id": "switched-limit-cycle-220",
            "result_poincare_state_end": [3.0, 9.0],
            "result_event_sequence": ["switch_off", "switch_on"],
            "result_floquet_multipliers_ri": [[1.2, 0.0]],
            "result_loss_cycle_energy_j": -3e-4,
            "accepted_step_sha256": "a" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "switched_limit_cycles_use_current_poincare_floquet_events_flux_charge_energy_step_waveform_and_result"
    ]


def test_v35_public_noise_input_output_referred_density_correlation_integration_gain_mismatch():
    payload = _v35()
    identity = payload["metrics"]["positive"][
        "noise_input_output_referred_density_correlation_integration_gain_circuit_result_identity"
    ]
    identity.update(
        {
            "input_noise_generation_id": "noise-referral-220",
            "result_frequency_hz": [2e3, 1e3],
            "result_output_referred_density_v_per_sqrt_hz": [1.0, 1.0],
            "result_source_correlation_matrix": [[1.0, 2.0], [-1.0, 1.0]],
            "result_spectral_density_unit": "V^2/Hz",
            "result_total_output_rms_noise_v": -1.0,
            "accepted_circuit_sha256": "b" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "noise_referral_uses_current_input_output_density_correlation_band_gain_circuit_and_result"
    ]


def test_v35_public_rejects_self_consistent_unstable_floquet_multiplier():
    payload = _v35()
    identity = payload["metrics"]["positive"][
        "switched_converter_poincare_limitcycle_floquet_event_flux_charge_cycle_energy_step_waveform_owner_result_identity"
    ]
    identity["floquet_multipliers_ri"] = identity["result_floquet_multipliers_ri"] = [[1.01, 0.0]]
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v35_public_rejects_self_consistent_wrong_noise_referral_gain():
    payload = _v35()
    identity = payload["metrics"]["positive"][
        "noise_input_output_referred_density_correlation_integration_gain_circuit_result_identity"
    ]
    identity["output_referred_density_v_per_sqrt_hz"] = [3e-9, 3e-9]
    identity["result_output_referred_density_v_per_sqrt_hz"] = [3e-9, 3e-9]
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
