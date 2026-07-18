import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v33 import _v33


_PROMOTED_CASE_IDS = (
    "v34_public_periodic_steady_state_charge_flux_cycle_energy_efficiency_phase_owner_mismatch",
    "v34_public_bias_small_signal_ac_transient_jacobian_pole_zero_transfer_consistency_mismatch",
)


def _v34():
    payload = _v33()
    positive = payload["metrics"]["positive"]
    generation = "periodic-steady-state-211"
    positive[
        "periodic_steady_state_charge_flux_cycle_energy_efficiency_phase_waveform_owner_result_identity"
    ] = {
        "periodic_generation_id": generation,
        **{
            key: generation
            for key in (
                "charge_periodic_generation_id", "flux_periodic_generation_id",
                "energy_periodic_generation_id", "efficiency_periodic_generation_id",
                "phase_periodic_generation_id", "waveform_periodic_generation_id",
                "owner_periodic_generation_id", "result_periodic_generation_id",
            )
        },
        "cycle_start_charge_c": [1e-3, -1e-3], "cycle_end_charge_c": [1e-3, -1e-3],
        "result_cycle_start_charge_c": [1e-3, -1e-3],
        "result_cycle_end_charge_c": [1e-3, -1e-3],
        "cycle_start_flux_wb": [2e-4], "cycle_end_flux_wb": [2e-4],
        "result_cycle_start_flux_wb": [2e-4], "result_cycle_end_flux_wb": [2e-4],
        "closure_tolerance": 1e-9, "result_closure_tolerance": 1e-9,
        "input_energy_j": 10.0, "result_input_energy_j": 10.0,
        "output_energy_j": 8.0, "result_output_energy_j": 8.0,
        "loss_energy_j": 2.0, "result_loss_energy_j": 2.0,
        "efficiency": 0.8, "result_efficiency": 0.8,
        "phase_window_rad": [0.0, 2 * math.pi],
        "result_phase_window_rad": [0.0, 2 * math.pi],
        "waveform_owner": "converter/pss-211", "accepted_waveform_owner": "converter/pss-211",
        "waveform_sha256": "1" * 64, "accepted_waveform_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "small-signal-211"
    positive[
        "bias_small_signal_jacobian_ac_transient_pole_zero_normalization_circuit_result_identity"
    ] = {
        "small_signal_generation_id": generation,
        **{
            key: generation
            for key in (
                "bias_small_signal_generation_id", "jacobian_small_signal_generation_id",
                "ac_small_signal_generation_id", "transient_small_signal_generation_id",
                "pole_zero_small_signal_generation_id", "normalization_small_signal_generation_id",
                "circuit_small_signal_generation_id", "result_small_signal_generation_id",
            )
        },
        "bias_state": [1.0], "result_bias_state": [1.0],
        "jacobian": [[-1.0]], "result_jacobian": [[-1.0]],
        "frequency_rad_s": [0.0, 1.0], "result_frequency_rad_s": [0.0, 1.0],
        "ac_transfer_ri": [[1.0, 0.0], [0.5, -0.5]],
        "result_ac_transfer_ri": [[1.0, 0.0], [0.5, -0.5]],
        "time_s": [0.0, 1.0], "result_time_s": [0.0, 1.0],
        "impulse_response": [1.0, math.exp(-1)],
        "result_impulse_response": [1.0, math.exp(-1)],
        "step_response": [0.0, 1.0 - math.exp(-1)],
        "result_step_response": [0.0, 1.0 - math.exp(-1)],
        "poles_ri": [[-1.0, 0.0]], "result_poles_ri": [[-1.0, 0.0]],
        "zeros_ri": [], "result_zeros_ri": [],
        "normalization": "monic_denominator_unit_dc_gain",
        "result_normalization": "monic_denominator_unit_dc_gain",
        "circuit_sha256": "3" * 64, "linearized_circuit_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v34_public_positive_periodic_and_small_signal_closure():
    assert ideal_transformer_identity_gate(_v34())["status"] == "ok"


def test_v34_public_periodic_steady_state_charge_flux_cycle_energy_efficiency_phase_owner_mismatch():
    payload = _v34()
    identity = payload["metrics"]["positive"][
        "periodic_steady_state_charge_flux_cycle_energy_efficiency_phase_waveform_owner_result_identity"
    ]
    identity.update({
        "charge_periodic_generation_id": "periodic-steady-state-210",
        "energy_periodic_generation_id": "periodic-steady-state-209",
        "result_periodic_generation_id": "periodic-steady-state-208",
        "result_cycle_end_charge_c": [2e-3, 0.0], "result_cycle_end_flux_wb": [5e-4],
        "result_closure_tolerance": 1e-2, "result_input_energy_j": 5.0,
        "result_output_energy_j": 8.0, "result_loss_energy_j": -3.0,
        "result_efficiency": 1.6, "result_phase_window_rad": [1.0, 2.0],
        "accepted_waveform_owner": "converter/old", "accepted_waveform_sha256": "8" * 64,
        "accepted_result_sha256": "9" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "periodic_steady_state_uses_current_charge_flux_energy_efficiency_phase_waveform_owner_and_result"
    ]


def test_v34_public_bias_small_signal_ac_transient_jacobian_pole_zero_transfer_consistency_mismatch():
    payload = _v34()
    identity = payload["metrics"]["positive"][
        "bias_small_signal_jacobian_ac_transient_pole_zero_normalization_circuit_result_identity"
    ]
    identity.update({
        "bias_small_signal_generation_id": "small-signal-210",
        "ac_small_signal_generation_id": "small-signal-209",
        "result_small_signal_generation_id": "small-signal-208",
        "result_bias_state": [2.0], "result_jacobian": [[1.0]],
        "result_frequency_rad_s": [1.0, 0.0], "result_ac_transfer_ri": [[2.0, 1.0]],
        "result_time_s": [1.0, 0.0], "result_impulse_response": [-1.0, 2.0],
        "result_step_response": [1.0, -1.0], "result_poles_ri": [[1.0, 0.0]],
        "result_zeros_ri": [[0.0, 0.0]], "result_normalization": "arbitrary",
        "linearized_circuit_sha256": "a" * 64, "accepted_result_sha256": "b" * 64,
    })
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "small_signal_uses_current_bias_jacobian_ac_transient_poles_zeros_normalization_circuit_and_result"
    ]


def test_v34_public_rejects_self_consistent_cycle_energy_gain():
    payload = _v34()
    identity = payload["metrics"]["positive"][
        "periodic_steady_state_charge_flux_cycle_energy_efficiency_phase_waveform_owner_result_identity"
    ]
    identity["output_energy_j"] = identity["result_output_energy_j"] = 12.0
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_unstable_small_signal_pole():
    payload = _v34()
    identity = payload["metrics"]["positive"][
        "bias_small_signal_jacobian_ac_transient_pole_zero_normalization_circuit_result_identity"
    ]
    identity["poles_ri"] = identity["result_poles_ri"] = [[1.0, 0.0]]
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
