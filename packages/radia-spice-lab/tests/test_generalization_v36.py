from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v35 import _v35

_PROMOTED_CASE_IDS = (
    "v36_public_pll_lock_phase_noise_jitter_loop_gain_waveform_owner_mismatch",
    "v36_public_electrothermal_runaway_thermal_impedance_power_temperature_fixedpoint_mismatch",
)


def _v36():
    payload = _v35()
    positive = payload["metrics"]["positive"]
    generation = "pll-closure-222"
    positive[
        "pll_lock_phase_noise_jitter_loop_gain_waveform_circuit_owner_result_identity"
    ] = {
        "pll_generation_id": generation,
        **{
            key: generation
            for key in (
                "lock_pll_generation_id",
                "phase_error_pll_generation_id",
                "noise_pll_generation_id",
                "jitter_pll_generation_id",
                "loop_gain_pll_generation_id",
                "waveform_pll_generation_id",
                "circuit_pll_generation_id",
                "result_pll_generation_id",
            )
        },
        "lock_state": "locked",
        "result_lock_state": "locked",
        "lock_time_s": 2e-4,
        "result_lock_time_s": 2e-4,
        "waveform_duration_s": 1e-3,
        "result_waveform_duration_s": 1e-3,
        "phase_error_start_rad": 0.4,
        "result_phase_error_start_rad": 0.4,
        "phase_error_end_rad": 2e-4,
        "result_phase_error_end_rad": 2e-4,
        "phase_error_tolerance_rad": 1e-3,
        "result_phase_error_tolerance_rad": 1e-3,
        "phase_noise_offset_hz": [1e3, 1e4, 1e5],
        "result_phase_noise_offset_hz": [1e3, 1e4, 1e5],
        "phase_noise_dbc_per_hz": [-80.0, -100.0, -120.0],
        "result_phase_noise_dbc_per_hz": [-80.0, -100.0, -120.0],
        "integrated_jitter_s": 2e-12,
        "result_integrated_jitter_s": 2e-12,
        "loop_gain_crossover_hz": 2e4,
        "result_loop_gain_crossover_hz": 2e4,
        "loop_gain_magnitude_at_crossover": 1.0,
        "result_loop_gain_magnitude_at_crossover": 1.0,
        "phase_margin_deg": 55.0,
        "result_phase_margin_deg": 55.0,
        "waveform_owner": "pll/waveform-222",
        "accepted_waveform_owner": "pll/waveform-222",
        "waveform_sha256": "1" * 64,
        "accepted_waveform_sha256": "1" * 64,
        "circuit_sha256": "2" * 64,
        "accepted_circuit_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }

    generation = "electrothermal-222"
    reference_power = 38.0
    coefficient = (40.0 / reference_power - 1.0) / 100.0
    positive[
        "electrothermal_thermal_impedance_power_temperature_device_fixedpoint_stability_circuit_result_identity"
    ] = {
        "electrothermal_generation_id": generation,
        **{
            key: generation
            for key in (
                "thermal_impedance_generation_id",
                "power_generation_id",
                "temperature_generation_id",
                "device_law_generation_id",
                "fixedpoint_generation_id",
                "stability_generation_id",
                "circuit_generation_id",
                "result_generation_id",
            )
        },
        "thermal_impedance_k_per_w": 2.5,
        "result_thermal_impedance_k_per_w": 2.5,
        "dissipated_power_w": 40.0,
        "result_dissipated_power_w": 40.0,
        "ambient_temperature_c": 25.0,
        "result_ambient_temperature_c": 25.0,
        "junction_temperature_c": 125.0,
        "result_junction_temperature_c": 125.0,
        "device_law": "linear_temperature_power",
        "result_device_law": "linear_temperature_power",
        "reference_power_w": reference_power,
        "result_reference_power_w": reference_power,
        "reference_temperature_c": 25.0,
        "result_reference_temperature_c": 25.0,
        "power_temperature_coefficient_per_k": coefficient,
        "result_power_temperature_coefficient_per_k": coefficient,
        "fixedpoint_residual_c": 0.0,
        "result_fixedpoint_residual_c": 0.0,
        "fixedpoint_tolerance_c": 1e-8,
        "result_fixedpoint_tolerance_c": 1e-8,
        "stability_slope": 0.05,
        "result_stability_slope": 0.05,
        "circuit_owner": "electrothermal/circuit-222",
        "accepted_circuit_owner": "electrothermal/circuit-222",
        "circuit_sha256": "4" * 64,
        "accepted_circuit_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return payload


def test_v36_public_positive_pll_and_electrothermal_closure():
    assert ideal_transformer_identity_gate(_v36())["status"] == "ok"


def test_v36_public_pll_lock_phase_noise_jitter_loop_gain_waveform_owner_mismatch():
    payload = _v36()
    identity = payload["metrics"]["positive"][
        "pll_lock_phase_noise_jitter_loop_gain_waveform_circuit_owner_result_identity"
    ]
    identity.update(
        {
            "lock_pll_generation_id": "pll-closure-221",
            "result_lock_state": "unlocked",
            "result_phase_error_end_rad": 0.2,
            "result_phase_noise_offset_hz": [1e5, 1e3],
            "result_integrated_jitter_s": -1.0,
            "result_loop_gain_magnitude_at_crossover": 2.0,
            "accepted_waveform_sha256": "a" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pll_uses_current_lock_phase_error_noise_jitter_loop_gain_waveform_circuit_and_result"
    ]


def test_v36_public_electrothermal_runaway_thermal_impedance_power_temperature_fixedpoint_mismatch():
    payload = _v36()
    identity = payload["metrics"]["positive"][
        "electrothermal_thermal_impedance_power_temperature_device_fixedpoint_stability_circuit_result_identity"
    ]
    identity.update(
        {
            "power_generation_id": "electrothermal-221",
            "result_thermal_impedance_k_per_w": -2.5,
            "result_junction_temperature_c": 25.0,
            "result_device_law": "constant_power",
            "result_fixedpoint_residual_c": 100.0,
            "result_stability_slope": 1.2,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = ideal_transformer_identity_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrothermal_uses_current_thermal_impedance_power_temperature_device_fixedpoint_stability_circuit_and_result"
    ]


def test_v36_public_rejects_self_consistent_unlocked_pll():
    payload = _v36()
    identity = payload["metrics"]["positive"][
        "pll_lock_phase_noise_jitter_loop_gain_waveform_circuit_owner_result_identity"
    ]
    identity["lock_state"] = identity["result_lock_state"] = "unlocked"
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_electrothermal_runaway_slope():
    payload = _v36()
    identity = payload["metrics"]["positive"][
        "electrothermal_thermal_impedance_power_temperature_device_fixedpoint_stability_circuit_result_identity"
    ]
    identity["stability_slope"] = identity["result_stability_slope"] = 1.2
    assert ideal_transformer_identity_gate(payload)["status"] == "needs_attention"
