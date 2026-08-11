from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from radia.ltspice.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v29 import _v29


_PROMOTED_CASE_IDS = (
    "v30_public_mosfet_soa_vds_id_pulse_width_duty_temperature_model_waveform_mismatch",
    "v30_public_monte_carlo_yield_distribution_tolerance_seed_failure_criteria_sample_owner_mismatch",
)


def _v30():
    summary = _v29()
    positive = summary["metrics"]["positive"]
    generation = "mosfet-soa-171"
    positive[
        "mosfet_soa_vds_id_pulse_width_duty_temperature_model_waveform_result_generation_identity"
    ] = {
        "soa_generation_id": generation,
        "voltage_soa_generation_id": generation,
        "current_soa_generation_id": generation,
        "pulse_soa_generation_id": generation,
        "duty_soa_generation_id": generation,
        "temperature_soa_generation_id": generation,
        "model_soa_generation_id": generation,
        "waveform_soa_generation_id": generation,
        "result_soa_generation_id": generation,
        "vds_v": 400.0,
        "result_vds_v": 400.0,
        "id_a": 20.0,
        "result_id_a": 20.0,
        "pulse_width_s": 1.0e-4,
        "result_pulse_width_s": 1.0e-4,
        "repetition_period_s": 1.0e-2,
        "result_repetition_period_s": 1.0e-2,
        "duty_cycle": 0.01,
        "result_duty_cycle": 0.01,
        "junction_temperature_c": 125.0,
        "result_junction_temperature_c": 125.0,
        "soa_limit_id_a": 25.0,
        "result_soa_limit_id_a": 25.0,
        "soa_margin_fraction": 0.2,
        "result_soa_margin_fraction": 0.2,
        "model_card_sha256": "1" * 64,
        "result_model_card_sha256": "1" * 64,
        "waveform_sha256": "2" * 64,
        "result_waveform_sha256": "2" * 64,
        "soa_result_sha256": "3" * 64,
        "accepted_soa_result_sha256": "3" * 64,
    }
    generation = "monte-carlo-yield-171"
    positive[
        "monte_carlo_yield_distribution_tolerance_seed_failure_sample_owner_result_generation_identity"
    ] = {
        "yield_generation_id": generation,
        "distribution_yield_generation_id": generation,
        "tolerance_yield_generation_id": generation,
        "seed_yield_generation_id": generation,
        "criterion_yield_generation_id": generation,
        "sample_yield_generation_id": generation,
        "owner_yield_generation_id": generation,
        "result_yield_generation_id": generation,
        "parameter_order": ["R1", "C1"],
        "result_parameter_order": ["R1", "C1"],
        "distribution_families": ["gaussian", "uniform"],
        "result_distribution_families": ["gaussian", "uniform"],
        "nominal_values": [1000.0, 1.0e-7],
        "result_nominal_values": [1000.0, 1.0e-7],
        "relative_tolerances": [0.01, 0.05],
        "result_relative_tolerances": [0.01, 0.05],
        "seed_schedule": [101, 103, 107, 109, 113],
        "result_seed_schedule": [101, 103, 107, 109, 113],
        "failure_criterion": "vout_min_lt_4p75",
        "result_failure_criterion": "vout_min_lt_4p75",
        "sample_ids": [0, 1, 2, 3, 4],
        "result_sample_ids": [0, 1, 2, 3, 4],
        "failed_sample_ids": [1],
        "result_failed_sample_ids": [1],
        "accepted_sample_ids": [0, 2, 3, 4],
        "result_accepted_sample_ids": [0, 2, 3, 4],
        "yield_fraction": 0.8,
        "result_yield_fraction": 0.8,
        "circuit_owner_sha256": "4" * 64,
        "result_circuit_owner_sha256": "4" * 64,
        "sample_table_sha256": "5" * 64,
        "result_sample_table_sha256": "5" * 64,
        "yield_result_sha256": "6" * 64,
        "accepted_yield_result_sha256": "6" * 64,
    }
    return summary


def test_v30_positive_mosfet_soa_and_monte_carlo_yield_identities():
    assert ideal_transformer_identity_gate(_v30())["status"] == "ok"


def test_v30_public_mosfet_soa_vds_id_pulse_width_duty_temperature_model_waveform_mismatch():
    summary = _v30()
    contract = summary["metrics"]["positive"][
        "mosfet_soa_vds_id_pulse_width_duty_temperature_model_waveform_result_generation_identity"
    ]
    contract.update(
        {
            "voltage_soa_generation_id": "mosfet-soa-170",
            "model_soa_generation_id": "mosfet-soa-169",
            "result_soa_generation_id": "mosfet-soa-168",
            "result_vds_v": 500.0,
            "result_id_a": 30.0,
            "result_pulse_width_s": 1.0e-3,
            "result_repetition_period_s": 5.0e-3,
            "result_duty_cycle": 0.2,
            "result_junction_temperature_c": 25.0,
            "result_soa_limit_id_a": 15.0,
            "result_soa_margin_fraction": -1.0,
            "result_model_card_sha256": "c" * 64,
            "result_waveform_sha256": "d" * 64,
            "accepted_soa_result_sha256": "e" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mosfet_soa_uses_current_voltage_current_pulse_duty_temperature_model_waveform_and_result"
    ]


def test_v30_public_monte_carlo_yield_distribution_tolerance_seed_failure_criteria_sample_owner_mismatch():
    summary = _v30()
    contract = summary["metrics"]["positive"][
        "monte_carlo_yield_distribution_tolerance_seed_failure_sample_owner_result_generation_identity"
    ]
    contract.update(
        {
            "distribution_yield_generation_id": "monte-carlo-yield-170",
            "sample_yield_generation_id": "monte-carlo-yield-169",
            "result_yield_generation_id": "monte-carlo-yield-168",
            "result_parameter_order": ["C1", "R1"],
            "result_distribution_families": ["uniform", "gaussian"],
            "result_nominal_values": [1.0e-7, 1000.0],
            "result_relative_tolerances": [0.1, 0.2],
            "result_seed_schedule": [101, 101, 101],
            "result_failure_criterion": "vout_min_lt_4p0",
            "result_sample_ids": [0, 1, 2],
            "result_failed_sample_ids": [0, 2],
            "result_accepted_sample_ids": [1, 4],
            "result_yield_fraction": 0.4,
            "result_circuit_owner_sha256": "f" * 64,
            "result_sample_table_sha256": "0" * 64,
            "accepted_yield_result_sha256": "1" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "monte_carlo_yield_uses_current_distributions_tolerances_seeds_failure_samples_owner_and_result"
    ]
