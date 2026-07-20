from copy import deepcopy
import math

from ltspice_converter.ltspice_v54_gates import LOOP, MONTE, validate_ltspice_v54_identity


PROMOTED_CASE_IDS = {
    "v54_public_montecarlo_seed_distribution_parameter_yield_sample_owner_mismatch",
    "v54_public_loopgain_injection_breakpoint_sign_crossover_phasemargin_owner_mismatch",
}


def _positive() -> dict[str, object]:
    monte_generation = "monte-public-v54"
    loop_generation = "loop-public-v54"
    distributions = {"R1": {"distribution": "normal", "mean": 1000.0, "stddev": 10.0}, "C1": {"distribution": "uniform", "lower": 0.95e-6, "upper": 1.05e-6}}
    samples = [{"sample_id": 0, "parameters": {"R1": 998.0, "C1": 0.99e-6}, "passed": True}, {"sample_id": 1, "parameters": {"R1": 1035.0, "C1": 1.04e-6}, "passed": False}]
    rows = [{"frequency_hz": 1.0e3, "real": -0.5, "imag": -1.5}, {"frequency_hz": 1.0e4, "real": -math.sqrt(0.5), "imag": -math.sqrt(0.5)}, {"frequency_hz": 1.0e5, "real": -0.1, "imag": -0.2}]
    return {
        MONTE: {
            "generation_id": monte_generation, **{field: monte_generation for field in ("seed_generation_id", "distribution_generation_id", "sample_generation_id", "yield_generation_id", "owner_generation_id", "result_generation_id")},
            "rng_seed": 314159, "result_rng_seed": 314159,
            "parameter_distributions": distributions, "result_parameter_distributions": distributions,
            "sampled_values": samples, "result_sampled_values": samples,
            "yield_count": 1, "result_yield_count": 1,
            "run_owner": "run:v54", "result_run_owner": "run:v54",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        LOOP: {
            "generation_id": loop_generation, **{field: loop_generation for field in ("injection_generation_id", "sign_generation_id", "crossover_generation_id", "margin_generation_id", "owner_generation_id", "result_generation_id")},
            "injection_point": "node:loop-break", "result_injection_point": "node:loop-break",
            "sign_convention": "negative_feedback_return_ratio", "result_sign_convention": "negative_feedback_return_ratio",
            "loop_gain_rows": rows, "result_loop_gain_rows": rows,
            "crossover_frequency_hz": 1.0e4, "result_crossover_frequency_hz": 1.0e4,
            "phase_margin_deg": 45.0, "result_phase_margin_deg": 45.0,
            "trace_owner": "trace:v54", "result_trace_owner": "trace:v54",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v54_positive_identity_is_accepted() -> None:
    assert validate_ltspice_v54_identity(_positive()) is True


def test_v54_frozen_public_counterfactuals_are_rejected() -> None:
    value = deepcopy(_positive())
    value[MONTE]["result_rng_seed"] = 42
    value[LOOP]["result_phase_margin_deg"] = -45.0
    assert validate_ltspice_v54_identity(value) is False


def test_v54_self_consistent_invalid_semantics_are_rejected() -> None:
    value = deepcopy(_positive())
    value[MONTE]["yield_count"] = value[MONTE]["result_yield_count"] = 2
    value[LOOP]["phase_margin_deg"] = value[LOOP]["result_phase_margin_deg"] = 90.0
    assert validate_ltspice_v54_identity(value) is False


def test_v54_malformed_values_reject_without_raising() -> None:
    value = deepcopy(_positive())
    value[MONTE]["sampled_values"] = [{"sample_id": [0]}]
    value[LOOP]["loop_gain_rows"] = [{"frequency_hz": [1.0e3]}]
    assert validate_ltspice_v54_identity(value) is False
