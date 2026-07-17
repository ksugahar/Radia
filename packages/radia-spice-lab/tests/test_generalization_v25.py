from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v24 import _v24


_PROMOTED_CASE_IDS = (
    "v25_public_transient_startup_initial_condition_uic_operating_point_waveform_generation_mismatch",
    "v25_public_stepped_monte_carlo_measure_aggregation_failed_row_seed_weight_mismatch",
)


def _v25():
    summary = _v24()
    positive = summary["metrics"]["positive"]
    positive[
        "transient_startup_initial_condition_uic_operating_point_waveform_generation_identity"
    ] = {
        "transient_generation_id": "startup-transient-121",
        "startup_mode_transient_generation_id": "startup-transient-121",
        "initial_condition_transient_generation_id": "startup-transient-121",
        "uic_transient_generation_id": "startup-transient-121",
        "operating_point_transient_generation_id": "startup-transient-121",
        "accepted_time_grid_transient_generation_id": "startup-transient-121",
        "waveform_transient_generation_id": "startup-transient-121",
        "result_transient_generation_id": "startup-transient-121",
        "startup_mode": "operating_point_then_transient",
        "result_startup_mode": "operating_point_then_transient",
        "initial_conditions": [["V(out)", 0.0], ["I(L1)", 0.0]],
        "result_initial_conditions": [["V(out)", 0.0], ["I(L1)", 0.0]],
        "uic_enabled": False,
        "result_uic_enabled": False,
        "operating_point_sha256": "3" * 64,
        "result_operating_point_sha256": "3" * 64,
        "accepted_time_s": [0.0, 1.0e-6, 2.0e-6, 4.0e-6],
        "result_time_s": [0.0, 1.0e-6, 2.0e-6, 4.0e-6],
        "waveform_trace_ids": ["V(out)", "I(L1)"],
        "result_waveform_trace_ids": ["V(out)", "I(L1)"],
        "waveform_table_sha256": "4" * 64,
        "result_waveform_table_sha256": "4" * 64,
    }
    positive[
        "stepped_monte_carlo_measure_aggregation_failed_row_seed_weight_generation_identity"
    ] = {
        "monte_carlo_generation_id": "stepped-mc-121",
        "measure_row_monte_carlo_generation_id": "stepped-mc-121",
        "seed_monte_carlo_generation_id": "stepped-mc-121",
        "filter_monte_carlo_generation_id": "stepped-mc-121",
        "weight_monte_carlo_generation_id": "stepped-mc-121",
        "aggregation_monte_carlo_generation_id": "stepped-mc-121",
        "result_monte_carlo_generation_id": "stepped-mc-121",
        "sample_ids": [1, 2, 3, 4],
        "result_sample_ids": [1, 2, 3, 4],
        "random_seeds": [101, 202, 303, 404],
        "result_random_seeds": [101, 202, 303, 404],
        "measure_statuses": ["passed", "passed", "passed", "passed"],
        "accepted_sample_ids": [1, 2, 3, 4],
        "result_accepted_sample_ids": [1, 2, 3, 4],
        "failed_sample_ids": [],
        "result_failed_sample_ids": [],
        "sample_values": [0.9, 1.0, 1.1, 1.2],
        "sample_weights": [1.0, 2.0, 2.0, 1.0],
        "result_sample_weights": [1.0, 2.0, 2.0, 1.0],
        "aggregation_rule": "weighted_mean",
        "result_aggregation_rule": "weighted_mean",
        "reported_weighted_mean": 1.05,
        "sample_table_sha256": "5" * 64,
        "result_sample_table_sha256": "5" * 64,
    }
    return summary


def test_v25_positive():
    assert ideal_transformer_identity_gate(_v25())["status"] == "ok"


def test_v25_public_transient_startup_initial_condition_uic_operating_point_waveform_generation_mismatch():
    summary = _v25()
    contract = summary["metrics"]["positive"][
        "transient_startup_initial_condition_uic_operating_point_waveform_generation_identity"
    ]
    contract.update(
        {
            "initial_condition_transient_generation_id": "startup-transient-120",
            "result_startup_mode": "uic_with_explicit_ic",
            "result_uic_enabled": True,
            "result_time_s": [0.0, 1.5e-6, 2.0e-6, 4.0e-6],
            "result_waveform_table_sha256": "d" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "transient_startup_uses_current_initial_conditions_uic_operating_point_grid_and_waveforms"
    ]


def test_v25_public_stepped_monte_carlo_measure_aggregation_failed_row_seed_weight_mismatch():
    summary = _v25()
    contract = summary["metrics"]["positive"][
        "stepped_monte_carlo_measure_aggregation_failed_row_seed_weight_generation_identity"
    ]
    contract.update(
        {
            "seed_monte_carlo_generation_id": "stepped-mc-120",
            "result_random_seeds": [101, 999, 303, 404],
            "measure_statuses": ["passed", "failed", "passed", "passed"],
            "failed_sample_ids": [2],
            "result_failed_sample_ids": [],
            "result_sample_weights": [1.0, 1.0, 1.0, 1.0],
            "result_aggregation_rule": "unweighted_mean",
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "stepped_monte_carlo_aggregation_uses_current_rows_seeds_filters_weights_and_rule"
    ]
