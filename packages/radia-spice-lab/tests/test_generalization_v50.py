from copy import deepcopy

from ltspice_converter.ltspice_v50_gates import NOISE, STEP, validate_ltspice_v50_identity


PROMOTED_CASE_IDS = (
    "v50_public_step_parameter_cartesian_nested_order_measure_row_owner_mismatch",
    "v50_public_noise_input_output_source_contribution_bandwidth_integration_owner_mismatch",
)


def _positive():
    step_generation = "step-v50-test"
    noise_generation = "noise-v50-test"
    names = ["RLOAD", "CLOAD"]
    grid = {"RLOAD": [10000.0, 20000.0], "CLOAD": [1.0e-9, 2.0e-9]}
    rows = [
        ["RLOAD", 10000.0, "CLOAD", 1.0e-9],
        ["RLOAD", 10000.0, "CLOAD", 2.0e-9],
        ["RLOAD", 20000.0, "CLOAD", 1.0e-9],
        ["RLOAD", 20000.0, "CLOAD", 2.0e-9],
    ]
    keys = [
        "RLOAD=10000;CLOAD=1e-09",
        "RLOAD=10000;CLOAD=2e-09",
        "RLOAD=20000;CLOAD=1e-09",
        "RLOAD=20000;CLOAD=2e-09",
    ]
    frequencies = [10.0, 100.0, 1000.0, 10000.0]
    sources = ["R1", "EAMP"]
    contributions = {
        "R1": [1.0e-18, 8.0e-19, 6.0e-19, 4.0e-19],
        "EAMP": [2.0e-18, 1.8e-18, 1.2e-18, 8.0e-19],
    }
    return {
        STEP: {
            "generation_id": step_generation,
            **{key: step_generation for key in ("parameter_generation_id", "cartesian_generation_id", "nesting_generation_id", "order_generation_id", "measure_generation_id", "owner_generation_id", "result_generation_id")},
            "parameter_names": names,
            "result_parameter_names": names,
            "parameter_value_grid": grid,
            "result_parameter_value_grid": grid,
            "nesting_order": names,
            "result_nesting_order": names,
            "cartesian_step_rows": rows,
            "result_cartesian_step_rows": rows,
            "measure_row_keys": keys,
            "result_measure_row_keys": keys,
            "measure_values": [0.1, 0.2, 0.3, 0.4],
            "result_measure_values": [0.1, 0.2, 0.3, 0.4],
            "sweep_owner": "sweep:test",
            "result_sweep_owner": "sweep:test",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        NOISE: {
            "generation_id": noise_generation,
            **{key: noise_generation for key in ("reference_generation_id", "source_generation_id", "frequency_generation_id", "bandwidth_generation_id", "integration_generation_id", "owner_generation_id", "result_generation_id")},
            "input_reference": "V(in)",
            "result_input_reference": "V(in)",
            "output_reference": "V(out)",
            "result_output_reference": "V(out)",
            "frequency_hz": frequencies,
            "result_frequency_hz": frequencies,
            "noise_source_order": sources,
            "result_noise_source_order": sources,
            "source_contribution_v2_per_hz": contributions,
            "result_source_contribution_v2_per_hz": contributions,
            "integration_band_hz": [10.0, 10000.0],
            "result_integration_band_hz": [10.0, 10000.0],
            "integrated_output_noise_v_rms": 2.5e-7,
            "result_integrated_output_noise_v_rms": 2.5e-7,
            "trace_owner": "trace:test",
            "result_trace_owner": "trace:test",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
    }


def test_v50_public_replay_identity():
    assert validate_ltspice_v50_identity(_positive()) is True


def test_v50_public_rejects_step_mutation():
    value = deepcopy(_positive())
    value[STEP].update({"result_nesting_order": ["CLOAD", "RLOAD"], "result_cartesian_step_rows": list(reversed(value[STEP]["cartesian_step_rows"])), "result_sweep_owner": "sweep:foreign"})
    assert validate_ltspice_v50_identity(value) is False


def test_v50_public_rejects_noise_mutation():
    value = deepcopy(_positive())
    value[NOISE].update({"result_frequency_hz": list(reversed(value[NOISE]["frequency_hz"])), "result_noise_source_order": ["EAMP", "R1"], "result_trace_owner": "trace:foreign"})
    assert validate_ltspice_v50_identity(value) is False


def test_v50_public_rejects_self_consistent_nonphysical_contracts():
    value = deepcopy(_positive())
    value[STEP]["cartesian_step_rows"] = value[STEP]["result_cartesian_step_rows"] = value[STEP]["cartesian_step_rows"][:-1]
    value[STEP]["measure_row_keys"] = value[STEP]["result_measure_row_keys"] = value[STEP]["measure_row_keys"][:-1]
    value[STEP]["measure_values"] = value[STEP]["result_measure_values"] = value[STEP]["measure_values"][:-1]
    value[NOISE]["frequency_hz"] = value[NOISE]["result_frequency_hz"] = list(reversed(value[NOISE]["frequency_hz"]))
    assert validate_ltspice_v50_identity(value) is False
