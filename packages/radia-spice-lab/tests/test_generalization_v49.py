from copy import deepcopy

from ltspice_converter.ltspice_v49_gates import MONTE_CARLO, SWITCH, validate_ltspice_v49_identity


PROMOTED_CASE_IDS = (
    "v49_public_monte_carlo_seed_distribution_sample_parameter_run_order_measure_owner_mismatch",
    "v49_public_switch_hysteresis_state_timestep_breakpoint_raw_trace_owner_mismatch",
)


def _positive():
    monte = "monte-v49"; switch = "switch-v49"
    samples = [-0.8, -0.1, 0.4, 1.2]; parameters = [9.2e3, 9.9e3, 10.4e3, 11.2e3]; order = [f"run:{index}" for index in range(4)]
    times = [0.0, 1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6]; states = ["off", "off", "on", "on", "off"]; breakpoints = [2.0e-6, 4.0e-6]
    trace = [[time, value] for time, value in zip(times, [0.0, 0.2, 1.0, 0.8, 0.0])]
    return {
        MONTE_CARLO: {
            "generation_id": monte, **{key: monte for key in ("seed_generation_id", "distribution_generation_id", "sample_generation_id", "parameter_generation_id", "run_generation_id", "measure_generation_id", "result_generation_id")},
            "random_seed": 4901, "result_random_seed": 4901, "distribution": {"kind": "normal", "mean": 0.0, "stddev": 1.0}, "result_distribution": {"kind": "normal", "mean": 0.0, "stddev": 1.0},
            "sample_values": samples, "result_sample_values": samples, "parameter_values": parameters, "result_parameter_values": parameters, "run_order": order, "result_run_order": order,
            "measure_values": [0.92, 0.99, 1.04, 1.12], "result_measure_values": [0.92, 0.99, 1.04, 1.12], "measure_owner": "measure:test", "result_measure_owner": "measure:test",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        SWITCH: {
            "generation_id": switch, **{key: switch for key in ("hysteresis_generation_id", "state_generation_id", "timestep_generation_id", "breakpoint_generation_id", "raw_generation_id", "trace_generation_id", "owner_generation_id", "result_generation_id")},
            "hysteresis_state": states, "result_hysteresis_state": states, "accepted_timesteps_s": times, "result_accepted_timesteps_s": times,
            "breakpoints_s": breakpoints, "result_breakpoints_s": breakpoints, "raw_trace_rows": trace, "result_raw_trace_rows": trace,
            "waveform_owner": "waveform:test", "result_waveform_owner": "waveform:test", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }


def test_v49_public_replay_identity(): assert validate_ltspice_v49_identity(_positive()) is True


def test_v49_public_rejects_monte_carlo_mutation():
    value = deepcopy(_positive()); value[MONTE_CARLO].update({"result_random_seed": 4902, "result_run_order": ["run:3", "run:2", "run:1", "run:0"], "result_measure_owner": "measure:old"})
    assert validate_ltspice_v49_identity(value) is False


def test_v49_public_rejects_switch_mutation():
    value = deepcopy(_positive()); value[SWITCH].update({"result_hysteresis_state": ["off", "on", "off", "on", "off"], "result_breakpoints_s": [1.0e-6, 3.0e-6], "result_waveform_owner": "waveform:old"})
    assert validate_ltspice_v49_identity(value) is False


def test_v49_public_rejects_self_consistent_nonphysical_contracts():
    value = deepcopy(_positive()); value[MONTE_CARLO]["distribution"] = value[MONTE_CARLO]["result_distribution"] = {"kind": "normal", "mean": 0.0, "stddev": 0.0}
    value[SWITCH]["breakpoints_s"] = value[SWITCH]["result_breakpoints_s"] = [1.0e-6, 3.0e-6]
    assert validate_ltspice_v49_identity(value) is False
