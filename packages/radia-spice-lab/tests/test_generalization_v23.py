from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_ideal_transformer_gate import _summary, _with_v22_switching_and_electrothermal_identity

def _v23():
    s = _with_v22_switching_and_electrothermal_identity(_summary()); p = s["metrics"]["positive"]
    p["subcircuit_monte_carlo_seed_model_include_raw_generation_identity"] = {
        "monte_carlo_generation_id":"mc-81","seed_monte_carlo_generation_id":"mc-81","model_include_monte_carlo_generation_id":"mc-81","parameter_override_monte_carlo_generation_id":"mc-81","raw_sample_monte_carlo_generation_id":"mc-81","statistic_monte_carlo_generation_id":"mc-81",
        "sample_ids":[1,2,3],"raw_sample_ids":[1,2,3],"random_seeds":[101,202,303],"raw_random_seeds":[101,202,303],"model_include_sha256":"1"*64,"raw_model_include_sha256":"1"*64,"parameter_override_sha256":"2"*64,"raw_parameter_override_sha256":"2"*64,"raw_sample_table_sha256":"3"*64,"statistic_sample_table_sha256":"3"*64}
    p["behavioral_switch_hysteresis_state_timestep_measure_generation_identity"] = {
        "transient_generation_id":"switch-hyst-81","hysteresis_state_transient_generation_id":"switch-hyst-81","event_history_transient_generation_id":"switch-hyst-81","accepted_timestep_transient_generation_id":"switch-hyst-81","measure_window_transient_generation_id":"switch-hyst-81","measure_result_transient_generation_id":"switch-hyst-81",
        "hysteresis_states":[0,1,1,0],"measure_hysteresis_states":[0,1,1,0],"event_times_s":[0,1e-6,3e-6,4e-6],"measure_event_times_s":[0,1e-6,3e-6,4e-6],"accepted_time_s":[0,1e-6,2e-6,3e-6,4e-6],"measure_time_s":[0,1e-6,2e-6,3e-6,4e-6],"measure_window_s":[1e-6,4e-6],"reported_measure_window_s":[1e-6,4e-6],"measure_table_sha256":"4"*64,"reported_measure_table_sha256":"4"*64}
    return s

def test_v23_positive(): assert ideal_transformer_identity_gate(_v23())["status"] == "ok"

def test_v23_public_subcircuit_monte_carlo_seed_model_include_raw_generation_mismatch():
    s=_v23(); c=s["metrics"]["positive"]["subcircuit_monte_carlo_seed_model_include_raw_generation_identity"]
    c.update({"seed_monte_carlo_generation_id":"mc-80","raw_sample_ids":[1,3,4],"raw_random_seeds":[101,999,404],"raw_model_include_sha256":"a"*64})
    r=ideal_transformer_identity_gate(s); assert r["status"]=="needs_attention"; assert not r["checks"]["monte_carlo_statistics_use_current_seeds_models_parameters_and_raw_samples"]

def test_v23_public_behavioral_switch_hysteresis_state_timestep_measure_generation_mismatch():
    s=_v23(); c=s["metrics"]["positive"]["behavioral_switch_hysteresis_state_timestep_measure_generation_identity"]
    c.update({"hysteresis_state_transient_generation_id":"switch-hyst-80","measure_hysteresis_states":[0,1,0,1],"measure_time_s":[0,2e-6,1e-6,3e-6,4e-6],"reported_measure_table_sha256":"d"*64})
    r=ideal_transformer_identity_gate(s); assert r["status"]=="needs_attention"; assert not r["checks"]["behavioral_switch_measures_use_current_hysteresis_events_timesteps_and_window"]
