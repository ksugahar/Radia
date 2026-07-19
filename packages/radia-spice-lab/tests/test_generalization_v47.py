from copy import deepcopy

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from ltspice_converter.ltspice_v47_gates import validate_ltspice_v47_identity


PROMOTED_CASE_IDS = (
    "v47_public_hierarchical_branch_current_sign_order_kcl_owner_mismatch",
    "v47_public_ac_tran_step_tuple_trace_measure_row_key_mismatch",
)


def _positive():
    return {
        "hierarchical_branch_v47_current_sign_order_kcl_owner_identity": {
            "generation_id": "hierarchy-v47",
            "hierarchy_generation_id": "hierarchy-v47",
            "current_generation_id": "hierarchy-v47",
            "kcl_generation_id": "hierarchy-v47",
            "result_generation_id": "hierarchy-v47",
            "hierarchy_paths": ["XTOP:XAMP:RIN", "XTOP:XAMP:CIN", "XTOP:RLOAD"],
            "result_hierarchy_paths": ["XTOP:XAMP:RIN", "XTOP:XAMP:CIN", "XTOP:RLOAD"],
            "branch_current_order": ["I(XTOP:XAMP:RIN)", "I(XTOP:XAMP:CIN)", "I(XTOP:RLOAD)"],
            "result_branch_current_order": ["I(XTOP:XAMP:RIN)", "I(XTOP:XAMP:CIN)", "I(XTOP:RLOAD)"],
            "current_sign_convention": "positive_into_pin1",
            "result_current_sign_convention": "positive_into_pin1",
            "kcl_owner": "node:XTOP:out/run:v47",
            "result_kcl_owner": "node:XTOP:out/run:v47",
            "result_sha256": "4" * 64,
            "accepted_result_sha256": "4" * 64,
        },
        "ac_tran_v47_step_tuple_trace_measure_row_identity": {
            "generation_id": "analysis-v47",
            "ac_generation_id": "analysis-v47",
            "tran_generation_id": "analysis-v47",
            "trace_generation_id": "analysis-v47",
            "measure_generation_id": "analysis-v47",
            "result_generation_id": "analysis-v47",
            "ac_step_tuples": [["RLOAD", 10.0], ["RLOAD", 20.0]],
            "result_ac_step_tuples": [["RLOAD", 10.0], ["RLOAD", 20.0]],
            "tran_step_tuples": [["VIN", 5.0], ["VIN", 12.0]],
            "result_tran_step_tuples": [["VIN", 5.0], ["VIN", 12.0]],
            "trace_row_keys": ["ac:RLOAD=10", "ac:RLOAD=20", "tran:VIN=5", "tran:VIN=12"],
            "result_trace_row_keys": ["ac:RLOAD=10", "ac:RLOAD=20", "tran:VIN=5", "tran:VIN=12"],
            "measure_row_keys": ["ac:RLOAD=10", "ac:RLOAD=20", "tran:VIN=5", "tran:VIN=12"],
            "result_measure_row_keys": ["ac:RLOAD=10", "ac:RLOAD=20", "tran:VIN=5", "tran:VIN=12"],
            "simulation_owner": "simulation:ac-tran-v47",
            "result_simulation_owner": "simulation:ac-tran-v47",
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
        },
    }


def test_v47_public_replay_identity():
    assert validate_ltspice_v47_identity(_positive()) is True


def test_v47_public_rejects_hierarchy_current_sign_and_kcl_owner_mutation():
    value = deepcopy(_positive())
    branch = value["hierarchical_branch_v47_current_sign_order_kcl_owner_identity"]
    branch["result_hierarchy_paths"] = ["XTOP:RLOAD", "XTOP:XAMP:RIN", "XTOP:XAMP:CIN"]
    branch["result_branch_current_order"] = ["I(XTOP:RLOAD)", "I(XTOP:XAMP:RIN)", "I(XTOP:XAMP:CIN)"]
    branch["result_current_sign_convention"] = "positive_out_of_pin1"
    branch["result_kcl_owner"] = "node:XTOP:old/run:v46"
    assert validate_ltspice_v47_identity(value) is False


def test_v47_public_rejects_mixed_analysis_tuple_and_row_mutation():
    value = deepcopy(_positive())
    rows = value["ac_tran_v47_step_tuple_trace_measure_row_identity"]
    rows["result_ac_step_tuples"] = [["RLOAD", 20.0], ["RLOAD", 10.0]]
    rows["result_tran_step_tuples"] = [["VIN", 12.0], ["VIN", 5.0]]
    rows["result_trace_row_keys"] = ["ac:RLOAD=20", "ac:RLOAD=10", "tran:VIN=12", "tran:VIN=5"]
    rows["result_measure_row_keys"] = ["ac:RLOAD=10", "tran:VIN=5"]
    assert validate_ltspice_v47_identity(value) is False
