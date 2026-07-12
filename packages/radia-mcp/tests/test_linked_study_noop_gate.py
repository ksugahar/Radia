import copy
import json

from radia_mcp.radia_ngsolve.linked_study_noop_gate import linked_study_silent_noop_gate
from radia_mcp.radia_ngsolve.server import linked_study_silent_noop_gate as mcp_gate


def _summary():
    return {
        "execution_mode": "native_hidden_linked_study_run",
        "source_digest_before": "a" * 64,
        "source_digest_after": "a" * 64,
        "work_copy_used": True,
        "owned_process_released": True,
        "owned_process_count_after": 0,
        "solver_success": False,
        "studies": [
            {
                "index": 0,
                "role": "field",
                "run_seconds": 0.0002,
                "has_result_before": False,
                "has_result_after": False,
                "result_file_count": 0,
                "result_table_count": 0,
            },
            {
                "index": 1,
                "role": "thermal",
                "run_seconds": 0.0001,
                "has_result_before": False,
                "has_result_after": False,
                "result_file_count": 0,
                "result_table_count": 0,
            },
        ],
    }


def test_linked_study_noop_is_classified_and_dispatched():
    payload = _summary()
    result = linked_study_silent_noop_gate(payload)
    assert result["status"] == "ok"
    assert result["classification"] == "verified_silent_noop"
    assert result["solver_result_accepted"] is False
    assert json.loads(mcp_gate(payload))["status"] == "ok"


def test_linked_study_noop_rejects_partial_result_and_process_leak():
    payload = copy.deepcopy(_summary())
    payload["studies"][1]["has_result_after"] = True
    payload["owned_process_count_after"] = 1
    result = linked_study_silent_noop_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_returned_without_results"] is False
    assert result["checks"]["no_owned_process_left"] is False
