import copy
import json

from radia_mcp.radia_ngsolve.server import source_free_static_null_solution_gate


def good():
    return {
        "analysis_kind": "magnetostatic",
        "source_count": 0,
        "condition_count": 0,
        "has_mesh": True,
        "has_result": True,
        "tables": [
            {"quantity_kind": kind, "rows": 1, "cols": 3, "values": [[0, 0, 0]]}
            for kind in ("joule_loss", "hysteresis_loss", "iron_loss")
        ],
    }


def test_accepts_exact_null_solution():
    result = json.loads(source_free_static_null_solution_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["maximum_absolute_observable"] == 0


def test_rejects_nonzero_loss_in_source_free_model():
    payload = copy.deepcopy(good())
    payload["tables"][0]["values"][0][2] = 1e-6
    result = json.loads(source_free_static_null_solution_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["null_solution_observables_zero"] is False


def test_rejects_driven_model_mislabeled_as_null():
    payload = good()
    payload["source_count"] = 1
    result = json.loads(source_free_static_null_solution_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"


def test_rejects_ragged_table_even_when_declared_shape_matches_first_row():
    payload = good()
    payload["tables"][0]["values"] = [[0, 0, 0], [0]]
    payload["tables"][0]["rows"] = 2
    result = json.loads(source_free_static_null_solution_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["table_shapes_consistent"] is False
