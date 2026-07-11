import copy, json
from radia_mcp.radia_ngsolve.global_local_optimization_gate import global_local_optimization_replay_gate
from radia_mcp.radia_ngsolve.server import global_local_optimization_replay_gate as mcp_gate


def _summary():
    short=[{"seed":s,"best_f":v,"history_monotone":True} for s,v in [(1,.3),(2,.2),(3,.4)]]
    long=[{"seed":s,"best_f":v,"history_monotone":True} for s,v in [(1,1e-6),(2,0),(3,2e-6)]]
    return {"short_runs":short,"long_runs":long,"analytic_minimum":0,"independent_global_best_f":0,
            "polished_best_f":0,"polished_gradient_norm":1e-10,"central_gradient_relative_error":1e-9,
            "complex_step_gradient_relative_error":1e-15,"source_objective_max_abs_error":0}


def test_global_local_gate_accepts_replayed_hybrid_and_dispatches():
    row=_summary(); assert global_local_optimization_replay_gate(json.dumps(row))["status"]=="ok"
    assert json.loads(mcp_gate(json.dumps(row)))["status"]=="ok"


def test_global_local_gate_rejects_budget_regression_and_bad_gradient():
    row=copy.deepcopy(_summary()); row["long_runs"][1]["best_f"]=.5; row["central_gradient_relative_error"]=.1
    result=global_local_optimization_replay_gate(json.dumps(row)); assert result["status"]=="needs_attention"
    assert result["checks"]["long_budget_not_worse"] is False
    assert result["checks"]["central_gradient_matches"] is False
