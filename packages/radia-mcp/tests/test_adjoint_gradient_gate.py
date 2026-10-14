import json

from radia_mcp.radia_ngsolve.adjoint_gradient_gate import adjoint_gradient_scaling_gate
from radia_mcp.radia_ngsolve.server import adjoint_gradient_scaling_gate as mcp_gate


def _rows():
    return [
        {"designVariableCount": count, "adjointSolves": 1,
         "gradientCheckRelativeError": error, "forwardAffineResidual": 1.0e-17,
         "plusAscentObjectiveRatio": 1.01, "minusAscentObjectiveRatio": 0.99,
         "fiftyStepObjectiveRatio": ratio, "fiftyStepMonotone": True}
        for count, error, ratio in [(4, 2e-10, 2.1), (8, 5e-10, 4.5), (16, 4e-9, 10.8)]
    ]


def test_adjoint_gradient_scaling_accepts_one_solve_and_signed_ascent():
    result = adjoint_gradient_scaling_gate(_rows())
    assert result["status"] == "ok"
    assert result["metrics"]["adjoint_solves"] == [1, 1, 1]


def test_adjoint_gradient_scaling_rejects_variable_cost_or_wrong_direction():
    rows = _rows()
    rows[2]["adjointSolves"] = 16
    rows[1]["plusAscentObjectiveRatio"] = 0.98
    result = adjoint_gradient_scaling_gate(rows)
    assert result["status"] == "needs_attention"
    assert result["checks"]["one_adjoint_solve_for_every_size"] is False
    assert result["checks"]["positive_direction_raises_objective"] is False


def test_adjoint_gradient_scaling_mcp_dispatches_json_rows():
    result = json.loads(mcp_gate(json.dumps(_rows())))
    assert result["status"] == "ok"
