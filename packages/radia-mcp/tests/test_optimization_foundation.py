"""Behavior and boundary tests for the shared optimization layer."""
import asyncio
import math

import pytest

from radia_mcp.optimization.diagnostics import optimization_gradient_check, optimization_stopping_audit
from radia_mcp.capability_packs.server import CapabilityServer


def audit(**changes):
    args = dict(x=[1.0], previous_x=[1.0], objective=1.0, previous_objective=1.0,
                gradient=[2.0], variable_scales=[1.0], objective_scale=1.0,
                problem_class="smooth_unconstrained")
    args.update(changes)
    return optimization_stopping_audit(**args)


def test_stagnation_is_not_stationarity_and_stationarity_is_not_optimality():
    assert audit()["status"] == "stalled_nonstationary"
    # Both x^2 and -x^2 have zero gradient at zero; no Hessian evidence supplied.
    stationary = audit(x=[0.0], previous_x=[0.0], objective=0, previous_objective=0, gradient=[0])
    assert stationary["status"] == "first_order_candidate"
    assert stationary["optimality_certified"] is False
    assert audit(gradient=[0], previous_x=[0])["status"] == "stationarity_only"
    assert audit(previous_x=[0])["status"] == "continue"


def test_units_and_objective_offset_do_not_change_scaled_evidence():
    base = audit(x=[2], previous_x=[1], objective=4, previous_objective=1, gradient=[4])
    converted = audit(x=[2000], previous_x=[1000], variable_scales=[1000],
                      objective=400, previous_objective=100, objective_scale=100, gradient=[0.4])
    shifted = audit(x=[2], previous_x=[1], objective=104, previous_objective=101, gradient=[4])
    assert base["metrics"] == converted["metrics"] == shifted["metrics"]
    assert base["checks"] == converted["checks"] == shifted["checks"]


@pytest.mark.parametrize("changes", [
    {"gradient": [math.nan]}, {"objective": math.inf}, {"variable_scales": [0]},
    {"variable_scales": [-1]}, {"gradient": []}, {"x": [1,2]},
    {"objective_scale": 0}, {"gradient_tolerance": -1}, {"objective": True},
    {"problem_class": "bounded"}, {"problem_class": "nonsmooth"},
])
def test_invalid_or_unsupported_evidence_fails_loudly(changes):
    with pytest.raises(ValueError):
        audit(**changes)


def test_independent_central_differences_and_wrong_gradient():
    x = [0.7, -0.4]
    scales = [2.0, 0.5]
    def f(v):
        return v[0] ** 4 + 3*v[1]**2
    references = []
    steps = [1e-4, 1e-5]
    for h in steps:
        row = []
        for i, scale in enumerate(scales):
            plus, minus = x.copy(), x.copy()
            plus[i] += h*scale
            minus[i] -= h*scale
            row.append((f(plus)-f(minus))/(2*h*scale))
        references.append(row)
    params = dict(reference_gradients=references, variable_scales=scales,
                  objective_scale=3, reference_method="central_difference", steps=steps)
    assert optimization_gradient_check([4*x[0]**3, 6*x[1]], **params)["status"] == "consistent"
    assert optimization_gradient_check([0, 6*x[1]], **params)["status"] == "needs_attention"
    params["steps"] = [1e-4, 1e-4]
    with pytest.raises(ValueError):
        optimization_gradient_check([1,2], **params)


def test_design_pack_owns_tools_without_changing_learning_profile():
    async def check():
        pack = await CapabilityServer("radia-design", "optimization").initialize()
        assert pack.status()["runtime_contract"]["complete"]
        assert pack.status()["owners"]["optimization_stopping_audit"] == "optimization"
        tool, source, _ = pack.tools["optimization_stopping_audit"]
        assert tool.annotations.readOnlyHint
        assert tool.outputSchema
        result = await source.call_tool("optimization_stopping_audit", dict(
            x=[1], previous_x=[1], objective=1, previous_objective=1, gradient=[2],
            variable_scales=[1], objective_scale=1, problem_class="smooth_unconstrained"))
        assert result[1]["status"] == "stalled_nonstationary"
        assert result[1]["optimality_certified"] is False
        learning = await CapabilityServer("radia-design", "learning").initialize()
        assert "optimization_stopping_audit" not in learning.tools
    asyncio.run(check())
