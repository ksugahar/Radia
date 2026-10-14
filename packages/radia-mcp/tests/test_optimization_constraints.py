"""Analytic counterexamples and public MCP contracts for constraint audits."""
import asyncio
import math

import pytest

from radia_mcp.optimization.constraints import (
    optimization_kkt_audit, optimization_projected_gradient_audit,
)
from radia_mcp.capability_packs.server import CapabilityServer


def kkt_args(**changes):
    # min x subject to -x<=0 at x=0: mu=1, L'=1-1=0.
    args = dict(gradient=[1.0], variable_scales=[1.0], objective_scale=1.0,
        inequality_values=[0.0], inequality_jacobian=[[-1.0]],
        inequality_multipliers=[1.0], inequality_scales=[1.0],
        equality_values=[], equality_jacobian=[], equality_multipliers=[],
        equality_scales=[], problem_class="smooth_constrained")
    return args | changes


def box_args(**changes):
    return dict(x=[0.0], gradient=[1.0], lower_bounds=[0.0], upper_bounds=[None],
        variable_scales=[1.0], objective_scale=1.0, problem_class="smooth_box_only") | changes


@pytest.mark.parametrize("changes,metric,status", [
    ({}, None, "first_order_candidate"),
    ({"inequality_values": [1]}, "inequality_violation_inf", "infeasible"),
    ({"inequality_multipliers": [-1], "gradient": [-1]}, "dual_violation_inf", "needs_attention"),
    ({"inequality_values": [-1]}, "complementarity_inf", "needs_attention"),
    ({"gradient": [2]}, "stationarity_inf", "needs_attention"),
    ({"equality_values": [1], "equality_jacobian": [[1]],
      "equality_multipliers": [0], "equality_scales": [1]}, "equality_violation_inf", "infeasible"),
])
def test_kkt_conditions_are_independent(changes, metric, status):
    result = optimization_kkt_audit(**kkt_args(**changes))
    assert result["status"] == status
    if metric:
        assert result["metrics"][metric] == 1
        assert not result["checks"][metric]
    assert result["constraint_qualification"] == "not_verified"
    assert result["optimality_certified"] is False


def test_mixed_constraints_with_unrestricted_equality_multiplier_and_scaling():
    # min x+2y; -x<=0, y=0, at (0,0), mu=1, nu=-2.
    args = kkt_args(gradient=[1,2], variable_scales=[1,1],
        inequality_jacobian=[[-1,0]], equality_values=[0], equality_jacobian=[[0,1]],
        equality_multipliers=[-2], equality_scales=[1])
    base = optimization_kkt_audit(**args)
    assert base["status"] == "first_order_candidate"
    # x'=1000x; f'=100f; g'=10g; h'=20h. Multipliers transform q/c.
    converted = args | dict(gradient=[.1,.2], variable_scales=[1000,1000], objective_scale=100,
        inequality_jacobian=[[-.01,0]], inequality_scales=[10], inequality_multipliers=[10],
        equality_jacobian=[[0,.02]], equality_scales=[20], equality_multipliers=[-10])
    assert optimization_kkt_audit(**converted)["metrics"] == base["metrics"]
    # Nonzero residuals must remain invariant too.
    args.update(inequality_values=[.25], equality_values=[.5], gradient=[3,2])
    converted.update(inequality_values=[2.5], equality_values=[10], gradient=[.3,.2])
    assert optimization_kkt_audit(**converted)["metrics"] == optimization_kkt_audit(**args)["metrics"]


def test_kkt_is_neither_sufficiency_nor_unconditional_necessity():
    # min -x^2 with h=y=0 at origin is a maximum along feasible x.
    maximum = kkt_args(gradient=[0,0], variable_scales=[1,1], inequality_values=[],
        inequality_jacobian=[], inequality_multipliers=[], inequality_scales=[],
        equality_values=[0], equality_jacobian=[[0,1]], equality_multipliers=[0], equality_scales=[1])
    assert optimization_kkt_audit(**maximum)["status"] == "first_order_candidate"
    assert not optimization_kkt_audit(**maximum)["optimality_certified"]
    # min x subject to x^2<=0 has unique feasible minimum, but no KKT multiplier.
    degenerate = optimization_kkt_audit(**kkt_args(inequality_jacobian=[[0]]))
    assert degenerate["status"] == "needs_attention"
    assert degenerate["constraint_qualification"] == "not_verified"


@pytest.mark.parametrize("changes,expected,status", [
    ({}, 0, "first_order_candidate"),
    ({"gradient": [-1]}, 1, "needs_attention"),
    ({"x": [1], "gradient": [-1], "upper_bounds": [1]}, 0, "first_order_candidate"),
    ({"x": [.5], "upper_bounds": [1]}, .5, "needs_attention"),
    ({"x": [0], "lower_bounds": [None], "gradient": [-2]}, 2, "needs_attention"),
    ({"upper_bounds": [0], "gradient": [-100]}, 0, "first_order_candidate"),
    ({"x": [-1], "gradient": [0]}, 1, "infeasible"),
])
def test_box_projection_and_feasibility(changes, expected, status):
    result = optimization_projected_gradient_audit(**box_args(**changes))
    assert result["status"] == status
    assert result["metrics"]["projected_gradient_inf"] == expected
    assert result["dimensionless_step"] == 1
    assert result["optimality_certified"] is False


def test_projection_scaling_and_large_coordinate_cancellation():
    base = optimization_projected_gradient_audit(**box_args(x=[.5], upper_bounds=[1]))
    converted = optimization_projected_gradient_audit(**box_args(x=[500], upper_bounds=[1000],
        variable_scales=[1000], objective_scale=100, gradient=[.1]))
    assert base["metrics"] == converted["metrics"]
    large = optimization_projected_gradient_audit(**box_args(x=[1e20], lower_bounds=[None]))
    assert large["metrics"]["projected_gradient_inf"] == 1


@pytest.mark.parametrize("changes", [
    {"inequality_jacobian": [[1,2]]}, {"inequality_multipliers": []},
    {"inequality_scales": [0]}, {"inequality_values": [math.nan]},
    {"gradient": [math.inf]}, {"equality_jacobian": [[1]]},
    {"objective_scale": 0}, {"variable_scales": [-1]},
    {"problem_class": "nonsmooth"}, {"tolerance": -1},
])
def test_kkt_rejects_malformed_evidence(changes):
    with pytest.raises(ValueError):
        optimization_kkt_audit(**kkt_args(**changes))


@pytest.mark.parametrize("changes", [
    {"lower_bounds": [2], "upper_bounds": [1]}, {"lower_bounds": []},
    {"upper_bounds": [math.inf]}, {"x": [math.nan]},
    {"gradient": []}, {"variable_scales": [0]}, {"tolerance": 0},
    {"problem_class": "smooth_constrained"},
])
def test_box_rejects_malformed_or_general_constraints(changes):
    with pytest.raises(ValueError):
        optimization_projected_gradient_audit(**box_args(**changes))


def test_composed_mcp_returns_diagnostic_status_and_errors():
    async def check():
        pack = await CapabilityServer("radia-design", "optimization").initialize()
        assert pack.status()["runtime_contract"]["complete"]
        for name, args in [("optimization_kkt_audit", kkt_args()),
                           ("optimization_projected_gradient_audit", box_args())]:
            tool, source, _ = pack.tools[name]
            assert tool.annotations.readOnlyHint and not tool.annotations.destructiveHint
            assert tool.outputSchema
            result = await source.call_tool(name, args)
            assert result[1]["status"] == "first_order_candidate"
            assert result[1]["optimality_certified"] is False
            with pytest.raises(Exception, match="Only"):
                await source.call_tool(name, args | {"problem_class": "nonsmooth"})
    asyncio.run(check())
