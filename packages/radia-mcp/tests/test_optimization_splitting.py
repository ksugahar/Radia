"""Independent shrinkage references and split-residual counterexamples."""
import asyncio
import math

import pytest

from radia_mcp.optimization.splitting import (
    optimization_admm_consensus_audit, optimization_proximal_gradient_audit,
)
from radia_mcp.capability_packs.server import CapabilityServer


def prox_args(**changes):
    return dict(x=[0.], smooth_gradient=[.5], proximal_point=[0.], step_size=1.,
                penalty="l1", weights=[1.], coordinate_system="dimensionless") | changes


def admm_args(**changes):
    return dict(x=[0.], z=[0.], previous_z=[0.], scaled_dual=[0.], rho=1.,
        coordinate_system="dimensionless", variant="two_block_consensus_fixed_rho") | changes


def soft(value, threshold):
    return math.copysign(max(abs(value)-threshold, 0.), value)


def test_proximal_stationarity_and_forged_stagnation():
    assert optimization_proximal_gradient_audit(**prox_args())["status"] == "first_order_candidate"
    bad = optimization_proximal_gradient_audit(**prox_args(smooth_gradient=[2.]))
    assert bad["checks"]["mapping_small"]
    assert not bad["checks"]["proximal_condition"]
    assert bad["status"] == "invalid_proximal_evidence"
    assert bad["metrics"]["proximal_stationarity_max_group_l2"] == 1.
    assert bad["optimality_certified"] is False
    # Avoid x-alpha*g cancellation masking a bad supplied prox at large x.
    large = optimization_proximal_gradient_audit(**prox_args(
        x=[1e20], proximal_point=[1e20], smooth_gradient=[1.], weights=[0.]))
    assert large["status"] == "invalid_proximal_evidence"


def test_l1_against_independent_soft_threshold_and_nonunit_step():
    x, grad, weights, alpha = [3.,-2.,.1], [.5,1.,2.], [1.,.5,2.], .25
    p = [soft(a-alpha*g, alpha*w) for a,g,w in zip(x,grad,weights)]
    result = optimization_proximal_gradient_audit(**prox_args(
        x=x, smooth_gradient=grad, proximal_point=p, weights=weights, step_size=alpha))
    assert result["status"] == "continue"
    assert result["checks"]["proximal_condition"]
    assert result["metrics"]["mapping_inf"] == pytest.approx(max(abs(a-b)/alpha for a,b in zip(x,p)))


@pytest.mark.parametrize("x,grad,p,weight,status", [
    ([3.,4.], [0.,0.], [2.4,3.2], 2., "continue"),
    ([3.,4.], [-1.2,-1.6], [3.,4.], 2., "first_order_candidate"),
    ([0.,0.], [1.,1.], [0.,0.], 2., "first_order_candidate"),
    ([0.,0.], [2.,2.], [0.,0.], 2., "invalid_proximal_evidence"),
])
def test_group_l2_proximal_ball_and_direction(x,grad,p,weight,status):
    out = optimization_proximal_gradient_audit(**prox_args(x=x, smooth_gradient=grad,
        proximal_point=p, step_size=.5, penalty="group_l2", groups=[[0,1]], weights=[weight]))
    assert out["status"] == status


@pytest.mark.parametrize("changes", [
    {"coordinate_system": "physical"}, {"penalty": "tv"}, {"weights": [-1.]},
    {"weights": []}, {"step_size": 0}, {"step_size": math.inf},
    {"smooth_gradient": [math.nan]}, {"proximal_point": [0.,1.]},
    {"groups": [[0]]}, {"penalty": "group_l2", "groups": [[0],[0]], "weights": [1.,1.]},
    {"penalty": "group_l2", "groups": [[1]]},
    {"penalty": "group_l2", "groups": [[True]]},
    {"penalty": "group_l2", "groups": [[]]}, {"proximal_tolerance": -1},
])
def test_proximal_rejects_unsupported_or_invalid_evidence(changes):
    with pytest.raises(ValueError):
        optimization_proximal_gradient_audit(**prox_args(**changes))


@pytest.mark.parametrize("changes,primal,dual", [
    ({"x": [1.]}, False, True),
    ({"x": [1.], "z": [1.]}, True, False),
    ({}, True, True),
])
def test_admm_requires_both_residuals(changes, primal, dual):
    out = optimization_admm_consensus_audit(**admm_args(**changes))
    assert out["checks"] == {"primal_small": primal, "dual_small": dual}
    assert out["status"] == ("residuals_small" if primal and dual else "continue")
    assert not out["optimality_certified"]
    assert out["subproblem_optimality"] == out["dual_update"] == "not_verified"


def test_admm_scaled_dual_and_dimension_tolerances():
    out = optimization_admm_consensus_audit(**admm_args(x=[3.,4.], z=[0.,0.],
        previous_z=[1.,2.], scaled_dual=[3.,4.], rho=2.,
        absolute_tolerance=.1, relative_tolerance=.01))
    assert out["metrics"]["primal_l2"] == 5.
    assert out["metrics"]["dual_l2"] == pytest.approx(math.sqrt(20))
    assert out["tolerances"]["primal"] == pytest.approx(math.sqrt(2)*.1+.05)
    assert out["tolerances"]["dual"] == pytest.approx(math.sqrt(2)*.1+.1)


def test_synthetic_diagonal_inverse_lasso_with_independent_exact_solution():
    # Tiny synthetic inverse map, not a claim about any physical coil geometry.
    a, b, weight, rho = [1.,2.,.5], [2.,-1.,.2], .3, 1.5
    z, u = [0.]*3, [0.]*3
    for _ in range(300):
        previous = z.copy()
        x = [(ai*bi+rho*(zi-ui))/(ai*ai+rho) for ai,bi,zi,ui in zip(a,b,z,u)]
        z = [soft(xi+ui, weight/rho) for xi,ui in zip(x,u)]
        u = [ui+xi-zi for ui,xi,zi in zip(u,x,z)]
    out = optimization_admm_consensus_audit(**admm_args(x=x,z=z,previous_z=previous,
        scaled_dual=u,rho=rho,absolute_tolerance=1e-10,relative_tolerance=1e-9))
    assert out["status"] == "residuals_small"
    exact = [soft(ai*bi, weight)/(ai*ai) for ai,bi in zip(a,b)]
    assert z == pytest.approx(exact, abs=1e-8)
    grad = [ai*(ai*zi-bi) for ai,bi,zi in zip(a,b,z)]
    p = [soft(zi-.2*gi, .2*weight) for zi,gi in zip(z,grad)]
    audit = optimization_proximal_gradient_audit(**prox_args(x=z,smooth_gradient=grad,
        proximal_point=p,step_size=.2,weights=[weight]*3))
    assert audit["status"] == "first_order_candidate"


@pytest.mark.parametrize("changes", [
    {"variant": "general_ax_bz"}, {"variant": "adaptive_rho"},
    {"coordinate_system": "physical"}, {"z": []}, {"scaled_dual": [0.,0.]},
    {"rho": 0}, {"rho": math.inf}, {"previous_z": [math.nan]},
    {"absolute_tolerance": 0}, {"relative_tolerance": -1},
])
def test_admm_rejects_invalid_or_unsupported_evidence(changes):
    with pytest.raises(ValueError):
        optimization_admm_consensus_audit(**admm_args(**changes))


def test_composed_tools_expose_structured_readonly_results():
    async def check():
        pack = await CapabilityServer("radia-design", "optimization").initialize()
        for name,args,status in [
            ("optimization_proximal_gradient_audit",prox_args(),"first_order_candidate"),
            ("optimization_admm_consensus_audit",admm_args(),"residuals_small"),
        ]:
            tool,source,_ = pack.tools[name]
            assert tool.outputSchema and tool.annotations.readOnlyHint
            assert not tool.annotations.destructiveHint
            result = await source.call_tool(name,args)
            assert result[1]["status"] == status
            assert result[1]["optimality_certified"] is False
        learning = await CapabilityServer("radia-design", "learning").initialize()
        assert "optimization_admm_consensus_audit" not in learning.tools
    asyncio.run(check())
