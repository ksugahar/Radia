"""Internal radia-design domain; no additional standalone server entry point."""
from mcp.server.fastmcp import FastMCP
from typing import Any
from mcp.types import ToolAnnotations
from .. import __version__
from ..common.mcp_contract import apply_tool_contract
from .diagnostics import optimization_gradient_check, optimization_stopping_audit
from .constraints import optimization_kkt_audit, optimization_projected_gradient_audit
from .splitting import optimization_proximal_gradient_audit, optimization_admm_consensus_audit

mcp = FastMCP("radia-optimization-domain", instructions=(
    "Solver-neutral optimization diagnostics on caller-supplied evidence. "
    "Use fixed physical scales and independent derivative checks. A small step "
    "is not stationarity; stationarity is not a minimum. PDE adjoints belong to "
    "topology_optimization. No solver or user code is executed here."
))
_read_only = ToolAnnotations(readOnlyHint=True, destructiveHint=False,
                             idempotentHint=True, openWorldHint=False)
mcp.add_tool(optimization_gradient_check, annotations=_read_only)
mcp.add_tool(optimization_stopping_audit, annotations=_read_only)
mcp.add_tool(optimization_kkt_audit, annotations=_read_only)
mcp.add_tool(optimization_projected_gradient_audit, annotations=_read_only)
mcp.add_tool(optimization_proximal_gradient_audit, annotations=_read_only)
mcp.add_tool(optimization_admm_consensus_audit, annotations=_read_only)


@mcp.tool(annotations=_read_only)
def optimization_guide() -> dict[str, Any]:
    """Explain common/domain ownership, fixed scaling, and evidence limitations."""
    return {
        "source": "Kanamori, Suzuki, Takeuchi, Sato: Continuous Optimization for Machine Learning (2016), sections 3.3 and 4.1",
        "scaling": "Choose fixed positive characteristic scales S_i and f_scale before the run. x=S*z, F=f/f_scale, grad_z(F)=S*grad_x(f)/f_scale. Do not choose a vanishing current objective as f_scale.",
        "stopping": "Inspect gradient, step and objective change separately. This conservative engineering policy is not a transcription of the book's machine-epsilon formula or a proof of convergence.",
        "gradient_check": "Compare the same physical gradient at the same point to independent central differences over multiple dimensionless steps; repeat at nonstationary points and inspect truncation/noise. Consistent samples cannot prove independence.",
        "routes": {"smooth_unconstrained": "optimization", "pde_adjoint": "topology_optimization", "black_box": "bayesian_opt/evolutionary", "linear_cg": "matrix_solvers; distinguish from nonlinear CG"},
        "constraints": "Use optimization_kkt_audit for smooth g<=0, h=0 with physical multipliers and row Jacobians; use optimization_projected_gradient_audit for box-only problems. Both need explicit scales and cannot certify optimality or constraint qualifications.",
        "constraint_source": "Kanamori et al. (2016), sections 3.2 and 10.1; fixed-scale residuals are engineering diagnostics, not verbatim source algorithms.",
        "helpers": {"least_squares": "radia_mcp.optimization.nonlinear_lsq", "regularized_inverse": "radia_mcp.optimization.linear_inverse", "linear_cg": "radia_mcp.matrix_solvers.krylov"},
        "helper_caution": "Legacy topology_optimization imports re-export the same helpers. LM converged is only an absolute gradient check; inspect termination_reason and independently check derivatives. Helpers are not additional MCP tools or production solver replacements.",
        "splitting": {
            "proximal": "optimization_proximal_gradient_audit checks supplied L1/disjoint group-L2 prox results and mappings in normalized coordinates. Differentiate only the smooth term. Overlap/TV need different operators.",
            "admm": "optimization_admm_consensus_audit checks both residuals for fixed-rho unrelaxed two-block x=z only; subproblem optimality and dual updates are not verified.",
            "scaling": "Normalize variables AND the objective, gradients, penalty weights and duals consistently before solving. Unequal within-group variable scaling generally changes an isotropic group-L2 penalty; do not reuse its old weights blindly.",
            "applications": "L1 can promote sparse coil excitations; disjoint group L2 can select grouped sources; TV penalizes spatial differences, not individual amplitudes. These are modeling choices, not guarantees of manufacturability or validated field designs.",
            "sources": ["Kanamori et al. (2016), sections 12.3.1 and 15.3.1", "Boyd et al. (2011), ADMM, section 3.3: https://web.stanford.edu/~boyd/papers/admm_distr_stats.html"],
        },
        "search_ownership": {
            "helpers": "radia_mcp.optimization.global_optimizers owns the existing DE and feasibility-selection helpers; old topology imports are identity aliases. Bayesian and evolutionary guidance remain separate. Prefer established ecosystem solvers for production.",
            "selection_limits": "Legacy selection assumes missing constraints mean unconstrained and does not validate finite inputs. Validate complete finite evidence and normalize constraint residuals before ranking; a least-violating fallback is not feasible.",
            "search_limits": "DE population spread is a termination heuristic, not stationarity or global optimality. Local refinement requires smoothness and independently checked derivatives; do not polish across discrete switches blindly.",
            "regression_gates": "Existing topology multistart and simplex gates remain reproduction-specific. Their ok means the specified regression evidence passed, not all solvers converged: multistart requires the old zero-start collapse and simplex requires detection of false convergence. Use common optimization audits for new runs.",
        },
        "deferred": ["general/relaxed/adaptive-rho ADMM", "overlapping-group/TV diagnostics", "independent physical inverse-design demonstration"],
    }


apply_tool_contract(mcp, server_name="radia-optimization-domain", version=__version__)
