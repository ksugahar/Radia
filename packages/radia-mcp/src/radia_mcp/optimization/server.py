"""Internal radia-design domain; no additional standalone server entry point."""
from mcp.server.fastmcp import FastMCP
from typing import Any
from mcp.types import ToolAnnotations
from .. import __version__
from ..common.mcp_contract import apply_tool_contract
from .diagnostics import optimization_gradient_check, optimization_stopping_audit

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


@mcp.tool(annotations=_read_only)
def optimization_guide() -> dict[str, Any]:
    """Explain common/domain ownership, fixed scaling, and evidence limitations."""
    return {
        "source": "Kanamori, Suzuki, Takeuchi, Sato: Continuous Optimization for Machine Learning (2016), sections 3.3 and 4.1",
        "scaling": "Choose fixed positive characteristic scales S_i and f_scale before the run. x=S*z, F=f/f_scale, grad_z(F)=S*grad_x(f)/f_scale. Do not choose a vanishing current objective as f_scale.",
        "stopping": "Inspect gradient, step and objective change separately. This conservative engineering policy is not a transcription of the book's machine-epsilon formula or a proof of convergence.",
        "gradient_check": "Compare the same physical gradient at the same point to independent central differences over multiple dimensionless steps; repeat at nonstationary points and inspect truncation/noise. Consistent samples cannot prove independence.",
        "routes": {"smooth_unconstrained": "optimization", "pde_adjoint": "topology_optimization", "black_box": "bayesian_opt/evolutionary", "linear_cg": "matrix_solvers; distinguish from nonlinear CG"},
        "deferred": ["KKT/projected-gradient diagnostics", "proximal/ADMM guidance", "remaining generic helper migration"],
    }


apply_tool_contract(mcp, server_name="radia-optimization-domain", version=__version__)
