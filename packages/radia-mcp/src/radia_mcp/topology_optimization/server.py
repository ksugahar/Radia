"""
Topology Optimization MCP Server (radia_mcp.topology_optimization)

Provides knowledge tools for shape and topology optimization in
nonlinear magnetostatics (electric motors, magnets, induction
heaters, accelerator pole faces).

Knowledge layer (theory) — pairs with `radia_mcp.radia_ngsolve`
(implementation) and `radia_mcp.mathematica` (symbolic adjoint
derivation).

Usage:
    mcp-server-topology-optimization              # Start MCP server
    mcp-server-topology-optimization --selftest   # Self-test

Distilled from:
  - Gangl-Langer-Laurain-Meftahi-Sturm 2015 (IPM motor)
  - Gangl PhD thesis Part I/II (nonlinear magnetostatics SD/TD)
  - Sturm 2015 (minimax Lagrangian for nonlinear PDE shape opt)
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from .shape_optimization_knowledge import get_shape_optimization_documentation
from .topology_derivative_knowledge import get_topology_derivative_documentation
from .applications_knowledge import get_applications_documentation
from .cae_ai_contract import cae_ai_artifact_gate as build_cae_ai_artifact_gate
from .simplex_stationarity_gate import simplex_stationarity_audit_gate as build_simplex_stationarity_audit_gate

mcp = FastMCP("mcp-server-topology-optimization")


@mcp.tool()
def topology_opt_shape_optimization(topic: str = "all") -> str:
    """Shape optimization for nonlinear magnetostatics.

    Topics:
      "all"
      "overview"        — Why nonlinear shape opt is hard, Gangl-Sturm 2015 fix
      "shape_derivative" — Velocity method, Hadamard formula, Lagrangian-based
                          derivation for nonlinear PDE constraints
      "algorithm"       — Gradient-based algorithm structure, NGSolve patterns
    """
    return get_shape_optimization_documentation(topic)


@mcp.tool()
def topology_opt_topology_derivative(topic: str = "all") -> str:
    """Topological derivative for changing topology (adding/removing material).

    Covers Sokolowski-Zochowski TD theory, the SD + TD hybrid approach,
    and nonlinear extensions (Beck-Sturm 2018, Amstutz-Gangl 2019).
    """
    return get_topology_derivative_documentation(topic)


@mcp.tool()
def topology_opt_applications(topic: str = "all") -> str:
    """Practical applications.

    Topics:
      "all"
      "motor"  — IPM brushless motor cogging-torque minimization
                 (Gangl 2015 case study, including the 4-layer
                 theory→Mathematica→NGSolve→Radia workflow)
      "field_synthesis" — the ANALYTIC / LINEAR-INVERSE branch (complement
                 to the density/gradient-shape branch the rest of this
                 server teaches): the VERIFIED permanent-magnet multipole
                 magnetization inverse (target field → radial M_r(r,θ) in
                 closed form, with the n=1 uniform-magnetization degeneracy),
                 the air-gap harmonic-content objective (SPM/Halbach), and a
                 pointer to the shipped stream-function / TSVD / target-field
                 COIL inverse in the `streamfunction` server.
      "linear_inverse" — the REGULARIZED LINEAR SOLVER behind field_synthesis:
                 truncated-SVD and Tikhonov filter factors phi=s^2/(s^2+lam^2)
                 + the L-curve, for the ill-conditioned forward map A x = b
                 (Sugahara Lab TSVD magnet/coil field-synthesis method),
                 verified against the SVD pseudo-inverse to ~1e-15.
      "outer_loop" — derivative-free OUTER-LOOP optimizers for the
                 non-closed-form / manufacturable-bounds case: Nelder-Mead
                 direct search + the fminsearchbnd bound-by-transformation
                 trick (Sugahara Lab MATLAB optimizer toolbox), verified on
                 Rosenbrock; pointer to the `evolutionary` and
                 `bayesian-opt` servers for population / global search.
    """
    return get_applications_documentation(topic)


@mcp.tool()
def topology_opt_cae_ai_artifact_gate(method_family: str, artifact_json: str) -> str:
    """Gate CAE-AI artifacts before they are promoted as engineering results.

    Supported method families are ``diffusion``, ``normalizing_flow``,
    ``reinforcement_learning``, and ``pseudoinverse``. Every family must save
    reproducibility metadata, named metrics with units, and an independent
    forward-solver verification block with observables and tolerances.

    Args:
        method_family: One of the four supported family names.
        artifact_json: JSON object containing the common and family-specific
            reproducibility fields returned in the gate report.
    """
    return build_cae_ai_artifact_gate(method_family, artifact_json)


@mcp.tool()
def topology_opt_simplex_stationarity_audit_gate(summary_json: str) -> str:
    """Audit derivative-free convergence using independent stationarity checks.

    A method can report convergence because its simplex or objective spread is
    small while stopping at a nonstationary point. This gate requires an
    independently evaluated gradient and trusted reference/control evidence.
    """
    return build_simplex_stationarity_audit_gate(summary_json)




register_status_tool(
    mcp,
    server_name='mcp-server-topology-optimization',
    description='Topology optimization: SIMP, level set, ON/OFF, MMA, Wakao autoencoder+LS SynRM',
    subpackage='radia_mcp.topology_optimization',
    related_servers=["motor", "bayesian-opt", "evolutionary"],
    optional_deps=["scipy"],
)


def main():
    if "--selftest" in sys.argv:
        print("Topology Optimization MCP server self-test:")
        for name, fn in [
            ("shape_optimization", topology_opt_shape_optimization),
            ("topology_derivative", topology_opt_topology_derivative),
            ("applications", topology_opt_applications),
        ]:
            print(f"  {name}('all'): {len(fn('all'))} chars")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
