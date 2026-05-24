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
    """
    return get_applications_documentation(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-topology-optimization',
    description='Topology optimization: SIMP, level set, ON/OFF, MMA, Wakao autoencoder+LS SynRM',
    subpackage='radia_mcp.topology_optimization',
    related_servers=["motor", "optuna", "evolutionary"],
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
