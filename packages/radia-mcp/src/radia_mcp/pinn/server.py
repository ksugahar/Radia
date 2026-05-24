"""
PINN MCP Server (radia_mcp.pinn) — Physics-Informed Neural Networks
and Gaussian Processes for Maxwell's equations.
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_pinn_documentation, TOPICS

mcp = FastMCP("mcp-server-pinn")


@mcp.tool()
def pinn(topic: str = "all") -> str:
    """Physics-Informed Neural Networks (PINN) + Gaussian Processes (PI-GP)
    for Maxwell's equations and computational EM.

    Topics:
      "all"
      "overview"      - PI-AI vs FEM, PINN vs PI-GP comparison
      "pi_gp"         - Physics-Informed GP foundation (Raissi 2017)
      "pi_gp_solvers" - PI-GP generalizes classical PDE solvers (Pförtner 2023)
      "pinn"          - PINN formulation, Maxwell applications
    """
    return get_pinn_documentation(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-pinn',
    description='Physics-Informed Neural Networks + Gaussian Processes for EM',
    subpackage='radia_mcp.pinn',
    related_servers=["gnn", "bayesian-opt", "fem"],
    optional_deps=["torch"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-pinn',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("PINN MCP server self-test:")
        print(f"  pinn('all'): {len(pinn('all'))} chars")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
