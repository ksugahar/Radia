"""
GNN MCP Server (radia_mcp.gnn)

Graph Neural Networks for PDE / EM problems, distilled from
public-safe curated corpus

Usage:
    mcp-server-gnn              # Start MCP server
    mcp-server-gnn --selftest   # Self-test
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool, register_topics_tool

from .knowledge import get_gnn_documentation, TOPICS

mcp = FastMCP("mcp-server-gnn")


@mcp.tool()
def gnn(topic: str = "all") -> str:
    """Graph Neural Networks for PDE / EM problems.

    Topics:
        "all"               - Everything
        "overview"          - GNN family tree; when GNN vs PINN/FEM
        "physics_embedded"  - Hard-constraint PDE GNN with mixed BCs
        "equivariant_gnn"   - E(n) / SO(3) equivariant for vector fields
                              (force, B-field) — SchNet / E(n)-GNN /
                              NequIP / MACE
    """
    return get_gnn_documentation(topic)




register_status_tool(
    mcp,
    server_name='mcp-server-gnn',
    description='Graph Neural Networks for PDE/EM. Physics-Embedded GNN, E(n)-GNN / NequIP / MACE',
    subpackage='radia_mcp.gnn',
    related_servers=["pinn", "fem"],
    optional_deps=["torch", "torch_geometric"],
)


register_topics_tool(
    mcp,
    server_name='mcp-server-gnn',
    topics=TOPICS,
)


def main():
    if "--selftest" in sys.argv:
        print("GNN MCP server self-test:")
        for t in ["overview", "physics_embedded", "equivariant_gnn",
                  "all"]:
            doc = gnn(t)
            if len(doc) < 200:
                print(f"  FAIL gnn('{t}'): only {len(doc)} chars")
                sys.exit(1)
            print(f"  gnn('{t}'): {len(doc)} chars OK")
        print("  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
