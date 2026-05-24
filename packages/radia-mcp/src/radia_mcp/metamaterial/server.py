"""
Metamaterial MCP Server (radia_mcp.metamaterial).

Auto-populated from the lab literature library.
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return ("Metamaterial bibliography index not yet generated. "
                "Run the auto-indexer.")

mcp = FastMCP("mcp-server-metamaterial")


@mcp.tool()
def metamaterial_bibliography(query: str = "") -> str:
    """Search the metamaterial bibliography catalog."""
    return get_bibliography_index(query)




register_status_tool(
    mcp,
    server_name='mcp-server-metamaterial',
    description='Metamaterials: homogenization, effective medium, periodic structures',
    subpackage='radia_mcp.metamaterial',
    related_servers=["wpt", "litz-transmission"],
)


def main():
    if "--selftest" in sys.argv:
        n = len(get_bibliography_index(""))
        print(f"Metamaterial MCP server self-test: bibliography {n} chars  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
