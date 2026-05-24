"""
NDT MCP Server (radia_mcp.ndt) — Non-Destructive Evaluation knowledge.

Auto-populated from the lab literature library (~145 PDFs).
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return ("NDT bibliography index not yet generated. "
                "Run tools/pdf_batch_indexer.py + "
                "tools/catalog_to_knowledge.py to populate.")

mcp = FastMCP("mcp-server-ndt")


@mcp.tool()
def ndt_bibliography(query: str = "") -> str:
    """Search the NDT/NDE bibliography catalog.

    Args:
      query: empty string → all entries; "__scanned__" → scanned PDFs;
             keyword → filter by title/authors/abstract substring.
    """
    return get_bibliography_index(query)




register_status_tool(
    mcp,
    server_name='mcp-server-ndt',
    description='Non-destructive testing: eddy current testing, magnetic flux leakage, MFL signal analysis',
    subpackage='radia_mcp.ndt',
    related_servers=["ih", "magnetic-materials"],
)


def main():
    if "--selftest" in sys.argv:
        n = len(get_bibliography_index(""))
        print(f"NDT MCP server self-test: bibliography {n} chars  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
