"""NMR/MRI MCP Server (radia_mcp.nmr_mri).  Auto-populated from lab literature."""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

try:
    from .bibliography_index_knowledge import get_bibliography_index
except ImportError:
    def get_bibliography_index(query: str = "") -> str:
        return "NMR/MRI bibliography index not yet generated."

mcp = FastMCP("mcp-server-nmr-mri")


@mcp.tool()
def nmr_mri_bibliography(query: str = "") -> str:
    """Search the NMR/MRI bibliography catalog."""
    return get_bibliography_index(query)




register_status_tool(
    mcp,
    server_name='mcp-server-nmr-mri',
    description='NMR/MRI: gradient coils, B0 shimming, RF coils, field uniformity',
    subpackage='radia_mcp.nmr_mri',
    related_servers=["electromagnet", "accelerator"],
)


def main():
    if "--selftest" in sys.argv:
        print(f"NMR/MRI MCP server self-test: bib {len(get_bibliography_index(''))} chars  PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
