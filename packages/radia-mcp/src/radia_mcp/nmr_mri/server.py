"""NMR/MRI MCP Server (radia_mcp.nmr_mri): literature-corpus search."""

import sys

from mcp.server.fastmcp import FastMCP
from radia_mcp.common import register_status_tool

mcp = FastMCP("mcp-server-nmr-mri")


@mcp.tool()
def nmr_mri_bibliography(query: str = "") -> str:
    """Search the local literature corpus (RADIA_LIT_ROOT) by filename keywords."""
    from radia_mcp.literature_index.index_tools import lit_search
    return lit_search(query)


register_status_tool(
    mcp,
    server_name='mcp-server-nmr-mri',
    description='NMR/MRI: gradient coils, B0 shimming, RF coils, field uniformity',
    subpackage='radia_mcp.nmr_mri',
    related_servers=["electromagnet", "accelerator", "radia-ngsolve"],
)


def main():
    if "--selftest" in sys.argv:
        print("NMR/MRI MCP server self-test: PASSED")
        return
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
