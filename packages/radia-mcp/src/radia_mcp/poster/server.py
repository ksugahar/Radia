"""MCP Server: radia_mcp.poster

Conference poster generation + lint: betterposter / Kelvin LaTeX templates, color-contrast (WCAG), font-size-by-viewing-distance, zone balance, QR audit, compile + print-readiness.

Promoted 2026-06-02 from mcp-server-document.poster (LAB-private)
to radia-mcp (public PyPI), alongside presentation / poster /
doc_convert / pdf / bibliography.

Usage:
    mcp-server-poster              # stdio
    mcp-server-poster --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-poster")

# Register all poster_* tools from tools.py (same loop the old
# mcp-server-document.poster.register used).
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-poster",
    description=(
        'Conference poster generation + lint: betterposter / Kelvin LaTeX templates, color-contrast (WCAG), font-size-by-viewing-distance, zone balance, QR audit, compile + print-readiness.'
    ),
    subpackage="radia_mcp.poster",
    related_servers=['presentation', 'figure'],
    optional_deps=['pymupdf', 'Pillow', 'matplotlib'],
)


def main():
    """Entry point for mcp-server-poster."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-poster self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
