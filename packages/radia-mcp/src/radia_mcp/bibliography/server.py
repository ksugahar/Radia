"""MCP Server: radia_mcp.bibliography

BibTeX / citation tooling: DOI / arXiv -> BibTeX, CrossRef search, parse / dedupe / canonicalize keys, lint, self-citation ratio, year distribution, journal-name normalization, health report.

Promoted 2026-06-02 from mcp-server-document.bibliography (LAB-private)
to radia-mcp (public PyPI), alongside presentation / poster /
doc_convert / pdf / bibliography.

Usage:
    mcp-server-bibliography              # stdio
    mcp-server-bibliography --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-bibliography")

# Register all bibliography_* tools from tools.py (same loop the old
# mcp-server-document.bibliography.register used).
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-bibliography",
    description=(
        'BibTeX / citation tooling: DOI / arXiv -> BibTeX, CrossRef search, parse / dedupe / canonicalize keys, lint, self-citation ratio, year distribution, journal-name normalization, health report.'
    ),
    subpackage="radia_mcp.bibliography",
    related_servers=['paper-writing', 'literature-index'],
    optional_deps=['requests'],
)


def main():
    """Entry point for mcp-server-bibliography."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-bibliography self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
