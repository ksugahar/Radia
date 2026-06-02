"""MCP Server: radia_mcp.presentation

Research-talk slide lint + PPTX tools: slide density / bullets / speaking-time-budget / figure / equation / title-verb checks, PPTX text extraction, and PPTX drawing helpers.

Promoted 2026-06-02 from mcp-server-document.presentation (LAB-private)
to radia-mcp (public PyPI), alongside presentation / poster /
doc_convert / pdf / bibliography.

Usage:
    mcp-server-presentation              # stdio
    mcp-server-presentation --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-presentation")

# Register all presentation_* tools from tools.py (same loop the old
# mcp-server-document.presentation.register used).
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-presentation",
    description=(
        'Research-talk slide lint + PPTX tools: slide density / bullets / speaking-time-budget / figure / equation / title-verb checks, PPTX text extraction, and PPTX drawing helpers.'
    ),
    subpackage="radia_mcp.presentation",
    related_servers=['poster', 'paper-writing', 'figure'],
    optional_deps=['python-pptx', 'pymupdf', 'Pillow'],
)


def main():
    """Entry point for mcp-server-presentation."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-presentation self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
