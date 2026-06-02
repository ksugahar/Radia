"""MCP Server: radia_mcp.pdf

PDF manipulation: merge / split / extract-pages / rotate, metadata get-set, bookmarks, watermark, compress, crop-whitespace, font-embed and image-DPI checks.

Promoted 2026-06-02 from mcp-server-document.pdf (LAB-private)
to radia-mcp (public PyPI), alongside presentation / poster /
doc_convert / pdf / bibliography.

Usage:
    mcp-server-pdf              # stdio
    mcp-server-pdf --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-pdf")

# Register all pdf_* tools from tools.py (same loop the old
# mcp-server-document.pdf.register used).
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-pdf",
    description=(
        'PDF manipulation: merge / split / extract-pages / rotate, metadata get-set, bookmarks, watermark, compress, crop-whitespace, font-embed and image-DPI checks.'
    ),
    subpackage="radia_mcp.pdf",
    related_servers=['doc-convert'],
    optional_deps=['pymupdf'],
)


def main():
    """Entry point for mcp-server-pdf."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-pdf self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
