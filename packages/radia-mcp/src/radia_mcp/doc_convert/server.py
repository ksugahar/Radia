"""MCP Server: radia_mcp.doc_convert

Document format conversion: PPTX<->PDF, PDF->JPG, slide / speaker-note extraction, figure-readability check, business-card OCR rename.

Promoted 2026-06-02 from mcp-server-document.doc_convert (LAB-private)
to radia-mcp (public PyPI), alongside presentation / poster /
doc_convert / pdf / bibliography.

Usage:
    mcp-server-doc-convert              # stdio
    mcp-server-doc-convert --selftest   # self-test
"""

import sys

from mcp.server.fastmcp import FastMCP

from ..common import register_status_tool
from . import register

mcp = FastMCP("mcp-server-doc-convert")

# Register all doc_convert_* tools from tools.py (same loop the old
# mcp-server-document.doc_convert.register used).
_n_tools = register(mcp)

register_status_tool(
    mcp,
    server_name="mcp-server-doc-convert",
    description=(
        'Document format conversion: PPTX<->PDF, PDF->JPG, slide / speaker-note extraction, figure-readability check, business-card OCR rename.'
    ),
    subpackage="radia_mcp.doc_convert",
    related_servers=['pdf', 'presentation'],
    optional_deps=['pymupdf', 'python-pptx', 'Pillow'],
)


def main():
    """Entry point for mcp-server-doc-convert."""
    if "--selftest" in sys.argv:
        print(f"mcp-server-doc-convert self-test: registered {_n_tools} tools")
        return
    mcp.run()


if __name__ == "__main__":
    main()
