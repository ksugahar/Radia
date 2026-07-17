"""presentation sub-module — 学会発表スライドの作文技術.

2026-07-17: MERGED into mcp-server-paper-writing (Sugahara decision:
slide decks cannot yet be authored end-to-end by AI, so the slide
lint / PPTX toolset does not warrant a standalone MCP server).  The
standalone server.py / mcp-server-presentation entry point were
retired; radia_mcp.paper_writing.server calls register() below so all
presentation_* tools are served by mcp-server-paper-writing.
"""
from . import tools as _tools


def register(mcp) -> int:
    """Register presentation_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("presentation_") and callable(getattr(_tools, name)):
            mcp.tool()(getattr(_tools, name))
            count += 1
    return count
