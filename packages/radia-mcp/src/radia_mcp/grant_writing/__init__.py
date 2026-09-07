"""grant_writing sub-module -- grant proposal lint and review helpers.

This package completes the public document-writing quartet:
paper_writing / figure / grant_writing / presentation.
"""
from __future__ import annotations

from . import tools as _tools
from ..common.server_hardening import ANN_READONLY


def register(mcp) -> int:
    """Register grant_writing_* tools with the given FastMCP instance."""
    count = 0
    for name in dir(_tools):
        if name.startswith("grant_writing_") and callable(getattr(_tools, name)):
            mcp.tool(annotations=ANN_READONLY)(getattr(_tools, name))
            count += 1
    return count
