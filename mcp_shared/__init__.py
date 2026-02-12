"""
MCP Shared Utilities

Shared infrastructure for Radia-NGSolve dual MCP server architecture.
Provides workspace management and inter-server communication.
"""

__version__ = "1.0.0"

from .workspace import SharedWorkspace

__all__ = ["SharedWorkspace"]
