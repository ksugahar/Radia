"""
Radia MCP Server

Model Context Protocol server for Radia electromagnetic simulation.
Provides magnetostatics simulation without NGSolve dependency.
"""

__version__ = "1.0.0"

from .server import RadiaMCPServer, main

__all__ = ["RadiaMCPServer", "main"]
