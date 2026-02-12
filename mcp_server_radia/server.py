#!/usr/bin/env python3
"""
Radia MCP Server

Model Context Protocol server for Radia electromagnetic simulation framework.
Provides magnetostatics simulation capabilities without NGSolve dependency.

Features:
- Geometry creation (magnets, coils)
- Material properties (linear, nonlinear, permanent magnets)
- Solver execution (LU, BiCGSTAB, HACApK)
- Field calculation (B, H, A, Phi)
- Workspace export for NGSolve coupling
"""

import sys
import json
import logging
from typing import Any, Dict, List, Optional, Sequence
from pathlib import Path

# Add Radia to path if needed
radia_path = Path(__file__).parent.parent / "build" / "Release"
if radia_path.exists() and str(radia_path) not in sys.path:
    sys.path.insert(0, str(radia_path))

try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import Tool, TextContent
except ImportError:
    print("Error: MCP SDK not installed. Install with: pip install mcp", file=sys.stderr)
    sys.exit(1)

# Import tool modules
from .tools import (
    geometry_tools,
    material_tools,
    solver_tools,
    field_tools,
    export_tools,
    ngsolve_integration_tools,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("radia-mcp-server")


class RadiaMCPServer:
    """Radia MCP Server implementation (NGSolve-independent)."""

    def __init__(self):
        """Initialize the Radia MCP server."""
        self.server = Server("radia-mcp-server")
        self.radia_state: Dict[str, Any] = {}

        # Register tool handlers
        self._register_handlers()

        logger.info("Radia MCP Server initialized")

    def _register_handlers(self):
        """Register all tool handlers."""

        @self.server.list_tools()
        async def list_tools() -> List[Tool]:
            """List all available Radia tools."""
            tools = []

            # Geometry tools
            tools.extend(geometry_tools.get_tools())

            # Material tools
            tools.extend(material_tools.get_tools())

            # Solver tools
            tools.extend(solver_tools.get_tools())

            # Field calculation tools
            tools.extend(field_tools.get_tools())

            # Export and workspace tools
            tools.extend(export_tools.get_tools())

            # NGSolve integration tools
            tools.extend(ngsolve_integration_tools.get_tools())

            logger.info(f"Listing {len(tools)} available tools")
            return tools

        @self.server.call_tool()
        async def call_tool(name: str, arguments: Dict[str, Any]) -> Sequence[TextContent]:
            """Execute a Radia tool."""
            logger.info(f"Calling tool: {name} with arguments: {arguments}")

            try:
                # Dispatch to appropriate tool handler
                if name.startswith("radia_geometry_"):
                    result = await geometry_tools.execute(name, arguments, self.radia_state)
                elif name.startswith("radia_material_"):
                    result = await material_tools.execute(name, arguments, self.radia_state)
                elif name.startswith("radia_solver_"):
                    result = await solver_tools.execute(name, arguments, self.radia_state)
                elif name.startswith("radia_field_"):
                    result = await field_tools.execute(name, arguments, self.radia_state)
                elif name.startswith("radia_export_") or name.startswith("radia_workspace_"):
                    result = await export_tools.execute(name, arguments, self.radia_state)
                elif name.startswith("radia_ngsolve_"):
                    result = await ngsolve_integration_tools.execute(name, arguments, self.radia_state)
                else:
                    raise ValueError(f"Unknown tool: {name}")

                # Format result
                return [TextContent(
                    type="text",
                    text=json.dumps(result, indent=2)
                )]

            except Exception as e:
                logger.error(f"Error executing tool {name}: {str(e)}", exc_info=True)
                return [TextContent(
                    type="text",
                    text=json.dumps({
                        "error": str(e),
                        "tool": name,
                        "arguments": arguments
                    }, indent=2)
                )]

    async def run(self):
        """Run the MCP server."""
        logger.info("Starting Radia MCP Server")
        async with stdio_server() as (read_stream, write_stream):
            await self.server.run(
                read_stream,
                write_stream,
                self.server.create_initialization_options()
            )


async def main():
    """Main entry point for the Radia MCP server."""
    server = RadiaMCPServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())
