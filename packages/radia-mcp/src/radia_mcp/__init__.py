"""radia-mcp: MCP servers for Radia CAE ecosystem.

Layout:
  radia_mcp.cubit       — standalone Cubit MCP server (Plan A: GUI + QTimer + file drop)
  radia_mcp.build123d   — standalone build123d MCP server (STEP → Cubit workflow)

Radia-core-coupled servers (radia_ngsolve, ih, electromagnet, peec, gmsh, elf)
live in the `radia` package, reachable via entry points when the `radia`
extra is installed (pip install radia-mcp[radia]).
"""

__version__ = "0.32.1"
