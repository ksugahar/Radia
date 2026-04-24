"""radia-mcp: MCP servers for Radia CAE ecosystem.

Layout:
  radia_mcp.cubit         — standalone Cubit MCP server (Plan A: GUI + QTimer + file drop)
  radia_mcp.build123d     — standalone build123d MCP server (STEP → Cubit workflow)
  radia_mcp.gmsh          — standalone gmsh MSH v4.1 inspect/validate/convert
  radia_mcp.elf           — historical ELF reference (no live deps)
  radia_mcp.interop       — Radia <-> third-party CAD interop
  radia_mcp.radia_ngsolve — Radia + NGSolve general FEM/BEM/Kelvin/PEEC/MSH knowledge
  radia_mcp.ih            — Induction Heating workflow (workpiece SIBC, ESIM,
                            Karl iteration, screening physics) — promoted from
                            s:/mcp-server/mcp-server-ih/ on 2026-04-24

The radia_ngsolve / ih servers reference radia from inside example
snippets in their knowledge modules, but importing them does not
require radia to be installed (knowledge is plain text + Python).
The `radia` extra (`pip install radia-mcp[radia]`) brings in the
runtime dependency for tools that actually exec radia code.
"""

__version__ = "0.32.2"
