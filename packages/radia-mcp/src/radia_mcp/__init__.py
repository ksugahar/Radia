"""radia-mcp: MCP servers for Radia CAE ecosystem.

Layout:
  radia_mcp.cubit         — standalone Cubit MCP server (Plan A: GUI + QTimer + file drop)
  radia_mcp.build123d     — standalone build123d MCP server (STEP → Cubit workflow)
  radia_mcp.gmsh          — standalone gmsh MSH v4.1 inspect/validate/convert
  radia_mcp.interop       — Radia <-> third-party CAD interop
  radia_mcp.radia_ngsolve — Radia + NGSolve general FEM/BEM/Kelvin/PEEC/MSH knowledge
  radia_mcp.ih            — Induction Heating workflow (workpiece SIBC, ESIM,
                            Karl iteration, screening physics) — promoted from
                            s:/mcp-server/mcp-server-ih/ on 2026-04-24
  radia_mcp.peec          — LAB PEEC workflow (Loop-Star, FastHenry, PyPEECBuilder,
                            PEECCircuitSolver, Bessel/Dowell/ESIM, PRIMA MOR,
                            SPICE extraction) — promoted from
                            s:/mcp-server/mcp-server-peec/ on 2026-04-24
  radia_mcp.electromagnet — Accelerator electromagnet (CoilBuilder, Hantila polarization,
                            B-input Play/Energy hysteresis, IMA sign selection,
                            multipole harmonics) — promoted from
                            s:/mcp-server/mcp-server-electromagnet/ on 2026-04-24

The radia_ngsolve / ih / peec / electromagnet servers reference radia
from inside example snippets in their knowledge modules, but importing
them does not require radia to be installed (knowledge is plain text +
Python). The `radia` extra (`pip install radia-mcp[radia]`) brings in
the runtime dependency for tools that actually exec radia code.
"""

__version__ = "0.45.3"
