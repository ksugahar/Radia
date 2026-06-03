"""radia_mcp.comsol_converter — COMSOL <-> IR <-> NGSolve model translation.

A neutral intermediate representation (IR) is the hub: each solver maps to
and from the IR instead of writing N x M direct translators. v1 scope is
low-frequency magnetics (magnetostatic A-formulation), covering the lab-core
TEAM problems (6 / 13 / 20).

Exposed as ``mcp-server-comsol-converter`` (tool prefix ``cc_``). The server
reuses existing radia-mcp knowledge (interop COMSOL Java API + lab tips, fem /
bem NGSolve formulation, team_benchmark reference solutions) and coordinates
the external COMSOL MCP (mph_reader / Java gen-run) + radia-ngsolve for the
actual solve. See ir.py for the IR and server.py for the tool surface.
"""
