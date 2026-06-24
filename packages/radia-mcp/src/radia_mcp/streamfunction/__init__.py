"""radia_mcp.streamfunction -- Stream-Function (SF) coil-design MCP server.

SF-focused server over the kernel-agnostic (ACA+)+TSVD least-norm solver,
FE-direct psi, regularisation / folded-Tikhonov Pareto front, single-stroke
chain, and sheet-metal levers.  The detailed knowledge lives in
``radia_mcp.streamfunction.knowledge.aca_tsvd`` (moved from
``radia_mcp.radia_ngsolve`` in 2026-06: SF is not a general NGSolve usage).
"""
