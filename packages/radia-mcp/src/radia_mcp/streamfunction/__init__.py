"""radia_mcp.streamfunction -- Stream-Function (SF) coil-design MCP server.

SF-focused front door over the kernel-agnostic (ACA+)+TSVD least-norm solver,
FE-direct psi, regularisation / folded-Tikhonov Pareto front, single-stroke
chain, and sheet-metal (板金) levers.  Detailed knowledge is shared with
``radia_mcp.radia_ngsolve`` ``aca_tsvd(<topic>)``.
"""
