"""Distribution copies generated from canonical Radia Python sources.

residual_gate.py is byte-identical to src/radia/_residual_gate.py. Edit only
that canonical source and run tools/sync_mcp_residual_gate.py. The copy lets
NGSolve-only MCP users enforce the same numerical contract without importing
Radia's native backend or creating a second numerical implementation.
"""
