"""radia_mcp.streamfunction.knowledge -- SF coil-design knowledge modules.

``aca_tsvd`` was moved here from ``radia_mcp.radia_ngsolve.knowledge`` (2026-06)
because the SF coil-design content is not a general NGSolve usage and was making
radia_ngsolve too large.  The streamfunction server now OWNS this knowledge.
"""
from .aca_tsvd import get_aca_tsvd_knowledge  # noqa: F401
