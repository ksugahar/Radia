"""radia_mcp.levitation: magnetic levitation FORCE-physics + stationary uses.

The force-physics counterpart to `radia_mcp.maglev_linear` (transport).
Covers the levitation mechanisms Radia / PEEC / NGSolve are built to
compute:
- induction (eddy-current) levitation -- jumping ring, EDS lift
- electromagnetic levitation MELTING (EML) -- containerless processing
- active magnetic bearings (AMB) -- force-current-displacement, flywheels
- superconducting levitation -- Meissner vs flux pinning, HTS bulk
- diamagnetic levitation -- pyrolytic graphite, water/frog
- Earnshaw's theorem + its loopholes (the WHY)
- how to compute the levitation force (Maxwell stress / virtual work /
  time-average Lorentz)

Pick by intent: move a vehicle -> maglev_linear; suspend / melt / bear /
compute a force -> levitation.
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
