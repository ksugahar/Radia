"""radia_mcp.fusion_reactor: Fusion reactor magnet knowledge.

Covers fusion plasma confinement EM analysis:
- ITER tokamak coil systems (TF, PF, CS)
- Stellarator (Heliotron, helical, modular)
- Eddy current in passive structures
- Mitsubishi heliotron + modular research (lab connection)

Distilled from public-safe curated corpus (52 files, 1.6 GB).
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
