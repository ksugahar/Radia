"""radia_mcp.litz_transmission: Litz wire + transmission line knowledge.

Covers high-frequency conductor analysis:
- Litz wire AC loss (skin/proximity, homogenization)
- Multiconductor transmission lines (impedance, coupling)
- Cross-link to Carstensen AC copper loss (radia_mcp.peec)
- Cross-link to Hollaus MSFEM (radia_mcp.motor)

Distilled from public-safe curated corpus (44 files, 743 MB) +
03_伝送線路/ (22 files, 190 MB).
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
