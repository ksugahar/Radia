"""radia_mcp.rna_mec: Reluctance Network Analysis + Magnetic Equivalent Circuit.

Covers RNA/MEC techniques + lab Sugahara group (田中/羽根) lineage:
- RNA fundamentals (nodal / mesh-based MEC)
- RNA + magnetic-moment model coupling for transformer core analysis
- Dynamic hysteresis in MEC (lab core)
- Topology optimization of magnetic actuators via RNA

Distilled from public-safe curated corpus (17 files, 45 MB).
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
