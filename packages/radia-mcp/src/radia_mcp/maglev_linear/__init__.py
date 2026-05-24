"""radia_mcp.maglev_linear: Magnetic levitation + linear drive knowledge.

Covers maglev (磁気浮上) and linear drives (リニアドライブ):
- Levitation principles (EMS attraction, EDS repulsion, superconducting)
- Linear induction motor (LIM), linear synchronous motor (LSM)
- Wheel motors, transit applications
- Lab specialty: bearingless motor (cross-link to motor MCP)

Distilled from W:/.../99_アプリケーション/07_磁気浮上/ (38 files, 634 MB)
+ 09_リニアドライブ/ (32 files, 555 MB).
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
