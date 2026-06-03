"""radia_mcp.maglev: Magnetic levitation knowledge.

Scope: magnetic levitation for transport / suspension / control.  The
former linear-drive (LIM/LSM, end-effect) material was removed; for the
levitation-FORCE physics see the sibling `radia_mcp.levitation`.

Covers:
- Levitation principles (EMS attraction, EDS repulsion, superconducting,
  Halbach/Inductrack, passive PM axial bearings)
- The lab's Radia-based maglev research line (CAE-AI Lab, Yano + Sugahara):
  * radia_iem_fem   -- Radia IEM (MMM/MSC) <-> reduced-potential FEM weak
    coupling for moving-magnet eddy-current levitation force
  * cln_mor_control -- Cauer Ladder Network model-order reduction for
    real-time control-coupled maglev (TEAM 28, ~1/500 speedup)
- Magnetic-wheel EDS (Kansai), Sumitomo Heavy industrial PM bearings.
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
