"""radia_mcp.maglev: Magnetic levitation knowledge (UNIFIED).

Covers both halves of magnetic levitation in one server (the former
separate radia_mcp.levitation server was consolidated into this one):

- maglev SYSTEMS: EMS attraction, EDS repulsion, superconducting,
  Halbach/Inductrack, magnetic wheels, passive PM axial bearings.
- levitation FORCE physics: induction/eddy-current lift, EML melting,
  active magnetic bearings (AMB), superconducting Meissner/pinning,
  diamagnetic, Earnshaw + loopholes, force computation, benchmarks.
- The lab's Radia-based maglev research line (CAE-AI Lab, Yano + Sugahara):
  * radia_iem_fem   -- Radia IEM (HDiv-VIM) <-> reduced-potential FEM weak
    coupling for moving-magnet eddy-current levitation force
  * cln_mor_control -- Cauer Ladder Network model-order reduction for
    real-time control-coupled maglev (TEAM 28, ~1/500 speedup)

The former linear-drive (LIM/LSM, end-effect) material was removed.
"""

from .knowledge import get_knowledge

__all__ = ["get_knowledge"]
