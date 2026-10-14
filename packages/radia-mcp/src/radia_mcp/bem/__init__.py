"""radia_mcp.bem: Moment Method / Boundary Element Method theory + decision layer.

Covers:
- MoM foundations (RWG, Galerkin/collocation, basis functions)
- Surface integral equations (EFIE / MFIE / CFIE / PMCHWT)
- Low-frequency stabilization (Loop-Star, Calderón, PMCHWT-LF)
- BEM for eddy current (SIBC, Ishibashi single-unknown formulation)
- FEM-BEM hybrid (open-boundary magnetostatics)
- H-matrix acceleration (ACA, BLR, HACApK — lab core)
- FMM (Greengard-Rokhlin, removed from Radia but documented)
- Application areas (PCB EMC, planar EM, shielding, transformer)

Distilled from public-safe curated corpus
(135 files, 1040 MB across 16 numbered top-level folders).

Theory/genealogy layer.  For code usage:
- Radia HDiv-VIM: see `radia_mcp.radia_ngsolve` (hdiv_vim)
- PEEC + filaments: see `radia_mcp.peec`
- ngsolve.bem: see `radia_mcp.radia_ngsolve` (ngsbem_inductance)
- IH workpiece SIBC: see `radia_mcp.ih` (sibc)
"""

from .overview_knowledge import get_overview_knowledge
from .mom_foundations_knowledge import get_mom_foundations_knowledge
from .surface_ie_knowledge import get_surface_ie_knowledge
from .low_freq_knowledge import get_low_freq_knowledge
from .h_matrix_knowledge import get_h_matrix_knowledge
from .fem_bem_hybrid_knowledge import get_fem_bem_hybrid_knowledge

__all__ = [
    "get_overview_knowledge",
    "get_mom_foundations_knowledge",
    "get_surface_ie_knowledge",
    "get_low_freq_knowledge",
    "get_h_matrix_knowledge",
    "get_fem_bem_hybrid_knowledge",
]
