"""radia_mcp.fem: Finite Element Method formulations theory layer.

Covers all major EM-FEM formulations + numerical techniques:
- Potential formulations: A-Omega, T-Omega, H-formulation, Reduced Potential, Darwin
- Element technology: edge (RT/Nedelec), high-order, hierarchical, XFEM, Isogeometric, DG
- Gauging: tree-cotree, Coulomb, A-V, Biro-Preis (links to matrix_solvers)
- Open boundary: Kelvin transform, Asymptotic BC, PML
- Domain-specific: axisymmetric (Henrotte), time-domain FEM, harmonic balance, high-frequency
- Coupling: circuit-FE coupling, multi-physics
- Scaling: large-scale FEM, error theory, a-posteriori estimators

Distilled from W:/03_文献・論文/00_電磁界解析/10_FEM_定式化/
(199 files, 4.9 GB across 19 numbered thematic subdirs).

This is the **theory/genealogy** layer.  For concrete code usage:
- NGSolve FEM API: see `radia_mcp.radia_ngsolve` (radia_usage, sparsesolv)
- Axisymmetric FEMM-canonical: see `radia_mcp.radia_ngsolve.axifem_documentation`
- Hollaus MSFEM lamination: see `radia_mcp.motor` (hollaus_eddy)
- Solver+preconditioner: see `radia_mcp.matrix_solvers`
- BEM/MoM theory: see `radia_mcp.bem`
"""

from .overview_knowledge import get_overview_knowledge
from .potential_formulations_knowledge import get_potential_formulations_knowledge
from .elements_knowledge import get_elements_knowledge
from .gauge_open_boundary_knowledge import get_gauge_open_boundary_knowledge
from .time_domain_axisym_knowledge import get_time_domain_axisym_knowledge
from .large_scale_special_knowledge import get_large_scale_special_knowledge
from .ngsolve_hierarchy_knowledge import get_ngsolve_hierarchy_knowledge

__all__ = [
    "get_overview_knowledge",
    "get_potential_formulations_knowledge",
    "get_elements_knowledge",
    "get_gauge_open_boundary_knowledge",
    "get_time_domain_axisym_knowledge",
    "get_large_scale_special_knowledge",
    "get_ngsolve_hierarchy_knowledge",
]
