"""radia_mcp.matrix_solvers: sparse linear solver theory + decision layer.

Covers:
- Direct sparse solvers (LU / PARDISO / MUMPS / SuperLU)
- Krylov subspace methods (CG / BiCGSTAB / GMRES / COCG / COCR / IDR(s))
- Preconditioners (Jacobi/SOR / IC/MIC/ILU / AMG / AMS Hiptmair-Xu)
- EM-specific (tree-cotree gauge / Biro-Preis A-V / shifted prec for air+conductor)
- Decision tree: which method for which problem class
- Lab stack reference: radia.sparsesolv_ngsolve concrete API

Distilled from public-safe curated corpus (18 files, ~236 MB).

This subpackage is the **theory/genealogy** layer.  For concrete code
usage of the lab stack (CompactAMS, COCR, etc.), see the `sparsesolv`
tool in `radia_mcp.radia_ngsolve`.
"""

from .overview_knowledge import get_overview_knowledge
from .direct_solvers_knowledge import get_direct_solvers_knowledge
from .krylov_knowledge import get_krylov_knowledge
from .preconditioners_knowledge import get_preconditioners_knowledge
from .em_specific_knowledge import get_em_specific_knowledge

__all__ = [
    "get_overview_knowledge",
    "get_direct_solvers_knowledge",
    "get_krylov_knowledge",
    "get_preconditioners_knowledge",
    "get_em_specific_knowledge",
]
