"""radia_mcp.team_benchmark: TEAM Workshop benchmark problems.

The TEAM Workshop (Testing Electromagnetic Analysis Methods) was
established in 1985 (1st workshop) as the standard EM-FEM/BEM
validation suite.  ~36 problems span:
- Magnetostatic (linear/nonlinear)
- Eddy current (3D, transient, with motion)
- Open boundary
- Force computation
- NDT (non-destructive testing)
- Inverse problems
- Hysteresis

Distilled from public-safe curated corpus
(45 files / 449 MB across 30 problem subfolders + _references/).

Source: TEAM Workshop reference papers (IGTE Symposium 1990,
Compumag, Cefc proceedings).

Use cases:
- "Which TEAM problem validates my HDiv-VIM implementation?"
- "What's the reference solution for problem 13 (nonlinear yoke)?"
- "Where do I find the geometry definition for problem 32 (vector hysteresis)?"
"""

from .catalog_knowledge import get_catalog_knowledge
from .magnetostatic_knowledge import get_magnetostatic_knowledge
from .eddy_current_knowledge import get_eddy_current_knowledge
from .force_motion_knowledge import get_force_motion_knowledge
from .ndt_inverse_knowledge import get_ndt_inverse_knowledge
from .special_knowledge import get_special_knowledge

__all__ = [
    "get_catalog_knowledge",
    "get_magnetostatic_knowledge",
    "get_eddy_current_knowledge",
    "get_force_motion_knowledge",
    "get_ndt_inverse_knowledge",
    "get_special_knowledge",
]
