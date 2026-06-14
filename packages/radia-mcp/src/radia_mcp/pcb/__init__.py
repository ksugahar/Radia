"""radia_mcp.pcb: Wireless Power Transfer knowledge layer.

Covers:
- Overview: WPT physics, decision tree, lab focus areas
- Coil design + compensation topology (SS / LCC / LCL / CLLLC)
- Efficiency, Q-factor, kQ product, k coupling
- Foreign Object Detection (FOD) — lab core
- Dynamic charging (EV running, robot, bearingless motor)
- Capacitive / microwave / rectenna alternatives

Distilled from W:/03_文献・論文/00_電磁界解析/99_アプリケーション/05_ワイヤレス給電/
(269 files, 2.1 GB across 13 numbered thematic subdirs including IEICE technical
reports 2010-2019, FOD literature, lab Sugahara group papers).

Cross-references:
- `radia_mcp.peec` — PEEC for coil inductance (kQ, coupling)
- `radia_mcp.ih` — analogous SIBC / eddy current physics
- `radia_mcp.bem` — BEM for open-boundary coil field
- `radia_mcp.matrix_solvers` — solvers for coupled WPT system
"""

from .overview_knowledge import get_overview_knowledge
from .coil_compensation_knowledge import get_coil_compensation_knowledge
from .efficiency_safety_knowledge import get_efficiency_safety_knowledge
from .fod_knowledge import get_fod_knowledge
from .applications_knowledge import get_applications_knowledge
from .alternatives_knowledge import get_alternatives_knowledge

__all__ = [
    "get_overview_knowledge",
    "get_coil_compensation_knowledge",
    "get_efficiency_safety_knowledge",
    "get_fod_knowledge",
    "get_applications_knowledge",
    "get_alternatives_knowledge",
]
