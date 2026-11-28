"""Import paths for the Radia application validation layer.

Production human interfaces are Simulink blocks. Notebook workbenches and the
old desktop ``radia_*.py`` PySide6 panels remain retired; Cubit's embedded
toolbar is tested in radia-mcp.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make local panel calc helpers importable for validation tests.
_REPO = Path(__file__).resolve().parents[2]
_RADIA = _REPO / "src" / "radia"
_PANELS = _RADIA / "panels"
_VALIDATION_PANELS = Path(__file__).resolve().parent
for p in (_RADIA, _PANELS, _VALIDATION_PANELS):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
