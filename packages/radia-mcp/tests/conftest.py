"""pytest fixtures + path setup for radia-mcp tests."""

import importlib.util
import os
import sys
from pathlib import Path

# Ensure tests resolve `radia_mcp` to THIS checkout's src/, not whatever
# `pip install -e` happens to point at on the editable-install machine.
_SRC = Path(__file__).resolve().parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# The radia-mcp matrix CI ("lightweight selftest") installs ONLY
# mcp + pytest + radia-mcp (--no-deps): NO ngsolve / netgen / scipy.
# The FEMM-parity + axisymmetric FEM cross-validation tests import
# ngsolve / netgen at MODULE level, so pytest's COLLECTION of them fails
# with ModuleNotFoundError on that runner and breaks the whole suite.
# When ngsolve is unavailable, skip collecting exactly those modules
# (detected by a source scan, so future heavy tests are covered
# automatically).  On LAB / the self-hosted runner (ngsolve present)
# nothing is skipped and every test runs.
#
# RADIA_MCP_FORCE_MINIMAL=1 forces the same skip path even when ngsolve
# IS installed, so `tools/ci_preflight.py` can faithfully reproduce the
# ubuntu minimal-dep matrix collection on a full-env LAB box BEFORE push
# (this is what would have caught the 2026-06-05 heavy-import CI break).
if (os.environ.get("RADIA_MCP_FORCE_MINIMAL") == "1"
        or importlib.util.find_spec("ngsolve") is None):
    collect_ignore = []
    _HEAVY_IMPORTS = ("import ngsolve", "from ngsolve",
                      "import netgen", "from netgen")
    for _f in sorted(Path(__file__).resolve().parent.glob("test_*.py")):
        try:
            _src = _f.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if any(_tok in _src for _tok in _HEAVY_IMPORTS):
            collect_ignore.append(_f.name)
