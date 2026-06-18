"""pytest fixtures + path setup for radia-mcp tests."""

import importlib.util
import os
import re
import sys
from pathlib import Path

# Ensure tests resolve `radia_mcp` to THIS checkout's src/, not whatever
# `pip install -e` happens to point at on the editable-install machine.
_SRC = Path(__file__).resolve().parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# The radia-mcp matrix CI ("lightweight selftest") installs ONLY mcp + pytest
# + radia-mcp (--no-deps): NO ngsolve / netgen / scipy / numpy / matplotlib /
# gmsh / chromadb / ...  A test that imports ANY such module AT MODULE LEVEL
# crashes pytest COLLECTION there (ModuleNotFoundError) and reddens the whole
# suite.  So skip collecting any test file whose module-level imports include a
# module that is not importable -- detected GENERICALLY (find_spec per imported
# module), so ANY future heavy dependency is covered without editing a list.
# On LAB / the self-hosted runner (everything installed) nothing is skipped and
# every test runs (no coverage lost: a present dependency is never "absent").
#
# RADIA_MCP_FORCE_MINIMAL=1 reproduces the matrix's minimal env on a full-env
# box (for tools/ci_preflight.py): treat everything OUTSIDE the minimal
# baseline (stdlib + mcp + pytest + radia_mcp) as absent.  This is what catches
# a heavy-import CI break (the 2026-06-05 ngsolve incident class) BEFORE push.
_FORCE_MINIMAL = os.environ.get("RADIA_MCP_FORCE_MINIMAL") == "1"
_MINIMAL_BASELINE = set(getattr(sys, "stdlib_module_names", ())) | {
    "mcp", "pytest", "_pytest", "pluggy", "iniconfig", "packaging",
    "anyio", "attr", "attrs", "typing_extensions", "radia_mcp", "__future__",
}


def _module_absent(top_module: str) -> bool:
    """Is this top-level module unavailable in the env that will collect the
    tests?  Baseline modules are always present; under FORCE_MINIMAL anything
    outside the baseline is treated absent (matrix sim); otherwise probe."""
    if not top_module or top_module in _MINIMAL_BASELINE:
        return False
    # A local sibling helper (tests/<mod>.py or tests/<mod>/) is resolvable at
    # test runtime (pytest puts the test dir on sys.path) and ships with the
    # repo, so it is present in EVERY env -- never "absent".
    _here = Path(__file__).resolve().parent
    if (_here / f"{top_module}.py").exists() or (_here / top_module).is_dir():
        return False
    if _FORCE_MINIMAL:
        return True
    try:
        return importlib.util.find_spec(top_module) is None
    except (ImportError, ValueError):
        return True


def _imported_modules(text: str) -> set:
    """Every top-level module name imported ANYWHERE in the file (module level
    OR inside a function -- e.g. a lazy `from ngsolve import ...` in a helper).
    Both forms make the test unrunnable when the module is absent: a
    module-level import crashes COLLECTION, an in-function import errors the
    TEST.  So scan the whole file (matches the prior whole-source behavior).
    Relative imports (`from . import x`) are intra-package and skipped;
    In the real GitHub-hosted minimal matrix, `pytest.importorskip("x")`
    collects and self-skips when x is absent.  Under
    RADIA_MCP_FORCE_MINIMAL=1 on LAB, x may actually be installed; include
    importorskip targets in that simulation so the local gate still behaves
    like the minimal matrix."""
    mods = set()
    for _line in text.splitlines():
        s = _line.strip()
        if s.startswith("from "):
            target = s[5:].split()[0]
            if not target.startswith("."):          # skip relative imports
                mods.add(target.split(".")[0])
        elif s.startswith("import "):
            for _part in s[len("import "):].split(","):
                tok = _part.strip().split()[0].split(".")[0]
                if tok:
                    mods.add(tok)
    if _FORCE_MINIMAL:
        for _m in re.finditer(r"(?:pytest\.)?importorskip\(\s*['\"]([^'\"]+)['\"]", text):
            mods.add(_m.group(1).split(".")[0])
    return mods


collect_ignore = []
for _f in sorted(Path(__file__).resolve().parent.glob("test_*.py")):
    try:
        _src = _f.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        continue
    if any(_module_absent(_m) for _m in _imported_modules(_src)):
        collect_ignore.append(_f.name)
