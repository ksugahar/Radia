"""pytest fixtures + path setup for radia-mcp tests."""

import ast
import importlib
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
_TEST_ROOT = Path(__file__).resolve().parent
_RADIA_MCP_ROOT = _SRC / "radia_mcp"

# The radia-mcp matrix CI ("lightweight selftest") installs ONLY mcp + pytest
# + radia-mcp (--no-deps): NO ngsolve / netgen / scipy / numpy / matplotlib /
# gmsh / chromadb / ...  A test that imports ANY such module AT MODULE LEVEL
# crashes pytest COLLECTION there (ModuleNotFoundError) and reddens the whole
# suite.  So skip collecting any test file whose imports include a module that
# is not importable, including imports under tests/ subdirectories and missing
# optional dependencies reached through a project module such as
# `radia_mcp.radia_ngsolve.solve -> ngsolve`.
# On LAB / a full-dependency runner, every remaining package test runs. Actual
# Netgen/NGSolve solves, convergence studies, and solver comparisons live in
# the repository-level `validation_test/radia_mcp/` suite instead.
#
# RADIA_MCP_FORCE_MINIMAL=1 reproduces the matrix's minimal env on a full-env
# box (for tools/ci_preflight.py): treat everything OUTSIDE the minimal
# baseline (stdlib + mcp + pytest + radia_mcp) as absent.  This is what catches
# a heavy-import CI break (the 2026-06-05 ngsolve incident class) BEFORE push.
_FORCE_MINIMAL = os.environ.get("RADIA_MCP_FORCE_MINIMAL") == "1"
_MINIMAL_BASELINE = set(getattr(sys, "stdlib_module_names", ())) | {
    "mcp", "pytest", "_pytest", "pluggy", "iniconfig", "packaging",
    "anyio", "attr", "attrs", "typing_extensions", "radia_mcp", "__future__",
    # Transitive dependencies installed by the workflow's `mcp>=1.0,<2` line.
    "annotated_types", "certifi", "click", "cffi", "cryptography", "dotenv",
    "h11", "httpcore", "httpx", "httpx_sse", "idna", "jsonschema", "jwt",
    "multipart", "pydantic", "pydantic_core", "pydantic_settings", "referencing",
    "rpds", "sse_starlette", "starlette", "typing_inspection", "uvicorn",
}
_PROJECT_IMPORT_CACHE = {}


def _top_module(module_name: str) -> str:
    return module_name.split(".")[0] if module_name else ""


def _module_absent(module_name: str) -> bool:
    """Is this top-level module unavailable in the env that will collect the
    tests?  Baseline modules are always present; under FORCE_MINIMAL anything
    outside the baseline is treated absent (matrix sim); otherwise probe."""
    top_module = _top_module(module_name)
    if not top_module or top_module in _MINIMAL_BASELINE:
        return False
    # A local sibling helper (tests/<mod>.py or tests/<mod>/) is resolvable at
    # test runtime (pytest puts the test dir on sys.path) and ships with the
    # repo, so it is present in EVERY env -- never "absent".
    if (_TEST_ROOT / f"{top_module}.py").exists() or (_TEST_ROOT / top_module).is_dir():
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
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return mods
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(alias.name for alias in node.names if alias.name)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            mods.add(node.module)
    if _FORCE_MINIMAL:
        for _m in re.finditer(r"(?:pytest\.)?importorskip\(\s*['\"]([^'\"]+)['\"]", text):
            mods.add(_m.group(1).split(".")[0])
    return mods


def _project_module_path(module_name: str) -> Path | None:
    if not module_name.startswith("radia_mcp."):
        return None
    parts = module_name.split(".")[1:]
    direct = _RADIA_MCP_ROOT.joinpath(*parts).with_suffix(".py")
    if direct.exists():
        return direct
    package_init = _RADIA_MCP_ROOT.joinpath(*parts) / "__init__.py"
    if package_init.exists():
        return package_init
    return None


def _project_import_has_absent_dependency(module_name: str, seen: set | None = None) -> bool:
    if not module_name.startswith("radia_mcp."):
        return False
    if module_name in _PROJECT_IMPORT_CACHE:
        return _PROJECT_IMPORT_CACHE[module_name]
    if seen is None:
        seen = set()
    if module_name in seen:
        return False
    seen.add(module_name)

    if _FORCE_MINIMAL:
        path = _project_module_path(module_name)
        if path is None:
            _PROJECT_IMPORT_CACHE[module_name] = False
            return False
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            _PROJECT_IMPORT_CACHE[module_name] = False
            return False
        for dep in _imported_modules(text):
            if _module_absent(dep):
                _PROJECT_IMPORT_CACHE[module_name] = True
                return True
            if dep.startswith("radia_mcp.") and _project_import_has_absent_dependency(dep, seen):
                _PROJECT_IMPORT_CACHE[module_name] = True
                return True
        _PROJECT_IMPORT_CACHE[module_name] = False
        return False

    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        missing = exc.name or ""
        # Skip only when the project module exists but one of its optional
        # third-party dependencies is absent.  A typo in the project import
        # itself must still surface as a real collection error.
        if missing and not module_name.startswith(missing):
            result = _module_absent(missing)
            _PROJECT_IMPORT_CACHE[module_name] = result
            return result
    _PROJECT_IMPORT_CACHE[module_name] = False
    return False


collect_ignore = []
for _f in sorted(_TEST_ROOT.rglob("test_*.py")):
    try:
        _src = _f.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        continue
    if any(_module_absent(_m) or _project_import_has_absent_dependency(_m)
           for _m in _imported_modules(_src)):
        collect_ignore.append(_f.relative_to(_TEST_ROOT).as_posix())
