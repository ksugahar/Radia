"""pytest fixtures + path setup for radia-mcp tests."""

import ast
import importlib
import importlib.util
import json
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

_CI_SELECTION_RAW = os.environ.get("RADIA_MCP_CI_SELECTION_JSON", "")
_CI_SELECTORS = tuple(json.loads(_CI_SELECTION_RAW)) if _CI_SELECTION_RAW else ()


def _selector_file(selector: str) -> str:
    return selector.replace("\\", "/").split("::", 1)[0]


def _relative_nodeid(nodeid: str) -> str:
    normalized = nodeid.replace("\\", "/")
    marker = "packages/radia-mcp/"
    if marker in normalized:
        return normalized.split(marker, 1)[1]
    return normalized


_CI_SELECTED_FILES = {_selector_file(selector) for selector in _CI_SELECTORS}
_CI_SELECT_ALL = not _CI_SELECTORS or "tests" in _CI_SELECTED_FILES

# The radia-mcp matrix CI ("lightweight selftest") installs mcp + pytest + numpy
# + radia-mcp[maintenance]: NO ngsolve / netgen / scipy / matplotlib /
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
# baseline (stdlib + mcp + pytest + numpy + radia_mcp) as absent. This catches
# a heavy-import CI break (the 2026-06-05 ngsolve incident class) BEFORE push.
_FORCE_MINIMAL = os.environ.get("RADIA_MCP_FORCE_MINIMAL") == "1"
_MINIMAL_BASELINE = set(getattr(sys, "stdlib_module_names", ())) | {
    "mcp", "pytest", "_pytest", "pluggy", "iniconfig", "packaging",
    "anyio", "attr", "attrs", "typing_extensions", "radia_mcp", "__future__",
    # Transitive dependencies installed from the package's MCP SDK requirement.
    "annotated_types", "certifi", "click", "cffi", "cryptography", "dotenv",
    "h11", "httpcore", "httpx", "httpx_sse", "idna", "jsonschema", "jwt",
    "multipart", "pydantic", "pydantic_core", "pydantic_settings", "referencing",
    "rpds", "sse_starlette", "starlette", "typing_inspection", "uvicorn",
    # Lightweight optional maintenance extra, exercised by the matrix.
    "tomlkit",
    # CSV/input contracts install NumPy in the test environment, not the wheel.
    "numpy",
}
_PROJECT_IMPORT_CACHE = {}

if _FORCE_MINIMAL:
    _real_find_spec = importlib.util.find_spec

    def _minimal_find_spec(name, package=None):
        """Make dynamic optional-dependency probes match GitHub minimal CI."""
        top_module = name.split(".", 1)[0]
        if top_module not in _MINIMAL_BASELINE and not (
            (_TEST_ROOT / f"{top_module}.py").exists()
            or (_TEST_ROOT / top_module).is_dir()
        ):
            return None
        return _real_find_spec(name, package)

    importlib.util.find_spec = _minimal_find_spec


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


_IMPORT_GUARDS = {"ImportError", "ModuleNotFoundError", "Exception", "BaseException"}


def _guards_import_errors(handler) -> bool:
    kinds = handler.type
    if kinds is None:
        return True
    names = kinds.elts if isinstance(kinds, ast.Tuple) else [kinds]
    return any(isinstance(n, ast.Name) and n.id in _IMPORT_GUARDS for n in names)


def _unguarded_imports(text: str) -> set:
    """Imports of a library module that fail when their target is absent.

    An import inside ``try: ... except ImportError`` is an optional backend by
    construction (the module reports it instead of crashing), so it does not
    make the tests that reach this module unrunnable; every other import,
    module level or lazy inside a function, does.
    """
    mods = set()
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return mods

    def visit(node):
        if isinstance(node, ast.Try) and any(_guards_import_errors(h) for h in node.handlers):
            for child in (*node.handlers, *node.orelse, *node.finalbody):
                visit(child)
            return
        if isinstance(node, ast.Import):
            mods.update(alias.name for alias in node.names if alias.name)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            mods.add(node.module)
        for child in ast.iter_child_nodes(node):
            visit(child)

    visit(tree)
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


def _project_import_has_absent_dependency(module_name: str, seen: set | None = None) -> str:
    """Return the absent module reached at import time through ``module_name`` ("" if none)."""
    if not module_name.startswith("radia_mcp."):
        return ""
    if module_name in _PROJECT_IMPORT_CACHE:
        return _PROJECT_IMPORT_CACHE[module_name]
    if seen is None:
        seen = set()
    if module_name in seen:
        return ""
    seen.add(module_name)

    # Never import a server merely to decide whether its tests are collectable.
    # FastMCP construction reads settings and may emit warnings or perform other
    # import-time work; with the repository's warnings-as-errors policy that can
    # abort collection before pytest reaches a single test.  Static traversal is
    # sufficient here because this gate only answers whether a declared import
    # is unavailable in the current (or FORCE_MINIMAL) environment.
    path = _project_module_path(module_name)
    if path is None:
        _PROJECT_IMPORT_CACHE[module_name] = ""
        return ""
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        _PROJECT_IMPORT_CACHE[module_name] = ""
        return ""
    for dep in _unguarded_imports(text):
        if _module_absent(dep):
            _PROJECT_IMPORT_CACHE[module_name] = _top_module(dep)
            return _top_module(dep)
        if dep.startswith("radia_mcp."):
            absent = _project_import_has_absent_dependency(dep, seen)
            if absent:
                _PROJECT_IMPORT_CACHE[module_name] = absent
                return absent
    _PROJECT_IMPORT_CACHE[module_name] = ""
    return ""


def _with_local_helper_imports(modules: set, seen: set | None = None) -> set:
    """Add the imports of sibling helper modules (``tests/<name>.py``).

    Shared payload builders live in such helpers; their dependencies decide
    whether the importing test file can be collected just as its own do.
    """
    seen = set() if seen is None else seen
    out = set(modules)
    for module in modules:
        top = _top_module(module)
        helper = _TEST_ROOT / f"{top}.py"
        if top in seen or not helper.is_file():
            continue
        seen.add(top)
        try:
            nested = _imported_modules(helper.read_text(encoding="utf-8", errors="ignore"))
        except OSError:
            continue
        out |= _with_local_helper_imports(nested, seen)
    return out


collect_ignore = []
# Files left out because an optional dependency is absent, with the reason.
# They are reported in the session header and summary, never dropped silently.
_DEPENDENCY_IGNORED: dict[str, str] = {}
if _CI_SELECT_ALL:
    _dependency_scan_files = sorted(_TEST_ROOT.rglob("test_*.py"))
else:
    # Impact CI discovers from tests/ so pytest honours collect_ignore before
    # importing a module.  Ignore every unselected file by name, but parse only
    # the selected files for optional dependencies; this keeps the SMB cost
    # bounded and prevents collection-time imports such as PIL/numpy in the
    # minimal matrix.
    _selected_files = {
        _TEST_ROOT / _selector_file(selector).removeprefix("tests/")
        for selector in _CI_SELECTORS
        if _selector_file(selector).endswith(".py")
    }
    _all_test_files = sorted(_TEST_ROOT.rglob("test_*.py"))
    collect_ignore.extend(
        _f.relative_to(_TEST_ROOT).as_posix()
        for _f in _all_test_files
        if _f not in _selected_files
    )
    _dependency_scan_files = sorted(_selected_files)

for _f in _dependency_scan_files:
    if not _f.is_file():
        continue
    _relative = _f.relative_to(_TEST_ROOT).as_posix()
    _package_relative = f"tests/{_relative}"
    try:
        _src = _f.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        continue
    for _m in sorted(_with_local_helper_imports(_imported_modules(_src))):
        _absent = _top_module(_m) if _module_absent(_m) else _project_import_has_absent_dependency(_m)
        if _absent:
            collect_ignore.append(_relative)
            _DEPENDENCY_IGNORED[_relative] = (
                _absent if _absent == _top_module(_m) else f"{_absent} (via {_m})")
            break


def pytest_report_header(config):
    if not _DEPENDENCY_IGNORED:
        return None
    return (f"radia-mcp: {len(_DEPENDENCY_IGNORED)} test file(s) not collected "
            "because an optional dependency is absent (listed in the summary)")


def pytest_terminal_summary(terminalreporter):
    if not _DEPENDENCY_IGNORED:
        return
    terminalreporter.section("radia-mcp files not collected (optional dependency absent)")
    for relative, reason in sorted(_DEPENDENCY_IGNORED.items()):
        terminalreporter.line(f"{relative}: missing {reason}")


def pytest_collection_modifyitems(config, items):
    """Keep node-specific selectors precise after file-level collection."""
    if _CI_SELECT_ALL:
        return

    whole_files = {
        selector for selector in _CI_SELECTORS if "::" not in selector
    }
    exact_nodes = {
        selector.replace("\\", "/")
        for selector in _CI_SELECTORS
        if "::" in selector
    }
    kept = []
    deselected = []
    for item in items:
        relative = _relative_nodeid(item.nodeid)
        file_part = relative.split("::", 1)[0]
        if file_part in whole_files or relative in exact_nodes:
            kept.append(item)
        else:
            deselected.append(item)
    if deselected:
        config.hook.pytest_deselected(items=deselected)
    items[:] = kept
