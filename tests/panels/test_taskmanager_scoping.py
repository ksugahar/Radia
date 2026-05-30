"""Regression test for the UnboundLocalError caused by a late
`from ngsolve import TaskManager` inside a function that ALREADY used
`with TaskManager():` earlier in the same function.

Background (keiko @ 100号機, radia 4.85.1, 2026-05-30):

    File "calc_verify_vol.py", line 38, in verify_vol
        with TaskManager():
    UnboundLocalError: cannot access local variable 'TaskManager' where
                       it is not associated with a value

Python compiles the entire function body before executing it.  If
ANY statement in the function does `from ngsolve import TaskManager`,
Python treats `TaskManager` as a LOCAL variable throughout the
function -- even at lines that appear lexically BEFORE the import.
So an earlier `with TaskManager():` raises UnboundLocalError on the
local that hasn't been assigned yet.

The fix is to put every late `from ngsolve import TaskManager` into
the function's TOP import block (next to `from ngsolve import Mesh,
Integrate, CF, ...`).  This test catches any future regression by
walking the AST of every `calc_*.py` in `src/radia/panels/` and
flagging any function where the FIRST `TaskManager` USAGE is BEFORE
the FIRST `TaskManager` IMPORT.

The test does NOT require Cubit / NGSolve to run -- it is pure
static analysis on the source file content.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

PANELS = (
    pathlib.Path(__file__).resolve().parents[2] / "src" / "radia" / "panels"
)


def _function_late_imports_taskmanager(fn: ast.FunctionDef):
    """Return (first_use_line, first_import_line) if `fn` uses
    `TaskManager` BEFORE importing it; else None."""
    first_use = None
    first_imp = None
    for node in ast.walk(fn):
        if (isinstance(node, ast.Name)
                and node.id == "TaskManager"
                and isinstance(node.ctx, ast.Load)):
            if first_use is None or node.lineno < first_use:
                first_use = node.lineno
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                name = alias.asname or alias.name
                if name == "TaskManager":
                    if first_imp is None or node.lineno < first_imp:
                        first_imp = node.lineno
    if first_use is not None and first_imp is not None and first_use < first_imp:
        return (first_use, first_imp)
    return None


@pytest.mark.parametrize(
    "calc_file",
    sorted(PANELS.glob("calc_*.py")),
    ids=lambda p: p.name,
)
def test_no_late_taskmanager_import(calc_file: pathlib.Path):
    """Each calc_*.py must not have any function that USES TaskManager
    before IMPORTING it (UnboundLocalError trap)."""
    tree = ast.parse(calc_file.read_text(encoding="utf-8"))
    offenders = []
    for fn in ast.walk(tree):
        if isinstance(fn, ast.FunctionDef):
            result = _function_late_imports_taskmanager(fn)
            if result:
                first_use, first_imp = result
                offenders.append(
                    f"{calc_file.name}::{fn.name}  uses TaskManager "
                    f"line {first_use}, imports it line {first_imp} "
                    f"(UnboundLocalError trap -- move import to top "
                    f"of function alongside other ngsolve imports)")
    assert not offenders, (
        "\n".join(offenders)
        + "\n\nFix: in every offending function, add `TaskManager` to "
          "the top-of-function `from ngsolve import ...` line and "
          "remove the later `from ngsolve import TaskManager`.")
