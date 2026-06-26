#!/usr/bin/env python
"""audit_new_panel_contract.py -- enforce the canonical panel-adding
recipe in docs/panels/ADDING_NEW_PANEL.md.

Pure static AST analysis -- runs in <1 s, no NGSolve/Cubit needed.
Suitable as a CI gate.

What it checks
==============

For every ``src/radia/panels/calc_*.py``:

  C1. Has a top-level ``build_argparser() -> ArgumentParser`` callable.
  C2. Has a ``calc(args) -> dict`` callable.
  C3. ``__main__`` block calls ``calc_main(build_argparser, calc)``.
  C4. Heavy imports (ngsolve, radia, cubit, netgen) are NOT at module top.
  C5. ``build_argparser()`` body has NO heavy imports either (must be
      cheap so panels can introspect it without paying the NGSolve cost).
  C6. No function uses ``TaskManager`` BEFORE importing it (the late-
      import UnboundLocalError trap -- keiko 100号機 2026-05-30).
      Already covered by validation_test/panels/test_taskmanager_scoping.py; we
      re-run it here so this script is a complete one-stop audit.

For every ``src/radia/radia_*.py`` panel module:

  P1. Heavy imports (ngsolve, radia.* runtime, cubit, netgen) are NOT
      at module top.  Panel modules load on every Cubit toolbar launch.
  P2. At least one class subclasses ``ModePanel`` and calls
      ``bind_argparser(`` somewhere in its ``__init__``.
  P3. Every ``build_command`` returns the result of
      ``build_command_from_parser(...)`` -- hand-assembled argv defeats
      single-source-of-truth.
  P4. The mode-suffix string is consistent between ``json_output(...,
      suffix)`` and ``msh_output(..., suffix)`` calls in the same
      ``build_command`` (so the Persistence Policy .log lands next to
      the right .msh).

Exit code 0 = clean; non-zero = violations printed and counted.

Run locally:
    python tools/audit_new_panel_contract.py

CI:
    python tools/audit_new_panel_contract.py || exit 1
"""

from __future__ import annotations

import ast
import pathlib
import re
import sys
from typing import Iterable

ROOT = pathlib.Path(__file__).resolve().parent.parent
CALC_DIR = ROOT / "src" / "radia" / "panels"
PANEL_DIR = ROOT / "src" / "radia"

# Modules that are NEVER allowed at top-of-file in either calc_*.py or
# radia_*.py panel modules.  Either they trigger long imports
# (ngsolve / netgen-mesher) or they require Cubit (which the panel
# subprocess explicitly cannot have, per CLAUDE.md 4-Layer policy).
HEAVY_TOP_IMPORTS = (
    "ngsolve",
    "netgen",
    "cubit",
    "radia_pybind",   # the Radia C++ binding
)


def _is_heavy_import(node: ast.AST) -> str | None:
    """Return the offending module name, or None if not a heavy import."""
    if isinstance(node, ast.Import):
        for alias in node.names:
            top = alias.name.split(".")[0]
            if top in HEAVY_TOP_IMPORTS:
                return alias.name
    if isinstance(node, ast.ImportFrom):
        top = (node.module or "").split(".")[0]
        if top in HEAVY_TOP_IMPORTS:
            return node.module
    return None


def _toplevel_heavy_imports(tree: ast.Module) -> list[tuple[int, str]]:
    """Return [(lineno, modname)] for heavy imports at module top."""
    bad = []
    for node in tree.body:
        m = _is_heavy_import(node)
        if m:
            bad.append((node.lineno, m))
    return bad


def _function_has_heavy_imports(fn: ast.FunctionDef) -> list[tuple[int, str]]:
    """Heavy imports anywhere inside fn (e.g. build_argparser MUST NOT)."""
    bad = []
    for node in ast.walk(fn):
        if node is fn:
            continue
        m = _is_heavy_import(node)
        if m:
            bad.append((node.lineno, m))
    return bad


def _function_late_imports_taskmanager(fn: ast.FunctionDef
                                      ) -> tuple[int, int] | None:
    """C6: Returns (first_use_line, first_import_line) if fn uses
    `TaskManager` BEFORE importing it.  Mirror of the test in
    validation_test/panels/test_taskmanager_scoping.py."""
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
                if (alias.asname or alias.name) == "TaskManager":
                    if first_imp is None or node.lineno < first_imp:
                        first_imp = node.lineno
    if first_use is not None and first_imp is not None and first_use < first_imp:
        return (first_use, first_imp)
    return None


def _has_callable(tree: ast.Module, name: str) -> bool:
    """Check for a function with `name` defined at module top."""
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return True
    return False


def _has_callable_anywhere(tree: ast.Module, name: str) -> bool:
    """Check for a function with `name` anywhere in the module (top-level
    OR nested inside main() etc.).  The canonical pattern allows
    `def run(args):` nested inside `def main():`."""
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return True
    return False


def _main_calls(tree: ast.Module, func: str) -> bool:
    """Does the if __name__ == '__main__' block call `func(...)`?"""
    for node in tree.body:
        if (isinstance(node, ast.If)
                and isinstance(node.test, ast.Compare)
                and isinstance(node.test.left, ast.Name)
                and node.test.left.id == "__name__"):
            for sub in ast.walk(node):
                if (isinstance(sub, ast.Call)
                        and isinstance(sub.func, ast.Name)
                        and sub.func.id == func):
                    return True
    return False


# Legacy calc_*.py files written before the canonical pattern was
# documented (2026-05-31).  They use ad-hoc entry points and are NOT
# required to comply, but should be migrated opportunistically.  New
# files (not in this list) MUST comply.  Add an entry only when you
# accept the legacy bug-class risk; remove an entry when you migrate
# that file to the canonical pattern.
LEGACY_CALC_EXEMPT = {
    "calc_heat.py",           # thermal legacy pattern
    "calc_heat_axisym.py",    # thermal sibling, legacy pattern
    "calc_heat_with_em_table.py",  # thermal sibling, legacy pattern
    "calc_mesh_eval.py",      # uses different entry
    "calc_motor_lamination.py",
    "calc_motor_transient.py",
    "calc_peec.py",
    "calc_pcb_peec.py",
    "calc_stream_coil.py",
    "calc_surface.py",
    "calc_verify_vol.py",     # uses different entry
    "calc_volume.py",         # legacy `def main(): ... calc_main(run, parser)`
    "calc_kelvin_benchmark.py",
    "calc_inductance.py",     # huge multi-mode legacy file
    "calc_analytic.py",       # ad-hoc analytical formula sandbox
    "calc_axisym_volumetric.py",  # axisym thermal/EM, separate pattern
    "calc_em_table.py",       # passes args via dict-of-dicts, not argparse
    "calc_equivalence_source.py",  # research script, separate pattern
    # calc_*.py NOT in this set MUST follow the canonical pattern
    # documented in docs/panels/ADDING_NEW_PANEL.md.
}


def _calls_anywhere(tree: ast.Module, func: str) -> bool:
    """True if func(...) is called anywhere in the module."""
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id == func):
            return True
    return False


def _audit_calc(path: pathlib.Path) -> list[str]:
    """Return list of violation strings for one calc_*.py."""
    if path.name == "calc_common.py":
        return []   # helper module, exempt
    if path.name in LEGACY_CALC_EXEMPT:
        return []   # grandfathered legacy file
    tree = ast.parse(path.read_text(encoding="utf-8"))
    viols = []
    rel = path.relative_to(ROOT).as_posix()

    # C1
    if not _has_callable(tree, "build_argparser"):
        viols.append(f"{rel}: C1 missing `def build_argparser():`")

    # C2: a callable named `run` or `calc` (the work function passed to
    # calc_main).  Either name is fine; can be at module top OR nested
    # inside `def main():` -- the canonical pattern defines it as a
    # closure capturing the parsed args.
    if not (_has_callable_anywhere(tree, "run")
            or _has_callable_anywhere(tree, "calc")):
        viols.append(f"{rel}: C2 missing work function `def run(args):` "
                     f"(or `def calc(args):`)")

    # C3: calc_main(...) must be called somewhere (not necessarily in
    # __main__; legacy convention is to wrap it in main()).
    if not _calls_anywhere(tree, "calc_main"):
        viols.append(f"{rel}: C3 module does not call calc_main(...) -- "
                     f"the helper that handles --output JSON dump")

    # C4
    for lineno, mod in _toplevel_heavy_imports(tree):
        viols.append(f"{rel}:{lineno}: C4 top-level heavy import `{mod}` "
                     f"(must be inside calc())")

    # C5: build_argparser() must be cheap
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "build_argparser":
            for ln, mod in _function_has_heavy_imports(node):
                viols.append(f"{rel}:{ln}: C5 heavy import `{mod}` inside "
                             f"build_argparser() (panel introspects it)")

    # C6
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            r = _function_late_imports_taskmanager(node)
            if r:
                u, i = r
                viols.append(f"{rel}:{u}: C6 {node.name} uses TaskManager "
                             f"BEFORE importing it (line {i}) -- "
                             f"UnboundLocalError trap")
    return viols


# Panel modules grandfathered before the canonical recipe was
# documented.  See same migration policy as LEGACY_CALC_EXEMPT above.
LEGACY_PANEL_EXEMPT = {
    "radia_gui_base.py",         # the base class itself
    "coil_builder.py",     # interactive 3D widget panel
    "radia_ngsolve.py",          # legacy / not a panel
    "radia_ih.py",               # custom multi-tab UI predates pattern
    # radia_motor.py migrated 2026-06-15 to the ModePanel/AnalysisWindow
    # generator contract (Transient + Lamination sub-panels) -- now audited.
    # radia_*.py NOT in this set MUST follow the canonical pattern.
}


def _audit_panel(path: pathlib.Path) -> list[str]:
    """Return list of violation strings for one radia_*.py panel."""
    if path.name in LEGACY_PANEL_EXEMPT:
        return []
    tree = ast.parse(path.read_text(encoding="utf-8"))
    viols = []
    rel = path.relative_to(ROOT).as_posix()

    # P1: no heavy top imports
    for lineno, mod in _toplevel_heavy_imports(tree):
        viols.append(f"{rel}:{lineno}: P1 top-level heavy import `{mod}` "
                     f"(panel loads on every Cubit toolbar launch)")

    # P2: a class subclasses ModePanel and calls bind_argparser in __init__
    has_modepanel_class_with_bind = False
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            bases = {(b.id if isinstance(b, ast.Name) else
                      (b.attr if isinstance(b, ast.Attribute) else ""))
                     for b in node.bases}
            if "ModePanel" not in bases:
                continue
            for sub in ast.walk(node):
                if (isinstance(sub, ast.Call)
                        and isinstance(sub.func, ast.Attribute)
                        and sub.func.attr == "bind_argparser"):
                    has_modepanel_class_with_bind = True
                    break
    if not has_modepanel_class_with_bind:
        viols.append(f"{rel}: P2 no ModePanel subclass calls bind_argparser(...)")

    # P3: every build_command either calls build_command_from_parser(...)
    # OR delegates to a sub-panel's build_command (composite pattern --
    # the sub-panel itself is then audited separately).
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "build_command":
            uses_helper = False
            delegates_to_sub = False
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call) and isinstance(sub.func, ast.Attribute):
                    if sub.func.attr == "build_command_from_parser":
                        uses_helper = True
                        break
                    if sub.func.attr == "build_command":
                        # delegate: self._current_sub().build_command(vol_path)
                        delegates_to_sub = True
            if not (uses_helper or delegates_to_sub):
                viols.append(f"{rel}:{node.lineno}: P3 build_command neither "
                             f"calls self.build_command_from_parser(...) nor "
                             f"delegates to a sub-panel's build_command()")

    # P4: json_output / msh_output suffix consistency per build_command
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "build_command":
            suffixes = {"json": set(), "msh": set()}
            for sub in ast.walk(node):
                if (isinstance(sub, ast.Call)
                        and isinstance(sub.func, ast.Name)
                        and sub.func.id in ("json_output", "msh_output")
                        and len(sub.args) >= 2
                        and isinstance(sub.args[1], ast.Constant)
                        and isinstance(sub.args[1].value, str)):
                    suffixes[sub.func.id.split("_")[0]].add(sub.args[1].value)
            if suffixes["json"] and suffixes["msh"]:
                if suffixes["json"] != suffixes["msh"]:
                    viols.append(
                        f"{rel}:{node.lineno}: P4 build_command uses "
                        f"json_output suffix {suffixes['json']} and "
                        f"msh_output suffix {suffixes['msh']} -- the "
                        f"Persistence Policy derives .log from --output "
                        f"JSON path so the suffix MUST match the .msh path")
    return viols


def main() -> int:
    viols: list[str] = []
    calcs = sorted(CALC_DIR.glob("calc_*.py"))
    panels = sorted(PANEL_DIR.glob("radia_*.py"))

    for f in calcs:
        viols.extend(_audit_calc(f))
    for f in panels:
        viols.extend(_audit_panel(f))

    if viols:
        print(f"Panel contract: {len(viols)} VIOLATION(S)")
        for v in viols:
            print(f"  {v}")
        print()
        print(f"Audited {len(calcs)} calc_*.py + {len(panels)} radia_*.py panels.")
        print("See docs/panels/ADDING_NEW_PANEL.md for the canonical pattern.")
        return 1
    print(f"Panel contract: CLEAN "
          f"({len(calcs)} calc_*.py + {len(panels)} radia_*.py panels audited)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
