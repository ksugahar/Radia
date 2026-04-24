#!/usr/bin/env python3
"""Static cross-checker between Radia GUI panels and their `calc_*.py`.

Catches three silent-bug classes at the panel/CLI boundary:

    (1) Panel emits `--foo` but no calc_*.py declares it  ->  argparse reject
    (2) calc_*.py declares `--bar` but no panel emits it  ->  silent default
    (3) Panel defines widget `X` that no `build_command` ever reads
                                                        ->  orphan widget

Run:

    python tests/panels/check_panel_cli.py
    python tests/panels/check_panel_cli.py --panel radia_ih.py
    python tests/panels/check_panel_cli.py --strict    # fail on (2) too

Exit code:
    0  clean (or only waived issues)
    1  REJECT (panel emits unknown flag)
    2  orphan widget or silent-default (when --strict)

Waivers:
    In a panel source, add a comment like::

        # CLI-DIFF: ignore --reg --shift-eps -- advanced solver knobs

    Listed flags are not counted as silent-default failures.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

ROOT = Path(__file__).resolve().parent.parent.parent
PANELS_DIR = ROOT / "src" / "radia" / "panels"
RADIA_DIR = ROOT / "src" / "radia"

# Helper functions that add CLI args to a parser indirectly.
# Map: qualified helper name -> list of flags it adds.
INDIRECT_HELPERS: Dict[str, List[str]] = {
    "add_material_args": [
        "--material", "--sigma", "--mu-r", "--bh-file",
        # optional include_hys=True adds --hys-file; we pick it up from
        # the call's keyword args below.
    ],
}


def _str_const(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


# ----------------------------------------------------------------------
# calc_*.py scanner
# ----------------------------------------------------------------------
def scan_calc_cli(path: Path) -> Set[str]:
    """Return the set of '--flag' strings accepted by argparse in *path*."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    flags: Set[str] = set()
    for node in ast.walk(tree):
        # parser.add_argument("--foo", ...)
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"
                and node.args):
            s = _str_const(node.args[0])
            if s and s.startswith("--"):
                flags.add(s)
        # add_material_args(parser, include_hys=True, ...) and similar
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id in INDIRECT_HELPERS):
            flags.update(INDIRECT_HELPERS[node.func.id])
            # include_hys=True adds --hys-file
            for kw in node.keywords:
                if (kw.arg == "include_hys"
                        and isinstance(kw.value, ast.Constant)
                        and kw.value.value is True):
                    flags.add("--hys-file")
    return flags


# ----------------------------------------------------------------------
# Panel _build_*_command() scanner
# ----------------------------------------------------------------------
_CALC_SCRIPT_RE = re.compile(r'calc_script\(\s*["\']([^"\']+)["\']\s*\)')


def _collect_flag_literals(func_node: ast.FunctionDef) -> Set[str]:
    """Walk the function body and collect every `"--flag"` string literal
    that appears in a list / list-extension / tuple — anywhere cmd-building
    code typically puts them.
    """
    flags: Set[str] = set()
    for sub in ast.walk(func_node):
        if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
            if sub.value.startswith("--"):
                flags.add(sub.value)
    return flags


def _find_calc_script_target(func_node: ast.FunctionDef) -> str | None:
    """Find `calc_script("xxx.py")` inside a function body."""
    src = ast.unparse(func_node) if hasattr(ast, "unparse") else ""
    m = _CALC_SCRIPT_RE.search(src)
    return m.group(1) if m else None


def scan_panel(path: Path) -> Dict:
    """Return a dict describing each build_command method + widget inventory."""
    text = path.read_text(encoding="utf-8")
    tree = ast.parse(text)

    methods: Dict[str, Dict] = {}  # method_name -> {calc, emits}
    widgets_defined: Set[str] = set()
    widgets_read: Set[str] = set()
    # var_to_widget_key: `coil_mat = self.add_combo("coil_material", ...)`
    # -> {"coil_mat": "coil_material"}
    var_to_widget_key: Dict[str, str] = {}
    used_vars: Set[str] = set()

    for node in ast.walk(tree):
        # class method def — both the per-method _build_X_command pattern
        # (IH with 3 methods) and the simpler overriding build_command
        # (single-method panels like radia_pcb.py).
        if (isinstance(node, ast.FunctionDef)
                and ((node.name.startswith("_build_")
                      and node.name.endswith("_command"))
                     or node.name == "build_command")):
            calc = _find_calc_script_target(node)
            emits = _collect_flag_literals(node)
            if calc is not None:
                methods[node.name] = {"calc": calc, "emits": emits}

        # widget definitions: self.add_line("key", ...), add_combo, add_spin, add_browse, add_check
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in ("add_line", "add_combo", "add_spin",
                                        "add_browse", "add_check")
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "self"
                and node.args):
            k = _str_const(node.args[0])
            if k:
                widgets_defined.add(k)

        # Track assignments like `coil_mat = self.add_combo("coil_material", ...)`
        # so that downstream `coil_mat.currentTextChanged.connect(...)` counts
        # as the widget being used.
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            tgt = node.targets[0]
            val = node.value
            if (isinstance(tgt, ast.Name)
                    and isinstance(val, ast.Call)
                    and isinstance(val.func, ast.Attribute)
                    and val.func.attr in ("add_line", "add_combo", "add_spin",
                                           "add_browse", "add_check")
                    and isinstance(val.func.value, ast.Name)
                    and val.func.value.id == "self"
                    and val.args):
                k = _str_const(val.args[0])
                if k:
                    var_to_widget_key[tgt.id] = k

        # Anything like `foo.method(...)` or `foo.attr` uses var `foo` —
        # if foo is a tracked widget-var, mark its key as read.
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
            used_vars.add(node.value.id)

        # widget reads: self.val("key") or self._widgets.get("key")
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "val"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "self"
                and node.args):
            k = _str_const(node.args[0])
            if k:
                widgets_read.add(k)
        # self._widgets.get("key") / self._widgets["key"]
        if (isinstance(node, ast.Subscript)
                and isinstance(node.value, ast.Attribute)
                and node.value.attr == "_widgets"):
            # subscript may be Index(Constant) (py38) or Constant (py39+)
            sub = node.slice
            if isinstance(sub, ast.Constant):
                k = _str_const(sub)
                if k:
                    widgets_read.add(k)
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "get"
                and isinstance(node.func.value, ast.Attribute)
                and node.func.value.attr == "_widgets"
                and node.args):
            k = _str_const(node.args[0])
            if k:
                widgets_read.add(k)

    # Also pick up self._set_row_visible("key", ...) as a read (visibility
    # toggle implies the widget is used by that method).
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "_set_row_visible"
                and node.args):
            k = _str_const(node.args[0])
            if k:
                widgets_read.add(k)

    # Promote: if `coil_mat = self.add_combo(...)` and `coil_mat` is used
    # anywhere (e.g. `.currentTextChanged.connect(...)`), the widget is
    # considered read.  This covers signal-driven widgets that never go
    # through `self.val()`.
    for var, key in var_to_widget_key.items():
        if var in used_vars:
            widgets_read.add(key)

    # Extract waivers: # CLI-DIFF: ignore --a --b -- reason
    waivers: Set[str] = set()
    waiver_re = re.compile(
        r"CLI-DIFF:\s*ignore\s+((?:--[\w-]+\s*)+)(?:--\s*.*)?", re.I)
    for line in text.splitlines():
        m = waiver_re.search(line)
        if m:
            for tok in m.group(1).split():
                if tok.startswith("--"):
                    waivers.add(tok)

    return {
        "methods": methods,
        "widgets_defined": widgets_defined,
        "widgets_read": widgets_read,
        "waivers": waivers,
    }


# ----------------------------------------------------------------------
# Main check
# ----------------------------------------------------------------------
def check(panels: List[Path], strict: bool = False) -> int:
    # Scan all calc_*.py first
    calc_index: Dict[str, Set[str]] = {}
    for p in sorted(PANELS_DIR.glob("calc_*.py")):
        calc_index[p.name] = scan_calc_cli(p)

    rc = 0
    for panel in panels:
        pinfo = scan_panel(panel)
        print(f"\n=== {panel.relative_to(ROOT)} ===")

        for mname, minfo in pinfo["methods"].items():
            calc = minfo["calc"]
            if calc is None:
                print(f"  method {mname}  (no calc_script target found)")
                continue
            accepted = calc_index.get(calc)
            if accepted is None:
                print(f"  method {mname} -> {calc}  MISSING CALC FILE")
                rc = max(rc, 1)
                continue
            emits = minfo["emits"]
            rejects = emits - accepted
            silent_defaults = accepted - emits - pinfo["waivers"]
            overlap = emits & accepted

            print(f"  method {mname:40s} -> {calc}")
            print(f"    OK ({len(overlap):2d})          : "
                  + " ".join(sorted(overlap)))
            if rejects:
                print(f"    REJECT ({len(rejects):2d})      : "
                      + " ".join(sorted(rejects)))
                rc = max(rc, 1)
            if silent_defaults:
                tag = "SILENT-DEF" if not strict else "SILENT-DEF (FAIL)"
                print(f"    {tag} ({len(silent_defaults):2d}) : "
                      + " ".join(sorted(silent_defaults)))
                if strict:
                    rc = max(rc, 2)

        orphans = pinfo["widgets_defined"] - pinfo["widgets_read"]
        if orphans:
            print(f"  ORPHAN WIDGETS ({len(orphans)}): "
                  + " ".join(sorted(orphans)))
            rc = max(rc, 2)

    print("\n---")
    if rc == 0:
        print("panel-cli-diff: CLEAN")
    else:
        print(f"panel-cli-diff: ISSUES (rc={rc})")
    return rc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", action="append", default=[],
                    help="Panel file(s) to check (default: all radia_*.py "
                         "under src/radia/)")
    ap.add_argument("--strict", action="store_true",
                    help="Treat silent-default + orphan widget as a "
                         "failure (exit rc=2)")
    args = ap.parse_args()

    if args.panel:
        panels = [RADIA_DIR / p for p in args.panel]
    else:
        panels = [p for p in sorted(RADIA_DIR.glob("radia_*.py"))
                  if p.name != "radia_gui_base.py"
                  and p.name != "radia_coil_builder.py"
                  and p.name != "radia_ngsolve.py"]

    sys.exit(check(panels, strict=args.strict))


if __name__ == "__main__":
    main()
