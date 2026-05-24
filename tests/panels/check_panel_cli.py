#!/usr/bin/env python3
"""Static cross-checker between Radia GUI panels and their `calc_*.py`.

Catches four silent-bug classes at the panel/CLI boundary:

    (1) Panel emits `--foo` but no calc_*.py declares it  ->  argparse reject
    (2) calc_*.py declares `--bar` but no panel emits it  ->  silent default
    (3) Panel defines widget `X` that no `build_command` ever reads
                                                        ->  orphan widget
    (4) Panel `_XXX_MAP` dict has a leaf value (CLI-flag-value-shaped
        string like "ams" / "dense-lu") that no calc_*.py argparse
        `choices=[]` accepts                            ->  map-value reject

Class (4) was added 2026-05-24 after the AMG->AMS typo (radia_em.py
mapped "AMG (Compact)" -> "amg" but calc_accel_magnet.py
choices=["auto", "pardiso", "bddc", "iccg", "ams"] rejected it,
producing a subprocess exit 2 every time a user clicked through Omega
formulation with AMG).  Class (1) couldn't catch this because the
panel emits the flag "--solver" itself (which IS accepted) -- the
bug is in the VALUE the panel sent, not the flag name.

Run:

    python tests/panels/check_panel_cli.py
    python tests/panels/check_panel_cli.py --panel radia_ih.py
    python tests/panels/check_panel_cli.py --strict    # fail on (2) too

Exit code:
    0  clean (or only waived issues)
    1  REJECT (panel emits unknown flag, or map-value not in choices)
    2  orphan widget or silent-default (when --strict)

Waivers:
    In a panel source, add a comment like::

        # CLI-DIFF: ignore --reg --shift-eps -- advanced solver knobs

    Listed flags are not counted as silent-default failures.

    For map-value waivers (e.g. a map value is consumed internally,
    not as a CLI flag argument), add::

        # CLI-DIFF: ignore-map-value some-literal-value -- reason

    Listed values are not counted as map-value-reject failures.
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


def scan_calc_choices(path: Path) -> Set[str]:
    """Return the union of every `choices=[...]` literal across all
    `add_argument(...)` calls in *path*.  Used by Check (4) to verify
    panel `_XXX_MAP` leaf values are accepted by at least one CLI flag.

    Each element must be a string Constant; non-string choices (int,
    enum) are skipped silently.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    choices: Set[str] = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"):
            for kw in node.keywords:
                if (kw.arg == "choices"
                        and isinstance(kw.value, (ast.List, ast.Tuple))):
                    for elt in kw.value.elts:
                        s = _str_const(elt)
                        if s is not None:
                            choices.add(s)
    return choices


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

    # Check (4) input: panel class-level _XXX_MAP dicts (any class attr
    # whose RHS is a Dict literal).  Recursively extract every string
    # leaf value -- these are CLI-flag-value candidates that must appear
    # in at least one calc_*.py argparse choices=[] (with the exception
    # of waived values via the ignore-map-value comment).
    map_leaf_values: Dict[str, Set[str]] = {}  # map_name -> {leaf strs}
    for cls in [n for n in ast.walk(tree) if isinstance(n, ast.ClassDef)]:
        for stmt in cls.body:
            # Class-level assignment: `_XXX_MAP = {...}` or annotated
            # form `_XXX_MAP: dict = {...}`.
            if isinstance(stmt, ast.Assign) and len(stmt.targets) == 1:
                tgt, val = stmt.targets[0], stmt.value
            elif isinstance(stmt, ast.AnnAssign) and stmt.value is not None:
                tgt, val = stmt.target, stmt.value
            else:
                continue
            if not (isinstance(tgt, ast.Name)
                    and tgt.id.endswith("_MAP")
                    and isinstance(val, ast.Dict)):
                continue
            leaves: Set[str] = set()

            def _collect(d: ast.Dict) -> None:
                for v in d.values:
                    if isinstance(v, ast.Constant) and isinstance(v.value, str):
                        leaves.add(v.value)
                    elif isinstance(v, ast.Dict):
                        _collect(v)

            _collect(val)
            if leaves:
                map_leaf_values[tgt.id] = leaves

    # Extract map-value waivers: # CLI-DIFF: ignore-map-value foo bar -- reason
    map_value_waivers: Set[str] = set()
    mv_waiver_re = re.compile(
        r"CLI-DIFF:\s*ignore-map-value\s+(\S+(?:\s+\S+)*?)(?:\s*--\s.*)?$",
        re.M)
    for line in text.splitlines():
        m = mv_waiver_re.search(line)
        if m:
            for tok in m.group(1).split():
                if tok and not tok.startswith("--"):
                    map_value_waivers.add(tok)

    return {
        "methods": methods,
        "widgets_defined": widgets_defined,
        "widgets_read": widgets_read,
        "waivers": waivers,
        "map_leaf_values": map_leaf_values,
        "map_value_waivers": map_value_waivers,
    }


# ----------------------------------------------------------------------
# Main check
# ----------------------------------------------------------------------
def check(panels: List[Path], strict: bool = False) -> int:
    # Scan all calc_*.py first
    calc_index: Dict[str, Set[str]] = {}
    calc_choices_union: Set[str] = set()
    for p in sorted(PANELS_DIR.glob("calc_*.py")):
        calc_index[p.name] = scan_calc_cli(p)
        calc_choices_union |= scan_calc_choices(p)

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

        # Check (4): map-value reject — leaf string values in any
        # _XXX_MAP dict that are NOT in any calc argparse choices=[].
        # Heuristic filter: only flag values that look CLI-flag-shaped
        # (lowercase letters / digits / hyphens, no spaces, len 2..20).
        # This skips human-readable labels (e.g. "Direct (PARDISO)" as a
        # key is unaffected; we only inspect VALUES).
        def _cli_value_shaped(s: str) -> bool:
            return (2 <= len(s) <= 20
                    and " " not in s and "(" not in s
                    and any(c.isalpha() for c in s)
                    and all(c.isalnum() or c in "-_." for c in s))

        for map_name, leaves in pinfo["map_leaf_values"].items():
            rejects: List[Tuple[str, str]] = []
            for v in sorted(leaves):
                if not _cli_value_shaped(v):
                    continue
                if v in calc_choices_union:
                    continue
                if v in pinfo["map_value_waivers"]:
                    continue
                rejects.append((map_name, v))
            if rejects:
                print(f"  MAP-VALUE REJECT ({len(rejects):2d}) : "
                      + ", ".join(f"{m}[{v!r}]" for m, v in rejects))
                rc = max(rc, 1)

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
