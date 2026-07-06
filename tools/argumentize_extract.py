"""
argumentize_extract.py -- retired helper for old desktop panel migration.

Current Radia panels are notebook workbenches backed by DesignSpec classes;
use `validation_test/panels/test_notebook_workbench.py` for the active
contract.  This script remains only for archaeology of removed `radia_*.py`
desktop panels.  It walks a `calc_*.py` and / or an archived `radia_*.py`
panel and reports candidates that would have needed promotion to argparse
arguments before the old ModePanel migration.

Severity:
  [A] panel widget WITHOUT an argparse counterpart -- must add to calc
      before generator migration, else the widget silently does nothing.
  [C] panel-side display-to-value dict / list (e.g. PCB_SOLVERS) --
      duplicate state, argparse should be authoritative; use
      `choices=` (string CLI) or panel `choice_map=` (numeric CLI).
  [E] dispatch branch inside build_command -- per-method conditional
      flag emission.  Resolve via per-method stacked sub-panels OR by
      harmonising argparse across the calc scripts.

Usage:
  python tools/argumentize_extract.py src/radia/panels/calc_em_table.py
  python tools/argumentize_extract.py --panel C:/temp/old_radia_em.py
  python tools/argumentize_extract.py --panel C:/temp/old_radia_em.py \
                                      --calc  src/radia/panels/calc_em_table.py
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path


# ============================================================
# AST extraction helpers
# ============================================================

def extract_argparse_args(path: Path) -> dict:
    """Walk a calc_*.py and pull every parser.add_argument("--flag", ...).

    Returns {dest: {"flag": str, "type": str|None, "default": str|None,
                    "choices": str|None, "required": bool, "line": int}}.
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    args: dict = {}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"):
            continue
        if not (node.args and isinstance(node.args[0], ast.Constant)):
            continue
        flag = node.args[0].value
        if not (isinstance(flag, str) and flag.startswith("--")):
            continue
        dest = flag[2:].replace("-", "_")
        info: dict = {"flag": flag, "line": node.lineno}
        for kw in node.keywords:
            if kw.arg in ("type", "default", "choices", "required",
                          "action", "help", "nargs"):
                try:
                    info[kw.arg] = ast.unparse(kw.value)
                except Exception:
                    info[kw.arg] = "?"
        # Honor explicit dest= override
        if "dest" in info:
            dest = info["dest"].strip("'\"")
        args[dest] = info
    return args


def panel_uses_bind_argparser(path: Path) -> bool:
    """True if the panel calls self.bind_argparser(...), meaning every
    argparse arg becomes a widget implicitly -- [A] / [i] diff against
    the calc no longer applies (the diff is structurally zero).
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "bind_argparser"):
            return True
    return False


def extract_panel_widgets(path: Path) -> list:
    """Walk a panel file and pull every self.add_*(key, label, ...) call.

    Returns [{"key": str, "method": str, "label": str, "line": int}].
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    widgets: list = []
    WIDGET_METHODS = {"add_line", "add_combo", "add_spin", "add_browse",
                      "add_check"}
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "self"
                and node.func.attr in WIDGET_METHODS):
            continue
        method = node.func.attr
        if not (node.args and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)):
            continue
        key = node.args[0].value
        label = ""
        if len(node.args) > 1 and isinstance(node.args[1], ast.Constant):
            label = node.args[1].value
        widgets.append({
            "key": key, "method": method, "label": label,
            "line": node.lineno,
        })
    return widgets


def extract_display_value_consts(path: Path) -> list:
    """Find module-level Display->Value mapping consts (dict or list-of-tuples).

    Returns [{"name": str, "kind": str, "n_items": int, "line": int}].
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    consts: list = []
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        # Convention: module-level UPPER_CASE constants are candidates.
        if not (target.id.isupper() or target.id.startswith("_")):
            continue
        value = node.value
        # dict {"display": value, ...}
        if isinstance(value, ast.Dict) and value.keys:
            if all(isinstance(k, ast.Constant) and isinstance(k.value, str)
                   for k in value.keys):
                consts.append({
                    "name": target.id, "kind": "dict",
                    "n_items": len(value.keys), "line": node.lineno,
                })
        # list [("display", value), ...]
        elif isinstance(value, ast.List) and value.elts:
            if all(isinstance(e, ast.Tuple) and len(e.elts) >= 2
                   for e in value.elts):
                consts.append({
                    "name": target.id, "kind": "list-of-tuples",
                    "n_items": len(value.elts), "line": node.lineno,
                })
    return consts


def extract_build_command_branches(path: Path) -> list:
    """Find `if <cond>: cmd += [flag, ...]` branches inside build_*
    methods.  These are per-method dispatch in the panel-side -- the
    generator-stacked-sub-panel pattern eliminates them.

    Returns [{"func": str, "cond": str, "line": int}].
    """
    tree = ast.parse(path.read_text(encoding="utf-8"))
    branches: list = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef)
                and node.name.startswith(("build_", "_build_"))):
            continue
        for sub in ast.walk(node):
            if not isinstance(sub, ast.If):
                continue
            # Look for `cmd += [...]` or `cmd.extend([...])` inside the if
            has_cmd_append = False
            for st in ast.walk(sub):
                if (isinstance(st, ast.AugAssign)
                        and isinstance(st.target, ast.Name)
                        and st.target.id == "cmd"):
                    has_cmd_append = True
                    break
            if has_cmd_append:
                try:
                    cond = ast.unparse(sub.test)
                except Exception:
                    cond = "?"
                branches.append({
                    "func": node.name, "cond": cond, "line": sub.lineno,
                })
    return branches


# ============================================================
# Report
# ============================================================

def report(panel_path: Path | None, calc_path: Path | None, quiet: bool):
    candidates: list = []  # [{severity, message, ...}]

    panel_widgets: list = []
    calc_args: dict = {}

    if panel_path:
        panel_widgets = extract_panel_widgets(panel_path)
        for c in extract_display_value_consts(panel_path):
            candidates.append({
                "severity": "C",
                "file": panel_path.name, "line": c["line"],
                "message": f"display->value mapping `{c['name']}` "
                           f"({c['kind']}, {c['n_items']} items) -- "
                           "argparse choices= or panel choice_map=",
            })
        for b in extract_build_command_branches(panel_path):
            candidates.append({
                "severity": "E",
                "file": panel_path.name, "line": b["line"],
                "message": f"dispatch branch in {b['func']}: "
                           f"if {b['cond'][:60]}: cmd += [...]",
            })

    if calc_path:
        calc_args = extract_argparse_args(calc_path)

    if panel_path and calc_path:
        # Generator-driven panels (bind_argparser) implicitly expose
        # every argparse arg as a widget -- [A] / [i] diff is moot.
        is_generator_driven = panel_uses_bind_argparser(panel_path)
        if is_generator_driven and not quiet:
            print(f"# (note) {panel_path.name} uses bind_argparser -- "
                  "[A] / [i] diff suppressed (generator covers them).")

        if not is_generator_driven:
            # [A] widget without argparse counterpart
            for w in panel_widgets:
                if w["key"] not in calc_args:
                    candidates.append({
                        "severity": "A",
                        "file": panel_path.name, "line": w["line"],
                        "message": (f"widget {w['key']!r} ({w['method']}) -- "
                                    f"no `--{w['key'].replace('_','-')}` in "
                                    f"{calc_path.name}"),
                    })
            # informational: arg without widget (panel-cli-diff covers
            # this at runtime; surface here for migration planning).
            for dest, info in calc_args.items():
                if dest == "output":
                    continue  # output auto-handled by the workbench
                if not any(w["key"] == dest for w in panel_widgets):
                    candidates.append({
                        "severity": "i",
                        "file": calc_path.name, "line": info["line"],
                        "message": (f"argparse {info['flag']} has no widget "
                                    f"in {panel_path.name}"),
                    })

    # Sort by severity (A > C > E > i) then by file:line
    severity_order = {"A": 0, "C": 1, "E": 2, "i": 3}
    candidates.sort(key=lambda c: (severity_order.get(c["severity"], 9),
                                    c["file"], c["line"]))

    # Print
    if not candidates:
        print("argumentize: no candidates found -- "
              "this (panel, calc) pair is already generator-ready.")
        return 0

    print(f"argumentize: {len(candidates)} candidate(s)")
    if not quiet:
        print("=" * 78)
        print("  [A] widget without argparse counterpart -- must add to calc")
        print("  [C] panel-side display->value map     -- duplicate state")
        print("  [E] dispatch branch in build_command  -- per-method panel split")
        print("  [i] argparse arg without widget       -- informational")
        print("=" * 78)
    for c in candidates:
        print(f"  [{c['severity']}] {c['file']}:{c['line']:>4}  {c['message']}")
    # Exit code: non-zero only if [A] candidates exist (= blockers)
    return 1 if any(c["severity"] == "A" for c in candidates) else 0


# ============================================================
# CLI
# ============================================================

def main():
    ap = argparse.ArgumentParser(
        description="Find variables that should become CLI arguments.",
        epilog="See .claude/skills/argumentize/SKILL.md for the full runbook.")
    ap.add_argument("path", nargs="?",
                    help="single calc_*.py or panel file (auto-detected).")
    ap.add_argument("--panel", type=Path,
                    help="panel file (radia_*.py)")
    ap.add_argument("--calc", type=Path,
                    help="calc file (calc_*.py)")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="omit the legend header.")
    args = ap.parse_args()

    panel_path: Path | None = args.panel
    calc_path: Path | None = args.calc
    if args.path and not (panel_path or calc_path):
        # Single-file mode: auto-detect
        p = Path(args.path)
        if p.name.startswith("calc_"):
            calc_path = p
        else:
            panel_path = p
    if not (panel_path or calc_path):
        ap.print_help()
        sys.exit(2)
    sys.exit(report(panel_path, calc_path, args.quiet))


if __name__ == "__main__":
    main()
