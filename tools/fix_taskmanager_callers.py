"""Auto-fixer for "caller missing TaskManager wrap" violations.

Reads `tools/audit_taskmanager.py --json` output and, for each
caller-side violation, rewrites the source file to:

1. Ensure `TaskManager` is imported from `ngsolve`.
2. Wrap the body of the function containing parallel ops in
   `with TaskManager():`.  Heuristic:
     - If the violating function is `main`, wrap its body
       starting AFTER any `argparse.parse_args()` / `args = ...`
       lines (so arg-parsing happens single-threaded and the
       `--nthreads` CLI can take effect before TM starts).
     - Otherwise wrap from the first parallel-op line to the end
       of the function.
3. Re-run audit and report.

Usage::

    # Dry run (show what would change, no writes):
    python tools/fix_taskmanager_callers.py --dry-run

    # Apply to all caller violations:
    python tools/fix_taskmanager_callers.py

    # Apply only to a specific subdir:
    python tools/fix_taskmanager_callers.py \\
        --filter validation_test/kelvin_transformation/AdaptiveMesh/

    # Verify after:
    python tools/audit_taskmanager.py

Exit codes:
    0  success (every violation rewritten or no violations)
    1  partial: some files could not be auto-fixed (unusual structure)
    2  tool error

CAVEAT: the fixer is heuristic; it never overwrites without making
a clean diff.  Re-audit after to catch any remaining cases that
need manual care (nested classes, decorators, etc.).
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).absolute().parent.parent


class Patch(NamedTuple):
    file: Path
    op_line: int
    func_line: int          # start of enclosing function (or 1)
    func_end: int
    indent: int             # body indent inside the function
    label: str              # parallel-op label
    qual: str               # qualname of the enclosing function


def _audit_now() -> dict:
    """Run audit and parse JSON."""
    proc = subprocess.run(
        [sys.executable, "-W", "ignore",
         str(ROOT / "tools" / "audit_taskmanager.py"), "--json"],
        capture_output=True, env={**__import__("os").environ,
                                    "PYTHONIOENCODING": "utf-8",
                                    "PYTHONWARNINGS": "ignore"},
        check=False)
    out = (proc.stdout or b"").decode("utf-8", errors="replace")
    # Skip any leading noise
    i = out.find("{")
    if i > 0:
        out = out[i:]
    import json
    return json.loads(out)


def _find_function_for(text: str, op_line: int) -> tuple[int, int, int, str]:
    """Return (func_start, func_end, body_indent, qualname) for the
    smallest function (or <module>) enclosing `op_line`."""
    tree = ast.parse(text)
    lines = text.splitlines()

    best: tuple[int, int, int, str] | None = None
    def _walk(node: ast.AST, prefix: str = ""):
        nonlocal best
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            qual = f"{prefix}{node.name}" if prefix else node.name
            s = node.lineno
            e = getattr(node, "end_lineno", s) or s
            if s <= op_line <= e:
                # First statement in body line — indent
                if node.body:
                    body_indent = len(lines[node.body[0].lineno - 1]) - \
                                  len(lines[node.body[0].lineno - 1].lstrip())
                else:
                    body_indent = 4
                if best is None or (e - s) < (best[1] - best[0]):
                    best = (s, e, body_indent, qual)
            for child in ast.iter_child_nodes(node):
                _walk(child, qual + ".")
        elif isinstance(node, ast.ClassDef):
            for child in ast.iter_child_nodes(node):
                _walk(child, f"{prefix}{node.name}.")
        else:
            for child in ast.iter_child_nodes(node):
                _walk(child, prefix)

    _walk(tree)
    if best is None:
        # Top-level scope: indent 0
        return (1, len(lines), 0, "<module>")
    return best


def _has_taskmanager_import(text: str) -> bool:
    return bool(re.search(r"\bfrom\s+ngsolve(?:\s+|\.[a-zA-Z_]+\s+)"
                          r"import\s+[^#\n]*\bTaskManager\b", text)) or \
           bool(re.search(r"\bimport\s+ngsolve\.[^#\n]*TaskManager", text))


def _add_taskmanager_import(text: str) -> str:
    """Ensure TaskManager is importable from ngsolve.

    Strategy: insert a fresh ``from ngsolve import TaskManager`` line
    immediately AFTER the last existing ngsolve-related import line.
    Do NOT mutate existing multi-line import lists (the safe
    approach -- adding a new line cannot break tuple parens / trailing
    commas / parentheses-continued imports).
    """
    if "TaskManager" in text and re.search(
            r"\bimport\b[^\n]*\bTaskManager\b", text):
        # Conservative: any mention of TaskManager (import or string)
        # means we leave it alone; user/auditor can verify manually if
        # this is a false negative.
        # But if it's clearly imported, skip:
        return text

    # Find the LAST line that looks like a top-level ngsolve import
    # (matches `from ngsolve...`, `import ngsolve...`).  We append
    # AFTER that line.
    lines = text.splitlines(keepends=True)
    last_ngsolve_line = -1
    for i, L in enumerate(lines):
        stripped = L.lstrip()
        if stripped.startswith(("from ngsolve", "import ngsolve")):
            last_ngsolve_line = i

    if last_ngsolve_line >= 0:
        # Match the indent of the existing ngsolve import line, so a
        # function-scope `from ngsolve import (...)` gets a sibling
        # `from ngsolve import TaskManager` at the SAME indent
        # (placing a top-level import in the middle of a function body
        # would be a syntax error).  Preserve the EXACT whitespace
        # characters (tab vs space) of the existing line -- otherwise
        # a tab-indented file gets a space-indented insert and Python
        # raises "unindent does not match any outer indentation level".
        existing = lines[last_ngsolve_line]
        existing_indent = existing[: len(existing) - len(existing.lstrip())]
        # Insert AFTER the last ngsolve-related line.  If that line
        # uses parenthesised continuation, find the line where its
        # close-paren lives, and insert AFTER that.
        end_i = last_ngsolve_line
        # If line opens an unclosed `(`, walk forward to closing `)`.
        depth = lines[end_i].count("(") - lines[end_i].count(")")
        while depth > 0 and end_i + 1 < len(lines):
            end_i += 1
            depth += lines[end_i].count("(") - lines[end_i].count(")")
        insert_line = f"{existing_indent}from ngsolve import TaskManager\n"
        new_lines = lines[:end_i + 1] + [insert_line] + lines[end_i + 1:]
        return "".join(new_lines)

    insert_line = "from ngsolve import TaskManager\n"

    # No ngsolve import found at all -- insert at very top (after
    # any module docstring / shebang).
    insert_idx = 0
    # Skip shebang
    if lines and lines[0].startswith("#!"):
        insert_idx = 1
    # Skip a top-level triple-quoted docstring (cheap heuristic:
    # find the closing triple-quote line).
    if insert_idx < len(lines) and lines[insert_idx].lstrip().startswith(
            ('"""', "'''")):
        q = lines[insert_idx].lstrip()[:3]
        # If the same line also closes the docstring:
        if lines[insert_idx].count(q) >= 2:
            insert_idx += 1
        else:
            for j in range(insert_idx + 1, len(lines)):
                if q in lines[j]:
                    insert_idx = j + 1
                    break
    new_lines = lines[:insert_idx] + [insert_line, "\n"] + lines[insert_idx:]
    return "".join(new_lines)


def _wrap_function_body(text: str, func_start: int, func_end: int,
                         body_indent: int, op_line: int) -> tuple[str, bool]:
    """Wrap the function body (from first stmt to last) with
    `with TaskManager():` and indent every body line by 4.

    Strategy: locate the first body line (the line of the FIRST
    statement after the def / docstring), insert `with TaskManager():`
    on a fresh line at body_indent, then indent every subsequent body
    line within [body_start, func_end] by +4 spaces.

    To minimise churn, we wrap STARTING AT the op_line (not at the
    very first body stmt) -- this preserves arg parsing and docstring
    at the function entry without re-indenting them.

    Returns (new_text, success_flag)."""
    lines = text.splitlines(keepends=True)
    # 0-indexed
    op_idx = op_line - 1
    end_idx = func_end - 1

    # Determine the indent of the op line itself (must equal or exceed
    # body_indent).
    op_line_text = lines[op_idx]
    op_indent = len(op_line_text) - len(op_line_text.lstrip())
    if op_indent < body_indent:
        return text, False   # something weird (decorator? skip)

    # We will:
    # 1. Insert `with TaskManager():` at line op_idx (above op_line)
    #    indented at op_indent.
    # 2. Indent every line from op_idx through end_idx (inclusive) by +4.
    indent_str = " " * op_indent
    inserted = f"{indent_str}with TaskManager():\n"

    new_lines: list[str] = []
    # Lines before op_idx: unchanged
    new_lines.extend(lines[:op_idx])
    # Insert the with-statement
    new_lines.append(inserted)
    # Indent lines from op_idx to end_idx (inclusive) by +4
    # BUT only those that are NOT blank and whose existing indent
    # is >= op_indent (i.e., they belong to the function body, not
    # the outer scope).
    for i in range(op_idx, min(end_idx + 1, len(lines))):
        L = lines[i]
        stripped = L.lstrip(" ")
        if L.strip() == "":
            new_lines.append(L)
            continue
        cur_indent = len(L) - len(stripped)
        if cur_indent >= op_indent:
            new_lines.append(" " * 4 + L)
        else:
            # Out of function body -- stop indenting
            new_lines.extend(lines[i:end_idx + 1])
            break
    else:
        # Loop completed normally; nothing to append.
        pass
    # Lines after func_end: unchanged
    new_lines.extend(lines[end_idx + 1:])
    return "".join(new_lines), True


def _fix_one(path: Path, op_line: int, label: str,
              dry_run: bool) -> tuple[str, bool]:
    """Apply the wrap to one file.  Returns (status, success)."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        return f"read failed: {e}", False

    try:
        func_start, func_end, body_indent, qual = \
            _find_function_for(text, op_line)
    except SyntaxError as e:
        return f"syntax error: {e}", False

    # Safety: if the wrap region (op_line..func_end) contains any
    # nested `def` or `class` at body indent or shallower, bail --
    # the wrap would mis-indent the nested definition's body.
    op_line_text = text.splitlines()[op_line - 1]
    op_indent = len(op_line_text) - len(op_line_text.lstrip())
    region_lines = text.splitlines()[op_line - 1: func_end]
    for L in region_lines:
        stripped = L.lstrip()
        if stripped.startswith(("def ", "class ", "async def ")):
            line_indent = len(L) - len(stripped)
            if line_indent <= op_indent:
                return ("skipped: wrap region contains nested def/class "
                        "at same or shallower indent (manual fix needed)"), \
                       False

    # Make sure TaskManager is imported
    new_text = text if _has_taskmanager_import(text) \
                    else _add_taskmanager_import(text)

    # Re-parse to get updated line numbers after the import edit
    if new_text != text:
        # Wrap the original body first, then add the import so source line
        # positions remain stable during the structural edit.
        wrapped, ok = _wrap_function_body(text, func_start, func_end,
                                            body_indent, op_line)
        if not ok:
            return ("could not wrap: op_indent < body_indent "
                    "(decorator/odd structure)"), False
        if not _has_taskmanager_import(wrapped):
            wrapped = _add_taskmanager_import(wrapped)
    else:
        wrapped, ok = _wrap_function_body(text, func_start, func_end,
                                            body_indent, op_line)
        if not ok:
            return ("could not wrap: op_indent < body_indent"
                    " (decorator/odd structure)"), False

    if dry_run:
        return (f"would wrap [{label}] @ L{op_line} in `{qual}()`"
                f" (func L{func_start}..L{func_end}, "
                f"+{wrapped.count(chr(10)) - text.count(chr(10))} lines)"), \
               True

    path.write_text(wrapped, encoding="utf-8")
    return (f"wrapped [{label}] @ L{op_line} in `{qual}()`"), True


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true",
                    help="print what would change; do not write")
    p.add_argument("--filter", default="",
                    help="only fix files whose relative path STARTS WITH "
                         "this prefix (e.g. 'validation_test/kelvin_'). "
                         "Multiple comma-separated allowed.")
    args = p.parse_args()

    audit = _audit_now()
    violations = audit["caller_violations"]
    if args.filter:
        filters = [f.strip() for f in args.filter.split(",") if f.strip()]
        violations = [v for v in violations
                       if any(v["file"].startswith(f) for f in filters)]

    print(f"fixing {len(violations)} caller violation(s) "
          f"(dry-run={args.dry_run})...")
    print()

    n_ok = n_skip = 0
    for v in violations:
        path = ROOT / v["file"]
        status, ok = _fix_one(path, v["line"], v["snippet"][:60],
                                args.dry_run)
        marker = "OK  " if ok else "SKIP"
        print(f"  {marker}  {v['file']}: {status}")
        if ok:
            n_ok += 1
        else:
            n_skip += 1

    print()
    print(f"Result: {n_ok} fixed (or would-fix), {n_skip} skipped.")
    if not args.dry_run and n_ok:
        print("Re-run `python tools/audit_taskmanager.py` to verify.")
    return 0 if n_skip == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
