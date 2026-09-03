"""Audit `with TaskManager():` placement against the
"caller-wraps, helper-does-not" policy (CLAUDE.md 2026-05-27).

Two checks repo-wide:

1. **Helper modules** under ``src/radia/`` (excluding executable panel and
   Simulink operator-assembly adapters) must have **ZERO**
   ``with TaskManager():`` blocks.  Helpers are
   composable building blocks; the caller wraps them.

2. **Caller modules** -- ``src/radia/panels/calc_*.py``,
   ``src/radia/simulink/*_operator_assembly.py``,
   ``validation_test/**.py``, ``tests/**.py``, and Python helpers under
   ``docs/`` -- that contain any
   NGSolve-parallel operation (``.Assemble()``,
   ``.Inverse(..., inverse=...)``, ``mesh.Curve(p)``,
   ``CGSolver(...)``, ``GMRESSolver(...)``,
   ``GridFunction.Set(coefficient, ...)``) MUST have
   at least one caller-owned ``with TaskManager():`` region.

This is intentionally a fast file-level minimum gate.  It detects a caller
that has no TaskManager region at all; it cannot prove that every runtime call
path is covered by the region.  Focused tests and review own that stronger
control-flow check.

Implementation: AST-based.  Comments, docstrings, string literals,
and f-strings containing the parallel-op pattern are NOT flagged.

Exit codes:
    0  All clean.
    1  Violations found (printed to stdout with file:line and the
       suggested fix).
    2  Tool error (path not found, parse error).

Usage::

    python tools/audit_taskmanager.py                 # repo-wide
    python tools/audit_taskmanager.py --quiet-ok      # only print violations
    python tools/audit_taskmanager.py --json          # machine output
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from fnmatch import fnmatch
from pathlib import Path
from typing import NamedTuple

ROOT = Path(__file__).absolute().parent.parent

# Paths classified as HELPERS (must NOT wrap internally).
HELPER_GLOBS = ["src/radia/**/*.py"]
HELPER_EXCLUDE = [
    "src/radia/panels/calc_*.py",
    "src/radia/simulink/*_operator_assembly.py",
    "**/__pycache__/**",
]

# Paths classified as CALLERS (MUST wrap when doing parallel ops).
CALLER_GLOBS = [
    "src/radia/panels/calc_*.py",
    "src/radia/simulink/*_operator_assembly.py",
    "validation_test/**/*.py",
    "tests/**/*.py",
    "docs/**/*.py",
]
CALLER_EXCLUDE = [
    "**/__pycache__/**",
    # Lint fixtures intentionally have bad code — don't gate on them.
    "tests/mcp_server/fixtures/*.py",
    # Lab notebooks: out of scope.
    "**/*.ipynb",
]
LANE_HELPERS = {
    "docs/electric_machine/planar_vim_motor_helpers.py",
    "tests/_ngsolve_2606.py",
    "tests/axifem/_vol_mesh.py",
    "validation_test/cubit/cubit_202512_helpers.py",
    "validation_test/feec/conftest.py",
    "validation_test/stream_function/regcoil_fusion_helpers.py",
}

_CALL_PREFIXES = (
    ".Assemble", ".Curve", ".Inverse", ".Set", ".GenerateMesh",
    "BilinearForm", "LinearForm", "CGSolver", "GMRESSolver",
    "BiCGStabSolver", "MinResSolver", "H1", "HCurl", "HDiv", "L2",
    "VectorH1", "VectorL2", "FacetFESpace", "NumberSpace", "SurfaceL2",
    "HDivDiv", "HDivSurface", "HCurlSurface", "TraceFESpace", "Periodic",
    "Discontinuous", "Compress", "Integrate", "Norm", "Mesh",
)
_CANDIDATE_TOKENS = ("TaskManager",) + tuple(
    f"{prefix}{spacing}("
    for prefix in _CALL_PREFIXES
    for spacing in ("", " ", "\t")
)


class Finding(NamedTuple):
    file: Path
    line: int
    kind: str         # "helper-wraps" | "caller-missing-wrap" | "parse-error"
    snippet: str
    suggested: str


def _matches_any(rel: str, patterns: list[str]) -> bool:
    return any(fnmatch(rel, p) for p in patterns)


def _gather_files() -> list[Path]:
    command = ["git", "grep", "--untracked", "-Il", "-F"]
    for token in _CANDIDATE_TOKENS:
        command.extend(("-e", token))
    command.extend(("--", "src/radia", "validation_test", "tests", "docs"))
    proc = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if proc.returncode in (0, 1):
        return sorted(
            ROOT / relative
            for relative in proc.stdout.splitlines()
            if relative.endswith(".py")
        )

    # Keep the audit usable in source archives that do not contain .git.
    seen: set[Path] = set()
    out: list[Path] = []
    for g in HELPER_GLOBS + CALLER_GLOBS:
        for p in ROOT.glob(g):
            if p.is_file() and p.suffix == ".py" and p not in seen:
                seen.add(p)
                out.append(p)
    return out


def _classify(path: Path) -> str | None:
    """Classify a file as 'caller', 'helper', or None.

    Note: fnmatch does NOT honour `**` recursion, so we use part-based
    matching instead of glob matching for the classification step.
    `_gather_files()` already used `Path.glob('**/...')` to find files,
    so we just have to classify what's already collected.
    """
    rel = path.relative_to(ROOT).as_posix()
    parts = path.relative_to(ROOT).parts
    name = path.name

    # Support modules are called by tests, validation drivers, or notebook
    # cells. Their execution-boundary callers own TaskManager.
    if rel in LANE_HELPERS:
        return "helper"

    # ---- Caller paths (checked first; more specific) ----
    is_panel_caller = (
        len(parts) >= 4
        and parts[0] == "src"
        and parts[1] == "radia"
        and parts[2] == "panels"
        and name.startswith("calc_")
        and name.endswith(".py")
    )
    is_simulink_assembly_caller = (
        len(parts) >= 4
        and parts[0] == "src"
        and parts[1] == "radia"
        and parts[2] == "simulink"
        and name.endswith("_operator_assembly.py")
    )
    is_lane_caller = (
        parts[0] in ("validation_test", "tests", "docs")
        and path.suffix == ".py"
    )
    if (
        is_panel_caller or is_simulink_assembly_caller or is_lane_caller
    ) and not _matches_any(rel, CALLER_EXCLUDE):
        return "caller"

    # ---- Helper paths ----
    if (
        parts[0] == "src"
        and parts[1] == "radia"
        and path.suffix == ".py"
        and not _matches_any(rel, HELPER_EXCLUDE)
    ):
        return "helper"

    return None


# ---------------------------------------------------------------------
# AST helpers
# ---------------------------------------------------------------------

def _is_taskmanager_with(node: ast.With | ast.AsyncWith) -> bool:
    """Return True if any context item in `with X:` calls TaskManager."""
    for item in node.items:
        ctx = item.context_expr
        if isinstance(ctx, ast.Call):
            fn = ctx.func
            if isinstance(fn, ast.Name) and fn.id == "TaskManager":
                return True
            if isinstance(fn, ast.Attribute) and fn.attr == "TaskManager":
                return True
    return False


# NGSolve callables that REQUIRE a TaskManager region around them.
# Per "TaskManager-Only Policy" (CLAUDE.md 2026-05-27), this list is
# intentionally broad: it includes ops that are currently mostly serial
# (FES construction, BilinearForm/LinearForm constructors) so that
# future NGSolve versions adding parallelism to them benefit the caller
# automatically, and so the canonical computation block sits inside
# one `with TaskManager():`.
_NGSOLVE_TM_REQUIRED_NAMES = frozenset({
    # Form assembly + solve
    "BilinearForm", "LinearForm",
    "CGSolver", "GMRESSolver", "BiCGStabSolver", "MinResSolver",
    # FE spaces (incl. transformations + surface variants for BEM)
    "H1", "HCurl", "HDiv", "L2", "VectorH1", "VectorL2",
    "FacetFESpace", "NumberSpace", "SurfaceL2", "HDivDiv",
    "HDivSurface", "HCurlSurface", "TraceFESpace",
    "Periodic", "Discontinuous", "Compress",
    # Integration / norms / inner products
    "Integrate", "Norm",
})

# Attribute-call ops (e.g. `mat.Inverse(...)`, `gf.Set(...)`) handled
# separately below because they need receiver disambiguation.


def _is_parallel_op(node: ast.Call) -> str | None:
    """Return a short label if `node` is a parallel NGSolve op, else None."""
    fn = node.func
    # foo.Assemble() — zero args
    if isinstance(fn, ast.Attribute):
        if fn.attr == "Assemble" and not node.args:
            return ".Assemble()"
        if fn.attr == "Curve" and node.args:
            # mesh.Curve(p) — at least one positional arg
            # Only flag if the receiver looks like a mesh
            rec = fn.value
            if isinstance(rec, ast.Name) and \
               ("mesh" in rec.id.lower() or rec.id == "ngmesh"):
                return "mesh.Curve()"
            if isinstance(rec, ast.Attribute) and \
               ("mesh" in rec.attr.lower()):
                return "mesh.Curve()"
        if fn.attr == "Inverse":
            # mat.Inverse(..., inverse="pardiso") — keyword `inverse=` present
            for kw in node.keywords:
                if kw.arg == "inverse":
                    return ".Inverse(inverse=...)"
        if fn.attr == "Set" and node.args:
            # gf.Set(cf, ...) — at least one positional arg + receiver
            #   that looks like a GridFunction
            rec = fn.value
            if isinstance(rec, ast.Name) and \
               ("gf" in rec.id.lower() or "grid" in rec.id.lower()):
                return "gf.Set(cf, ...)"
            if isinstance(rec, ast.Attribute) and \
               ("gf" in rec.attr.lower() or "grid" in rec.attr.lower()):
                return "gf.Set(cf, ...)"
        # Mesh(geo.GenerateMesh(...)) — flag the inner GenerateMesh
        if fn.attr == "GenerateMesh":
            return ".GenerateMesh()"
    # CGSolver(...) / Integrate(...) / BilinearForm(...) / H1(...) etc.
    if isinstance(fn, ast.Name) and fn.id in _NGSOLVE_TM_REQUIRED_NAMES:
        return f"{fn.id}(...)"
    if isinstance(fn, ast.Attribute) and fn.attr in _NGSOLVE_TM_REQUIRED_NAMES:
        return f".{fn.attr}(...)"
    # ngsolve.Mesh(...) constructor with a Netgen mesh / geometry inside
    if (isinstance(fn, ast.Name) and fn.id == "Mesh" and node.args
            and isinstance(node.args[0], ast.Call)):
        # Only flag `Mesh(<call>)` to avoid catching `Mesh(filename)` literals
        return "Mesh(...)"
    if (isinstance(fn, ast.Attribute) and fn.attr == "Mesh" and node.args
            and isinstance(node.args[0], ast.Call)):
        return "Mesh(...)"
    return None


def _find_function_ranges(tree: ast.AST) -> list[tuple[int, int, str]]:
    """Return list of (start_line, end_line, qualname) for every
    FunctionDef / AsyncFunctionDef in tree."""
    out: list[tuple[int, int, str]] = []

    def _walk(node: ast.AST, prefix: str = ""):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            qual = f"{prefix}{node.name}" if prefix else node.name
            start = node.lineno
            end = getattr(node, "end_lineno", start) or start
            out.append((start, end, qual))
            for child in ast.iter_child_nodes(node):
                _walk(child, qual + ".")
        elif isinstance(node, ast.ClassDef):
            for child in ast.iter_child_nodes(node):
                _walk(child, f"{prefix}{node.name}.")
        else:
            for child in ast.iter_child_nodes(node):
                _walk(child, prefix)

    _walk(tree)
    return out


def _find_tm_wraps(tree: ast.AST) -> list[tuple[int, int]]:
    """Return list of (start_line, end_line) of every `with TaskManager():`
    block in tree."""
    out: list[tuple[int, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.With, ast.AsyncWith)) and \
                _is_taskmanager_with(node):
            start = node.lineno
            end = getattr(node, "end_lineno", start) or start
            out.append((start, end))
    return out


def _find_parallel_ops(tree: ast.AST) -> list[tuple[int, str]]:
    """Return list of (line, op_label) for every parallel op call."""
    out: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            label = _is_parallel_op(node)
            if label:
                out.append((node.lineno, label))
    return out


# ---------------------------------------------------------------------
# Audit routines
# ---------------------------------------------------------------------

def _audit_helper(path: Path) -> list[Finding]:
    """Helper rule: no `with TaskManager():`.  AST-based: docstrings,
    comments, string literals are ignored."""
    findings: list[Finding] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        if "TaskManager" not in text:
            return findings
        tree = ast.parse(text, filename=str(path))
    except SyntaxError as exc:
        return [Finding(path, exc.lineno or 0, "parse-error",
                         str(exc)[:80], "fix syntax")]
    for start, _end in _find_tm_wraps(tree):
        snippet = text.splitlines()[start - 1].strip()[:80]
        findings.append(Finding(
            file=path, line=start, kind="helper-wraps",
            snippet=snippet,
            suggested="DELETE this `with TaskManager():`; the caller wraps.",
        ))
    return findings


def _audit_caller(path: Path) -> list[Finding]:
    """Caller rule (PER-FILE binary):

    If a caller file contains ANY parallel-op call (`.Assemble()`,
    `mesh.Curve(...)`, `mat.Inverse(inverse=...)`, `gf.Set(cf, ...)`,
    `CGSolver(...)`, `GMRESSolver(...)`), the file MUST contain at
    least one `with TaskManager():` block somewhere.

    This is a file-level minimum gate.  The policy is "callers wrap";
    we do NOT infer runtime control flow or police which function contains
    the wrap.  A wrap in an entry function can correctly cover calls into
    local helper functions.

    A single violation is emitted per file (pointing at the first
    parallel op line) when the file has parallel ops but no TM wrap
    at all.
    """
    findings: list[Finding] = []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        if not any(token in text for token in _CANDIDATE_TOKENS):
            return findings
        tree = ast.parse(text, filename=str(path))
    except SyntaxError as exc:
        return [Finding(path, exc.lineno or 0, "parse-error",
                         str(exc)[:80], "fix syntax")]
    parallel_ops = _find_parallel_ops(tree)
    if not parallel_ops:
        return findings   # no parallel work; nothing to require

    tm_wraps = _find_tm_wraps(tree)
    if tm_wraps:
        return findings   # file has at least one wrap; policy satisfied

    # File has parallel ops but zero TM wraps -> single per-file finding
    first_line, first_label = parallel_ops[0]
    snippet = text.splitlines()[first_line - 1].strip()[:80]
    findings.append(Finding(
        file=path, line=first_line, kind="caller-missing-wrap",
        snippet=f"[{first_label}] {snippet}  (+{len(parallel_ops)-1} more)",
        suggested=("Add `with TaskManager():` around the FEM block. "
                   "File has parallel ops but NO TM wrap; per CLAUDE.md "
                   "'Caller Wraps, Helper Does NOT' the caller wraps."),
    ))
    return findings


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quiet-ok", action="store_true",
                    help="suppress 'clean' output")
    p.add_argument("--json", action="store_true",
                    help="emit machine-readable JSON to stdout")
    args = p.parse_args()

    files = _gather_files()
    helpers = [f for f in files if _classify(f) == "helper"]
    callers = [f for f in files if _classify(f) == "caller"]

    helper_findings: list[Finding] = []
    for f in helpers:
        helper_findings.extend(_audit_helper(f))
    caller_findings: list[Finding] = []
    for f in callers:
        caller_findings.extend(_audit_caller(f))

    total = len(helper_findings) + len(caller_findings)

    if args.json:
        out = {
            "n_helpers_scanned": len(helpers),
            "n_callers_scanned": len(callers),
            "helper_violations": [
                {"file": v.file.relative_to(ROOT).as_posix(),
                 "line": v.line, "snippet": v.snippet,
                 "fix": v.suggested}
                for v in helper_findings],
            "caller_violations": [
                {"file": v.file.relative_to(ROOT).as_posix(),
                 "line": v.line, "snippet": v.snippet,
                 "fix": v.suggested}
                for v in caller_findings],
            "total": total,
            "ok": total == 0,
        }
        print(json.dumps(out, indent=2))
        return 0 if total == 0 else 1

    if not args.quiet_ok:
        print(f"audit_taskmanager: scanned "
              f"{len(helpers)} helpers + {len(callers)} callers")
        print("  policy: helper has 0 wraps; every parallel-op-bearing "
              "caller file has a caller-owned region")
        print()

    if helper_findings:
        print(f"== {len(helper_findings)} HELPER VIOLATIONS "
              f"(remove the internal wrap) ==")
        for v in helper_findings:
            rel = v.file.relative_to(ROOT).as_posix()
            print(f"  {rel}:{v.line}  {v.snippet}")
            print(f"      -> {v.suggested}")
        print()

    if caller_findings:
        # Group by file for readability.
        by_file: dict[str, list[Finding]] = {}
        for v in caller_findings:
            by_file.setdefault(
                v.file.relative_to(ROOT).as_posix(), []).append(v)
        print(f"== {len(caller_findings)} CALLER VIOLATIONS in "
              f"{len(by_file)} file(s) ==")
        for fname, vs in sorted(by_file.items()):
            print(f"  {fname}:")
            for v in vs:
                print(f"    L{v.line}: {v.snippet}")
                print(f"        -> {v.suggested}")
        print()

    if total == 0:
        if not args.quiet_ok:
            print("OK: all helpers and callers comply with the "
                  "caller-wraps-helper-does-not policy.")
        return 0
    print(f"FAIL: {total} violation(s).  "
          f"See CLAUDE.md 'TaskManager Wrap Policy: "
          f"Caller Wraps, Helper Does NOT'.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
