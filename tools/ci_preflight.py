#!/usr/bin/env python
"""tools/ci_preflight.py -- run the CI gates LOCALLY, BEFORE you push.

Make "commit -> CI check -> push" a single command so CI does not
discover a failure you could have caught locally in ~2 minutes.

WHY THIS EXISTS (empirical, 2026-06-05 analysis of the last 80 failed
GitHub Actions runs, 2026-05-25..06-04):

    30x  radia-mcp matrix  "Pytest (meta health + ... + versioning)"
    16x  CI                "Run basic tests"
     9x  radia-mcp matrix  "TOOLS.md drift gate"
     7x  Policy Lint       "Policy 4: No CblasColMajor in core"
     3x  radia-mcp matrix  "Meta health"

EVERY one of these is detectable locally before push.  This script runs
the same gates the three CI workflows run (policy-lint.yml,
radia-mcp-matrix.yml, build-test.yml "Run basic tests"), fast-first.

Gates:
  1. policy-lint        -- the 7 static policies (CblasColMajor etc.)
  2. version-consistency-- pyproject == __init__ for radia / radia-mcp / cme
  3. tools-md-drift     -- regenerate docs/TOOLS.md and diff (WIP-aware)
  4. radia-mcp-matrix   -- compile + meta_health + the FULL radia-mcp
                           pytest run UNDER minimal-dep simulation
                           (RADIA_MCP_FORCE_MINIMAL=1), which is what the
                           ubuntu matrix actually runs -- this is the gate
                           that catches heavy-import collection errors
                           (ngsolve/netgen imported at module level)
  5. toplevel-collect   -- `pytest tests/ --collect-only` with the CI
                           ignore-set: catches import errors that break
                           the self-hosted "Run basic tests" at collection
  6. toplevel-run       -- (only with --full) actually run top-level
                           tests/ with the CI ignore-set + reruns

Usage:
    python tools/ci_preflight.py            # gates 1-5 (~2 min)
    python tools/ci_preflight.py --full     # also gate 6 (slow: FEM solves)
    python tools/ci_preflight.py --fix      # auto-regenerate TOOLS.md on drift
    python tools/ci_preflight.py --only policy,tools-md

Exit 0 = all green, safe to push.  Non-zero = a gate CI would fail is red.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

# --- UNC-safe repo root -------------------------------------------------
# NEVER use Path.resolve() here: on the LAB box the repo lives on an S:
# drive mapped to \\192.168.11.100\..., and resolve() canonicalises to the
# UNC form, which Python's open() then rejects with OSError [Errno 22]
# (this is exactly what broke `release_triple.py preflight`).  os.path
# joins keep the drive letter the script was invoked with.
_THIS = os.path.abspath(__file__)
REPO = os.path.dirname(os.path.dirname(_THIS))
# LAB: the repo is an S: drive mapped onto a \\192.168.11.100\... UNC share.
# When invoked from the pre-push hook, git's cwd (and abspath(__file__)) come
# back in the UNC form, and Python file reads on that form fail.  Remap to the
# S: drive (what every manual run uses + what open() handles) by anchoring on
# the repo's known path suffix -- robust even when REPO == the UNC base
# exactly (the earlier `len(prefix)` slice produced a bare "S:").
_norm = REPO.replace("\\", "/")
_anchor = "/Radia/01_GitHub"
if "192.168.11.100" in _norm and _anchor in _norm:
    REPO = "S:" + _norm[_norm.index(_anchor):]
MCP = os.path.join(REPO, "packages", "radia-mcp")

# CI ignore-set for the self-hosted "Run basic tests" step
# (mirror .github/workflows/build-test.yml).
TOPLEVEL_IGNORES = [
    "tests/cubit", "tests/panels",
    "tests/test_far_field_accuracy.py",
    "tests/test_rad_ngsolve_function.py",
    "tests/test_tetrahedral_solver.py",
    "tests/test_batch_evaluation.py",
    "tests/test_curlA_equals_B.py",
    "tests/test_mesh_import.py",
    "tests/test_scalar_bie_sibc.py",
]

# Policy 4 allowlist (genuine LAPACK / HACApK column-major interop).
CBLAS_COLMAJOR_ALLOW = (
    "rad_mmm_matrices.cpp", "rad_relaxation_methods.cpp",
    "rad_hacapk.cpp", "rad_stream_function.cpp",
)

GREEN, RED, YEL, DIM, RST = "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"


def _sh(cmd, cwd=REPO, env=None, timeout=None):
    """Run a command, return (returncode, stdout+stderr)."""
    e = dict(os.environ)
    if env:
        e.update(env)
    p = subprocess.run(cmd, cwd=cwd, env=e, capture_output=True,
                       text=True, encoding="utf-8", errors="replace",
                       timeout=timeout)
    return p.returncode, (p.stdout or "") + (p.stderr or "")


def _git_grep(pattern, pathspecs, extra=()):
    """git grep -n in the working tree; return list of matching lines."""
    cmd = ["git", "grep", "-n", *extra, pattern, "--", *pathspecs]
    rc, out = _sh(cmd)
    # git grep: 0 = matches, 1 = no matches, >1 = error
    if rc not in (0, 1):
        return None  # error (e.g. pathspec doesn't exist)
    return [ln for ln in out.splitlines() if ln.strip()]


# ======================================================================
# Gate 1: policy lint  (mirror .github/workflows/policy-lint.yml)
# ======================================================================
def gate_policy_lint():
    fails = []

    # Policy 1: no FldUnits() in examples/*.py
    hits = _git_grep("FldUnits", ["examples/*.py", "examples/**/*.py"])
    if hits:
        fails.append("Policy 1 (FldUnits removed): " + hits[0])

    # Policy 2: no tracked binaries
    rc, out = _sh(["git", "ls-files", "*.pyd", "*.dll", "*.so",
                   "*.lib", "*.exe", "*.obj"])
    bins = [b for b in out.splitlines() if b.strip()]
    if bins:
        fails.append(f"Policy 2 (no tracked binaries): {bins[:3]}")

    # Policy 3: no Helmholtz WAVE kernel in src/core (hodge/decomp excepted)
    hits = _git_grep("helmholtz", ["src/core/*.cpp", "src/core/*.h",
                                   "src/core/**/*.cpp", "src/core/**/*.h"],
                     extra=["-i"])
    if hits:
        bad = [h for h in hits
               if "helmholtz-hodge" not in h.lower()
               and "helmholtz hodge" not in h.lower()
               and "helmholtz-decomposition" not in h.lower()
               and "helmholtz decomposition" not in h.lower()]
        if bad:
            fails.append("Policy 3 (no Helmholtz wave kernel): " + bad[0])

    # Policy 4: no CblasColMajor in src/core (LAPACK/HACApK allowlisted)
    hits = _git_grep("CblasColMajor", ["src/core/*.cpp", "src/core/*.h",
                                       "src/core/**/*.cpp", "src/core/**/*.h"])
    if hits:
        bad = [h for h in hits
               if not any(a in h for a in CBLAS_COLMAJOR_ALLOW)]
        if bad:
            fails.append("Policy 4 (use CblasRowMajor): " + bad[0])

    # Policy 5: no generated files at repo root.  TRACKED only -- CI runs on
    # a fresh checkout, so untracked local scratch (.msh/.vtu left by a demo
    # run) is NOT what the workflow's `find` sees; checking os.listdir here
    # would be a false positive for "will CI fail".
    rc, out = _sh(["git", "ls-files", "*.msh", "*.vtu", "*.vtk",
                   "*.vol", "*.vts"])
    root_gen = [f for f in out.splitlines() if f.strip() and "/" not in f.strip()]
    if root_gen:
        fails.append(f"Policy 5 (no tracked generated files at root): {root_gen}")

    # Policy 6: forbid the legacy import path (use src/radia instead).
    # Build the needle from fragments so this gate never flags its OWN
    # source: CI's Policy 6 greps for that literal token together with a
    # sys.path line, and a contiguous occurrence here would self-match.
    _legacy = "src/" + "python"
    hits = _git_grep(_legacy, ["*.py", "**/*.py"])
    if hits:
        bad = [h for h in hits if "sys.path" in h]
        if bad:
            fails.append(f"Policy 6 (use src/radia not {_legacy}): " + bad[0])

    # Policy 7: every examples/*/ has README.md
    ex = os.path.join(REPO, "examples")
    missing = []
    if os.path.isdir(ex):
        for d in sorted(os.listdir(ex)):
            dp = os.path.join(ex, d)
            if os.path.isdir(dp) and not os.path.exists(os.path.join(dp, "README.md")):
                missing.append(f"examples/{d}/")
    if missing:
        fails.append(f"Policy 7 (README.md per example dir): {missing[:3]}")

    return (not fails), "; ".join(fails) if fails else "7 policies pass"


# ======================================================================
# Gate 2: version consistency  (pyproject == __init__)
# ======================================================================
def _ver_in(path, key="version"):
    import re
    try:
        with open(path, encoding="utf-8") as f:
            txt = f.read()
    except OSError as e:
        return None, f"unreadable: {e}"
    m = re.search(r'(?:^version|^__version__)\s*=\s*"([^"]+)"', txt, re.M)
    return (m.group(1) if m else None), None


def gate_version_consistency():
    checks = [
        ("radia", "pyproject.toml", "src/radia/__init__.py"),
        ("radia-mcp", "packages/radia-mcp/pyproject.toml",
         "packages/radia-mcp/src/radia_mcp/__init__.py"),
        ("cubit-mesh-export", "packages/cubit-mesh-export/pyproject.toml",
         "packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py"),
    ]
    fails = []
    for name, pp, init in checks:
        v1, e1 = _ver_in(os.path.join(REPO, pp))
        v2, e2 = _ver_in(os.path.join(REPO, init))
        if e1 or e2:
            fails.append(f"{name}: {e1 or e2}")
        elif v1 != v2:
            fails.append(f"{name}: pyproject={v1} != __init__={v2}")
    return (not fails), "; ".join(fails) if fails else "radia/radia-mcp/cme versions in lockstep"


# ======================================================================
# Gate 3: TOOLS.md drift  (WIP-aware)
# ======================================================================
def gate_tools_md(fix=False):
    gen = os.path.join("packages", "radia-mcp", "scripts", "gen_tools_doc.py")
    doc = "packages/radia-mcp/docs/TOOLS.md"

    # WIP-awareness: regenerating with uncommitted radia_mcp/src changes
    # produces a TOOLS.md that reflects WIP tools.  If you then commit only
    # TOOLS.md (not the WIP code), CI's drift gate goes red against the
    # committed code.  Warn loudly.
    rc, out = _sh(["git", "status", "--porcelain", "--",
                   "packages/radia-mcp/src"])
    wip = [ln for ln in out.splitlines() if ln.strip()]
    warn = ""
    if wip:
        warn = (f" [WARN: {len(wip)} uncommitted radia_mcp/src file(s); "
                "commit code + regenerated TOOLS.md TOGETHER or CI will drift]")

    rc, out = _sh([sys.executable, gen])
    if rc != 0:
        return False, f"gen_tools_doc.py failed: {out[-200:]}"
    rc, _ = _sh(["git", "diff", "--exit-code", "--", doc])
    if rc == 0:
        return True, "docs/TOOLS.md matches generated inventory" + warn
    # drifted
    if fix:
        return True, "docs/TOOLS.md was STALE -- regenerated (stage it before commit)" + warn
    _sh(["git", "checkout", "--", doc])  # restore (don't leave dirty)
    return False, ("docs/TOOLS.md is STALE vs code -- run "
                   "`python tools/ci_preflight.py --fix` (or "
                   "`python packages/radia-mcp/scripts/gen_tools_doc.py`) "
                   "and commit it" + warn)


# ======================================================================
# Gate 4: radia-mcp matrix  (compile + health + minimal-dep pytest)
# ======================================================================
def gate_radia_mcp_matrix():
    # 4a compile
    rc, out = _sh([sys.executable, "-m", "compileall", "-q",
                   "packages/radia-mcp/src"])
    if rc != 0:
        return False, f"compileall failed: {out[-200:]}"

    # 4b meta health (all subpackages import + register_status_tool)
    probe = ("import sys; sys.path.insert(0, r'%s');"
             "from radia_mcp.meta.server import radia_mcp_health as h;"
             "r=h(); print('HEALTH', r['n_servers_healthy'], r['n_servers_total']);"
             "sys.exit(0 if r['all_healthy'] else 7)"
             % os.path.join(MCP, "src"))
    rc, out = _sh([sys.executable, "-c", probe])
    health = next((l for l in out.splitlines() if l.startswith("HEALTH")), "")
    if rc != 0:
        return False, f"meta_health FAILED ({health}): {out[-300:]}"

    # 4c the FULL radia-mcp pytest, UNDER minimal-dep simulation -- this is
    #     exactly the ubuntu matrix "Pytest" step (the 30x top failure).
    rc, out = _sh([sys.executable, "-m", "pytest", "tests/", "-q",
                   "-p", "no:cacheprovider", "--no-header"],
                  cwd=MCP, env={"RADIA_MCP_FORCE_MINIMAL": "1"})
    tail = out.strip().splitlines()[-1] if out.strip() else "(no output)"
    if rc != 0:
        # surface collection errors / failures
        errs = [l for l in out.splitlines()
                if ("ERROR" in l or "error" in l or "FAILED" in l
                    or "ModuleNotFoundError" in l)][:6]
        return False, ("radia-mcp pytest (minimal-dep sim) FAILED: " + tail
                       + ("\n      " + "\n      ".join(errs) if errs else ""))
    return True, f"compile + meta_health({health.replace('HEALTH ','')}) + minimal-dep pytest: {tail}"


# ======================================================================
# Gate 5: top-level collect-only  (CI "Run basic tests" import check)
# ======================================================================
def gate_toplevel_collect():
    cmd = [sys.executable, "-m", "pytest", "tests/", "--collect-only", "-q",
           "-p", "no:cacheprovider"]
    for ig in TOPLEVEL_IGNORES:
        cmd.append(f"--ignore={ig}")
    rc, out = _sh(cmd)
    if rc != 0:
        errs = [l for l in out.splitlines()
                if "error" in l.lower() or "Interrupted" in l][:6]
        return False, "top-level collection FAILED:\n      " + "\n      ".join(errs)
    last = next((l for l in reversed(out.splitlines()) if "collected" in l), out.strip()[-120:])
    return True, last.strip()


# ======================================================================
# Gate 6: top-level run  (--full only; slow FEM solves)
# ======================================================================
def gate_toplevel_run():
    cmd = [sys.executable, "-m", "pytest", "tests/", "-q", "--no-header",
           "-p", "no:cacheprovider", "--reruns", "1", "--reruns-delay", "1"]
    for ig in TOPLEVEL_IGNORES:
        cmd.append(f"--ignore={ig}")
    rc, out = _sh(cmd, timeout=3600)
    tail = out.strip().splitlines()[-1] if out.strip() else "(no output)"
    if rc != 0:
        fails = [l for l in out.splitlines() if "FAILED" in l][:10]
        return False, "top-level pytest FAILED: " + tail + (
            "\n      " + "\n      ".join(fails) if fails else "")
    return True, tail


ALL_GATES = [
    ("policy",          "Policy Lint (7 static policies)",        gate_policy_lint),
    ("version",         "Version consistency (pyproject==init)",  gate_version_consistency),
    ("tools-md",        "TOOLS.md drift gate",                    gate_tools_md),
    ("radia-mcp",       "radia-mcp matrix (minimal-dep pytest)",  gate_radia_mcp_matrix),
    ("toplevel-collect","Top-level collect-only (import check)",  gate_toplevel_collect),
]
FULL_GATES = [
    ("toplevel-run",    "Top-level pytest (full, slow)",          gate_toplevel_run),
]


def _changed_since(ref):
    """Files changed between <ref> and HEAD; None if range can't be resolved."""
    rc, out = _sh(["git", "diff", "--name-only", f"{ref}..HEAD"])
    if rc != 0:
        return None
    return [ln.strip().replace("\\", "/") for ln in out.splitlines() if ln.strip()]


def _gates_for_changes(changed):
    """Path-aware gate selection (used by the pre-push hook so most pushes
    stay fast): policy + version ALWAYS; the radia-mcp gates only when
    packages/radia-mcp changed; the top-level collect gate only when tests/
    or src/ changed."""
    sel = {"policy", "version"}
    if any(f.startswith("packages/radia-mcp/") for f in changed):
        sel |= {"tools-md", "radia-mcp"}
    if any(f.startswith("tests/") or f.startswith("src/") for f in changed):
        sel |= {"toplevel-collect"}
    return [g for g in ALL_GATES if g[0] in sel]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", action="store_true",
                    help="also run the full top-level pytest (slow FEM solves)")
    ap.add_argument("--fix", action="store_true",
                    help="auto-regenerate TOOLS.md if drifted (stage it yourself)")
    ap.add_argument("--only", default="",
                    help="comma-separated gate keys to run (default: all fast gates)")
    ap.add_argument("--since", default="",
                    help="select gates by what changed since REF (e.g. "
                         "origin/main); used by the pre-push hook to stay fast")
    args = ap.parse_args(argv)

    gates = list(ALL_GATES)
    if args.full:
        gates += FULL_GATES
    if args.since:
        changed = _changed_since(args.since)
        if changed is None:
            print(f"{YEL}--since {args.since}: range unresolved -> running ALL "
                  f"gates (safe default){RST}")
        elif not changed:
            print(f"{DIM}--since {args.since}: nothing changed -> no gates{RST}")
            gates = []
        else:
            gates = _gates_for_changes(changed)
            print(f"{DIM}--since {args.since}: {len(changed)} file(s) changed "
                  f"-> {len(gates)} gate(s){RST}")
    if args.only:
        keys = {k.strip() for k in args.only.split(",")}
        gates = [g for g in (ALL_GATES + FULL_GATES) if g[0] in keys]

    print(f"{DIM}repo: {REPO}{RST}")
    print(f"CI preflight: {len(gates)} gate(s)\n")
    results = []
    for key, label, fn in gates:
        sys.stdout.write(f"  {label:42} ... ")
        sys.stdout.flush()
        try:
            ok, detail = (fn(fix=args.fix) if key == "tools-md" else fn())
        except Exception as e:  # noqa: BLE001
            ok, detail = False, f"gate crashed: {type(e).__name__}: {e}"
        results.append((key, label, ok, detail))
        mark = f"{GREEN}PASS{RST}" if ok else f"{RED}FAIL{RST}"
        print(mark)
        print(f"      {DIM if ok else ''}{detail}{RST}")

    nfail = sum(1 for *_, ok, _ in results if not ok)
    print()
    if nfail == 0:
        print(f"{GREEN}[OK] all {len(results)} gate(s) green -- safe to push.{RST}")
        return 0
    print(f"{RED}[FAIL] {nfail}/{len(results)} gate(s) red -- DO NOT push; "
          f"fix locally first (this is what 'commit -> CI check -> push' prevents).{RST}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
