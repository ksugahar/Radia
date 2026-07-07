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
  4. radia-mcp-matrix   -- compile + meta_health + the FAST radia-mcp
                           pytest run UNDER minimal-dep simulation
                           (RADIA_MCP_FORCE_MINIMAL=1), which is what the
                           ubuntu matrix actually runs -- this is the gate
                           that catches heavy-import collection errors
                           (ngsolve/netgen imported at module level)
  5. toplevel-collect   -- `pytest tests/ --collect-only` for the lightweight
                           debug/CI suite
  6. toplevel-run       -- (only with --full) actually run the lightweight
                           tests/ suite with reruns
  7. validation-collect -- (only with --validation) collect the heavy
                           validation_test/ suite
  8. validation-run     -- (only with --validation --full) run the heavy
                           validation_test/ suite

Usage:
    python tools/ci_preflight.py            # gates 1-5 (~2 min)
    python tools/ci_preflight.py --full     # also run lightweight tests/
    python tools/ci_preflight.py --validation  # also collect validation_test/
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
# (this is exactly what broke the old release preflight).  os.path
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

# The lightweight tests/ gate should not require directory or file ignores.
# Validation-class tests belong in validation_test/.
TOPLEVEL_IGNORES = [
]

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


# ======================================================================
# Gate 1: policy lint  (delegates to tools/policy_lint.py)
# ======================================================================
def gate_policy_lint():
    # SINGLE SOURCE OF TRUTH: tools/policy_lint.py is also what
    # .github/workflows/policy-lint.yml runs, so the local gate and CI can
    # never drift (previously the 7 policies were re-implemented inline here).
    rc, out = _sh([sys.executable, os.path.join(REPO, "tools", "policy_lint.py"),
                   "--quiet"])
    if rc == 0:
        return True, "7 policies pass"
    fails = [ln[6:].strip() for ln in out.splitlines() if ln.startswith("FAIL")]
    return False, "; ".join(fails) if fails else (out.strip()[-200:] or "policy lint failed")


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

    rc, out = _sh(["git", "status", "--porcelain", "--", "packages/radia-mcp/src"])
    wip = [ln for ln in out.splitlines() if ln.strip()]

    if not wip:
        # Fast path: working tree == committed; regenerate in place + diff.
        rc, out = _sh([sys.executable, os.path.join(REPO, gen)])
        if rc != 0:
            return False, f"gen_tools_doc.py failed: {out[-200:]}"
        rc, _ = _sh(["git", "diff", "--exit-code", "--", doc])
        if rc == 0:
            return True, "docs/TOOLS.md matches generated inventory"
        if fix:
            return True, "docs/TOOLS.md was STALE -- regenerated (stage it before commit)"
        _sh(["git", "checkout", "--", doc])  # restore (don't leave dirty)
        return False, ("docs/TOOLS.md is STALE vs code -- run "
                       "`python tools/ci_preflight.py --fix` and commit it")

    # WIP present: CI checks out the COMMITTED code, so check the COMMITTED
    # state (HEAD), NOT the working tree -- otherwise uncommitted radia_mcp/src
    # changes (this or another concurrent session) false-positive this gate.
    # Extract HEAD's radia-mcp to a temp dir, regenerate from that committed
    # code, and compare to HEAD's committed TOOLS.md (exactly what CI does).
    import io
    import shutil
    import tarfile
    import tempfile
    tmp = tempfile.mkdtemp(prefix="radia_toolsmd_")
    try:
        ar = subprocess.run(["git", "archive", "HEAD", "packages/radia-mcp"],
                            cwd=REPO, capture_output=True)
        if ar.returncode != 0:
            return False, "git archive HEAD failed: " + ar.stderr[-200:].decode(errors="replace")
        with tarfile.open(fileobj=io.BytesIO(ar.stdout)) as tf:
            tf.extractall(tmp, filter="data")  # filter: Py3.12+ safe-extract
        gen_tmp = os.path.join(tmp, "packages", "radia-mcp", "scripts", "gen_tools_doc.py")
        rc, out = _sh([sys.executable, gen_tmp])  # regen from COMMITTED code
        if rc != 0:
            return False, f"gen_tools_doc.py (committed) failed: {out[-200:]}"
        with open(os.path.join(tmp, "packages", "radia-mcp", "docs", "TOOLS.md"),
                  encoding="utf-8") as f:
            regen = f.read()
        show = subprocess.run(["git", "show", f"HEAD:{doc}"], cwd=REPO,
                              capture_output=True, text=True,
                              encoding="utf-8", errors="replace")
        committed = show.stdout
        note = f" (checked HEAD; {len(wip)} working-tree WIP file(s) ignored)"
        if regen.replace("\r\n", "\n") == committed.replace("\r\n", "\n"):
            return True, "committed docs/TOOLS.md matches committed code" + note
        if fix:
            with open(os.path.join(REPO, doc), "w", encoding="utf-8", newline="") as f:
                f.write(regen)
            return True, "committed docs/TOOLS.md was STALE -- regenerated from HEAD (stage it)" + note
        return False, ("committed docs/TOOLS.md is STALE vs COMMITTED code (HEAD): a "
                       "tool was added/removed without regenerating -- run "
                       "`python tools/ci_preflight.py --fix` and commit it" + note)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


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

    # 4c the FAST radia-mcp pytest, UNDER minimal-dep simulation -- this is
    #     exactly the ubuntu matrix "Pytest" step. Slow corpus / subprocess
    #     selftest gates are explicit validation or dedicated CI steps.
    rc, out = _sh([sys.executable, "-m", "pytest", "tests/", "-q",
                   "-p", "no:cacheprovider", "--no-header",
                   "-m", "not xval and not slow",
                   "--reruns", "1", "--reruns-delay", "1"],  # match matrix flaky-retry
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
# Gate 6: top-level run  (--full only; lightweight tests/)
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


# ======================================================================
# Gate 7: validation collect  (manual heavy validation_test/ import check)
# ======================================================================
def gate_validation_collect():
    cmd = [sys.executable, "-m", "pytest", "validation_test/",
           "--collect-only", "-q", "-p", "no:cacheprovider"]
    rc, out = _sh(cmd)
    if rc != 0:
        errs = [l for l in out.splitlines()
                if "error" in l.lower() or "Interrupted" in l][:8]
        return False, "validation collection FAILED:\n      " + "\n      ".join(errs)
    last = next((l for l in reversed(out.splitlines()) if "collected" in l), out.strip()[-120:])
    return True, last.strip()


# ======================================================================
# Gate 8: validation run  (manual heavy validation_test/ run)
# ======================================================================
def gate_validation_run():
    cmd = [sys.executable, "-m", "pytest", "validation_test/", "-q",
           "--no-header", "-p", "no:cacheprovider",
           "--reruns", "1", "--reruns-delay", "1"]
    rc, out = _sh(cmd, timeout=7200)
    tail = out.strip().splitlines()[-1] if out.strip() else "(no output)"
    if rc != 0:
        fails = [l for l in out.splitlines() if "FAILED" in l][:12]
        return False, "validation pytest FAILED: " + tail + (
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
    ("toplevel-run",    "Top-level pytest (full lightweight)",    gate_toplevel_run),
]
VALIDATION_GATES = [
    ("validation-collect", "Validation collect-only (heavy suite)", gate_validation_collect),
]
VALIDATION_FULL_GATES = [
    ("validation-run",     "Validation pytest (full, very slow)",   gate_validation_run),
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
                    help="also run the lightweight tests/ pytest")
    ap.add_argument("--validation", action="store_true",
                    help="also collect validation_test/; with --full, run it")
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
    if args.validation:
        gates += VALIDATION_GATES
        if args.full:
            gates += VALIDATION_FULL_GATES
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
        gates = [g for g in (ALL_GATES + FULL_GATES + VALIDATION_GATES + VALIDATION_FULL_GATES)
                 if g[0] in keys]

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
