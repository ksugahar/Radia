#!/usr/bin/env python
"""Run impact-scoped CI gates locally before pushing.

The default compares HEAD and the working tree with ``origin/main``. It runs
only the compact gates owned by affected package paths. Full tests and solver
validation remain explicit operations rather than automatic push costs.

Gates:
  1. policy-lint        -- the 8 static policies (CblasColMajor etc.)
  2. publish-boundary   -- radia-mcp provenance/internal-path lint + selftest
  3. version-consistency-- pyproject == __init__ for radia / radia-mcp / cme
  4. tools-md-drift     -- regenerate docs/TOOLS.md and diff (WIP-aware)
  5. radia-mcp          -- compile + health + impact-selected package tests
                           and affected server selftests under the same
                           minimal-dependency simulation as GitHub CI
  6. toplevel-collect   -- explicit `pytest tests/ --collect-only` diagnostic
  7. toplevel-run       -- (only with --full) run the lightweight tests/
  8. validation-collect -- (only with --validation) collect the heavy
                           validation_test/ suite
  9. validation-run     -- (only with --validation --full) run the heavy
                           validation_test/ suite

Usage:
    python tools/ci_preflight.py            # affected compact gates
    python tools/ci_preflight.py --full     # also run lightweight tests/
    python tools/ci_preflight.py --validation  # also collect validation_test/
    python tools/ci_preflight.py --fix      # auto-regenerate TOOLS.md on drift
    python tools/ci_preflight.py --only policy,tools-md

Exit 0 = all green, safe to push.  Non-zero = a gate CI would fail is red.
"""
from __future__ import annotations

import argparse
import json
import os
import runpy
import subprocess
import sys

# --- UNC-safe repo root -------------------------------------------------
# NEVER use Path.resolve() here: on the LAB box the repo lives on an S:
# drive mapped to the Radia NAS share, and resolve() canonicalises to the
# UNC form, which Python's open() then rejects with OSError [Errno 22]
# (this is exactly what broke the old release preflight).  os.path
# joins keep the drive letter the script was invoked with.
_THIS = os.path.abspath(__file__)
REPO = os.path.dirname(os.path.dirname(_THIS))
def _remap_lab_unc(path: str) -> str:
    """Map the LAB Radia UNC tree to its stable ``S:`` drive spelling."""
    normalized = path.replace("\\", "/")
    lowered = normalized.lower()
    if not any(host in lowered for host in ("192.168.11.100", "192.168.121.100")):
        return path
    anchor_index = lowered.find("/radia/")
    if anchor_index < 0 and lowered.endswith("/radia"):
        anchor_index = len(normalized) - len("/Radia")
    if anchor_index < 0:
        return path
    return "S:" + normalized[anchor_index:].replace("/", "\\")


# A pre-push hook may report the mapped share as UNC. Keep both the historical
# and current NAS addresses readable, including isolated worktrees.
REPO = _remap_lab_unc(REPO)
MCP = os.path.join(REPO, "packages", "radia-mcp")

# The lightweight tests/ gate should not require directory or file ignores.
# Validation-class tests belong in validation_test/.
TOPLEVEL_IGNORES = [
]

GREEN, RED, YEL, DIM, RST = "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"

_PREFLIGHT_CHANGED_FILES = None


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
    # never drift (previously the policies were re-implemented inline here).
    rc, out = _sh([sys.executable, os.path.join(REPO, "tools", "policy_lint.py"),
                   "--quiet"])
    if rc == 0:
        return True, "8 policies pass"
    fails = [ln[6:].strip() for ln in out.splitlines() if ln.startswith("FAIL")]
    return False, "; ".join(fails) if fails else (out.strip()[-200:] or "policy lint failed")


# ======================================================================
# Gate 2: radia-mcp publish boundary (same commands as GitHub CI)
# ======================================================================
def gate_publish_boundary_lint():
    script = os.path.join(MCP, "tools", "policy_lint.py")
    rc, out = _sh([sys.executable, script, "--selftest"], timeout=120)
    if rc != 0:
        tail = out.strip().splitlines()[-8:]
        return False, "publish-boundary selftest failed: " + " | ".join(tail)

    rc, out = _sh([sys.executable, script], timeout=180)
    if rc == 0:
        return True, "radia-mcp provenance and internal-path scan passes"
    findings = [
        line.strip()
        for line in out.splitlines()
        if "finding(s)" in line or "lab-local absolute" in line
    ]
    detail = " | ".join(findings[-6:]) or out.strip()[-500:]
    return False, f"radia-mcp publish-boundary lint failed: {detail}"


# ======================================================================
# Gate 3: version consistency  (pyproject == __init__)
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
# Gate 4: TOOLS.md drift  (WIP-aware)
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
# Gate 5: radia-mcp impact lane
# ======================================================================
def gate_radia_mcp_matrix():
    selector = os.path.join(MCP, "tools", "select_ci_tests.py")
    selector_cmd = [sys.executable, selector]
    if _PREFLIGHT_CHANGED_FILES is None:
        selector_cmd.append("--full")
    else:
        selector_cmd.extend(
            [
                "--changed-files-json",
                json.dumps(_PREFLIGHT_CHANGED_FILES),
                "--base",
                "origin/main",
            ]
        )
    rc, out = _sh(selector_cmd)
    if rc != 0:
        return False, f"impact selector failed: {out[-400:]}"
    try:
        plan = json.loads(out)
    except json.JSONDecodeError as exc:
        return False, f"impact selector emitted invalid JSON: {exc}"

    rc, out = _sh([sys.executable, "-m", "compileall", "-q",
                   "packages/radia-mcp/src"])
    if rc != 0:
        return False, f"compileall failed: {out[-200:]}"

    # The compact health probe catches shared import/registration regressions.
    probe = ("import sys; sys.path.insert(0, r'%s');"
             "from radia_mcp.meta.server import radia_mcp_health as h;"
             "r=h(); print('HEALTH', r['n_servers_healthy'], r['n_servers_total']);"
             "sys.exit(0 if r['all_healthy'] else 7)"
             % os.path.join(MCP, "src"))
    rc, out = _sh([sys.executable, "-c", probe])
    health = next((l for l in out.splitlines() if l.startswith("HEALTH")), "")
    if rc != 0:
        return False, f"meta_health FAILED ({health}): {out[-300:]}"

    selected_tests = plan["package_tests"]
    test_env = {
        "RADIA_MCP_FORCE_MINIMAL": "1",
        "PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1",
        # Tests that spawn a child Python must resolve the same checkout as
        # the parent pytest process, even in an isolated preflight venv.
        "PYTHONPATH": os.pathsep.join(
            part for part in (
                os.path.join(MCP, "src"),
                os.environ.get("PYTHONPATH", ""),
            ) if part
        ),
    }
    if selected_tests != ["tests"]:
        test_env["RADIA_MCP_CI_SELECTION_JSON"] = json.dumps(selected_tests)
    print(
        "  [radia-mcp] running "
        f"{len(selected_tests)} impact selector(s) (timeout: 5 min)...",
        flush=True,
    )
    try:
        # Discover from the package test root even for impact-scoped runs.
        # conftest.py ignores unselected and unavailable-optional-dependency
        # files before module import, then deselects unrequested node IDs.
        # Passing files explicitly bypasses pytest's collect_ignore contract.
        pytest_targets = ["tests/"]
        rc, out = _sh([sys.executable, "-m", "pytest", *pytest_targets, "-q",
                       "-p", "no:cacheprovider", "--no-header",
                       "-m", "not xval and not slow"],
                      cwd=MCP, env=test_env, timeout=300)
    except subprocess.TimeoutExpired as exc:
        partial = (exc.stdout or exc.stderr or "")
        return False, ("radia-mcp pytest timed out after 300 seconds; "
                       "the selected lane is too slow or hung" +
                       (f" (last output: {partial[-300:]})" if partial else ""))
    tail = out.strip().splitlines()[-1] if out.strip() else "(no output)"
    if rc != 0:
        # surface collection errors / failures
        errs = [l for l in out.splitlines()
                if ("ERROR" in l or "error" in l or "FAILED" in l
                    or "ModuleNotFoundError" in l)][:6]
        return False, ("radia-mcp impact pytest FAILED: " + tail
                       + ("\n      " + "\n      ".join(errs) if errs else ""))

    catalog = runpy.run_path(
        os.path.join(MCP, "src", "radia_mcp", "meta", "catalog.py")
    )["CATALOG"]
    for short in plan["server_selftests"]:
        info = catalog[short]
        command = (
            "import sys; sys.path.insert(0, r'%s'); "
            "from %s.server import main; "
            "sys.argv = [%r, '--selftest']; main()"
            % (os.path.join(MCP, "src"), info["subpackage"], info["entry_point"])
        )
        try:
            rc, server_out = _sh(
                [sys.executable, "-c", command],
                timeout=120 if short == "radia-ngsolve" else 60,
            )
        except subprocess.TimeoutExpired:
            return False, f"{short} selftest timed out"
        if rc != 0:
            return False, f"{short} selftest failed: {server_out[-400:]}"

    return True, (
        f"compile + health({health.replace('HEALTH ', '')}) + {tail}; "
        f"{len(plan['server_selftests'])} server selftest(s)"
    )


# ======================================================================
# Gate 6: top-level collect-only  (CI "Run basic tests" import check)
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
# Gate 7: top-level run  (--full only; lightweight tests/)
# ======================================================================
def gate_toplevel_run():
    cmd = [sys.executable, "-m", "pytest", "tests/", "-q", "--no-header",
           "-p", "no:cacheprovider"]
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
# Gate 8: validation collect  (manual heavy validation_test/ import check)
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
# Gate 9: validation run  (manual heavy validation_test/ run)
# ======================================================================
def gate_validation_run():
    cmd = [sys.executable, "-m", "pytest", "validation_test/", "-q",
           "--no-header", "-p", "no:cacheprovider",
           "-m", "not compute_host"]
    rc, out = _sh(cmd, timeout=7200)
    tail = out.strip().splitlines()[-1] if out.strip() else "(no output)"
    if rc != 0:
        fails = [l for l in out.splitlines() if "FAILED" in l][:12]
        return False, "validation pytest FAILED: " + tail + (
            "\n      " + "\n      ".join(fails) if fails else "")
    return True, tail


ALL_GATES = [
    ("policy",          "Policy Lint (8 static policies)",        gate_policy_lint),
    ("publish-boundary","radia-mcp publish-boundary lint",        gate_publish_boundary_lint),
    ("version",         "Version consistency (pyproject==init)",  gate_version_consistency),
    ("tools-md",        "TOOLS.md drift gate",                    gate_tools_md),
    ("radia-mcp",       "radia-mcp impact lane",                  gate_radia_mcp_matrix),
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
    """Committed and working-tree files changed from ``ref``."""
    def git_names(command):
        completed = subprocess.run(
            command,
            cwd=REPO,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if completed.returncode:
            return None
        # Git may emit core.autocrlf warnings on stderr. They are diagnostics,
        # never candidate paths for impact selection.
        return completed.stdout

    out = git_names(["git", "diff", "--name-only", f"{ref}...HEAD"])
    if out is None:
        return None
    changed = {ln.strip().replace("\\", "/") for ln in out.splitlines() if ln.strip()}
    for command in (
        ["git", "diff", "--name-only", "HEAD"],
        ["git", "diff", "--name-only", "--cached"],
        ["git", "ls-files", "--others", "--exclude-standard"],
    ):
        out = git_names(command)
        if out is None:
            return None
        changed.update(
            ln.strip().replace("\\", "/") for ln in out.splitlines() if ln.strip()
        )
    return sorted(changed)


def _gates_for_changes(changed):
    """Path-aware gate selection (used by the pre-push hook so most pushes
    stay fast): policy + version ALWAYS; the radia-mcp gates only when
    packages/radia-mcp changed. Repository contracts are owned by the fixed
    fast mdx lane, so broad top-level collection is an explicit diagnostic."""
    sel = {"policy", "version"}
    mcp_changes = [
        f for f in changed if f.startswith("packages/radia-mcp/")
    ]
    mcp_docs_only = bool(mcp_changes) and all(
        f.startswith("packages/radia-mcp/docs/")
        or f in {
            "packages/radia-mcp/README.md",
            "packages/radia-mcp/CHANGELOG.md",
        }
        for f in mcp_changes
    )
    if mcp_changes and not mcp_docs_only:
        sel |= {"publish-boundary", "radia-mcp"}
    if any(
        f in {
            "packages/radia-mcp/docs/TOOLS.md",
            "packages/radia-mcp/scripts/gen_tools_doc.py",
        }
        for f in changed
    ):
        sel.add("tools-md")
    return [g for g in ALL_GATES if g[0] in sel]


def main(argv=None):
    global _PREFLIGHT_CHANGED_FILES
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
    scope_ref = args.since or "origin/main"
    changed = _changed_since(scope_ref)
    if changed is None:
        print(f"{YEL}--since {scope_ref}: range unresolved -> running ALL "
              f"gates (safe default){RST}")
    elif not changed:
        print(f"{DIM}--since {scope_ref}: nothing changed -> no gates{RST}")
        gates = []
        _PREFLIGHT_CHANGED_FILES = []
    else:
        gates = _gates_for_changes(changed)
        _PREFLIGHT_CHANGED_FILES = changed
        print(f"{DIM}--since {scope_ref}: {len(changed)} file(s) changed "
              f"-> {len(gates)} gate(s){RST}")
    if args.full:
        gates += FULL_GATES
    if args.validation:
        gates += VALIDATION_GATES
        if args.full:
            gates += VALIDATION_FULL_GATES
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
