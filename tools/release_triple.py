#!/usr/bin/env python
"""release_triple.py — orchestrator for the 3-package release flow.

Walks Phase 0 -> 9 of the release-triple skill in order, gating each
phase on the success of the previous one. Refuses to skip steps that
have caused real outages (2026-04-14 incident series).

Usage:
    python tools/release_triple.py preflight
        Read-only: report current state and consistency. Use anytime.

    python tools/release_triple.py phase0
        Mandatory clean rebuild of the Cubit plugin (~3-4 min).

    python tools/release_triple.py phase8 [--target lab|100|all]
        Run Phase 8a..8d on each target: kill Cubit, install from NAS,
        cubit-plugin-install, --verify-only, cubit-smoke-test. Refuses
        to start if Phase 0 has not been done since the last source
        change in src/cubit_plugin/.

    python tools/release_triple.py phase8e
        Upgrade mdx from PyPI. Refuses to run if pip index versions
        radia / cubit-mesh-export / radia-mcp don't match the local
        repo (i.e. PyPI hasn't propagated yet).

    python tools/release_triple.py phase9
        Cross-machine consistency probe. Final gate.

    python tools/release_triple.py all
        phase8 -> phase8e -> phase9 with all preconditions enforced.

Exit codes:
    0   success
    2   precondition failure (skip detected)
    3   action failure (external command)
    4   verification mismatch
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path

# Force UTF-8 stdout/stderr so en/em dashes and CJK in messages do not
# crash the script on ja-JP cp932 consoles.
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

REPO = Path(__file__).resolve().parent.parent
NAS_REPO_LAB = "S:/Radia/01_GitHub"
NAS_REPO_100 = "//192.168.11.100/work/00_CAE/Radia/01_GitHub"
SSH_100 = "192.168.11.100"
SSH_MDX = "mdx"


# ============================================================
# tiny io helpers
# ============================================================

def _color(c, s):
    return f"\033[{c}m{s}\033[0m"


def info(msg):  print("  " + msg)
def ok(msg):    print("  " + _color("32;1", "[OK]    ") + msg)
def warn(msg):  print("  " + _color("33;1", "[WARN]  ") + msg)
def fail(msg):  print("  " + _color("31;1", "[FAIL]  ") + msg)
def step(msg):  print("\n" + _color("36;1", f"=== {msg} ==="))


def run(cmd, *, check=True, capture=False, shell=False, **kw):
    """Subprocess wrapper that prints what it runs."""
    print("  $ " + (cmd if shell else " ".join(str(c) for c in cmd)))
    p = subprocess.run(cmd, shell=shell, capture_output=capture,
                        text=True, **kw)
    if check and p.returncode != 0:
        fail(f"command failed (exit {p.returncode})")
        if capture and p.stderr:
            print(p.stderr.strip())
        sys.exit(3)
    return p


# ============================================================
# state inspectors
# ============================================================

def _read_repo_versions():
    """Parse versions out of pyproject.toml + __init__.py (no toml dep)."""
    out = {}
    import re
    for label, path in [
        ("radia",             REPO / "pyproject.toml"),
        ("radia.__version__", REPO / "src/radia/__init__.py"),
        ("cubit-mesh-export", REPO / "packages/cubit-mesh-export/pyproject.toml"),
        ("cme.__version__",   REPO / "packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py"),
        ("radia-mcp",         REPO / "packages/radia-mcp/pyproject.toml"),
    ]:
        text = path.read_text(encoding="utf-8")
        m = re.search(r'(?:^version|^__version__)\s*=\s*"([^"]+)"', text, re.M)
        out[label] = m.group(1) if m else None
    return out


def _newest_mtime(root: Path, suffixes):
    latest = 0.0
    for r, dirs, files in os.walk(root):
        # skip build dirs
        dirs[:] = [d for d in dirs
                    if d not in ("build-pyd", "build-ccm", "build", "compact_netgen")]
        for f in files:
            if Path(f).suffix.lower() in suffixes:
                p = Path(r) / f
                try:
                    mt = p.stat().st_mtime
                    if mt > latest:
                        latest = mt
                except OSError:
                    pass
    return latest


def _bundled_plugin_mtime():
    """Newest mtime of bundled .ccm/.ccl in cubit-mesh-export package."""
    pkg = REPO / "packages/cubit-mesh-export/src/cubit_mesh_export"
    times = []
    for name in ("radia_cubit.ccm", "radia_cubit.ccl"):
        p = pkg / name
        if p.is_file():
            times.append(p.stat().st_mtime)
    return max(times) if times else 0.0


# ============================================================
# Phases
# ============================================================

def cmd_preflight(args):
    """Read-only state report. Always safe to run."""
    step("Phase preflight: state report")

    # Versions in repo
    v = _read_repo_versions()
    info(f"Repo versions:")
    info(f"  radia              pyproject={v['radia']}  __version__={v['radia.__version__']}")
    info(f"  cubit-mesh-export  pyproject={v['cubit-mesh-export']}  __version__={v['cme.__version__']}")
    info(f"  radia-mcp          pyproject={v['radia-mcp']}")

    pp_radia = (v["radia"] == v["radia.__version__"])
    pp_cme   = (v["cubit-mesh-export"] == v["cme.__version__"])
    if pp_radia: ok("radia pyproject == __init__")
    else:        fail("radia pyproject != __init__ — fix before any release")
    if pp_cme: ok("cubit-mesh-export pyproject == __init__")
    else:      fail("cubit-mesh-export pyproject != __init__ — fix before any release")

    # Cubit plugin freshness
    src_dir = REPO / "src/cubit_plugin"
    src_mtime = _newest_mtime(src_dir, {".cpp", ".cc", ".cxx", ".c", ".h",
                                          ".hpp", ".hh", ".hxx", ".cmake", ".txt"})
    bin_mtime = _bundled_plugin_mtime()
    if src_mtime == 0:
        warn("could not measure src/cubit_plugin/ mtime")
    elif bin_mtime == 0:
        fail("bundled .ccm/.ccl missing — Phase 0 not done")
    elif bin_mtime + 1 < src_mtime:
        from datetime import datetime
        fail(f"bundled .ccm/.ccl ({datetime.fromtimestamp(bin_mtime)}) older than "
              f"src/cubit_plugin/ ({datetime.fromtimestamp(src_mtime)}). "
              "Run `python tools/release_triple.py phase0`.")
        return 2
    else:
        ok("bundled plugin .ccm/.ccl >= src/cubit_plugin/ mtime")

    return 0


def cmd_phase0(args):
    """Clean rebuild of Cubit plugin (.ccm + .ccl)."""
    step("Phase 0: clean rebuild of Cubit plugin (~3-4 min)")
    build_pyd = REPO / "src/cubit_plugin/build-pyd"
    build_ccm = REPO / "src/cubit_plugin/build-ccm"
    if build_pyd.exists():
        run(["rm", "-rf", str(build_pyd)])
    if build_ccm.exists():
        run(["rm", "-rf", str(build_ccm)])

    # Build via the same ps1 we used in the 2026-04-14 manual run.
    ps1 = REPO / "tools/_build_cubit_plugin.ps1"
    if not ps1.is_file():
        fail(f"missing helper script: {ps1}")
        return 3
    run(["pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(ps1)])

    # Propagate to both source-of-truth dirs.
    for src_name, dst_dirs in [
        ("build-ccm/radia_cubit.ccm", ["src/radia", "packages/cubit-mesh-export/src/cubit_mesh_export"]),
        ("build-ccm/radia_cubit.ccl", ["src/radia", "packages/cubit-mesh-export/src/cubit_mesh_export"]),
    ]:
        src = REPO / "src/cubit_plugin" / src_name
        if not src.is_file():
            fail(f"build did not produce {src}")
            return 3
        for d in dst_dirs:
            dst = REPO / d / src.name
            run(["cp", str(src), str(dst)])

    ok("Phase 0 complete; .ccm + .ccl propagated to src/radia + cubit-mesh-export pkg")
    return 0


def _kill_cubit_local():
    info("force-kill any local Cubit process")
    run(["pwsh", "-NoProfile", "-Command",
         "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
         "$_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit' "
         "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
         "Start-Sleep -Seconds 2"], check=False)


def _kill_mcp_local():
    info("force-kill any local mcp-server-*.exe (otherwise radia-mcp install fails)")
    run(["pwsh", "-NoProfile", "-Command",
         "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
         "$_.Name -like 'mcp-server*' "
         "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
         "Start-Sleep -Seconds 2"], check=False)


def _deploy_lab():
    step("Phase 8 (LAB): kill, install from NAS, plugin install, verify, smoke")
    _kill_cubit_local()
    _kill_mcp_local()
    repo = NAS_REPO_LAB
    for sub in ("", "/packages/cubit-mesh-export", "/packages/radia-mcp"):
        run(["pip", "install", "--force-reinstall", "--no-deps",
             "--no-cache-dir", repo + sub])
    run(["cubit-plugin-install"])
    run(["cubit-plugin-install", "--verify-only"])
    run(["cubit-smoke-test"])
    ok("Phase 8 complete on LAB")


def _deploy_100():
    step("Phase 8 (100号機): kill + install + plugin install + verify + smoke (over SSH)")
    repo = NAS_REPO_100
    ps_block = (
        "$ErrorActionPreference = 'Continue'; "
        "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
        "$_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit' "
        "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
        "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
        "$_.Name -like 'mcp-server*' "
        "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
        "Start-Sleep -Seconds 2; "
        f"pip install --force-reinstall --no-deps --no-cache-dir '{repo}'; "
        f"pip install --force-reinstall --no-deps --no-cache-dir '{repo}/packages/cubit-mesh-export'; "
        f"pip install --force-reinstall --no-deps --no-cache-dir '{repo}/packages/radia-mcp'; "
        "cubit-plugin-install --all-users; "
        "cubit-plugin-install --verify-only; "
        "cubit-smoke-test"
    )
    run(["ssh", SSH_100, "pwsh", "-ExecutionPolicy", "Bypass",
         "-Command", ps_block])
    ok("Phase 8 complete on 100号機")


def cmd_phase8(args):
    """Deploy + verify + smoke on LAB and/or 100号機."""
    # precondition: Phase 0 freshness
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed; refusing Phase 8")
        return rc

    targets = args.target.split(",") if args.target else ["lab", "100"]
    for t in targets:
        t = t.strip().lower()
        if t == "lab":  _deploy_lab()
        elif t in ("100", "100号機", "100goki"): _deploy_100()
        elif t == "all": _deploy_lab(); _deploy_100()
        else:
            fail(f"unknown target: {t!r}")
            return 2
    return 0


def cmd_phase8e(args):
    """Upgrade mdx from PyPI (only after PyPI propagation)."""
    step("Phase 8e: mdx PyPI install")
    # check repo vs PyPI versions match — refuse if PyPI hasn't caught up
    v = _read_repo_versions()
    info("checking PyPI propagation...")
    for pkg, want in [("radia", v["radia"]),
                       ("cubit-mesh-export", v["cubit-mesh-export"]),
                       ("radia-mcp", v["radia-mcp"])]:
        p = run(["python", "-m", "pip", "index", "versions", pkg],
                capture=True, check=False)
        first = p.stdout.splitlines()[0] if p.stdout else ""
        if want and want in first:
            ok(f"PyPI {pkg} live at {want}")
        else:
            fail(f"PyPI {pkg} not yet at {want} (got {first!r}). "
                 "Wait for CI / PyPI propagation, then retry.")
            return 2

    run(["ssh", SSH_MDX, "pwsh", "-Command",
         "pip install --upgrade --no-deps radia cubit-mesh-export radia-mcp"])
    ok("Phase 8e complete on mdx")
    return 0


CROSS_MACHINE_PROBE = '''import hashlib, os
import importlib.metadata as md

def hsh_text(p):
    h = hashlib.sha256()
    try:
        d = open(p, "rb").read().replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
        h.update(d); return h.hexdigest()[:12]
    except Exception: return "MISSING"

def hsh_bin(p):
    h = hashlib.sha256()
    try:
        with open(p, "rb") as f:
            for c in iter(lambda: f.read(65536), b""): h.update(c)
        return h.hexdigest()[:12]
    except Exception: return "MISSING"

def ver(n):
    try: return md.version(n)
    except Exception: return "MISSING"

import radia, cubit_mesh_export
rad = os.path.dirname(radia.__file__)
cme = os.path.dirname(cubit_mesh_export.__file__)
print(f"VER radia              = {radia.__version__}")
print(f"VER cubit-mesh-export  = {cubit_mesh_export.__version__}")
print(f"VER radia-mcp          = {ver('radia-mcp')}")
print(f"COMPAT cme  -> radia   = [{cubit_mesh_export.COMPAT_RADIA_MIN}, {cubit_mesh_export.COMPAT_RADIA_MAX}]")
print(f"COMPAT rad  -> cme     = [{radia.COMPAT_CUBIT_MESH_EXPORT_MIN}, {radia.COMPAT_CUBIT_MESH_EXPORT_MAX}]")
for r in ["panels/register_toolbar.py",
          "panels/calc_peec_bem.py",
          "panels/calc_peec_inductance.py",
          "panels/calc_fem_kelvin.py",
          "panels/calc_fem_coilmesh.py"]:
    print(f"SHA radia/{r:35s} = {hsh_text(os.path.join(rad,r))}")
for r in ["radia_cubit.ccm","radia_cubit.ccl"]:
    print(f"SHA cme/{r:35s} = {hsh_bin(os.path.join(cme,r))}")
'''


def _probe(host_label, cmd_prefix):
    """Run the probe on a target (cmd_prefix is the python invocation)."""
    p = subprocess.run(cmd_prefix, input=CROSS_MACHINE_PROBE,
                        capture_output=True, text=True, shell=False)
    if p.returncode != 0:
        fail(f"probe failed on {host_label}: {p.stderr.strip()}")
        return None
    return p.stdout


def cmd_phase9(args):
    """Cross-machine consistency probe."""
    step("Phase 9: cross-machine consistency (LAB / 100号機 / mdx)")
    out_lab = _probe("LAB", ["python", "-"])
    out_100 = _probe("100号機", ["ssh", SSH_100, "python", "-"])
    out_mdx = _probe("mdx", ["ssh", SSH_MDX, "python", "-"])
    if not (out_lab and out_100 and out_mdx):
        fail("could not collect probe data from all 3 machines")
        return 4

    rows_lab = out_lab.strip().splitlines()
    rows_100 = out_100.strip().splitlines()
    rows_mdx = out_mdx.strip().splitlines()
    n = min(len(rows_lab), len(rows_100), len(rows_mdx))
    drift = 0

    print(f"\n  {'field':<44} | {'LAB':<18} | {'100号機':<18} | {'mdx':<18}")
    print("  " + "-" * 110)
    for i in range(n):
        a, b, c = rows_lab[i], rows_100[i], rows_mdx[i]

        def split(s):
            k, _, v = s.partition("=")
            k = k.strip()
            # strip leading "VER " / "SHA " / "COMPAT " from key
            for prefix in ("VER ", "SHA ", "COMPAT "):
                if k.startswith(prefix):
                    k = k[len(prefix):]
            return k.strip(), v.strip()

        ka, va = split(a)
        kb, vb = split(b)
        kc, vc = split(c)
        match = (va == vb == vc)
        marker = ok.__name__.upper() if match else "DRIFT"
        marker = _color("32;1", "OK") if match else _color("31;1", "DRIFT")
        if not match:
            drift += 1
        print(f"  {ka:<44} | {va:<18} | {vb:<18} | {vc:<18}  [{marker}]")

    if drift:
        print("")
        fail(f"{drift} field(s) drift across machines — release NOT done.")
        return 4
    print("")
    ok("all fields match across LAB / 100号機 / mdx — release verified.")
    return 0


def cmd_all(args):
    """Run the full deploy + verify chain (phase8 LAB+100, phase8e mdx, phase9)."""
    rc = cmd_phase8(argparse.Namespace(target="lab,100"))
    if rc != 0: return rc
    rc = cmd_phase8e(args)
    if rc != 0:
        warn("Phase 8e failed (PyPI not yet live or mdx unreachable). "
              "Phase 9 will be skipped — re-run later when PyPI propagates.")
        return rc
    return cmd_phase9(args)


def cmd_done(args):
    """Definition-of-done check: preflight (repo) + phase9 (3 machines).

    Read-only. Exit 0 means the release is consistent across LAB / 100号機 /
    mdx and the repo is in a clean release-ready state. Exit non-zero
    means do NOT tell the user "release done" yet.
    """
    step("Definition-of-done check (preflight + phase9)")
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed — repo state not release-ready.")
        return rc
    rc = cmd_phase9(args)
    if rc != 0:
        fail("phase9 drift detected — at least one machine is out of sync.")
        return rc
    print("")
    ok("DEFINITION OF DONE met. Release is consistent across LAB / 100号機 / mdx.")
    return 0


# ============================================================
# CLI
# ============================================================

def main():
    p = argparse.ArgumentParser(prog="release_triple",
                                 description="Enforce the release-triple flow.")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("preflight",
                    help="read-only state report (always safe)")
    sub.add_parser("phase0",
                    help="clean rebuild of the Cubit plugin")
    s8 = sub.add_parser("phase8",
                         help="deploy + verify + smoke on LAB / 100号機")
    s8.add_argument("--target", default="lab,100",
                     help="comma list: lab, 100, all (default lab,100)")
    sub.add_parser("phase8e",
                    help="upgrade mdx from PyPI (after PyPI propagation)")
    sub.add_parser("phase9",
                    help="cross-machine consistency probe")
    sub.add_parser("all",
                    help="phase8 -> phase8e -> phase9 in one shot")
    sub.add_parser("done",
                    help="definition-of-done: preflight + phase9 (read-only)")

    args = p.parse_args()
    handler = {
        "preflight": cmd_preflight,
        "phase0":    cmd_phase0,
        "phase8":    cmd_phase8,
        "phase8e":   cmd_phase8e,
        "phase9":    cmd_phase9,
        "all":       cmd_all,
        "done":      cmd_done,
    }[args.cmd]
    raise SystemExit(handler(args))


if __name__ == "__main__":
    main()
