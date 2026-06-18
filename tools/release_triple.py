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
import shutil
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


def remove_tree(path: Path):
    """Remove a build directory without relying on POSIX rm being on PATH."""
    print(f"  $ remove-tree {path}")
    shutil.rmtree(path)


def copy_file(src: Path, dst: Path):
    """Copy a file without relying on POSIX cp being on PATH."""
    print(f"  $ copy-file {src} {dst}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


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
    """Newest mtime of bundled .ccm in cubit-mesh-export package.

    Note (radia 4.80.0): the .ccl was removed (Qt5 GUI deleted; PySide6
    toolbar at src/radia/panels/radia_export_menu.py replaces it).  The
    freshness gate now tracks only the .ccm.
    """
    pkg = REPO / "packages/cubit-mesh-export/src/cubit_mesh_export"
    times = []
    for name in ("cubit_mesh_export.ccm",):
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
        fail("bundled .ccm missing — Phase 0 not done")
    elif bin_mtime + 1 < src_mtime:
        from datetime import datetime
        fail(f"bundled .ccm ({datetime.fromtimestamp(bin_mtime)}) older than "
              f"src/cubit_plugin/ ({datetime.fromtimestamp(src_mtime)}). "
              "Run `python tools/release_triple.py phase0`.")
        return 2
    else:
        ok("bundled plugin .ccm >= src/cubit_plugin/ mtime")

    return 0


def cmd_phase0(args):
    """Clean rebuild of Cubit plugin (.ccm + .pyd).

    Note (radia 4.80.0): the .ccl target was removed (Qt5 GUI deleted;
    PySide6 toolbar at src/radia/panels/radia_export_menu.py replaces
    it).  Phase 0 now builds only the Qt-independent .ccm + .pyd.
    """
    step("Phase 0: clean rebuild of Cubit plugin (~2-3 min)")
    build_pyd = REPO / "src/cubit_plugin/build-pyd"
    build_ccm = REPO / "src/cubit_plugin/build-ccm"
    if build_pyd.exists():
        remove_tree(build_pyd)
    if build_ccm.exists():
        remove_tree(build_ccm)

    # Build via the same ps1 we used in the 2026-04-14 manual run.
    ps1 = REPO / "tools/_build_cubit_plugin.ps1"
    if not ps1.is_file():
        fail(f"missing helper script: {ps1}")
        return 3
    run(["pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(ps1)])

    # Propagate to the cubit-mesh-export package ONLY (Tier-2, 2026-06-01):
    # cme is the sole shipper of the Cubit plugin binary; radia no longer
    # bundles cubit_mesh_export.ccm, so radia + cme release fully independently.
    # Only .ccm here (the .ccl target was removed in radia 4.80.0; the .pyd
    # is propagated by the full Build.ps1, not this fast .ccm-only phase0).
    for src_name, dst_dirs in [
        ("build-ccm/cubit_mesh_export.ccm", ["packages/cubit-mesh-export/src/cubit_mesh_export"]),
    ]:
        src = REPO / "src/cubit_plugin" / src_name
        if not src.is_file():
            fail(f"build did not produce {src}")
            return 3
        for d in dst_dirs:
            dst = REPO / d / src.name
            copy_file(src, dst)

    ok("Phase 0 complete; .ccm propagated to cubit-mesh-export pkg (radia no longer bundles it)")
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


def _check_pypi_propagation(versions):
    """Refuse to deploy if PyPI hasn't propagated to repo's current versions.

    Returns 0 on success, 2 if any package is stale.  Used by every PyPI
    install target (100号機 + mdx) to prevent installing the OLD version
    while CI is still publishing the new one.
    """
    info("checking PyPI propagation...")
    for pkg, want in [("radia", versions["radia"]),
                       ("cubit-mesh-export", versions["cubit-mesh-export"]),
                       ("radia-mcp", versions["radia-mcp"])]:
        p = run(["python", "-m", "pip", "index", "versions", pkg],
                capture=True, check=False)
        first = p.stdout.splitlines()[0] if p.stdout else ""
        if want and want in first:
            ok(f"PyPI {pkg} live at {want}")
        else:
            fail(f"PyPI {pkg} not yet at {want} (got {first!r}). "
                 "Wait for CI / PyPI propagation, then retry.")
            return 2
    return 0


def _deploy_pypi(ssh_host, label):
    """PyPI-install recipe for downstream Cubit-equipped machines.

    Used identically for 100号機 and mdx per the 2-tier distribution
    policy (CLAUDE.md: "LAB = NAS editable, 100号機 + mdx = PyPI").
    Before this refactor, _deploy_100 used `pip install <NAS path>` which
    SHIPPED LAB WORKTREE BITS (including unreleased dev bumps) onto the
    21-user shared lab box -- a quiet policy violation that defeated the
    point of the "PyPI is the verified release channel" gate.

    Recipe:
      1. PyPI propagation check (refuse if stale)
      2. force-kill Cubit + mcp-server-*.exe (otherwise pip install blocks
         on locked Scripts/mcp-server-*.exe)
      3. pip install --upgrade --no-cache-dir from PyPI, pinned to the
         repo's current versions
      4. cubit-plugin-install --all-users (regular-file deploy of the
         freshly-installed wheel's plugin)
      5. cubit-plugin-install --verify-only (sha256 sanity)
      6. cubit-smoke-test (Cubit 2025.3 -batch run on ih_bem_sample.jou)
    """
    step(f"Phase 8 ({label}): kill + PyPI install + plugin install + verify + smoke (over SSH)")
    v = _read_repo_versions()
    rc = _check_pypi_propagation(v)
    if rc != 0:
        return rc

    v_radia = v["radia"]
    v_cme   = v["cubit-mesh-export"]
    v_mcp   = v["radia-mcp"]

    ps_block = (
        "$ErrorActionPreference = 'Continue'; "
        "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
        "$_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit' "
        "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
        "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
        "$_.Name -like 'mcp-server*' "
        "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
        "Start-Sleep -Seconds 2; "
        # --force-reinstall is mandatory: `pip install --upgrade X==Y` on a
        # machine already at X==Y is a NO-OP and leaves the on-disk files
        # untouched. That bit us 2026-05-26: 100号機 was already at
        # cubit-mesh-export==0.10.1 from a prior NAS-source install, so
        # `--upgrade 0.10.1` did nothing and the worktree binary
        # (1fd45675...) stayed on 100号機 even though PyPI 0.10.1 had a
        # different binary (ef49da18...). Force-reinstall guarantees the
        # PyPI wheel's bytes overwrite whatever is on disk, which is the
        # whole point of "PyPI is the canonical channel" in the 2-tier policy.
        f"pip install --upgrade --force-reinstall --no-deps --no-cache-dir "
        f"'radia[cubit,gui]=={v_radia}' "
        f"'radia-mcp=={v_mcp}' 'cubit-mesh-export=={v_cme}'; "
        "cubit-plugin-install --all-users; "
        "cubit-plugin-install --verify-only; "
        "cubit-smoke-test"
    )
    run(["ssh", ssh_host, "pwsh", "-ExecutionPolicy", "Bypass",
         "-Command", ps_block])
    ok(f"Phase 8 complete on {label}")
    return 0


def _deploy_100():
    return _deploy_pypi(SSH_100, "100号機")


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
        if t == "lab":
            _deploy_lab()
        elif t in ("100", "100号機", "100goki"):
            rc = _deploy_100()
            if rc != 0:
                return rc
        elif t == "all":
            _deploy_lab()
            rc = _deploy_100()
            if rc != 0:
                return rc
        else:
            fail(f"unknown target: {t!r}")
            return 2
    return 0


def cmd_phase8e(args):
    """Upgrade mdx from PyPI (only after PyPI propagation).  Same recipe
    as Phase 8 for 100号機; the only difference is the SSH host.  Both
    are PyPI verification points per CLAUDE.md 2-tier policy."""
    return _deploy_pypi(SSH_MDX, "mdx")


CROSS_MACHINE_PROBE = '''import hashlib, os
import importlib.metadata as md

def hsh_text(p):
    h = hashlib.sha256()
    try:
        d = open(p, "rb").read().replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
        h.update(d); return h.hexdigest()[:12]
    except Exception: return "MISSING"

def ver(n):
    try: return md.version(n)
    except Exception: return "MISSING"

import radia, cubit_mesh_export
rad = os.path.dirname(radia.__file__)
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
'''


# LAB probe (2026-05-28 fix).  LAB is the editable DEV checkout, so
# os.path.dirname(radia.__file__) is the *working tree* -- full of
# uncommitted dev WIP that has nothing to do with the released wheel.
# Hashing that guaranteed perpetual false-positive drift (the whole reason
# this variant exists).  Instead, hash each tracked panel file as it exists
# at the RELEASE TAG (v{radia.__version__}) via `git show`: byte-identical to
# what the consumers' wheel was built from, immune to (a) uncommitted WIP and
# (b) post-release commits on main.  Versions/COMPAT come from installed
# metadata (re-synced to the release in Phase 8), identical to the consumer
# probe so the 10 rows line up 1:1 for the row-by-row comparison.
CROSS_MACHINE_PROBE_LAB = '''import hashlib, os, subprocess
import importlib.metadata as md

def ver(n):
    try: return md.version(n)
    except Exception: return "MISSING"

import radia, cubit_mesh_export
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(radia.__file__))))
tag = "v" + radia.__version__

def hsh_git(relpath):
    try:
        d = subprocess.run(["git", "-C", root, "show", tag + ":" + relpath],
                           capture_output=True).stdout
        NL = bytes([10]); CR = bytes([13])
        d = d.replace(CR + NL, NL).replace(CR, NL)
        if not d: return "MISSING"
        h = hashlib.sha256(); h.update(d); return h.hexdigest()[:12]
    except Exception: return "MISSING"

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
    print(f"SHA radia/{r:35s} = {hsh_git('src/radia/' + r)}")
'''


def _probe(host_label, cmd_prefix, probe_src=CROSS_MACHINE_PROBE):
    """Run the probe on a target (cmd_prefix is the python invocation).

    probe_src defaults to the consumer probe (hashes the installed wheel
    files).  LAB passes CROSS_MACHINE_PROBE_LAB (hashes tracked files at the
    release tag via git) -- see those probe strings for the rationale.
    """
    p = subprocess.run(cmd_prefix, input=probe_src,
                        capture_output=True, text=True, shell=False)
    if p.returncode != 0:
        fail(f"probe failed on {host_label}: {p.stderr.strip()}")
        return None
    return p.stdout


def cmd_phase9(args):
    """Cross-machine consistency probe."""
    step("Phase 9: cross-machine consistency (LAB / 100号機 / mdx)")
    out_lab = _probe("LAB", ["python", "-"], CROSS_MACHINE_PROBE_LAB)
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


# ============================================================
# LAB editable-install verifier (POLICY 2026-05-27, CLAUDE.md
# "release 後の LAB editable 再確認")
# ============================================================

# (package name, expected Editable project location prefix on LAB)
# Path-prefix match (case-insensitive, slash-normalised) so a UNC vs
# drive-letter representation of the same NAS path is accepted.
LAB_EDITABLE_PKGS = [
    ("radia",               "S:/Radia/01_GitHub"),
    ("cubit-mesh-export",   "S:/Radia/01_GitHub/packages/cubit-mesh-export"),
    ("radia-mcp",           "S:/Radia/01_GitHub/packages/radia-mcp"),
    # LAB-private; tolerated as missing if pip show says not installed.
    ("mcp-server-document", "S:/mcp-server"),
]


def _norm_path(p):
    """Lower-case + forward-slashes + strip trailing slash so a UNC and
    a drive-letter form of the same NAS path compare equal."""
    p = (p or "").replace("\\", "/").rstrip("/").lower()
    # Treat the NAS UNC (//192.168.11.100/work/00_cae/radia/01_github) as
    # equivalent to S:/radia/01_github -- both resolve to the same files.
    return p.replace("//192.168.11.100/work/00_cae/radia/01_github",
                     "s:/radia/01_github")


def _pip_show(pkg):
    """Return parsed dict from `pip show <pkg>`, or None if not
    installed.  Keys are lower-cased; values are stripped strings."""
    p = run(["pip", "show", pkg], check=False, capture=True)
    if p.returncode != 0:
        return None
    d = {}
    for line in (p.stdout or "").splitlines():
        if ":" in line:
            k, _, v = line.partition(":")
            d[k.strip().lower()] = v.strip()
    return d


def _verify_lab_editable():
    """Check the 4 LAB-editable packages still point at NAS source.

    Returns (n_ok, n_drift, n_missing, details) tuple.  Drift is the
    one we care about for the "release-done" gate -- it means LAB
    cannot dev-loop on that package because edits won't reflect.

    Missing for `mcp-server-document` is OK (LAB-private, may not be
    installed on every developer's box).  Missing for `radia` /
    `cubit-mesh-export` / `radia-mcp` is itself a drift state (LAB
    should always have these editable).
    """
    step("LAB editable verify (POLICY 2026-05-27)")
    n_ok = n_drift = n_missing = 0
    details = []

    for pkg, want_prefix in LAB_EDITABLE_PKGS:
        d = _pip_show(pkg)
        if d is None:
            if pkg == "mcp-server-document":
                info(f"{pkg:25s}  not installed  (LAB-private; tolerated)")
                n_missing += 1
            else:
                fail(f"{pkg:25s}  NOT INSTALLED  -- expected editable @ "
                     f"{want_prefix}")
                n_drift += 1
                details.append((pkg, "not_installed", want_prefix))
            continue

        version = d.get("version", "?")
        editable = d.get("editable project location") \
                   or d.get("editable-project-location") \
                   or d.get("location")
        if not editable:
            fail(f"{pkg:25s}  v{version}  -- no Location field in "
                 f"pip show output")
            n_drift += 1
            details.append((pkg, "no_location", "?"))
            continue

        if d.get("editable project location"):
            kind = "editable"
        else:
            kind = "non-editable (Location only)"

        got_norm = _norm_path(editable)
        want_norm = _norm_path(want_prefix)
        if got_norm == want_norm:
            ok(f"{pkg:25s}  v{version}  {kind}  -> {editable}")
            n_ok += 1
        else:
            fail(f"{pkg:25s}  v{version}  DRIFT  {kind}\n"
                 f"        got:      {editable}\n"
                 f"        expected: {want_prefix}")
            n_drift += 1
            details.append((pkg, "drift", editable))

    print("")
    if n_drift == 0:
        ok(f"All {n_ok} LAB-editable package(s) point at LAB source"
           + (f" ({n_missing} LAB-private skipped)" if n_missing else ""))
    else:
        fail(f"{n_drift} LAB-editable package(s) DRIFTED.  Fix:")
        for pkg, why, _got in details:
            print(f"        # {pkg} ({why})")
            print(f"        Get-Process | Where-Object {{ $_.Name -like "
                  f"'mcp-server*' }} | Stop-Process -Force")
            print(f"        pip uninstall -y {pkg}")
            print(f"        pip install -e {dict(LAB_EDITABLE_PKGS)[pkg]} "
                  f"--no-deps --no-cache-dir")
        print("")
        print("        See CLAUDE.md \"POLICY (2026-05-27): release 後の "
              "LAB editable 再確認\" for the full recovery procedure.")
    return n_drift


def cmd_verify_editable(args):
    """Standalone editable verifier (no preflight, no phase9).

    Use this between deploys, or any time pip operations have run on
    LAB and you want to confirm dev-loop integrity is intact.
    """
    return _verify_lab_editable()


# ============================================================
# Phase 5.5 gate: CI-green BEFORE tagging (gh-free)
# ============================================================
# The self-hosted [windows-radia] runner runs ON LAB, so we read CI
# state from the local runner workspace instead of `gh` (abandoned on
# LAB 2026-05-28; not on PATH).  This is the enforcement behind the
# "push main -> confirm CI green -> THEN tag" policy that stops the
# v4.80.0 -> v4.80.5 version-burning (tag CI failing on a broken commit).
CI_WORKSPACE = r"C:\actions-runner\_work\Radia\Radia"


def _ci_worker_running():
    """True iff the GitHub Actions Runner.Worker (a running CI job) exists."""
    p = subprocess.run(
        ["tasklist", "/FI", "IMAGENAME eq Runner.Worker.exe", "/NH"],
        capture_output=True, text=True)
    return "Runner.Worker" in (p.stdout or "")


def _git_repo_owner_name():
    """Extract 'owner/name' from `git config remote.origin.url`."""
    import re
    url = subprocess.check_output(
        ["git", "config", "--get", "remote.origin.url"],
        cwd=str(REPO), text=True).strip()
    m = re.search(r"github\.com[:/]([^/\s]+)/([^/\s.]+?)(?:\.git)?$", url)
    if not m:
        raise RuntimeError(f"cannot parse GitHub repo from origin {url!r}")
    return f"{m.group(1)}/{m.group(2)}"


def _check_github_hosted_workflows(sha, *, timeout_sec=1800, poll_sec=20):
    """Poll the GitHub check-runs API for `sha` until all complete, then
    verify every conclusion is green.

    Closes the documented cmd_ci_verify caveat: this catches
    policy-lint.yml and radia-mcp-matrix.yml which run on
    github-hosted ubuntu-latest and leave nothing in CI_WORKSPACE.
    Public-repo check-runs endpoint does not require authentication.

    Historical incident (2026-05-30): policy-lint Policy 4 was silently
    red since commit 6c50c4cc (rad_stream_function.cpp added without
    updating the CblasColMajor exception list).  cmd_ci_verify reported
    GREEN because the junit XMLs from self-hosted Windows CI were fine.
    The user saw RED in the GitHub UI for weeks.  This helper closes
    that blind spot.

    Returns (ok: bool, message: str).
    """
    import urllib.request, urllib.error, json, time as _time

    repo = _git_repo_owner_name()
    url = (f"https://api.github.com/repos/{repo}/commits/{sha}/check-runs"
           "?per_page=100")
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "radia-release_triple",
    }

    def _fetch():
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read())

    started = _time.time()
    deadline = started + timeout_sec
    last_pending = []
    while True:
        try:
            data = _fetch()
        except urllib.error.HTTPError as e:
            return False, f"GitHub API HTTPError {e.code}: {e.reason}"
        except urllib.error.URLError as e:
            return False, f"GitHub API network error: {e.reason}"
        except Exception as e:
            return False, f"GitHub API error: {type(e).__name__}: {e}"

        runs = data.get("check_runs", [])
        # Push-triggered workflows may take ~30 s to register; allow up
        # to 90 s before declaring "no runs found" (the push never
        # triggered any github-hosted workflow, e.g. paths filter
        # excluded them, or workflows are disabled).
        if not runs:
            if _time.time() - started > 90:
                return False, ("no check-runs registered for "
                               f"{sha[:8]} after 90 s (paths filter?)")
            _time.sleep(poll_sec)
            continue

        pending = [r for r in runs if r["status"] != "completed"]
        if not pending:
            break
        last_pending = pending
        if _time.time() > deadline:
            return False, ("timeout: github-hosted runs still pending: "
                           + ", ".join(r["name"] for r in last_pending))
        _time.sleep(poll_sec)

    # Every run completed. Green conclusions: success / skipped / neutral.
    failures = [r for r in runs
                if r["conclusion"] not in ("success", "skipped", "neutral")]
    if failures:
        msg = "; ".join(f"{r['name']}: {r['conclusion']}" for r in failures)
        return False, "github-hosted CI RED -- " + msg
    names = sorted(set(r["name"] for r in runs))
    return True, (f"all {len(runs)} github-hosted check-runs GREEN "
                  f"({', '.join(names)})")


def cmd_ci_verify(args):
    """gh-free CI-green gate.  Run AFTER `git push origin main`, BEFORE tags.

    1. Wait for the self-hosted runner job (Runner.Worker) to appear, then
       finish.  If no run starts within 10 min, fail (did the push trigger
       CI?) -- do NOT tag on an unverified commit.
    2. Assert every junit XML in the CI workspace is FRESH (written by THIS
       run, not a stale prior one) AND failures=errors=0.

    Exit 0 = green (safe to create + push the release tags).  Non-zero =
    red / unverified (fix-forward on main, re-run, do NOT tag).

    CAVEAT (documented, not authoritative): this covers the build + test
    steps -- a build failure leaves the XMLs stale/absent (reads as
    not-green), but a failure in a post-test step is not caught.  gh's
    check-runs API would be authoritative; gh is unavailable on LAB.
    """
    import time
    import datetime as _dt
    from xml.etree import ElementTree as ET

    step("Phase 5.5: CI verify (gh-free) -- wait for runner, then check test XMLs")
    t0 = time.time()
    saw = False
    while True:
        now = time.time()
        if _ci_worker_running():
            if not saw:
                print("  Runner.Worker active -- CI job running; waiting for it to finish...")
            saw = True
        elif saw:
            print("  Runner.Worker exited -- CI job complete.")
            break
        elif now - t0 > 600:
            fail("no CI run started within 10 min -- did `git push origin main` "
                 "trigger a run? Not verified; do NOT tag.")
            return 4
        if now - t0 > 40 * 60:
            fail("CI run did not finish within 40 min; not verified.")
            return 4
        time.sleep(20)

    import glob
    xmls = sorted(glob.glob(os.path.join(CI_WORKSPACE, "*results*.xml")))
    if not xmls:
        fail(f"no *results*.xml in {CI_WORKSPACE} -- CI produced no test output "
             "(the build step likely failed). NOT green.")
        return 4

    bad = 0
    print(f"\n  {'xml':<30} | {'tests':>6} | {'fail':>4} | {'err':>4} | mtime")
    print("  " + "-" * 72)
    for x in xmls:
        try:
            root = ET.parse(x).getroot()
            suites = root.findall("testsuite") or [root]
            tests = sum(int(s.get("tests", 0) or 0) for s in suites)
            fails = sum(int(s.get("failures", 0) or 0) for s in suites)
            errs = sum(int(s.get("errors", 0) or 0) for s in suites)
        except Exception as e:
            fail(f"could not parse {os.path.basename(x)}: {e}")
            bad += 1
            continue
        mt = os.path.getmtime(x)
        fresh = mt >= (t0 - 180)   # written by this run (allow 3 min pre-slack)
        mts = _dt.datetime.fromtimestamp(mt).strftime("%H:%M:%S")
        is_ok = (fails == 0 and errs == 0 and fresh)
        if not is_ok:
            bad += 1
        print(f"  {os.path.basename(x):<30} | {tests:>6} | {fails:>4} | {errs:>4} | "
              f"{mts}{'  STALE' if not fresh else ''}")
    if bad:
        print("")
        fail(f"{bad} test XML(s) failing or stale -- CI NOT green. "
             "Fix-forward on main and re-run; do NOT tag.")
        return 4
    print("")
    ok("self-hosted (build-test) CI GREEN "
       "(all test XMLs fresh, failures=errors=0).")

    # ALSO verify github-hosted workflows (policy-lint, radia-mcp-matrix).
    # See _check_github_hosted_workflows for the why (2026-05-30 incident).
    step("CI verify (github-hosted): policy-lint / radia-mcp-matrix")
    head_sha = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=str(REPO), text=True).strip()
    print(f"  HEAD = {head_sha[:8]}")
    gh_ok, gh_msg = _check_github_hosted_workflows(head_sha)
    print("  " + gh_msg)
    if not gh_ok:
        fail("github-hosted CI is RED -- inspect at "
             f"github.com/{_git_repo_owner_name()}/actions, fix-forward.")
        return 5

    print("")
    ok("CI is fully GREEN (self-hosted + github-hosted). "
       "Safe to create + push the release tags (Phase 6).")
    return 0


def _run_pyside6_health_audit():
    """Layer A pyside6-health audit (gated in `done` since 2026-05-30).

    Runs tools/audit_pyside6_only.py which checks:
      1. zero real PyQt5 / PySide2 / PyQt6 imports in tracked *.py
      2. core GUI modules import PySide6
      3. radia pyproject declares PySide6, no PyQt5 dependency
      4. cubit_mesh_export.ccm has no Qt5 DLL dependency
      5. headless panel smoke (IH/EM/PCB construct + ExportDialog
         builds all 6 formats) in an isolated offscreen subprocess

    Returns the audit script's exit code (0 = CLEAN).
    """
    step("pyside6-health audit (Layer A: static + headless smoke)")
    audit = REPO / "tools" / "audit_pyside6_only.py"
    if not audit.exists():
        fail(f"audit script missing: {audit}")
        return 1
    result = subprocess.run(
        [sys.executable, str(audit)],
        cwd=str(REPO),
        capture_output=True, text=True,
    )
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip())
    if result.returncode == 0:
        ok("pyside6-health: CLEAN (PySide6-only + panels construct).")
    return result.returncode


def cmd_done(args):
    """Definition-of-done check: preflight + LAB-editable verify + phase9 + pyside6-health.

    Read-only. Exit 0 means the release is consistent across LAB / 100号機 /
    mdx, the repo is release-ready, the LAB dev loop is intact, AND the
    panel GUI is PySide6-only with every panel window constructible
    headlessly. Exit non-zero means do NOT tell the user "release done" yet.
    """
    step("Definition-of-done check "
         "(preflight + LAB editable + phase9 + pyside6-health)")
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed — repo state not release-ready.")
        return rc

    drift = _verify_lab_editable()
    if drift > 0:
        fail(f"{drift} LAB-editable package(s) drifted.  "
             "Run the printed recovery commands, then re-run "
             "`release_triple done`.")
        return 1

    rc = cmd_phase9(args)
    if rc != 0:
        fail("phase9 drift detected — at least one machine is out of sync.")
        return rc

    rc = _run_pyside6_health_audit()
    if rc != 0:
        fail("pyside6-health audit failed -- a Qt5/PyQt5 reference or "
             "panel construction regression slipped in. Fix per the audit "
             "output, then re-run `release_triple done`.")
        return rc

    print("")
    ok("DEFINITION OF DONE met. Release is consistent across LAB / 100号機 / "
       "mdx, LAB dev loop is intact, and the panel GUI is PySide6-only.")
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
    sub.add_parser("verify-editable",
                    help="LAB editable-install pointers check (read-only)")
    sub.add_parser("ci-verify",
                    help="Phase 5.5: gh-free CI-green gate (run after push main, before tag)")
    sub.add_parser("done",
                    help="definition-of-done: preflight + LAB editable + "
                         "phase9 (read-only)")

    args = p.parse_args()
    handler = {
        "preflight":        cmd_preflight,
        "phase0":           cmd_phase0,
        "phase8":           cmd_phase8,
        "phase8e":          cmd_phase8e,
        "phase9":           cmd_phase9,
        "all":              cmd_all,
        "verify-editable":  cmd_verify_editable,
        "ci-verify":        cmd_ci_verify,
        "done":             cmd_done,
    }[args.cmd]
    raise SystemExit(handler(args))


if __name__ == "__main__":
    main()
