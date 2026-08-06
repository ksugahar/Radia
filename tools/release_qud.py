#!/usr/bin/env python
"""release_qud.py — orchestrator for the 3-package / 4-machine release flow.

Walks Phase 0 -> 9 of the release-qud skill in order, gating each
phase on the success of the previous one. Refuses to skip steps that
have caused real outages (2026-04-14 incident series).

Usage:
    python tools/release_qud.py preflight
        Read-only: report current state and consistency. Use anytime.

    python tools/release_qud.py phase0
        Mandatory clean rebuild of the Cubit plugin (~3-4 min).

    python tools/release_qud.py phase8 [--target lab|100|hibino|all]
        Run Phase 8a..8d on each target: kill Cubit, install by the
        target's tier (LAB/100 editable, hibino PyPI), cubit-plugin-install,
        --verify-only, cubit-smoke-test. Refuses
        to start if Phase 0 has not been done since the last source
        change in src/cubit_plugin/.

    python tools/release_qud.py phase8e
        Upgrade mdx from PyPI. Refuses to run if pip index versions
        radia / cubit-mesh-export don't match the local repo
        (i.e. PyPI hasn't propagated yet). radia-mcp is intentionally
        not installed on mdx -- and is actively uninstalled if a prior
        release left it behind (mdx is a compute consumer, no MCP).

    python tools/release_qud.py phase9
        Cross-machine consistency probe. Final gate.

    python tools/release_qud.py simulink-candidate --package <zip> --target all
        Extract and execute the exact Simulink package on all four MATLAB machines.

    python tools/release_qud.py all
        phase8 -> phase8e -> phase9 with all preconditions enforced.

        When the canonical LAB worktree contains parallel WIP, set
        RADIA_RELEASE_EDITABLE_REPO_LAB and RADIA_RELEASE_EDITABLE_REPO_100
        to the LAB and 100-machine views of one clean NAS release worktree.
        Both editable deployments then fail before install unless that
        worktree is tracked-clean and matches the invoking release SHA.

    python tools/release_qud.py done --simulink-package <zip>
        Require both the normal release gate and the matching four-machine
        Simulink candidate state.

Exit codes:
    0   success
    2   precondition failure (skip detected)
    3   action failure (external command)
    4   verification mismatch
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import shutil
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import url2pathname

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
NAS_REPO_100 = r"W:\00_CAE\Radia\01_GitHub"
EDITABLE_REPO_LAB_ENV = "RADIA_RELEASE_EDITABLE_REPO_LAB"
EDITABLE_REPO_100_ENV = "RADIA_RELEASE_EDITABLE_REPO_100"
SSH_100 = "192.168.11.100"
SSH_MDX = "mdx"
SSH_HIBINO = "hibino"
PY_HIBINO = "py -3.12"
MATLAB_EXE = r"C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
SIMULINK_GATE_ROOT = Path(r"C:\temp\radia-release-qud")
SIMULINK_TARGETS = {
    "lab": ("LAB", None, "python"),
    "100": ("100号機", SSH_100, "python"),
    "mdx": ("mdx", SSH_MDX, "python"),
    "hibino": ("hibino", SSH_HIBINO, "py -3.12"),
}


def _editable_repo_lab():
    """Return the LAB editable source, allowing an exact release worktree."""
    return os.environ.get(EDITABLE_REPO_LAB_ENV, NAS_REPO_LAB).strip().rstrip("/\\")


def _editable_repo_100():
    """Return the 100-machine view of the same release worktree."""
    return os.environ.get(EDITABLE_REPO_100_ENV, NAS_REPO_100).strip().rstrip("/\\")


def _release_head():
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=str(REPO), text=True
    ).strip().lower()


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

    Note: the retired Qt5 .ccl target is gone. The Cubit-embedded PySide
    toolbar is Python package data and is checked by deploy probes, not by this
    compiled-plugin freshness gate.
    """
    pkg = REPO / "packages/cubit-mesh-export/src/cubit_mesh_export"
    times = []
    for name in ("cubit_mesh_export.ccm",):
        p = pkg / name
        if p.is_file():
            times.append(p.stat().st_mtime)
    return max(times) if times else 0.0


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _simulink_manifest(package: Path) -> dict:
    verifier = REPO / "tools/verify_simulink_release.py"
    result = subprocess.run(
        [sys.executable, str(verifier), str(package), "--manifest-only"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    with zipfile.ZipFile(package) as bundle:
        return json.loads(bundle.read("manifest.json"))


def _simulink_state_path(package_sha256: str) -> Path:
    return SIMULINK_GATE_ROOT / f"simulink-{package_sha256}.json"


def _write_simulink_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def _run_simulink_candidate_target(
        key: str, package: Path, package_sha256: str,
        success_marker: str) -> tuple[bool, str]:
    label, host, python_command = SIMULINK_TARGETS[key]
    verifier = REPO / "tools/verify_simulink_release.py"
    if host is None:
        command = [
            sys.executable, str(verifier), str(package),
            "--matlab", MATLAB_EXE,
        ]
        result = subprocess.run(command, capture_output=True, text=True)
    else:
        remote_root_posix = f"C:/temp/radia-release-qud/{package_sha256[:16]}"
        remote_root_windows = remote_root_posix.replace("/", "\\")
        prepare = f"New-Item -ItemType Directory -Force -Path '{remote_root_windows}' | Out-Null\n"
        created = subprocess.run(
            ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-Command", "-"],
            input=prepare, capture_output=True, text=True,
        )
        if created.returncode != 0:
            return False, created.stderr.strip() or created.stdout.strip()
        remote_package = f"{remote_root_posix}/{package.name}"
        remote_verifier = f"{remote_root_posix}/{verifier.name}"
        for source, destination in (
                (package, remote_package), (verifier, remote_verifier)):
            copied = subprocess.run(
                ["scp", str(source), f"{host}:{destination}"],
                capture_output=True, text=True,
            )
            if copied.returncode != 0:
                return False, copied.stderr.strip() or copied.stdout.strip()
        invocation = (
            f"& {python_command} '{remote_verifier}' '{remote_package}' "
            f"--matlab '{MATLAB_EXE}'\nexit $LASTEXITCODE\n"
        )
        result = subprocess.run(
            ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-Command", "-"],
            input=invocation, capture_output=True, text=True,
        )
    output = ((result.stdout or "") + (result.stderr or "")).strip()
    if result.returncode != 0:
        return False, output
    if success_marker not in output or '"status": "passed"' not in output:
        return False, f"{label} did not emit both release success markers:\n{output}"
    return True, output


def cmd_simulink_candidate(args):
    """Verify one extracted Simulink archive on the requested MATLAB machines."""
    package = Path(args.package).resolve()
    step("Simulink candidate gate (LAB / 100号機 / mdx / hibino)")
    try:
        manifest = _simulink_manifest(package)
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as error:
        fail(f"invalid Simulink candidate: {error}")
        return 2
    package_sha256 = _sha256_file(package)
    success_marker = (
        "RADIA_SIMULINK_RELEASE_OK"
        if manifest.get("schema") in {
            "radia.simulink.library-release-manifest.v1",
            "radia.simulink.library-release-manifest.v2",
        }
        else "RADIA_IH_RELEASE_OK"
    )
    state_path = _simulink_state_path(package_sha256)
    state = {
        "schema": "radia.release-qud.simulink-candidate.v1",
        "package": str(package),
        "package_sha256": package_sha256,
        "version": manifest.get("version"),
        "commit": manifest.get("commit"),
        "targets": {},
    }
    if state_path.is_file():
        previous = json.loads(state_path.read_text(encoding="utf-8"))
        if previous.get("package_sha256") == package_sha256:
            state["targets"] = previous.get("targets", {})

    requested = [part.strip().lower() for part in args.target.split(",")]
    if "all" in requested:
        requested = list(SIMULINK_TARGETS)
    unknown = sorted(set(requested) - set(SIMULINK_TARGETS))
    if unknown:
        fail(f"unknown Simulink target(s): {', '.join(unknown)}")
        return 2

    failed = 0
    for key in requested:
        label = SIMULINK_TARGETS[key][0]
        info(f"verifying extracted package on {label}")
        passed, output = _run_simulink_candidate_target(
            key, package, package_sha256, success_marker)
        state["targets"][key] = {
            "label": label,
            "status": "passed" if passed else "failed",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "output_tail": output[-4000:],
        }
        _write_simulink_state(state_path, state)
        if passed:
            ok(f"Simulink package passed on {label}")
        else:
            failed += 1
            fail(f"Simulink package failed on {label}")
            if output:
                print(output[-4000:])
    if failed:
        return 4
    ok(f"Simulink candidate state: {state_path}")
    return 0


def _verify_simulink_candidate_state(package_arg: str) -> int:
    package = Path(package_arg).resolve()
    try:
        manifest = _simulink_manifest(package)
        package_sha256 = _sha256_file(package)
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as error:
        fail(f"invalid Simulink candidate: {error}")
        return 2
    state_path = _simulink_state_path(package_sha256)
    if not state_path.is_file():
        fail("Simulink candidate has no release-qud state. Run "
             "`release_qud simulink-candidate --package <zip> --target all`.")
        return 4
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if state.get("package_sha256") != package_sha256 or \
            state.get("commit") != manifest.get("commit"):
        fail("Simulink candidate state does not match the supplied archive")
        return 4
    missing = [key for key in SIMULINK_TARGETS
               if state.get("targets", {}).get(key, {}).get("status") != "passed"]
    if missing:
        fail(f"Simulink candidate has not passed: {', '.join(missing)}")
        return 4
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO,
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    if manifest.get("commit") != head:
        fail("Simulink manifest commit differs from HEAD; rebuild the archive")
        return 4
    ok("supplied Simulink candidate passed LAB / 100号機 / mdx / hibino")
    return 0


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
              "Run `python tools/release_qud.py phase0`.")
        return 2
    else:
        ok("bundled plugin .ccm >= src/cubit_plugin/ mtime")

    return 0


def cmd_phase0(args):
    """Clean rebuild of Cubit plugin (.ccm + .pyd).

    Note: the retired Qt5 .ccl target is gone. Phase 0 builds only the
    C++/APREPRO plugin (.ccm + .pyd); the Cubit-embedded PySide toolbar is
    shipped as Python package data.
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


_CONSOLE_SCRIPT_WRAPPER_RE = (
    r'(?i)^\s*"?[^\"]*pythonw?\.exe"?\s+"?'
    r'[^\"]*(?:mcp-server-|radia[-_])[^\"\s]*\.exe(?:[\"\s]|$)'
)


def _kill_cubit_local():
    info("force-kill any local Cubit process")
    run(["pwsh", "-NoProfile", "-Command",
         "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
         "$_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit' "
         "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
         "Start-Sleep -Seconds 2"], check=False)


def _kill_mcp_local():
    info("force-kill local mcp-server-*.exe and Radia panel launchers")
    run(["pwsh", "-NoProfile", "-Command",
         "Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | "
         "Where-Object { $_.ProcessId -ne $PID -and ("
         "$_.Name -like 'mcp-server*' -or $_.Name -like 'radia-*' -or "
         "$_.Name -like 'radia_*' -or "
         "((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and "
         f"$_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}') "
         ") } | ForEach-Object { Stop-Process -Id $_.ProcessId -Force "
         "-ErrorAction SilentlyContinue }; "
         "Start-Sleep -Seconds 2"], check=False)


def _deploy_lab():
    step("Phase 8 (LAB): kill, install from NAS, plugin install, verify, smoke")
    repo = _editable_repo_lab()
    rc = _verify_local_release_source(repo, _release_head())
    if rc != 0:
        return rc
    _kill_cubit_local()
    _kill_mcp_local()
    for sub in ("", "/packages/cubit-mesh-export", "/packages/radia-mcp"):
        run(["pip", "install", "-e", repo + sub, "--no-deps",
             "--no-cache-dir"])
    run(["cubit-plugin-install"])
    run(["cubit-plugin-install", "--verify-only"])
    run(["cubit-smoke-test"])
    ok("Phase 8 complete on LAB")
    return 0


def _deploy_editable_remote(ssh_host, label, repo):
    """Editable-install recipe for machines that should read NAS source.

    LAB and 100号機 are the editable tier.  PyPI propagation is not a
    precondition for this tier; hibino/mdx remain the wheel-consumer
    verification tier.
    """
    step(f"Phase 8 ({label}): kill + NAS editable install + plugin install + verify + smoke (over SSH)")
    expected_sha = _release_head()
    ps_block = f"""
$ErrorActionPreference = 'Continue'
$sourceHead = (& git -C "{repo}" rev-parse HEAD).Trim().ToLowerInvariant()
if ($LASTEXITCODE -ne 0 -or $sourceHead -ne "{expected_sha}") {{
  Write-Error "Release source SHA mismatch: expected {expected_sha}, got $sourceHead"
  exit 41
}}
$sourceDirty = (& git -C "{repo}" status --porcelain --untracked-files=no) -join "`n"
if ($LASTEXITCODE -ne 0 -or $sourceDirty) {{
  Write-Error "Release source has tracked changes: $sourceDirty"
  exit 42
}}
Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | Where-Object {{
  $_.ProcessId -ne $PID -and (
    $_.Name -eq 'coreform_cubit.exe' -or $_.Name -eq 'cubit.exe' -or
    $_.Name -like 'mcp-server*' -or
    $_.Name -like 'radia-*' -or $_.Name -like 'radia_*' -or
    ((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and
      $_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}')
  )
}} | ForEach-Object {{
  Write-Host "Stopping $($_.Name) pid=$($_.ProcessId)"
  Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
}}
Start-Sleep -Seconds 2
pip install -e "{repo}" -e "{repo}\\packages\\cubit-mesh-export" -e "{repo}\\packages\\radia-mcp" --no-deps --no-cache-dir
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-plugin-install --all-users
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-plugin-install --verify-only
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-smoke-test
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
"""
    encoded = base64.b64encode(ps_block.encode("utf-16le")).decode("ascii")
    run(["ssh", ssh_host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
         "-EncodedCommand", encoded])
    ok(f"Phase 8 complete on {label}")
    return 0


def _check_pypi_propagation(versions, *, include_mcp=True):
    """Refuse to deploy if PyPI hasn't propagated to repo's current versions.

    Returns 0 on success, 2 if any package is stale.  Used by every PyPI
    install target (hibino + mdx) to prevent installing the OLD version
    while CI is still publishing the new one.
    """
    info("checking PyPI propagation...")
    packages = [("radia", versions["radia"]),
                ("cubit-mesh-export", versions["cubit-mesh-export"])]
    if include_mcp:
        packages.append(("radia-mcp", versions["radia-mcp"]))
    for pkg, want in packages:
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


def _deploy_pypi(ssh_host, label, *, include_mcp=True, python_cmd="python", cubit_optional=False):
    """PyPI-install recipe for downstream Cubit-equipped machines.

    Used for hibino and mdx.  hibino gets radia-mcp; mdx is a compute
    consumer and intentionally skips the MCP server package -- and
    actively uninstalls radia-mcp if a prior release left it behind.

    Recipe:
      1. PyPI propagation check (refuse if stale)
      2. force-kill Cubit + mcp-server-*.exe (otherwise pip install blocks
         on locked Scripts/mcp-server-*.exe)
      3. pip install --upgrade --no-cache-dir from PyPI, pinned to the
         repo's current versions
      4. cubit-plugin-install --all-users (regular-file deploy of the
         freshly-installed wheel's plugin; skipped on optional-Cubit
         targets when Cubit is not installed)
      5. cubit-plugin-install --verify-only (sha256 sanity)
      6. cubit-smoke-test (Cubit 2025.3 -batch run on ih_bem_sample.jou)
    """
    step(f"Phase 8 ({label}): kill + PyPI install + plugin install + verify + smoke (over SSH)")
    v = _read_repo_versions()
    rc = _check_pypi_propagation(v, include_mcp=include_mcp)
    if rc != 0:
        return rc

    v_radia = v["radia"]
    v_cme   = v["cubit-mesh-export"]
    v_mcp   = v["radia-mcp"]
    mcp_pin = f' "radia-mcp=={v_mcp}"' if include_mcp else ""
    # mdx is a compute consumer: radia-mcp must NOT be present there.  Older
    # (pre-policy) releases left radia-mcp installed, so on any
    # include_mcp=False target actively uninstall it rather than merely
    # skipping the install -- "absent" is the enforced invariant, not a
    # passive side effect.  `pip uninstall -y` on an already-absent package
    # is a no-op (exit 0), so this is safe to run unconditionally there.
    mcp_uninstall_ps = (
        "" if include_mcp
        else f"{python_cmd} -m pip uninstall -y radia-mcp\n"
    )
    cubit_optional_ps = "$true" if cubit_optional else "$false"

    ps_block = f"""
$ErrorActionPreference = 'Continue'
Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | Where-Object {{
  $_.ProcessId -ne $PID -and (
    $_.Name -eq 'coreform_cubit.exe' -or $_.Name -eq 'cubit.exe' -or
    $_.Name -like 'mcp-server*' -or
    $_.Name -like 'radia-*' -or $_.Name -like 'radia_*' -or
    ((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and
      $_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}')
  )
}} | ForEach-Object {{
  Write-Host "Stopping $($_.Name) pid=$($_.ProcessId)"
  Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
}}
Start-Sleep -Seconds 2
{mcp_uninstall_ps}        # --force-reinstall is mandatory: `pip install --upgrade X==Y` on a
        # machine already at X==Y is a NO-OP and leaves the on-disk files
        # untouched. That bit us 2026-05-26: 100号機 was already at
        # cubit-mesh-export==0.10.1 from a prior NAS-source install, so
        # `--upgrade 0.10.1` did nothing and the worktree binary
        # (1fd45675...) stayed on 100号機 even though PyPI 0.10.1 had a
        # different binary (ef49da18...). Force-reinstall guarantees the
        # PyPI wheel's bytes overwrite whatever is on disk, which is the
        # whole point of "PyPI is the canonical channel" in the 2-tier policy.
{python_cmd} -m pip install --upgrade --force-reinstall --no-deps --no-cache-dir "radia[cubit]=={v_radia}" "cubit-mesh-export=={v_cme}"{mcp_pin}
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
$pyExe = (& {python_cmd} -c "import sys; print(sys.executable)").Trim()
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
$scriptDir = Join-Path (Split-Path -Parent $pyExe) "Scripts"
$env:PATH = $scriptDir + ";" + $env:PATH
$cubitOptional = {cubit_optional_ps}
$cubitExeCandidates = @()
if ($env:CUBIT_PATH) {{ $cubitExeCandidates += (Join-Path $env:CUBIT_PATH "coreform_cubit.exe") }}
$cubitExeCandidates += "C:\\Program Files\\Coreform Cubit 2025.12\\bin\\coreform_cubit.exe"
$cubitFound = $false
foreach ($candidate in $cubitExeCandidates) {{
  if (Test-Path $candidate) {{ $cubitFound = $true; break }}
}}
if (-not $cubitFound -and $cubitOptional) {{
  Write-Host "[WARN] Coreform Cubit 2025.12+ not found; skipping Cubit plugin install on optional-Cubit target."
  exit 0
}}
cubit-plugin-install --all-users
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-plugin-install --verify-only
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-smoke-test
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
"""
    encoded = base64.b64encode(ps_block.encode("utf-16le")).decode("ascii")
    run(["ssh", ssh_host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
         "-EncodedCommand", encoded])
    ok(f"Phase 8 complete on {label}")
    return 0


def _deploy_100():
    return _deploy_editable_remote(
        SSH_100, "100号機", _editable_repo_100()
    )


def _deploy_hibino():
    return _deploy_pypi(SSH_HIBINO, "hibino", python_cmd=PY_HIBINO, cubit_optional=True)


def cmd_phase8(args):
    """Deploy + verify + smoke on LAB, 100号機, and/or hibino."""
    # precondition: Phase 0 freshness
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed; refusing Phase 8")
        return rc

    targets = args.target.split(",") if args.target else ["lab", "100"]
    for t in targets:
        t = t.strip().lower()
        if t == "lab":
            rc = _deploy_lab()
            if rc != 0:
                return rc
        elif t in ("100", "100号機", "100goki"):
            rc = _deploy_100()
            if rc != 0:
                return rc
        elif t == "hibino":
            rc = _deploy_hibino()
            if rc != 0:
                return rc
        elif t == "all":
            rc = _deploy_lab()
            if rc != 0:
                return rc
            rc = _deploy_100()
            if rc != 0:
                return rc
            rc = _deploy_hibino()
            if rc != 0:
                return rc
        else:
            fail(f"unknown target: {t!r}")
            return 2
    return 0


def cmd_phase8e(args):
    """Upgrade mdx from PyPI (radia + cubit-mesh-export only)."""
    return _deploy_pypi(SSH_MDX, "mdx", include_mcp=False)


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
          "panels/radia_export_menu.py",
          "simulink/application.py",
          "panels/calc_inductance.py",
          "panels/calc_fem_kelvin.py",
          "panels/calc_fem_coilmesh.py"]:
    print(f"SHA radia/{r:35s} = {hsh_text(os.path.join(rad,r))}")
'''


CROSS_MACHINE_PROBE_NO_MCP = CROSS_MACHINE_PROBE.replace(
    "print(f\"VER radia-mcp          = {ver('radia-mcp')}\")",
    'print("VER radia-mcp          = N/A")',
)


# Editable-tier probe (2026-05-28 fix).  LAB/100号機 are editable DEV checkouts, so
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
          "panels/radia_export_menu.py",
          "simulink/application.py",
          "panels/calc_inductance.py",
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
    targets = [
        ("LAB", ["python", "-"], CROSS_MACHINE_PROBE_LAB),
        ("100号機", ["ssh", SSH_100, "python", "-"], CROSS_MACHINE_PROBE_LAB),
        ("mdx", ["ssh", SSH_MDX, "python", "-"], CROSS_MACHINE_PROBE_NO_MCP),
        ("hibino", ["ssh", SSH_HIBINO, "py", "-3.12", "-"], CROSS_MACHINE_PROBE),
    ]
    step("Phase 9: cross-machine consistency (LAB / 100号機 / mdx / hibino)")
    outputs = []
    for label, cmd_prefix, probe_src in targets:
        out = _probe(label, cmd_prefix, probe_src)
        if not out:
            fail(f"could not collect probe data from {label}")
            return 4
        outputs.append((label, out.strip().splitlines()))

    if not outputs:
        fail("could not collect probe data from any machine")
        return 4

    n = min(len(rows) for _, rows in outputs)
    drift = 0

    def split(s):
        k, _, v = s.partition("=")
        k = k.strip()
        # strip leading "VER " / "SHA " / "COMPAT " from key
        for prefix in ("VER ", "SHA ", "COMPAT "):
            if k.startswith(prefix):
                k = k[len(prefix):]
        return k.strip(), v.strip()

    labels = [label for label, _ in outputs]
    header = f"\n  {'field':<44} | " + " | ".join(f"{label:<18}" for label in labels)
    print(header)
    print("  " + "-" * max(110, len(header) - 2))
    for i in range(n):
        parsed = [split(rows[i]) for _, rows in outputs]
        key = parsed[0][0]
        values = [value for _, value in parsed]
        comparable = [value for value in values if value != "N/A"]
        match = bool(comparable) and all(value == comparable[0] for value in comparable)
        marker = ok.__name__.upper() if match else "DRIFT"
        marker = _color("32;1", "OK") if match else _color("31;1", "DRIFT")
        if not match:
            drift += 1
        row = f"  {key:<44} | " + " | ".join(f"{value:<18}" for value in values)
        print(f"{row}  [{marker}]")

    if drift:
        print("")
        fail(f"{drift} field(s) drift across machines — release NOT done.")
        return 4
    print("")
    ok("all fields match across LAB / 100号機 / mdx / hibino — release verified.")
    return 0


def cmd_all(args):
    """Run the full deploy + verify chain (phase8 LAB+100+hibino, phase8e mdx, phase9)."""
    rc = cmd_phase8(argparse.Namespace(target="lab,100,hibino"))
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
def _lab_editable_packages():
    root = _editable_repo_lab()
    return [
        ("radia", root),
        ("cubit-mesh-export", root + "/packages/cubit-mesh-export"),
        ("radia-mcp", root + "/packages/radia-mcp"),
        # LAB-private; tolerated as missing if pip show says not installed.
        ("mcp-server-document", "S:/mcp-server"),
    ]


def _remote_100_editable_packages():
    root = _editable_repo_100()
    return [
        ("radia", root),
        ("cubit-mesh-export", root + r"\packages\cubit-mesh-export"),
        ("radia-mcp", root + r"\packages\radia-mcp"),
    ]


def _verify_local_release_source(repo, expected_sha):
    """Fail before install unless the editable source is the exact clean SHA."""
    path = Path(repo)
    if not path.is_dir():
        fail(f"LAB editable release source does not exist: {repo}")
        return 2
    head = run(
        ["git", "-C", repo, "rev-parse", "HEAD"],
        capture=True, check=False,
    )
    got_sha = (head.stdout or "").strip().lower()
    if head.returncode != 0 or got_sha != expected_sha.lower():
        fail(
            "LAB editable release source SHA mismatch: "
            f"expected {expected_sha}, got {got_sha or '<unavailable>'}"
        )
        return 4
    dirty = run(
        ["git", "-C", repo, "status", "--porcelain", "--untracked-files=no"],
        capture=True, check=False,
    )
    tracked_changes = (dirty.stdout or "").strip()
    if dirty.returncode != 0 or tracked_changes:
        fail(
            "LAB editable release source has tracked changes: "
            f"{tracked_changes or '<status unavailable>'}"
        )
        return 4
    ok(f"LAB editable release source is exact and clean ({got_sha[:10]})")
    return 0


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
    try:
        dist = importlib_metadata.distribution(pkg)
    except importlib_metadata.PackageNotFoundError:
        return None

    d = {
        "version": dist.version,
        "location": str(dist.locate_file("")),
    }
    direct_url = dist.read_text("direct_url.json")
    if direct_url:
        try:
            data = json.loads(direct_url)
        except json.JSONDecodeError:
            data = {}
        if data.get("dir_info", {}).get("editable") and data.get("url"):
            parsed = urlparse(data["url"])
            if parsed.scheme == "file":
                path = url2pathname(parsed.path)
                if parsed.netloc:
                    path = f"//{parsed.netloc}{path}"
                elif os.name == "nt" and path.startswith("/") and len(path) > 2 and path[2] == ":":
                    path = path[1:]
                d["editable project location"] = path
            else:
                d["editable project location"] = data["url"]
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

    packages = _lab_editable_packages()
    for pkg, want_prefix in packages:
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
            print(f"        pip install -e {dict(packages)[pkg]} "
                  f"--no-deps --no-cache-dir")
        print("")
        print("        See CLAUDE.md \"POLICY (2026-05-27): release 後の "
              "LAB editable 再確認\" for the full recovery procedure.")
    return n_drift


REMOTE_EDITABLE_VERIFY = r'''
import importlib.metadata as md
import json
import sys
from urllib.parse import urlparse
from urllib.request import url2pathname

EXPECT = __EXPECT__

def norm(p):
    p = (p or "").replace("\\", "/").rstrip("/").lower()
    p = p.replace("//192.168.11.100/work/00_cae/radia/01_github",
                  "s:/radia/01_github")
    return p

def editable_location(dist):
    direct_url = dist.read_text("direct_url.json")
    if not direct_url:
        return None
    try:
        data = json.loads(direct_url)
    except json.JSONDecodeError:
        return None
    if not data.get("dir_info", {}).get("editable"):
        return None
    url = data.get("url")
    if not url:
        return None
    parsed = urlparse(url)
    if parsed.scheme == "file":
        path = url2pathname(parsed.path)
        if parsed.netloc:
            path = f"//{parsed.netloc}{path}"
        elif sys.platform.startswith("win") and path.startswith("/") and len(path) > 2 and path[2] == ":":
            path = path[1:]
        return path
    return url

bad = 0
for pkg, want in EXPECT:
    try:
        dist = md.distribution(pkg)
    except md.PackageNotFoundError:
        print(f"FAIL {pkg}: not installed; expected editable @ {want}")
        bad += 1
        continue
    got = editable_location(dist)
    if not got:
        print(f"FAIL {pkg}: installed v{dist.version}, but not editable")
        bad += 1
        continue
    if norm(got) != norm(want):
        print(f"FAIL {pkg}: editable drift")
        print(f"  got:      {got}")
        print(f"  expected: {want}")
        bad += 1
    else:
        print(f"OK   {pkg}: v{dist.version} editable -> {got}")

sys.exit(1 if bad else 0)
'''


def _verify_remote_editable(ssh_host, label, expected):
    """Check remote editable installs by reading package direct_url.json."""
    step(f"{label} editable verify")
    script = REMOTE_EDITABLE_VERIFY.replace("__EXPECT__", repr(expected))
    p = subprocess.run(["ssh", ssh_host, "python", "-"], input=script,
                       capture_output=True, text=True)
    if p.stdout:
        print(p.stdout.rstrip())
    if p.stderr:
        print(p.stderr.rstrip())
    if p.returncode == 0:
        ok(f"{label} editable packages point at NAS source")
        return 0
    fail(f"{label} editable package drift detected")
    return 1


def _verify_100_editable():
    return _verify_remote_editable(
        SSH_100, "100号機", _remote_100_editable_packages()
    )


def cmd_verify_editable(args):
    """Standalone editable verifier (no preflight, no phase9).

    Use this between deploys, or any time pip operations have run on
    LAB/100号機 and you want to confirm editable-loop integrity is intact.
    """
    drift = _verify_lab_editable()
    drift += _verify_100_editable()
    return 1 if drift else 0


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
        "User-Agent": "radia-release_qud",
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


def _run_retired_standalone_pyside_guard():
    """Guard that retired non-Cubit PySide panels have not been reintroduced."""
    step("retired standalone PySide panel guard")
    retired = [
        REPO / "src/radia/radia_gui_base.py",
        REPO / "src/radia/radia_ih.py",
        REPO / "src/radia/radia_em.py",
        REPO / "src/radia/radia_pcb.py",
        REPO / "src/radia/radia_motor.py",
        REPO / "src/radia/radia_streamfunction.py",
        REPO / "src/radia/_heat_panel.py",
    ]
    present = [str(p.relative_to(REPO)) for p in retired if p.exists()]
    pyproject = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    forbidden_tokens = [
        "PySide6>=6.5",
        "radia-ih = \"radia.radia_ih:main\"",
        "radia-em = \"radia.radia_em:main\"",
        "radia-pcb = \"radia.radia_pcb:main\"",
        "radia-streamfunction = \"radia.radia_streamfunction:main\"",
    ]
    present += [f"pyproject token: {t}" for t in forbidden_tokens if t in pyproject]
    if present:
        for item in present:
            fail(f"retired standalone PySide panel still present: {item}")
        return 1
    ok("retired standalone PySide panel files and entry points are absent")
    return 0


def cmd_done(args):
    """Definition-of-done check: preflight + editable verify + phase9 + retired panel guard.

    Read-only. Exit 0 means the release is consistent across LAB / 100号機 /
    mdx / hibino, the repo is release-ready, the editable tier is intact, AND
    the retired non-Cubit PySide panel surface has not been reintroduced.
    """
    step("Definition-of-done check "
         "(preflight + editable tier + phase9 + retired standalone panel guard)")
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed — repo state not release-ready.")
        return rc

    drift = _verify_lab_editable()
    drift += _verify_100_editable()
    if drift > 0:
        fail(f"{drift} editable-tier check(s) drifted.  "
             "Run the printed recovery commands, then re-run "
             "`release_qud done`.")
        return 1

    rc = cmd_phase9(args)
    if rc != 0:
        fail("phase9 drift detected — at least one machine is out of sync.")
        return rc

    rc = _run_retired_standalone_pyside_guard()
    if rc != 0:
        fail("retired standalone PySide panel surface reappeared. Remove it, then re-run "
             "`release_qud done`.")
        return rc

    if getattr(args, "simulink_package", None):
        rc = _verify_simulink_candidate_state(args.simulink_package)
        if rc != 0:
            fail("Simulink candidate did not satisfy the four-machine gate.")
            return rc

    print("")
    suffix = (" The supplied Simulink candidate also passed all four MATLAB "
              "machines." if getattr(args, "simulink_package", None) else "")
    ok("DEFINITION OF DONE met. Release is consistent across LAB / 100号機 / "
       "mdx / hibino, the editable tier is intact, and the retired standalone "
       "PySide panel surface is absent." + suffix)
    return 0


# ============================================================
# CLI
# ============================================================

def main():
    p = argparse.ArgumentParser(prog="release_qud",
                                 description="Enforce the release-qud flow.")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("preflight",
                    help="read-only state report (always safe)")
    sub.add_parser("phase0",
                    help="clean rebuild of the Cubit plugin")
    s8 = sub.add_parser("phase8",
                         help="deploy + verify + smoke on LAB / 100号機 / hibino")
    s8.add_argument("--target", default="lab,100",
                     help="comma list: lab, 100, hibino, all (default lab,100)")
    sub.add_parser("phase8e",
                    help="upgrade mdx from PyPI (after PyPI propagation)")
    sub.add_parser("phase9",
                    help="cross-machine consistency probe")
    ss = sub.add_parser(
        "simulink-candidate",
        help="verify one extracted Simulink package on four MATLAB machines")
    ss.add_argument("--package", required=True,
                    help="path to an IH preview or full Radia Simulink ZIP")
    ss.add_argument("--target", default="all",
                    help="comma list: lab, 100, mdx, hibino, all")
    sub.add_parser("all",
                    help="phase8 -> phase8e -> phase9 in one shot")
    sub.add_parser("verify-editable",
                    help="LAB/100号機 editable-install pointers check (read-only)")
    sub.add_parser("ci-verify",
                    help="Phase 5.5: gh-free CI-green gate (run after push main, before tag)")
    done = sub.add_parser(
        "done",
        help="definition-of-done: preflight + editable-tier + phase9 + guards")
    done.add_argument(
        "--simulink-package",
        help="also require a matching four-machine Simulink candidate pass")

    args = p.parse_args()
    handler = {
        "preflight":        cmd_preflight,
        "phase0":           cmd_phase0,
        "phase8":           cmd_phase8,
        "phase8e":          cmd_phase8e,
        "phase9":           cmd_phase9,
        "simulink-candidate": cmd_simulink_candidate,
        "all":              cmd_all,
        "verify-editable":  cmd_verify_editable,
        "ci-verify":        cmd_ci_verify,
        "done":             cmd_done,
    }[args.cmd]
    raise SystemExit(handler(args))


if __name__ == "__main__":
    main()
