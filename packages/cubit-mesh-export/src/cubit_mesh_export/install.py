"""
Cubit plugin installer for cubit-mesh-export.

Deploys the Cubit plugin binaries (.ccm, .pyd) and Netgen DLLs to the
Coreform Cubit 2025.12+ installation directory.

After ``pip install cubit-mesh-export``, run::

    cubit-plugin-install                # install plugin (current user Cubit profile)
    cubit-plugin-install --all-users    # install for all user profiles (admin)
    cubit-plugin-install --check-only   # preflight only, no writes

This is the SINGLE entry point for Cubit plugin and Radia panel
deployment. If the ``radia`` package is installed, this command also
registers the PySide6 Radia-NGSolve toolbar and verifies that Cubit's
startup files point at it.

Safety policy (2026-04-14 -- post-incident hardening):

* All errors are **loud**. ``OSError`` during unlink / copy exits non-zero
  with the offending path printed -- never silently swallowed.
* Cubit **must not be running** on the target machine before install. We
  detect live ``cubit*`` / ``coreform*`` processes and abort.
* Every destination file is **verified** after copy (exists, size matches
  source, mtime newer than preflight start). A mismatch aborts.
* ``--check-only`` runs the preflight (Cubit dir + lock + process check)
  without mutating anything. Use this to gate deploys on shared machines.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import importlib.metadata as importlib_metadata
import os
import shutil
import sys
import time
from pathlib import Path

MIN_CUBIT_VERSION = (2025, 12)
MIN_CUBIT_VERSION_TEXT = "2025.12"
_REQUIRED_PACKAGE_BINARIES = (
    "cubit_mesh_export.ccm",
    "cubit_mesh_curver.pyd",
)


# Windows console defaults to cp932 for ja-JP locales; our user-facing
# messages use en-dashes and em-dashes. Force utf-8 so print() never
# crashes the installer with UnicodeEncodeError (2026-04-14 hardening).
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


# ============================================================
# Preflight helpers
# ============================================================

def _find_cubit_dir():
    """Find a supported Coreform Cubit installation directory.

    Search order:
      1. CUBIT_PATH environment variable
      2. Common install locations (newest supported version first)
    """
    def _cubit_version_key(p):
        name = os.path.basename(os.path.normpath(p))
        if name.lower() == "bin":
            name = os.path.basename(os.path.dirname(os.path.normpath(p)))
        v = name.replace("Coreform Cubit ", "", 1).strip()
        try:
            return tuple(int(x) for x in v.split("."))
        except (ValueError, AttributeError):
            return (0,)

    def _supported(root: Path) -> bool:
        return (
            (root / "bin" / "plugins").is_dir()
            and _cubit_version_key(str(root)) >= MIN_CUBIT_VERSION
        )

    cubit_path = os.environ.get("CUBIT_PATH")
    if cubit_path:
        p = Path(cubit_path)
        candidates = [p.parent if (p / "plugins").is_dir() else p]
        if (p / "bin" / "plugins").is_dir():
            candidates.insert(0, p)
        for candidate in candidates:
            if _supported(candidate):
                return candidate
        return None

    if sys.platform == "win32":
        for base in [os.environ.get("ProgramFiles", "")]:
            if not base:
                continue
            candidates = sorted(
                glob.glob(os.path.join(base, "Coreform Cubit *")),
                key=_cubit_version_key,
                reverse=True,
            )
            for c in candidates:
                root = Path(c)
                if _supported(root):
                    return root
    return None


def _package_dir():
    """Directory containing this package's bundled binaries."""
    return Path(__file__).resolve().parent


def _parse_version(v: str):
    """Parse 'X.Y.Z' -> (X, Y, Z) for comparison. Extras after '+' dropped."""
    core = v.split("+", 1)[0].split("-", 1)[0]
    parts = core.split(".")
    out = []
    for p in parts[:3]:
        try:
            out.append(int(p))
        except ValueError:
            out.append(0)
    while len(out) < 3:
        out.append(0)
    return tuple(out)


def _check_radia_compat():
    """If radia is importable, verify version compatibility.

    Compares radia.__version__ against COMPAT_RADIA_MIN/MAX declared in
    this package's __init__, and reciprocally checks
    radia.COMPAT_CUBIT_MESH_EXPORT_MIN/MAX.

    Returns:
        (ok: bool, message: str)

    ok=True when radia is not installed (nothing to check) or versions
    match. ok=False with a descriptive message on mismatch.
    """
    try:
        import cubit_mesh_export as _cme
        import radia as _rad
    except ImportError:
        return True, "radia not installed in this Python env -- skipping"

    cme_ver = _parse_version(getattr(_cme, "__version__", "0.0.0"))
    rad_ver = _parse_version(getattr(_rad, "__version__", "0.0.0"))

    cme_min = _parse_version(getattr(_cme, "COMPAT_RADIA_MIN", "0.0.0"))
    cme_max = _parse_version(getattr(_cme, "COMPAT_RADIA_MAX", "999.999.999"))
    rad_min = _parse_version(getattr(_rad, "COMPAT_CUBIT_MESH_EXPORT_MIN", "0.0.0"))
    rad_max = _parse_version(getattr(_rad, "COMPAT_CUBIT_MESH_EXPORT_MAX", "999.999.999"))

    problems = []
    if not (cme_min <= rad_ver <= cme_max):
        problems.append(
            f"cubit-mesh-export {_cme.__version__} expects radia in "
            f"[{_cme.COMPAT_RADIA_MIN}, {_cme.COMPAT_RADIA_MAX}], got {_rad.__version__}")
    if not (rad_min <= cme_ver <= rad_max):
        problems.append(
            f"radia {_rad.__version__} expects cubit-mesh-export in "
            f"[{_rad.COMPAT_CUBIT_MESH_EXPORT_MIN}, {_rad.COMPAT_CUBIT_MESH_EXPORT_MAX}], "
            f"got {_cme.__version__}")
    if problems:
        return False, "; ".join(problems)
    return True, (f"radia {_rad.__version__} <-> cubit-mesh-export "
                   f"{_cme.__version__} compatible")


def _find_netgen_dlls():
    """Find nglib.dll and ngcore.dll from pip-installed netgen."""
    try:
        import netgen
        ng_dir = Path(netgen.__file__).parent
        nglib = ng_dir / "nglib.dll"
        ngcore = ng_dir / "ngcore.dll"
        if nglib.is_file() and ngcore.is_file():
            return nglib, ngcore
    except ImportError:
        pass
    return None, None


def _radia_distribution_installed() -> bool:
    try:
        importlib_metadata.version("radia")
        return True
    except importlib_metadata.PackageNotFoundError:
        return False


def _is_cubit_process_name(name: str) -> bool:
    """Return True iff *name* is an actual Cubit GUI binary (holds plugin FDs).

    Accept:
      - coreform_cubit.exe / coreform_cubit            (Cubit 2024+)
      - cubit.exe / cubit                              (legacy)

    Reject (these contain 'cubit' but do not hold the plugin):
      - cubit-plugin-install.exe (this installer itself)
      - mcp-server-cubit.exe (Radia MCP server)
      - any python.exe running a cubit helper
    """
    if not name:
        return False
    n = name.lower()
    # Strip the .exe suffix so the compare is uniform.
    if n.endswith(".exe"):
        n = n[:-4]
    return n in ("coreform_cubit", "cubit")


def _running_cubit_processes():
    """Return a list of (pid, name) for running Cubit GUI processes.

    A non-empty list means the plugin files are likely locked -- installing
    over them is not safe. The 2026-04-14 incident on 100号機 was caused by
    a running Cubit holding cubit_mesh_export.ccm open, which made copy silently
    produce a half-updated plugin set.
    """
    hits = []
    try:
        import psutil  # optional; we fall back to os-native on Windows
        for p in psutil.process_iter(["pid", "name"]):
            n = p.info.get("name") or ""
            if _is_cubit_process_name(n):
                hits.append((p.info["pid"], n))
        return hits
    except ImportError:
        pass

    if sys.platform == "win32":
        try:
            import subprocess
            out = subprocess.run(
                ["tasklist", "/FO", "CSV", "/NH"],
                capture_output=True, text=True, timeout=10,
                encoding="cp932", errors="replace")
            for line in out.stdout.splitlines():
                parts = [p.strip('"') for p in line.split('","')]
                if len(parts) < 2:
                    continue
                name = parts[0].strip('"')
                if _is_cubit_process_name(name):
                    try:
                        pid = int(parts[1])
                    except ValueError:
                        pid = 0
                    hits.append((pid, name))
        except Exception:
            # Best-effort -- if we cannot detect, fall through to the file-lock
            # check which is more reliable anyway.
            pass
    return hits


def _file_is_locked(path: Path) -> str:
    """Return '' if the file is writable (or missing), else a reason string.

    On Windows, if a process has the file mapped (e.g. a running Cubit with
    cubit_mesh_export.ccm loaded), attempting to open it with exclusive
    read-write will raise PermissionError. We use that as the lock probe.
    """
    if not path.exists():
        return ""
    try:
        fd = os.open(str(path), os.O_RDWR)
        os.close(fd)
        return ""
    except PermissionError as e:
        return f"PermissionError: {e}"
    except OSError as e:
        return f"OSError: {e}"


def _critical_plugin_files(cubit_dir: Path):
    """Files whose state gates a safe install."""
    plugins = cubit_dir / "bin" / "plugins"
    bin_ = cubit_dir / "bin"
    return [
        plugins / "cubit_mesh_export.ccm",
        plugins / "cubit_mesh_curver.cp312-win_amd64.pyd",
        plugins / "nglib.dll",
        plugins / "ngcore.dll",
        bin_ / "cubit_mesh_export.ccl",      # historical stale Qt5 location
        plugins / "cubit_mesh_export.ccl",  # historical stale Qt5 location
    ]


def _expected_deployments(pkg_dir: Path, cubit_dir: Path):
    """Return list[(pkg_src, deployed_dst)] for each binary cubit-plugin-install
    is responsible for."""
    plugins_dir = cubit_dir / "bin" / "plugins"
    bin_dir = cubit_dir / "bin"
    pairs = []
    if (pkg_dir / "cubit_mesh_export.ccm").is_file():
        pairs.append((pkg_dir / "cubit_mesh_export.ccm",
                       plugins_dir / "cubit_mesh_export.ccm"))
    if (pkg_dir / "cubit_mesh_curver.pyd").is_file():
        pairs.append((pkg_dir / "cubit_mesh_curver.pyd",
                       plugins_dir / "cubit_mesh_curver.cp312-win_amd64.pyd"))
    nglib, ngcore = _find_netgen_dlls()
    if nglib:
        pairs.append((nglib, plugins_dir / nglib.name))
        pairs.append((ngcore, plugins_dir / ngcore.name))
    # Cubit-side Python helpers (Kelvin transformation, etc.).  Deployed
    # to bin/plugins/cubit_helpers/ so the .ccm code can locate them via
    # the plugin directory and add them to the embedded Python sys.path
    # before invoking `play "<dir>/auto_kelvin_entry.py"`.
    helpers_src = pkg_dir / "cubit_helpers"
    if helpers_src.is_dir():
        helpers_dst = plugins_dir / "cubit_helpers"
        for py in sorted(helpers_src.glob("*.py")):
            pairs.append((py, helpers_dst / py.name))
    return pairs


def verify_deployment(pkg_dir: Path, cubit_dir: Path, *, verbose: bool = True):
    """Confirm that every expected file is deployed and matches the package.

    Returns:
        (ok: bool, issues: list[str])

    Checks:
      - every expected destination file exists
      - destination file size + sha256 match the package source
      - no stale ``cubit_mesh_export.ccl`` lingering from the retired Qt5
        plugin path

    Independent from ``preflight()`` (which is read-only safety-before-write).
    Use this after ``cubit-plugin-install`` finishes, OR at any later
    time to re-confirm the install state.
    """
    issues = []

    for name in _REQUIRED_PACKAGE_BINARIES:
        p = pkg_dir / name
        if not p.is_file():
            issues.append(
                f"required package binary missing: {p}. Rebuild before deploy."
            )

    nglib, ngcore = _find_netgen_dlls()
    if not nglib or not ngcore:
        issues.append(
            "required Netgen DLLs missing from the Python environment "
            "(nglib.dll / ngcore.dll)"
        )

    pairs = _expected_deployments(pkg_dir, cubit_dir)
    if verbose:
        print("  Verification:")

    for src, dst in pairs:
        if not dst.is_file():
            msg = f"missing deployed file: {dst}"
            issues.append(msg)
            if verbose: print(f"    [MISS] {msg}")
            continue
        src_size = src.stat().st_size
        dst_size = dst.stat().st_size
        if src_size != dst_size:
            msg = (f"size mismatch: {dst} = {dst_size} bytes, "
                   f"package {src} = {src_size} bytes")
            issues.append(msg)
            if verbose: print(f"    [SIZE] {msg}")
            continue
        if _sha256(src) != _sha256(dst):
            msg = f"sha256 mismatch: {dst} != {src}"
            issues.append(msg)
            if verbose: print(f"    [HASH] {msg}")
            continue
        if verbose:
            print(f"    [OK] {dst.name}  ({dst_size} bytes, sha256 match)")

    # Retired Qt5 .ccl guard.  Since radia 4.80.0 Cubit UI is PySide6-only;
    # any remaining .ccl can load stale menus or confuse first installs.
    for stale_ccl in [
        cubit_dir / "bin" / "cubit_mesh_export.ccl",
        cubit_dir / "bin" / "plugins" / "cubit_mesh_export.ccl",
        cubit_dir / "bin" / "radia_cubit.ccl",
        cubit_dir / "bin" / "plugins" / "radia_cubit.ccl",
    ]:
        if stale_ccl.is_file():
            msg = (f"stale retired Qt5 .ccl remains: {stale_ccl}. "
                   "Re-run cubit-plugin-install to clean it.")
            issues.append(msg)
            if verbose: print(f"    [STALE] {msg}")

    # Compat report (read-only).
    compat_ok, compat_msg = _check_radia_compat()
    if verbose:
        print(f"    compat: {compat_msg}")
    if not compat_ok:
        issues.append(f"Version incompatibility: {compat_msg}")

    if verbose:
        if issues:
            print(f"  [FAIL] {len(issues)} issue(s) found.")
        else:
            print("  [OK] every expected binary present and matches package source.")

    return (not issues), issues


def preflight(cubit_dir: Path, *, verbose: bool = True):
    """Run the read-only safety checks.

    Returns:
        (ok: bool, problems: list[str])

    ok=True means: Cubit is not running, no plugin file is locked.
    """
    problems = []

    procs = _running_cubit_processes()
    if procs:
        for pid, name in procs:
            problems.append(
                f"Cubit process running: pid={pid} name={name!r}. "
                "Ask the owner to save and close Cubit before deploying.")

    for f in _critical_plugin_files(cubit_dir):
        reason = _file_is_locked(f)
        if reason:
            problems.append(f"File locked: {f} ({reason})")

    compat_ok, compat_msg = _check_radia_compat()
    if not compat_ok:
        problems.append(f"Version incompatibility: {compat_msg}")

    if verbose:
        print("  Preflight:")
        print(f"    compat: {compat_msg}")
        if problems:
            for p in problems:
                print(f"    [!!] {p}")
        else:
            print("    [OK] no Cubit processes, no locked plugin files")

    return (not problems), problems


# ============================================================
# Cleanup + copy
# ============================================================

_CLEAN_PATTERNS = [
    "cubit_mesh_export.ccm",
    "cubit_mesh_export.ccl",
    "cubit_mesh_curver*.pyd",
    "nglib.dll",
    "ngcore.dll",
    # Tier-2 rename cleanup (2026-06-01): also remove the OLD radia_cubit.*
    # plugin from machines deployed before the cubit_mesh_export rename, so
    # Cubit does not load BOTH the old and new .ccm and double-register the
    # `export` APREPRO commands.  Safe to leave indefinitely (a no-op
    # once every machine is migrated).
    "radia_cubit.ccm",
    "radia_cubit.ccl",
    "radia_cubit_mesh*.pyd",
]


def _clean_old_plugins(cubit_dir: Path):
    """Remove all radia / cubit-mesh-export plugin files from Cubit.

    Returns:
        (removed: int, errors: list[str])

    Errors are collected rather than raised so the caller can surface
    every locked/missing-privilege case in one go, then abort.
    """
    plugins_dir = cubit_dir / "bin" / "plugins"
    bin_dir = cubit_dir / "bin"
    removed = 0
    errors = []

    # Clean from BOTH plugins/ and bin/ -- .ccl has historically landed in
    # either location. A stale copy in the wrong place silently shadows a
    # fresh one (2026-04-14 incident).
    for pattern in _CLEAN_PATTERNS:
        for dir_ in (plugins_dir, bin_dir):
            for f in dir_.glob(pattern):
                try:
                    f.unlink()
                    removed += 1
                except OSError as e:
                    errors.append(f"Could not delete {f}: {e}")

    # Clean Python helpers (cubit_helpers/*.py) so a removed helper
    # doesn't linger.  Also drop __pycache__ from the previous install.
    helpers_dst = plugins_dir / "cubit_helpers"
    if helpers_dst.is_dir():
        for f in helpers_dst.glob("*.py"):
            try:
                f.unlink()
                removed += 1
            except OSError as e:
                errors.append(f"Could not delete {f}: {e}")
        pyc_dir = helpers_dst / "__pycache__"
        if pyc_dir.is_dir():
            for f in pyc_dir.glob("*"):
                try:
                    f.unlink()
                    removed += 1
                except OSError as e:
                    errors.append(f"Could not delete {f}: {e}")
            try:
                pyc_dir.rmdir()
            except OSError:
                pass  # __pycache__ may have new entries; harmless

    for name in ["cubit_radia.bat"]:
        f = bin_dir / name
        if f.is_file():
            try:
                f.unlink()
                removed += 1
            except OSError as e:
                errors.append(f"Could not delete {f}: {e}")

    return removed, errors


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _copy_verified(src: Path, dst: Path, start_time: float):
    """Copy src -> dst and verify that the result is correct.

    Post-conditions:
      - dst exists
      - dst size == src size
      - dst sha256 == src sha256 (catches silent truncation / antivirus
        interference on Windows)
      - dst mtime >= start_time (catches the case where an AV reverted
        the file or a symlink pointed somewhere unexpected)

    Raises:
        RuntimeError on any post-condition failure.
    """
    if not src.is_file():
        raise RuntimeError(f"Source file missing: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)

    if not dst.is_file():
        raise RuntimeError(f"Copy succeeded but destination missing: {dst}")

    src_size = src.stat().st_size
    dst_size = dst.stat().st_size
    if src_size != dst_size:
        raise RuntimeError(
            f"Size mismatch after copy: {dst} has {dst_size} bytes, "
            f"expected {src_size} (from {src})")

    src_hash = _sha256(src)
    dst_hash = _sha256(dst)
    if src_hash != dst_hash:
        raise RuntimeError(
            f"SHA256 mismatch after copy: {dst} != {src}. "
            "Antivirus or filesystem interference suspected.")

    dst_mtime = dst.stat().st_mtime
    if dst_mtime < start_time - 1:
        raise RuntimeError(
            f"Stale mtime on {dst}: {dst_mtime} < preflight {start_time}. "
            "The file was not actually rewritten.")


# ============================================================
# Main install flow
# ============================================================

def install_plugin(*, all_users: bool = False, check_only: bool = False,
                    verify_only: bool = False):
    """Deploy Cubit plugin binaries to Cubit installation.

    Args:
        all_users: if True, also register the Radia panel on every user's
            ~/.cubit. Requires admin privileges on Windows.
        check_only: if True, run preflight (find Cubit, process / lock
            scan) and exit without writing anything. Returns True iff
            preflight is clean.

    Returns:
        True on full success. False (or raises RuntimeError) on failure.

    Exit codes (via ``main``):
        0 -- install ok / preflight clean
        1 -- Cubit not found
        2 -- preflight failed (Cubit running, file locked, missing privs)
        3 -- clean or copy step failed
        4 -- post-install verification failed
    """
    pkg_dir = _package_dir()
    cubit_dir = _find_cubit_dir()
    start_time = time.time()

    print("=" * 60)
    print("  cubit-mesh-export: Plugin Install")
    if check_only:
        print("  (--check-only: no writes will be performed)")
    elif verify_only:
        print("  (--verify-only: confirming current deployment, no writes)")
    print("=" * 60)
    print()

    if not cubit_dir:
        print(f"  [FAIL] Coreform Cubit {MIN_CUBIT_VERSION_TEXT}+ not found.")
        print("         Set CUBIT_PATH to the Cubit 2025.12 bin directory.")
        return False

    plugins_dir = cubit_dir / "bin" / "plugins"
    print(f"  Cubit:   {cubit_dir}")
    print(f"  Package: {pkg_dir}")
    print()

    if verify_only:
        ok, issues = verify_deployment(pkg_dir, cubit_dir, verbose=True)
        try:
            from radia.install_panels import verify_panel_installation
        except ImportError as e:
            verify_panel_installation = None
            if _radia_distribution_installed():
                issues = list(issues) + [
                    f"radia is installed but panel verifier could not import: {e}"
                ]
                ok = False
        if verify_panel_installation is not None:
            print()
            panel_ok, panel_issues = verify_panel_installation(
                all_users=all_users, verbose=True)
            if not panel_ok:
                issues = list(issues) + list(panel_issues)
                ok = False
        if not ok:
            raise SystemExit(4)
        return True

    ok, problems = preflight(cubit_dir, verbose=True)
    print()
    if not ok:
        print("  [FAIL] Preflight refused to proceed:")
        for p in problems:
            print(f"    - {p}")
        print()
        print("  RESOLUTION:")
        print("    1. Ask all users with Cubit open to save and close it.")
        print("    2. Re-run cubit-plugin-install (or --check-only to verify).")
        raise SystemExit(2)

    if check_only:
        print("  [OK] preflight clean -- install would succeed. Exiting.")
        return True

    # Clean old files -- abort if any unlink failed.
    removed, clean_errs = _clean_old_plugins(cubit_dir)
    if clean_errs:
        print("  [FAIL] clean-up step hit errors:")
        for e in clean_errs:
            print(f"    - {e}")
        print()
        print("  Likely cause: Cubit slipped back in while the preflight was")
        print("  running, or a file is held by AV / indexer. Close Cubit,")
        print("  pause AV on the plugins dir, and retry.")
        raise SystemExit(3)
    if removed:
        print(f"  Cleaned {removed} old plugin file(s)")
        print()

    # Copy artefacts with per-file verification.
    copy_jobs = []
    ccm_src = pkg_dir / "cubit_mesh_export.ccm"
    if ccm_src.is_file():
        copy_jobs.append((ccm_src, plugins_dir / "cubit_mesh_export.ccm"))
    else:
        print(f"  [FAIL] cubit_mesh_export.ccm not found in {pkg_dir}")

    pyd_src = pkg_dir / "cubit_mesh_curver.pyd"
    if pyd_src.is_file():
        copy_jobs.append(
            (pyd_src, plugins_dir / "cubit_mesh_curver.cp312-win_amd64.pyd"))
    else:
        print(f"  [FAIL] cubit_mesh_curver.pyd not found in {pkg_dir}")

    nglib, ngcore = _find_netgen_dlls()
    if nglib:
        copy_jobs.append((nglib, plugins_dir / nglib.name))
        copy_jobs.append((ngcore, plugins_dir / ngcore.name))
    else:
        print("  [FAIL] Netgen DLLs not found (nglib.dll / ngcore.dll)")

    # Cubit-side Python helpers (Kelvin transformation, etc.).
    helpers_src = pkg_dir / "cubit_helpers"
    if helpers_src.is_dir():
        helpers_dst = plugins_dir / "cubit_helpers"
        for py in sorted(helpers_src.glob("*.py")):
            copy_jobs.append((py, helpers_dst / py.name))
    else:
        print(f"  [FAIL] cubit_helpers/ not found in {pkg_dir}")

    copy_errors = []
    for required in _REQUIRED_PACKAGE_BINARIES:
        if not (pkg_dir / required).is_file():
            copy_errors.append(f"required package binary missing: {pkg_dir / required}")
    if not nglib or not ngcore:
        copy_errors.append("required Netgen DLLs missing")
    if not helpers_src.is_dir():
        copy_errors.append(f"required cubit_helpers/ missing: {helpers_src}")

    for src, dst in copy_jobs:
        try:
            _copy_verified(src, dst, start_time)
            print(f"  [OK] {src.name} -> {dst}")
        except RuntimeError as e:
            copy_errors.append(str(e))
            print(f"  [FAIL] {src.name} -> {dst}: {e}")

    if copy_errors:
        print()
        print("  [FAIL] one or more files did not deploy correctly:")
        for e in copy_errors:
            print(f"    - {e}")
        raise SystemExit(3)

    # Final post-install sanity -- delegated to verify_deployment for a
    # single source of truth (also exposed via --verify-only).
    print()
    ok_v, _ = verify_deployment(pkg_dir, cubit_dir, verbose=True)
    if not ok_v:
        raise SystemExit(4)

    print()
    print("  Plugin installed. All destinations verified.")
    print("=" * 60)

    # Install Radia-NGSolve panels if the radia package is available.
    try:
        from radia.install_panels import install_panels
        print()
        if not install_panels(all_users=all_users):
            raise SystemExit(4)
    except ImportError as e:
        if _radia_distribution_installed():
            print(f"  [FAIL] radia is installed but panel installer could not import: {e}")
            raise SystemExit(4)
        pass  # radia not installed, panels not needed

    return True


def main():
    """Console script entry point for cubit-plugin-install."""
    parser = argparse.ArgumentParser(
        prog="cubit-plugin-install",
        description="Install the Radia / cubit-mesh-export Cubit plugin "
                    "(.ccm / .pyd + Netgen DLLs) and Radia panels.")
    parser.add_argument("--all-users", action="store_true",
                        help="also register Radia toolbar in every "
                             "user's ~/.cubit (requires admin)")
    parser.add_argument("--check-only", action="store_true",
                        help="run preflight only; do not modify anything")
    parser.add_argument("--verify-only", action="store_true",
                        help="check that the currently-deployed binaries "
                             "match this package's source-of-truth (size + "
                             "sha256). No writes.")
    args = parser.parse_args()
    if args.check_only and args.verify_only:
        parser.error("--check-only and --verify-only are mutually exclusive")

    try:
        ok = install_plugin(all_users=args.all_users,
                            check_only=args.check_only,
                            verify_only=args.verify_only)
    except SystemExit:
        raise
    except Exception as e:
        print(f"  [FATAL] {type(e).__name__}: {e}")
        raise SystemExit(5)

    if not ok:
        # install_plugin already printed a diagnostic.
        raise SystemExit(1)


if __name__ == "__main__":
    main()
