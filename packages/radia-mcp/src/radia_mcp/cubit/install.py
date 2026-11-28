"""
Cubit plugin installer for radia-mcp.

Deploys the ``export vol`` command plugin (``cubit_to_vol_export.dll``,
built from ``packages/radia-mcp/native/cubit_to_vol_export`` with ``build.py``)
into the Coreform Cubit 2026.8 plugin directory, then registers the Netgen
Export toolbar and menu (:mod:`radia_mcp.cubit.gui`) through
:mod:`radia_mcp.cubit.toolbar_install`.

The radia-mcp wheel ships no native plugin.  The built plugin directory is
``--plugin-dir``, else ``$RADIA_CUBIT_PLUGIN_DIR`` (the directory the export
worker also uses), and must hold ``cubit_to_vol_export.dll``::

    python -m radia_mcp.cubit.install --plugin-dir <dir>              # current user
    python -m radia_mcp.cubit.install --plugin-dir <dir> --all-users  # every profile (admin)
    python -m radia_mcp.cubit.install --check-only                    # preflight, no writes
    python -m radia_mcp.cubit.install --plugin-dir <dir> --verify-only

The plugin links the fork ``nglib``/``ngcore`` it was built against
(``build.py --netgen-dir``).  Those DLLs are not copied: Cubit must be able to
load them for the GUI menu's in-process ``export vol``.  The headless export
route (:mod:`radia_mcp.cubit.export_runner`) loads them itself.

Safety policy:

* All errors are **loud**. ``OSError`` during unlink / copy exits non-zero
  with the offending path printed -- never silently swallowed.
* Cubit **must not be running** on the target machine before install. We
  detect live Cubit GUI processes and abort.
* Every destination file is **verified** after copy (exists, size and
  sha256 match the source, mtime newer than preflight start).
* ``--check-only`` runs the preflight (Cubit dir + lock + process check)
  without mutating anything.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
import time
from pathlib import Path

from radia_mcp.cubit.layout import (
    COMMAND_PLUGIN, SUPPORTED_CUBIT_RELEASE, command_plugin_path, install_root,
    is_supported_cubit_root, plugin_directory,
)


# Windows consoles may default to cp932; never crash the installer on print().
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
    """The supported Cubit install root, from the one radia-mcp Cubit locator."""
    from radia_mcp.cubit.toolbar_install import find_cubit_bin
    cubit_bin = find_cubit_bin()
    if not cubit_bin:
        return None
    root = install_root(Path(cubit_bin))
    if plugin_directory(root).is_dir() and is_supported_cubit_root(root):
        return root
    return None


def resolve_plugin_dir(plugin_dir: str | os.PathLike[str] | None = None) -> Path:
    """The built plugin directory: ``plugin_dir``, else ``$RADIA_CUBIT_PLUGIN_DIR``."""
    from radia_mcp.cubit.export_runner import resolve_plugin_dir as _resolve
    return _resolve(plugin_dir)


def _is_cubit_process_name(name: str) -> bool:
    """Return True iff *name* is a Cubit GUI binary that can hold the plugin open."""
    if not name:
        return False
    n = name.lower()
    if n.endswith(".exe"):
        n = n[:-4]
    return n in {"coreform_cubit", "cubitx", "clarox"}


def _running_cubit_processes():
    """Return a list of (pid, name) for running Cubit GUI processes."""
    hits = []
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
            # Best effort; the file-lock check below is the reliable guard.
            pass
    return hits


def _file_is_locked(path: Path) -> str:
    """Return '' if the file is writable (or missing), else a reason string."""
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
    return [command_plugin_path(cubit_dir)]


def verify_deployment(plugin_dir: Path, cubit_dir: Path, *, verbose: bool = True):
    """Confirm the deployed plugin matches the built one (size and sha256).

    Returns:
        (ok: bool, issues: list[str])
    """
    issues = []
    src = Path(plugin_dir) / COMMAND_PLUGIN
    dst = command_plugin_path(cubit_dir)
    if verbose:
        print("  Verification:")
    if not src.is_file():
        issues.append(f"built plugin missing: {src}")
    elif not dst.is_file():
        issues.append(f"missing deployed file: {dst}")
    elif src.stat().st_size != dst.stat().st_size:
        issues.append(f"size mismatch: {dst} = {dst.stat().st_size} bytes, "
                      f"built {src} = {src.stat().st_size} bytes")
    elif _sha256(src) != _sha256(dst):
        issues.append(f"sha256 mismatch: {dst} != {src}")
    elif verbose:
        print(f"    [OK] {dst.name}  ({dst.stat().st_size} bytes, sha256 match)")
    if verbose:
        for issue in issues:
            print(f"    [FAIL] {issue}")
        if not issues:
            print("  [OK] deployed plugin matches the built plugin.")
    return (not issues), issues


def preflight(cubit_dir: Path, *, verbose: bool = True):
    """Run the read-only safety checks: no Cubit running, no plugin file locked."""
    problems = []
    for pid, name in _running_cubit_processes():
        problems.append(
            f"Cubit process running: pid={pid} name={name!r}. "
            "Ask the owner to save and close Cubit before deploying.")
    for f in _critical_plugin_files(cubit_dir):
        reason = _file_is_locked(f)
        if reason:
            problems.append(f"File locked: {f} ({reason})")
    if verbose:
        print("  Preflight:")
        if problems:
            for p in problems:
                print(f"    [!!] {p}")
        else:
            print("    [OK] no Cubit processes, no locked plugin files")
    return (not problems), problems


# ============================================================
# Cleanup + copy
# ============================================================

def _clean_old_plugins(cubit_dir: Path):
    """Remove the deployed plugin.

    Returns:
        (removed: int, errors: list[str])
    """
    removed = 0
    errors = []
    for path in _critical_plugin_files(cubit_dir):
        if path.is_file():
            try:
                path.unlink()
                removed += 1
            except OSError as e:
                errors.append(f"Could not delete {path}: {e}")
    return removed, errors


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _copy_verified(src: Path, dst: Path, start_time: float):
    """Copy src -> dst; raise RuntimeError unless size, sha256 and mtime verify."""
    if not src.is_file():
        raise RuntimeError(f"Source file missing: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    if not dst.is_file():
        raise RuntimeError(f"Copy succeeded but destination missing: {dst}")
    if src.stat().st_size != dst.stat().st_size:
        raise RuntimeError(
            f"Size mismatch after copy: {dst} has {dst.stat().st_size} bytes, "
            f"expected {src.stat().st_size} (from {src})")
    if _sha256(src) != _sha256(dst):
        raise RuntimeError(
            f"SHA256 mismatch after copy: {dst} != {src}. "
            "Antivirus or filesystem interference suspected.")
    if dst.stat().st_mtime < start_time - 1:
        raise RuntimeError(
            f"Stale mtime on {dst}: {dst.stat().st_mtime} < preflight {start_time}. "
            "The file was not actually rewritten.")


# ============================================================
# Main install flow
# ============================================================

def install_plugin(*, all_users: bool = False, check_only: bool = False,
                   verify_only: bool = False,
                   plugin_dir: str | os.PathLike[str] | None = None):
    """Deploy the command plugin and register the toolbar.

    Args:
        all_users: also register the toolbar in every user's ``~/.cubit``
            (requires admin privileges on Windows).
        check_only: run the preflight and exit without writing anything.
        verify_only: verify the deployed plugin and toolbar without writing.
        plugin_dir: built plugin directory (default ``$RADIA_CUBIT_PLUGIN_DIR``);
            not needed for ``check_only``.

    Returns:
        True on full success, False when Cubit is not found.

    Exit codes (via ``main``):
        0 -- install ok / preflight clean
        1 -- Cubit not found
        2 -- preflight failed (Cubit running, file locked)
        3 -- clean or copy step failed
        4 -- post-install verification failed
        5 -- any other error (missing plugin directory included)
    """
    cubit_dir = _find_cubit_dir()
    start_time = time.time()

    print("=" * 60)
    print("  radia-mcp: Cubit plugin install")
    if check_only:
        print("  (--check-only: no writes will be performed)")
    elif verify_only:
        print("  (--verify-only: confirming current deployment, no writes)")
    print("=" * 60)
    print()

    if not cubit_dir:
        print(f"  [FAIL] Coreform Cubit {SUPPORTED_CUBIT_RELEASE} not found.")
        print("         Set CUBIT_BIN_DIR to the selected Cubit bin directory.")
        return False
    print(f"  Cubit:   {cubit_dir}")

    built = None if check_only else resolve_plugin_dir(plugin_dir)
    if built is not None:
        print(f"  Plugin:  {built / COMMAND_PLUGIN}")
    print()

    if verify_only:
        ok, issues = verify_deployment(built, cubit_dir, verbose=True)
        from radia_mcp.cubit.toolbar_install import verify_panel_installation
        print()
        panel_ok, panel_issues = verify_panel_installation(all_users=all_users, verbose=True)
        if not (ok and panel_ok):
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
        print("    2. Re-run python -m radia_mcp.cubit.install (or --check-only).")
        raise SystemExit(2)

    if check_only:
        print("  [OK] preflight clean -- install would succeed. Exiting.")
        return True

    removed, clean_errs = _clean_old_plugins(cubit_dir)
    if clean_errs:
        print("  [FAIL] clean-up step hit errors:")
        for e in clean_errs:
            print(f"    - {e}")
        raise SystemExit(3)
    if removed:
        print(f"  Cleaned {removed} old plugin file(s)")
        print()

    src, dst = built / COMMAND_PLUGIN, command_plugin_path(cubit_dir)
    try:
        _copy_verified(src, dst, start_time)
    except RuntimeError as e:
        print(f"  [FAIL] {src.name} -> {dst}: {e}")
        raise SystemExit(3)
    print(f"  [OK] {src.name} -> {dst}")

    print()
    ok_v, _ = verify_deployment(built, cubit_dir, verbose=True)
    if not ok_v:
        raise SystemExit(4)

    print()
    print("  Plugin installed and verified.")
    print("=" * 60)

    from radia_mcp.cubit.toolbar_install import install_panels
    print()
    if not install_panels(all_users=all_users):
        raise SystemExit(4)
    return True


def main(argv: list[str] | None = None):
    """``python -m radia_mcp.cubit.install`` entry point."""
    parser = argparse.ArgumentParser(
        prog="python -m radia_mcp.cubit.install", allow_abbrev=False,
        description="Install the radia-mcp Cubit command plugin "
                    f"({COMMAND_PLUGIN}) and the Vol Export toolbar.")
    parser.add_argument("--plugin-dir", default=None,
                        help="directory holding the built plugin (build.py --out); "
                             "default $RADIA_CUBIT_PLUGIN_DIR")
    parser.add_argument("--all-users", action="store_true",
                        help="also register the toolbar in every user's ~/.cubit (requires admin)")
    parser.add_argument("--check-only", action="store_true",
                        help="run preflight only; do not modify anything")
    parser.add_argument("--verify-only", action="store_true",
                        help="check that the deployed plugin matches the built one "
                             "(size + sha256) and the toolbar is registered; no writes")
    args = parser.parse_args(argv)
    if args.check_only and args.verify_only:
        parser.error("--check-only and --verify-only are mutually exclusive")

    try:
        ok = install_plugin(all_users=args.all_users, check_only=args.check_only,
                            verify_only=args.verify_only, plugin_dir=args.plugin_dir)
    except SystemExit:
        raise
    except Exception as e:
        print(f"  [FATAL] {type(e).__name__}: {e}")
        raise SystemExit(5)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
