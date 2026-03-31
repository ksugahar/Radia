#!/usr/bin/env python3
"""
Radia Full Installation

Install or update the complete Radia environment:

    python install_full.py                # current user
    python install_full.py --all-users    # all existing user profiles (admin)

Steps:
  1. pip install --upgrade radia  -- from PyPI (NGSolve, MKL, MCP servers)
  2. Install Cubit panels         -- if Coreform Cubit is detected
  3. Install Cubit plugin (.ccm)  -- copy plugin + Netgen DLLs to Cubit plugins dir

High-order mesh curving is handled by the Cubit plugin (C++ ACIS kernel),
not by netgen fork. Standard upstream netgen-mesher is sufficient.

Requirements:
    Python 3.12+ on Windows (x64)
"""

import platform
import shutil
import subprocess
import sys
from pathlib import Path


# Cubit plugin
CUBIT_PLUGIN_NAME = "radia_cubit.ccm"
CUBIT_INSTALL_DIRS = [
    Path(r"C:\Program Files\Coreform Cubit 2025.3"),
    Path(r"C:\Program Files\Coreform Cubit 2024.8"),
]


def _pip(*args):
    cmd = [sys.executable, "-m", "pip", "install"] + list(args)
    print(f"  $ {' '.join(cmd)}")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(f"pip install failed: {' '.join(args)}")


def _find_cubit_dir():
    """Find Cubit installation directory."""
    for d in CUBIT_INSTALL_DIRS:
        if (d / "bin" / "plugins").is_dir():
            return d
    return None


def _find_plugin_ccm():
    """Find radia_cubit.ccm: next to this script, in dist/, or in build dir."""
    here = Path(__file__).resolve().parent
    candidates = [
        here / CUBIT_PLUGIN_NAME,
        here / "dist" / CUBIT_PLUGIN_NAME,
        here / "build-cubit-plugin" / "Release" / CUBIT_PLUGIN_NAME,
        here / "src" / "cubit_plugin" / "build-test" / CUBIT_PLUGIN_NAME,
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


def _find_netgen_dlls():
    """Find nglib.dll and ngcore.dll from pip-installed netgen package."""
    try:
        r = subprocess.run(
            [sys.executable, "-c",
             "import netgen; print(netgen.__file__)"],
            capture_output=True, text=True, timeout=10)
        if r.returncode == 0 and r.stdout.strip():
            ng_dir = Path(r.stdout.strip()).parent
            nglib = ng_dir / "nglib.dll"
            ngcore = ng_dir / "ngcore.dll"
            if nglib.is_file() and ngcore.is_file():
                return nglib, ngcore
    except Exception:
        pass
    return None, None


def main():
    all_users = "--all-users" in sys.argv

    print("=" * 60)
    print("  Radia Full Installation")
    print("=" * 60)
    print(f"  Python: {sys.version}")
    print(f"  Platform: {sys.platform} {platform.machine()}")
    if all_users:
        print("  Mode: --all-users (all existing user profiles)")
    print()

    # Step 1: radia from PyPI (always upgrade)
    print("[1/3] Installing/updating radia from PyPI...")
    _pip("--upgrade", "radia")
    print()

    # Step 2: Cubit panels (always update)
    print("[2/3] Installing Cubit panels...")
    try:
        cmd = [sys.executable, "-m", "radia.install_panels"]
        if all_users:
            cmd.append("--all-users")
        r = subprocess.run(cmd, timeout=30)
        if r.returncode != 0:
            print("  Skipped (Cubit not found)")
    except Exception as e:
        print(f"  Skipped ({e})")
    print()

    # Step 3: Cubit plugin (.ccm + Netgen DLLs for high-order curving)
    print("[3/3] Cubit plugin (radia_cubit.ccm + Netgen DLLs)...")
    cubit_dir = _find_cubit_dir()
    if cubit_dir:
        plugins_dir = cubit_dir / "bin" / "plugins"
        ccm = _find_plugin_ccm()
        if ccm:
            dst = plugins_dir / CUBIT_PLUGIN_NAME
            shutil.copy2(ccm, dst)
            print(f"  Copied {ccm.name} -> {dst}")
        else:
            print(f"  WARNING: {CUBIT_PLUGIN_NAME} not found (build it first)")

        nglib, ngcore = _find_netgen_dlls()
        if nglib:
            for dll in (nglib, ngcore):
                dst = plugins_dir / dll.name
                shutil.copy2(dll, dst)
                print(f"  Copied {dll.name} -> {dst}")
        else:
            print("  WARNING: Netgen DLLs not found (high-order curving disabled)")
    else:
        print("  Skipped (Cubit not found)")
    print()

    # Summary
    print("=" * 60)
    print("  Done!")
    print("=" * 60)
    print()
    print("  Verify:")
    print('    python -c "import radia; print(radia.__version__)"')
    print()


if __name__ == "__main__":
    main()
