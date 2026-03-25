#!/usr/bin/env python3
"""
Radia Full Installation

The single installation command for the complete Radia environment:

    python install_full.py

Steps:
  1. pip install radia (from PyPI) -- includes NGSolve 6.2.2602, MKL, MCP servers
  2. Install ksugahar/netgen fork  -- CallbackGeometry + SetGeomInfo (PR#232 pending)
  3. Install Cubit panels          -- if Coreform Cubit is detected

Requirements:
    Python 3.12+ on Windows (x64)
"""

import argparse
import json
import platform
import subprocess
import sys
import urllib.request


NETGEN_FORK_REPO = "ksugahar/netgen"
NETGEN_FORK_TAG = "v6.2.2602.post1-setgeominfo"


def _py_tag():
    v = sys.version_info
    return f"cp{v.major}{v.minor}"


def _plat_tag():
    if sys.platform == "win32" and platform.machine().lower() in ("amd64", "x86_64"):
        return "win_amd64"
    raise RuntimeError(f"Unsupported platform: {sys.platform} {platform.machine()}")


def _fetch_netgen_wheel_url():
    py, plat = _py_tag(), _plat_tag()
    api_url = (f"https://api.github.com/repos/{NETGEN_FORK_REPO}"
               f"/releases/tags/{NETGEN_FORK_TAG}")
    req = urllib.request.Request(
        api_url, headers={"Accept": "application/vnd.github.v3+json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())
    for asset in data.get("assets", []):
        name = asset["name"]
        if name.endswith(".whl") and py in name and plat in name:
            return asset["browser_download_url"], name
    available = [a["name"] for a in data.get("assets", [])
                 if a["name"].endswith(".whl")]
    raise RuntimeError(
        f"No matching wheel for {py}-{plat}.\n"
        f"Available: {available}\n"
        f"https://github.com/{NETGEN_FORK_REPO}/releases/tag/{NETGEN_FORK_TAG}")


def _pip(*args):
    cmd = [sys.executable, "-m", "pip", "install"] + list(args)
    print(f"  $ {' '.join(cmd)}")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        raise RuntimeError(f"pip install failed: {' '.join(args)}")


def main():
    parser = argparse.ArgumentParser(description="Radia full installation")
    parser.add_argument("--upgrade", action="store_true",
                        help="Upgrade radia to latest PyPI version")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be installed")
    args = parser.parse_args()

    print("=" * 60)
    print("  Radia Full Installation")
    print("=" * 60)
    print(f"  Python: {sys.version}")
    print(f"  Platform: {sys.platform} {platform.machine()}")
    print()

    # Step 1: radia from PyPI (ngsolve==6.2.2602.post2 pinned)
    print("[1/3] Installing radia from PyPI...")
    if args.dry_run:
        print("  [DRY RUN] pip install radia")
    else:
        pip_args = ["radia"]
        if args.upgrade:
            pip_args.insert(0, "--upgrade")
        _pip(*pip_args)
    print()

    # Step 2: netgen fork (replaces standard netgen-mesher)
    print("[2/3] Installing ksugahar/netgen fork...")
    print("  (CallbackGeometry + SetGeomInfo for Cubit mesh curving)")
    try:
        url, name = _fetch_netgen_wheel_url()
        print(f"  Wheel: {name}")
        if args.dry_run:
            print(f"  [DRY RUN] pip install {url} --force-reinstall")
        else:
            _pip(url, "--force-reinstall")
    except Exception as e:
        print(f"  WARNING: {e}")
    print()

    # Step 3: Cubit panels
    print("[3/3] Installing Cubit panels...")
    try:
        cmd = [sys.executable, "-m", "radia.install_panels"]
        print(f"  $ {' '.join(cmd)}")
        if args.dry_run:
            print("  [DRY RUN]")
        else:
            r = subprocess.run(cmd)
            if r.returncode != 0:
                print("  Cubit panels: skipped (Cubit not found or not configured)")
    except Exception as e:
        print(f"  Cubit panels: skipped ({e})")
    print()

    # Summary
    print("=" * 60)
    print("  Installation complete!")
    print("=" * 60)
    print()
    print("  Verify:")
    print('    python -c "import radia; print(radia.__version__)"')
    print('    python -c "from netgen.meshing import Mesh; print(hasattr(Mesh, \'SetGeomInfo\'))"')
    print()


if __name__ == "__main__":
    main()
