#!/usr/bin/env python3
"""
Radia Full Installation

Install or update the complete Radia environment:

    python install_full.py

Steps:
  1. pip install --upgrade radia  -- from PyPI (NGSolve 6.2.2602, MKL, MCP servers)
  2. Install ksugahar/netgen fork -- CallbackGeometry + SetGeomInfo (if not already installed)
  3. Install Cubit panels         -- if Coreform Cubit is detected

Requirements:
    Python 3.12+ on Windows (x64)
"""

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


def _has_netgen_fork():
    """Check if netgen fork is already installed."""
    try:
        r = subprocess.run(
            [sys.executable, "-c",
             "from netgen.meshing import Mesh; print(hasattr(Mesh, 'SetGeomInfo'))"],
            capture_output=True, text=True, timeout=30)
        return "True" in r.stdout
    except Exception:
        return False


def main():
    print("=" * 60)
    print("  Radia Full Installation")
    print("=" * 60)
    print(f"  Python: {sys.version}")
    print(f"  Platform: {sys.platform} {platform.machine()}")
    print()

    # Step 1: radia from PyPI (always upgrade)
    print("[1/3] Installing/updating radia from PyPI...")
    _pip("--upgrade", "radia")
    print()

    # Step 2: netgen fork (skip if already installed)
    print("[2/3] Netgen fork...")
    if _has_netgen_fork():
        print("  Already installed (SetGeomInfo found). Skipping.")
    else:
        print("  Installing ksugahar/netgen fork...")
        try:
            url, name = _fetch_netgen_wheel_url()
            print(f"  Wheel: {name}")
            _pip(url, "--force-reinstall")
        except Exception as e:
            print(f"  WARNING: {e}")
    print()

    # Step 3: Cubit panels (always update)
    print("[3/3] Installing Cubit panels...")
    try:
        r = subprocess.run(
            [sys.executable, "-m", "radia.install_panels", "--all-users"],
            timeout=30)
        if r.returncode != 0:
            print("  Skipped (Cubit not found)")
    except Exception as e:
        print(f"  Skipped ({e})")
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
