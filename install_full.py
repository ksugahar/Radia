#!/usr/bin/env python3
"""
Radia Full Installation Script

Installs the complete Radia environment including:
  1. radia (from PyPI) - includes NGSolve, sparsesolv, MKL, MCP servers
  2. ksugahar/netgen fork (from GitHub Releases) - CallbackGeometry, SetGeomInfo

Usage:
    python install_full.py          # Install everything
    python install_full.py --skip-netgen-fork   # Skip netgen fork (basic mode)
    python install_full.py --dry-run            # Show what would be installed

After installation:
    pip install radia       <-- already done by this script
    import radia            <-- Radia magnetostatics
    from ngsolve.la import CompactAMSPreconditioner  <-- Compact HX preconditioner
    mcp-server-radia        <-- MCP server for Claude Code / Cursor

Requirements:
    Python 3.12+ on Windows (x64)
"""

import argparse
import json
import platform
import subprocess
import sys
import urllib.request
from pathlib import Path


NETGEN_FORK_REPO = "ksugahar/netgen"
NETGEN_FORK_TAG = "v6.2.2602.post1-setgeominfo"


def get_python_tag():
    """Get Python version tag for wheel matching (e.g., cp312)."""
    v = sys.version_info
    return f"cp{v.major}{v.minor}"


def get_platform_tag():
    """Get platform tag for wheel matching."""
    if sys.platform == "win32":
        if platform.machine().lower() in ("amd64", "x86_64"):
            return "win_amd64"
    elif sys.platform == "linux":
        return "manylinux"
    raise RuntimeError(f"Unsupported platform: {sys.platform} {platform.machine()}")


def fetch_netgen_wheel_url():
    """Find the matching netgen fork wheel URL from GitHub Releases."""
    py_tag = get_python_tag()
    plat_tag = get_platform_tag()

    api_url = f"https://api.github.com/repos/{NETGEN_FORK_REPO}/releases/tags/{NETGEN_FORK_TAG}"
    print(f"  Fetching release info from {NETGEN_FORK_REPO} tag={NETGEN_FORK_TAG}")

    req = urllib.request.Request(api_url, headers={"Accept": "application/vnd.github.v3+json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())

    for asset in data.get("assets", []):
        name = asset["name"]
        if name.endswith(".whl") and py_tag in name and plat_tag in name:
            return asset["browser_download_url"], name

    available = [a["name"] for a in data.get("assets", []) if a["name"].endswith(".whl")]
    raise RuntimeError(
        f"No matching wheel for {py_tag}-{plat_tag}.\n"
        f"Available: {available}\n"
        f"Check: https://github.com/{NETGEN_FORK_REPO}/releases/tag/{NETGEN_FORK_TAG}"
    )


def pip_install(*args):
    """Run pip install with given arguments."""
    cmd = [sys.executable, "-m", "pip", "install"] + list(args)
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        raise RuntimeError(f"pip install failed: {' '.join(args)}")


def main():
    parser = argparse.ArgumentParser(description="Install Radia with all dependencies")
    parser.add_argument("--skip-netgen-fork", action="store_true",
                        help="Skip netgen fork installation (basic mode)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be installed without actually installing")
    parser.add_argument("--upgrade", action="store_true",
                        help="Upgrade existing packages")
    args = parser.parse_args()

    print("=" * 60)
    print("  Radia Full Installation")
    print("=" * 60)
    print(f"  Python: {sys.version}")
    print(f"  Platform: {sys.platform} {platform.machine()}")
    print()

    # Step 1: Install radia from PyPI (pulls NGSolve, sparsesolv, MKL)
    print("[1/2] Installing radia from PyPI...")
    print("  (includes NGSolve, ngsolve-sparsesolv, MKL, MCP servers)")
    if args.dry_run:
        print("  [DRY RUN] pip install radia")
    else:
        pip_args = ["radia"]
        if args.upgrade:
            pip_args.insert(0, "--upgrade")
        pip_install(*pip_args)
    print()

    # Step 2: Install netgen fork (replaces standard netgen-mesher)
    if args.skip_netgen_fork:
        print("[2/2] Skipping netgen fork (--skip-netgen-fork)")
        print("  Note: export_NGSolveCurvedMesh() will not work without the fork.")
        print("  Basic Radia functionality (magnets, solvers, PEEC) works fine.")
    else:
        print("[2/2] Installing ksugahar/netgen fork...")
        print("  (CallbackGeometry + SetGeomInfo for high-order Cubit mesh curving)")
        try:
            url, name = fetch_netgen_wheel_url()
            print(f"  Wheel: {name}")
            if args.dry_run:
                print(f"  [DRY RUN] pip install {url} --force-reinstall")
            else:
                pip_install(url, "--force-reinstall")
        except Exception as e:
            print(f"  WARNING: Could not install netgen fork: {e}")
            print("  You can install it manually later:")
            print(f"    pip install <wheel> --force-reinstall")
            print(f"    from https://github.com/{NETGEN_FORK_REPO}/releases")
    print()

    # Summary
    print("=" * 60)
    print("  Installation complete!")
    print("=" * 60)
    print()
    print("  Verify:")
    print("    python -c \"import radia; print(radia.__version__)\"")
    print("    python -c \"from ngsolve.la import CompactAMSPreconditioner; print('OK')\"")
    print()
    print("  MCP servers (for Claude Code / Cursor):")
    print("    mcp-server-radia")
    print("    mcp-server-ngsolve")
    print("    mcp-server-cubit")
    print()


if __name__ == "__main__":
    main()
