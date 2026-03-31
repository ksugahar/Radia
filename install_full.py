#!/usr/bin/env python3
"""
Radia Full Installation

Install or update the complete Radia environment:

    python install_full.py                # current user
    python install_full.py --all-users    # all existing user profiles (admin)

Steps:
  1. pip install --upgrade radia  -- from PyPI (NGSolve 6.2.2602, MKL, MCP servers)
  2. pip install netgen-mesher    -- upstream (CallbackGeometry included since 6.2.2602)
  3. Install Cubit panels         -- if Coreform Cubit is detected
  4. Install Cubit plugin         -- copy .ccm + .pyd + Netgen DLLs to Cubit plugins dir

Requirements:
    Python 3.12+ on Windows (x64)
"""

import platform
import shutil
import subprocess
import sys
from pathlib import Path


NETGEN_VERSION = "6.2.2602"
NETGEN_FORK_REPO = "ksugahar/netgen"
NETGEN_FORK_TAG = "v6.2.2602.post1-setgeominfo"

# Cubit plugin
CUBIT_PLUGIN_NAME = "radia_cubit.ccm"
CUBIT_INSTALL_DIRS = [
    Path(r"C:\Program Files\Coreform Cubit 2025.3"),
    Path(r"C:\Program Files\Coreform Cubit 2024.8"),
]


def _has_netgen_fork():
    """Check if netgen fork (with SetGeomInfo) is already installed."""
    try:
        r = subprocess.run(
            [sys.executable, "-c",
             "from netgen.meshing import Mesh; print(hasattr(Mesh, 'SetGeomInfo'))"],
            capture_output=True, text=True, timeout=30)
        return "True" in r.stdout
    except Exception:
        return False


def _py_tag():
    v = sys.version_info
    return f"cp{v.major}{v.minor}"


def _plat_tag():
    if sys.platform == "win32" and platform.machine().lower() in ("amd64", "x86_64"):
        return "win_amd64"
    raise RuntimeError(f"Unsupported platform: {sys.platform} {platform.machine()}")


def _fetch_netgen_wheel_url():
    """Fetch netgen fork wheel URL from GitHub Releases."""
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


def _find_plugin_pyd():
    """Find radia_cubit_mesh .pyd: next to this script, in dist/, or in build dir."""
    here = Path(__file__).resolve().parent
    py, plat = _py_tag(), _plat_tag()
    pyd_name = f"radia_cubit_mesh.{py}-{plat}.pyd"
    candidates = [
        here / pyd_name,
        here / "dist" / pyd_name,
        here / "src" / "cubit_plugin" / "build-pyd" / pyd_name,
    ]
    for p in candidates:
        if p.is_file():
            return p
    return None


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
    print("[1/4] Installing/updating radia from PyPI...")
    _pip("--upgrade", "radia")
    print()

    # Step 2: netgen fork (CallbackGeometry + SetGeomInfo)
    print("[2/4] Netgen fork (ksugahar/netgen)...")
    if _has_netgen_fork():
        print("  Already installed (SetGeomInfo found). Skipping.")
    else:
        print("  Installing from GitHub Releases...")
        try:
            url, name = _fetch_netgen_wheel_url()
            print(f"  Wheel: {name}")
            _pip(url, "--force-reinstall")
        except Exception as e:
            print(f"  WARNING: Fork install failed: {e}")
            print(f"  Falling back to upstream netgen-mesher=={NETGEN_VERSION}")
            try:
                _pip("netgen-mesher==" + NETGEN_VERSION)
            except Exception:
                pass
    print()

    # Step 3: Cubit panels (always update)
    print("[3/4] Installing Cubit panels...")
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

    # Step 4: Cubit plugin (.ccm + .pyd + Netgen DLLs)
    print("[4/4] Cubit plugin (.ccm + .pyd + Netgen DLLs)...")
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

        pyd = _find_plugin_pyd()
        if pyd:
            dst = plugins_dir / pyd.name
            shutil.copy2(pyd, dst)
            print(f"  Copied {pyd.name} -> {dst}")
        else:
            print("  WARNING: radia_cubit_mesh .pyd not found (build it first)")

        nglib, ngcore = _find_netgen_dlls()
        if nglib:
            for dll in (nglib, ngcore):
                dst = plugins_dir / dll.name
                shutil.copy2(dll, dst)
                print(f"  Copied {dll.name} -> {dst}")
        else:
            print("  WARNING: Netgen DLLs not found (high-order disabled)")
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
