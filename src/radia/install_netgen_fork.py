"""
Install ksugahar/netgen fork (CallbackGeometry + SetGeomInfo).

This replaces the standard netgen-mesher with the fork version that adds:
- CallbackGeometry: ACIS direct curving (no STEP/OCC needed)
- SetGeomInfo: high-order mesh curving for externally imported meshes

Required for: export_NGSolveCurvedMesh(), Cubit panels, mesh.Curve(order)

Usage:
    radia-install-netgen-fork          # from command line (after pip install radia)
    python -m radia.install_netgen_fork  # alternative

The fork is pinned to match the ngsolve version in radia's dependencies.
PR: https://github.com/NGSolve/netgen/pull/232
"""

import json
import platform
import subprocess
import sys
import urllib.request


NETGEN_FORK_REPO = "ksugahar/netgen"
NETGEN_FORK_TAG = "v6.2.2602.post1-setgeominfo"


def _get_python_tag():
    v = sys.version_info
    return f"cp{v.major}{v.minor}"


def _get_platform_tag():
    if sys.platform == "win32":
        if platform.machine().lower() in ("amd64", "x86_64"):
            return "win_amd64"
    elif sys.platform == "linux":
        return "manylinux"
    raise RuntimeError(f"Unsupported platform: {sys.platform} {platform.machine()}")


def _fetch_wheel_url():
    py_tag = _get_python_tag()
    plat_tag = _get_platform_tag()

    api_url = (f"https://api.github.com/repos/{NETGEN_FORK_REPO}"
               f"/releases/tags/{NETGEN_FORK_TAG}")
    req = urllib.request.Request(
        api_url, headers={"Accept": "application/vnd.github.v3+json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read())

    for asset in data.get("assets", []):
        name = asset["name"]
        if name.endswith(".whl") and py_tag in name and plat_tag in name:
            return asset["browser_download_url"], name

    available = [a["name"] for a in data.get("assets", [])
                 if a["name"].endswith(".whl")]
    raise RuntimeError(
        f"No matching wheel for {py_tag}-{plat_tag}.\n"
        f"Available: {available}\n"
        f"Check: https://github.com/{NETGEN_FORK_REPO}/releases/tag/{NETGEN_FORK_TAG}"
    )


def main():
    print(f"Installing ksugahar/netgen fork ({NETGEN_FORK_TAG})...")
    print(f"  Python: {_get_python_tag()}, Platform: {_get_platform_tag()}")

    # Check current state
    try:
        from netgen.meshing import Mesh
        has_fork = hasattr(Mesh, 'SetGeomInfo')
        if has_fork:
            print("  Netgen fork is already installed (SetGeomInfo found).")
            print("  Use --force to reinstall.")
            if "--force" not in sys.argv:
                return
    except ImportError:
        pass

    url, name = _fetch_wheel_url()
    print(f"  Wheel: {name}")
    print(f"  URL: {url}")

    cmd = [sys.executable, "-m", "pip", "install", url, "--force-reinstall"]
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print("  ERROR: pip install failed.")
        sys.exit(1)

    print()
    print("  Netgen fork installed successfully.")
    print("  Verify: python -c \"from netgen.meshing import Mesh; print(hasattr(Mesh, 'SetGeomInfo'))\"")


if __name__ == "__main__":
    main()
