"""
Cubit plugin installer for cubit-mesh-export.

Deploys the Cubit plugin binaries (.ccm, .ccl, .pyd) and Netgen DLLs
to the Coreform Cubit installation directory.

After `pip install cubit-mesh-export`, run:
    cubit-plugin-install               # install plugin
    cubit-plugin-install --all-users   # install for all user profiles

This is the SINGLE entry point for Cubit plugin deployment.
Radia-NGSolve panels are installed separately via `radia-setup`.
"""

import glob
import os
import shutil
import sys
from pathlib import Path


def _find_cubit_dir():
    """Find Coreform Cubit installation directory.

    Search order:
      1. CUBIT_PATH environment variable
      2. Common install locations (newest version first)
    """
    cubit_path = os.environ.get("CUBIT_PATH")
    if cubit_path and os.path.isdir(os.path.join(cubit_path, "plugins")):
        return Path(cubit_path).parent  # CUBIT_PATH points to bin/
    if cubit_path and os.path.isdir(os.path.join(cubit_path, "bin", "plugins")):
        return Path(cubit_path)

    if sys.platform == "win32":
        for base in [os.environ.get("ProgramFiles", "")]:
            if not base:
                continue
            candidates = sorted(glob.glob(os.path.join(base, "Coreform Cubit *")),
                                reverse=True)
            for c in candidates:
                if os.path.isdir(os.path.join(c, "bin", "plugins")):
                    return Path(c)
    return None


def _package_dir():
    """Directory containing this package's bundled binaries."""
    return Path(__file__).resolve().parent


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


def _clean_old_plugins(cubit_dir):
    """Remove ALL radia/cubit-mesh-export plugin files from Cubit."""
    plugins_dir = cubit_dir / "bin" / "plugins"
    bin_dir = cubit_dir / "bin"
    removed = 0

    for pattern in ["radia_cubit.ccm", "radia_cubit_mesh*.pyd",
                    "nglib.dll", "ngcore.dll"]:
        for f in plugins_dir.glob(pattern):
            try:
                f.unlink()
                removed += 1
            except OSError:
                pass

    for name in ["radia_cubit.ccl", "cubit_radia.bat"]:
        f = bin_dir / name
        if f.is_file():
            try:
                f.unlink()
                removed += 1
            except OSError:
                pass

    return removed


def install_plugin():
    """Deploy Cubit plugin binaries to Cubit installation.

    Copies .ccm, .ccl, .pyd, and Netgen DLLs to the Cubit directory.

    Returns:
        True if successful, False if Cubit not found.
    """
    pkg_dir = _package_dir()
    cubit_dir = _find_cubit_dir()

    print("=" * 60)
    print("  cubit-mesh-export: Plugin Install")
    print("=" * 60)
    print()

    if not cubit_dir:
        print("  Cubit not found. Set CUBIT_PATH if installed elsewhere.")
        return False

    plugins_dir = cubit_dir / "bin" / "plugins"
    print(f"  Cubit:   {cubit_dir}")
    print(f"  Package: {pkg_dir}")
    print()

    # Clean old files
    n = _clean_old_plugins(cubit_dir)
    if n:
        print(f"  Cleaned {n} old plugin file(s)")
        print()

    # Copy .ccm (APREPRO export commands)
    ccm_src = pkg_dir / "radia_cubit.ccm"
    if ccm_src.is_file():
        dst = plugins_dir / "radia_cubit.ccm"
        shutil.copy2(ccm_src, dst)
        print(f"  [OK] radia_cubit.ccm -> {dst}")
    else:
        print(f"  [--] radia_cubit.ccm not found in {pkg_dir}")

    # Copy .ccl (Qt5 GUI -> bin/, NOT plugins/)
    ccl_src = pkg_dir / "radia_cubit.ccl"
    if ccl_src.is_file():
        dst = cubit_dir / "bin" / "radia_cubit.ccl"
        shutil.copy2(ccl_src, dst)
        print(f"  [OK] radia_cubit.ccl -> {dst}")

    # Copy .pyd (high-order mesh curving)
    pyd_src = pkg_dir / "radia_cubit_mesh.pyd"
    if pyd_src.is_file():
        dst = plugins_dir / "radia_cubit_mesh.cp312-win_amd64.pyd"
        shutil.copy2(pyd_src, dst)
        print(f"  [OK] radia_cubit_mesh.pyd -> {dst}")
    else:
        print(f"  [--] radia_cubit_mesh.pyd not found in {pkg_dir}")

    # Copy Netgen DLLs
    nglib, ngcore = _find_netgen_dlls()
    if nglib:
        for dll in (nglib, ngcore):
            dst = plugins_dir / dll.name
            shutil.copy2(dll, dst)
            print(f"  [OK] {dll.name} -> {dst}")
    else:
        print("  [--] Netgen DLLs not found (high-order curving disabled)")

    print()
    print("  Plugin installed. Restart Cubit to load.")
    print("=" * 60)
    return True


def main():
    """Console script entry point for cubit-plugin-install."""
    install_plugin()


if __name__ == "__main__":
    main()
