"""
Cubit integration setup: install plugin + panels in one command.

After `pip install radia`, run:
    radia-setup                    # current user
    radia-setup --all-users        # all user profiles (admin)

This command:
  1. Copies radia_cubit.ccm to Cubit plugins/ (APREPRO export commands)
  2. Copies radia_cubit.ccl to Cubit bin/ (Qt5 GUI Export Mesh menu)
  3. Copies radia_cubit_mesh.pyd to Cubit plugins/ (high-order mesh curving)
  4. Copies Netgen DLLs (nglib.dll, ngcore.dll) to Cubit plugins/
  5. Registers toolbar panels in ~/.cubit startup
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


def _find_radia_package_dir():
    """Find the radia package directory (where .pyd and .ccm are installed)."""
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


def _clean_old_plugins(plugins_dir):
    """Remove old radia plugin files from Cubit plugins/ before installing new ones."""
    old_patterns = [
        "radia_cubit.ccm",
        "radia_cubit_mesh*.pyd",
        "nglib.dll",
        "ngcore.dll",
    ]
    removed = 0
    for pattern in old_patterns:
        for f in plugins_dir.glob(pattern):
            try:
                f.unlink()
                removed += 1
            except OSError:
                pass
    return removed


def _clean_pip_artifacts():
    """Remove pip's stale temp directories (~ prefixed) in site-packages."""
    site_packages = Path(sys.prefix) / "Lib" / "site-packages"
    if not site_packages.is_dir():
        return 0
    removed = 0
    for entry in site_packages.iterdir():
        if entry.is_dir() and entry.name.startswith("~"):
            try:
                shutil.rmtree(entry)
                removed += 1
            except OSError:
                pass
    return removed


def _clean_panels_json():
    """Remove stale entries from radia_panels.json settings."""
    settings_path = Path.home() / ".cubit" / "radia_panels.json"
    if not settings_path.is_file():
        return
    try:
        import json
        with open(settings_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            return
        # Remove keys for dialogs that no longer exist
        # (all Qt dialogs removed; standalone app uses radia_app.json)
        valid_dialogs = set()
        stale = [k for k in data if k not in valid_dialogs]
        if stale:
            for k in stale:
                del data[k]
            with open(settings_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
    except Exception:
        pass


def setup_cubit(all_users=False):
    """Install Cubit plugin + panels.

    Cleans old installations first, then installs fresh.

    Args:
        all_users: Install panels for all user profiles.

    Returns:
        True if successful.
    """
    print("=" * 60)
    print("  Radia - Cubit Setup")
    print("=" * 60)
    print()

    radia_dir = _find_radia_package_dir()
    cubit_dir = _find_cubit_dir()

    # Clean pip artifacts
    n = _clean_pip_artifacts()
    if n:
        print(f"  Cleaned {n} stale pip artifact(s)")

    # Clean panels settings
    _clean_panels_json()

    if not cubit_dir:
        print("  Cubit not found. Skipping plugin installation.")
        print("  Set CUBIT_PATH environment variable if installed elsewhere.")
        print()
        _install_panels(all_users)
        return True

    plugins_dir = cubit_dir / "bin" / "plugins"
    print(f"  Cubit:  {cubit_dir}")
    print(f"  Radia:  {radia_dir}")
    print()

    # 0. Clean old plugin files
    n = _clean_old_plugins(plugins_dir)
    if n:
        print(f"  Cleaned {n} old plugin file(s)")
        print()

    # 1. Copy .ccm
    ccm_src = radia_dir / "radia_cubit.ccm"
    if ccm_src.is_file():
        dst = plugins_dir / "radia_cubit.ccm"
        shutil.copy2(ccm_src, dst)
        print(f"  [OK] radia_cubit.ccm -> {dst}")
    else:
        print(f"  [--] radia_cubit.ccm not found in {radia_dir}")

    # 2. Copy .ccl (GUI component -> bin/, NOT plugins/)
    ccl_src = radia_dir / "radia_cubit.ccl"
    if ccl_src.is_file():
        dst = cubit_dir / "bin" / "radia_cubit.ccl"
        shutil.copy2(ccl_src, dst)
        print(f"  [OK] radia_cubit.ccl -> {dst}")

    # 3. Copy .pyd
    pyd_src = radia_dir / "radia_cubit_mesh.pyd"
    if pyd_src.is_file():
        dst = plugins_dir / "radia_cubit_mesh.cp312-win_amd64.pyd"
        shutil.copy2(pyd_src, dst)
        print(f"  [OK] radia_cubit_mesh.pyd -> {dst}")
    else:
        print(f"  [--] radia_cubit_mesh.pyd not found in {radia_dir}")

    # 4. Copy Netgen DLLs
    nglib, ngcore = _find_netgen_dlls()
    if nglib:
        for dll in (nglib, ngcore):
            dst = plugins_dir / dll.name
            shutil.copy2(dll, dst)
            print(f"  [OK] {dll.name} -> {dst}")
    else:
        print("  [--] Netgen DLLs not found (high-order curving disabled)")

    print()

    # 4. Install panels (install_panels handles old .cubit block removal)
    _install_panels(all_users)

    # Summary
    print()
    print("=" * 60)
    print("  Setup complete! Restart Cubit to load the plugin.")
    print("=" * 60)
    print()
    return True


def _install_panels(all_users):
    """Install Cubit toolbar panels."""
    print("  Installing Cubit panels...")
    try:
        from radia.install_panels import install_panels
        install_panels(all_users=all_users)
    except Exception as e:
        print(f"  [--] Panel install failed: {e}")


def main():
    """Console script entry point for radia-setup."""
    all_users = "--all-users" in sys.argv
    setup_cubit(all_users=all_users)


if __name__ == "__main__":
    main()
