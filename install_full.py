#!/usr/bin/env python3
"""
Radia Full Installation

Install or update the complete Radia environment:

    python install_full.py                    # from PyPI
    python install_full.py --wheel path.whl   # from local wheel
    python install_full.py --all-users        # all user profiles (admin)

Steps:
  1. pip install radia  -- from PyPI or local wheel
  2. radia-setup        -- deploy Cubit plugin + panels

Requirements:
    Python 3.12+ on Windows (x64)
"""

import glob
import os
import subprocess
import sys


def main():
    all_users = "--all-users" in sys.argv

    # Find --wheel argument
    wheel_path = None
    for i, arg in enumerate(sys.argv):
        if arg == "--wheel" and i + 1 < len(sys.argv):
            wheel_path = sys.argv[i + 1]

    # Auto-detect wheel next to this script
    if not wheel_path:
        here = os.path.dirname(os.path.abspath(__file__))
        candidates = sorted(glob.glob(os.path.join(here, "dist", "radia-*.whl")),
                            reverse=True)
        if not candidates:
            candidates = sorted(glob.glob(os.path.join(here, "radia-*.whl")),
                                reverse=True)
        if candidates:
            wheel_path = candidates[0]

    print("=" * 60)
    print("  Radia Full Installation")
    print("=" * 60)
    print(f"  Python: {sys.version.split()[0]}")
    if wheel_path:
        print(f"  Wheel:  {wheel_path}")
    else:
        print("  Source: PyPI")
    print()

    # Step 1: pip install
    if wheel_path:
        print(f"[1/1] pip install --force-reinstall {os.path.basename(wheel_path)}")
        r = subprocess.run([sys.executable, "-m", "pip", "install",
                            "--force-reinstall", wheel_path])
    else:
        print("[1/1] pip install --upgrade radia")
        r = subprocess.run([sys.executable, "-m", "pip", "install",
                            "--upgrade", "radia"])
    if r.returncode != 0:
        print("ERROR: pip install failed")
        sys.exit(1)
    print()

    # Step 2 (repo clones only): install the developer git hooks so that
    # pushes to main run tools/ci_preflight.py via the pre-push gate.  A pure
    # PyPI install has no clone (no .git), so this is skipped there.
    here = os.path.dirname(os.path.abspath(__file__))
    hook_installer = os.path.join(here, "tools", "install_git_hooks.py")
    if os.path.exists(os.path.join(here, ".git")) and os.path.isfile(hook_installer):
        print("[hooks] install pre-push CI-preflight gate")
        subprocess.run([sys.executable, hook_installer])   # best-effort
        print()

    print("Verify:")
    print('  python -c "import radia; print(radia.__version__)"')


if __name__ == "__main__":
    main()
