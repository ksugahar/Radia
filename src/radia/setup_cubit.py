"""
Radia-NGSolve setup: installs Cubit plugin + Radia-NGSolve panels.

After installing ``radia[cubit,gui]``, run:
    cubit-plugin-install                    # current user
    cubit-plugin-install --all-users        # all user profiles (admin)

This command:
  1. Calls cubit-plugin-install (from cubit-mesh-export package)
  2. Installs Radia-NGSolve toolbar panels in Cubit's startup files
"""

import sys


def setup_cubit(all_users=False):
    """Install Cubit plugin + Radia-NGSolve panels."""

    # 1. Install Cubit plugin (from cubit-mesh-export)
    try:
        from cubit_mesh_export.install import install_plugin
        install_plugin(all_users=all_users)
    except ImportError:
        print("cubit-mesh-export not installed. Run: pip install cubit-mesh-export")
        raise SystemExit(1)

    print()
    print("  Setup complete. Restart Cubit 2025.12 to load.")


def main():
    """Legacy console script entry point; prefer cubit-plugin-install."""
    all_users = "--all-users" in sys.argv
    setup_cubit(all_users=all_users)


if __name__ == "__main__":
    main()
