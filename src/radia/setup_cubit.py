"""
Radia-NGSolve setup: installs Cubit plugin + Radia-NGSolve panels.

After `pip install radia`, run:
    radia-setup                    # current user
    radia-setup --all-users        # all user profiles (admin)

This command:
  1. Calls cubit-plugin-install (from cubit-mesh-export package)
  2. Installs Radia-NGSolve toolbar panels in ~/.cubit startup
"""

import sys


def setup_cubit(all_users=False):
    """Install Cubit plugin + Radia-NGSolve panels."""

    # 1. Install Cubit plugin (from cubit-mesh-export)
    try:
        from cubit_mesh_export.install import install_plugin
        install_plugin()
    except ImportError:
        print("cubit-mesh-export not installed. Run: pip install cubit-mesh-export")
        print("Skipping plugin installation.")
    except Exception as e:
        print(f"Plugin install failed: {e}")

    print()

    # 2. Install Radia-NGSolve panels
    print("  Installing Radia-NGSolve panels...")
    try:
        from radia.install_panels import install_panels
        install_panels(all_users=all_users)
    except Exception as e:
        print(f"  [--] Panel install failed: {e}")

    print()
    print("  Setup complete! Restart Cubit to load.")


def main():
    """Console script entry point for radia-setup."""
    all_users = "--all-users" in sys.argv
    setup_cubit(all_users=all_users)


if __name__ == "__main__":
    main()
