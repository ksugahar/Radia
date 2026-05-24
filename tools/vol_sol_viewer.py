"""vol_sol_viewer.py -- standalone shim.

The canonical implementation lives in `radia.tools.vol_sol_viewer`
and is shipped to PyPI as `radia-vol-viewer` / `radia-vol-viewer-gui`
console scripts.  This top-level `tools/vol_sol_viewer.py` is kept
only so existing Windows file associations that point at the script
path (`pythonw.exe tools/vol_sol_viewer.py %1`) keep working.

For new installs prefer:

    radia-vol-viewer --register

which writes the association to `radia-vol-viewer-gui.exe` (no
console window).  See `src/radia/tools/vol_sol_viewer.py` and the
[project.gui-scripts] section of pyproject.toml.
"""

import sys

if __name__ == "__main__":
    try:
        from radia.tools.vol_sol_viewer import main
    except ImportError as e:
        print(f"Error importing radia.tools.vol_sol_viewer: {e}", file=sys.stderr)
        print("Install radia first: pip install -e <Radia repo root>", file=sys.stderr)
        sys.exit(1)
    sys.exit(main())
