"""Launch .vol/.msh via python -B -m radia_mcp.mesh_viewer; see --help for post views."""
from .launcher import main

if __name__ == "__main__":
    raise SystemExit(main())
