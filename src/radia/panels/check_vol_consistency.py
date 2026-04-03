"""
Check .vol mesh consistency against CAD reference values.

This module re-exports from cubit_mesh_export (canonical source).
Install: pip install cubit-mesh-export  (or pip install radia[cubit])

Usage:
    python check_vol_consistency.py model.vol
    python check_vol_consistency.py model.vol --json model.vol.json
    python check_vol_consistency.py model.vol --threshold 0.5

Exit code: 0 = all checks pass, 1 = warnings found, 2 = error.
"""

try:
    from cubit_mesh_export.check import (
        check_consistency,
        print_table,
        main,
    )
except ImportError:
    raise ImportError(
        "cubit_mesh_export is required for mesh consistency checking.\n"
        "Install with: pip install cubit-mesh-export  (or pip install radia[cubit])"
    ) from None

__all__ = ["check_consistency", "print_table", "main"]

if __name__ == "__main__":
    main()
