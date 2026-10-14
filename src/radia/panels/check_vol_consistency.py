"""
Check .vol mesh consistency against CAD reference values.

This module re-exports from cubit_mesh_export (canonical source).
Install: pip install cubit-mesh-export  (or pip install radia[cubit])

Usage:
    python check_vol_consistency.py model.vol
    python check_vol_consistency.py model.vol --strict-labels
    python check_vol_consistency.py model.vol --contract labels.json

Exit code: 0 = all checks pass, 1 = validation finding, 2 = input error.
"""

try:
    from cubit_mesh_export.check import (
        LABEL_CONTRACT_SCHEMA,
        REPORT_SCHEMA,
        check_consistency,
        check_label_contract,
        check_mesh_quality,
        print_table,
        main,
    )
except ImportError:
    raise ImportError(
        "cubit_mesh_export is required for mesh consistency checking.\n"
        "Install with: pip install cubit-mesh-export  (or pip install radia[cubit])"
    ) from None

__all__ = [
    "LABEL_CONTRACT_SCHEMA",
    "REPORT_SCHEMA",
    "check_consistency",
    "check_label_contract",
    "check_mesh_quality",
    "print_table",
    "main",
]

if __name__ == "__main__":
    raise SystemExit(main())
