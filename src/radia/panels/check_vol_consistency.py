"""
Check .vol mesh consistency against CAD reference values.

This module re-exports from radia.vol_check (canonical source, no Cubit needed).

Usage:
    python check_vol_consistency.py model.vol
    python check_vol_consistency.py model.vol --strict-labels
    python check_vol_consistency.py model.vol --contract labels.json

Exit code: 0 = all checks pass, 1 = validation finding, 2 = input error.
"""

from radia.vol_check import (
    LABEL_CONTRACT_SCHEMA,
    REPORT_SCHEMA,
    check_consistency,
    check_label_contract,
    check_mesh_quality,
    print_table,
    main,
)

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
