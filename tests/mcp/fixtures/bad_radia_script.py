"""Intentionally bad Radia script for lint testing.

This script contains one violation per Radia-specific rule.
Used by --selftest when examples/ directory is not available.

Expected findings: 8 (one per Radia rule that can fire on a standalone file)
"""
import sys
import os
import radia as rad

# Rule: hardcoded-absolute-path (HIGH)
sys.path.insert(0, r'S:\Radia\src')

# Rule: removed-fldunits (HIGH)
rad.FldUnits("mm")

# Rule: removed-fldbatch (HIGH)
result = rad.FldBatch(mag, pts)

# Rule: removed-solver-api (HIGH)
rad.SetHACApKParams(1e-4, 10, 2.0)

# Rule: build-release-path (LOW)
sys.path.insert(0, 'build/Release')

# Rule: objbckg-needs-callable (CRITICAL)
bkg = rad.ObjBckg([0, 0, 0.1])

if __name__ == "__main__":
    pm = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 954930])
    B = rad.Fld(pm, "b", [0.02, 0, 0])
    print(B)
    # Rule: missing-utidelall (HIGH) -- no cleanup call
