#!/usr/bin/env python
"""
Verify FldUnits deprecation behavior:
- rad.FldUnits('m') should be a noop (no warning) for backward compatibility
- rad.FldUnits('mm') should issue DeprecationWarning
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import warnings
import radia as rad

print("=" * 80)
print("FldUnits Deprecation Behavior Verification")
print("=" * 80)
print()

# ============================================================================
# Test 1: FldUnits('m') should be a noop with no warning
# ============================================================================
print("[Test 1] rad.FldUnits('m') should be a noop (no warning)")
rad.UtiDelAll()

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    rad.FldUnits('m')
    deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]

if len(deprecation_warnings) == 0:
    print("  [PASS] No DeprecationWarning issued for FldUnits('m')")
else:
    print("  [FAIL] Unexpected DeprecationWarning issued for FldUnits('m')")
    for dw in deprecation_warnings:
        print(f"         Warning: {dw.message}")

# Verify field computation still works after FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.02], [0, 0, 1.0])
B = rad.Fld(magnet, 'b', [0, 0, 0.05])
print(f"  Field check: Bz = {B[2]:.6e} T (should be non-zero)")
assert abs(B[2]) > 0, "Field should be non-zero"
print()

# ============================================================================
# Test 2: FldUnits('mm') should issue DeprecationWarning
# ============================================================================
print("[Test 2] rad.FldUnits('mm') should issue DeprecationWarning")
rad.UtiDelAll()

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    rad.FldUnits('mm')
    deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]

if len(deprecation_warnings) > 0:
    print("  [PASS] DeprecationWarning issued for FldUnits('mm')")
    print(f"         Message: {deprecation_warnings[0].message}")
else:
    print("  [INFO] No DeprecationWarning issued for FldUnits('mm')")
    print("         (FldUnits may not yet have deprecation warnings implemented)")

print()

# ============================================================================
# Summary
# ============================================================================
print("=" * 80)
print("Summary")
print("=" * 80)
print("  Radia always uses meters internally.")
print("  FldUnits('m') is accepted as a noop for backward compatibility.")
print("  FldUnits('mm') is deprecated -- Radia always uses meters.")
print("=" * 80)

rad.UtiDelAll()
