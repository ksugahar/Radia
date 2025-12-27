# Verify curl(A) = B

Verification script for the Maxwell relation B = curl(A) using Radia and NGSolve.

## Purpose

This script verifies that:
1. Vector potential A is correctly computed by Radia for ObjPolyhdr hexahedral permanent magnets
2. The Maxwell relation B = curl(A) holds when proper unit conversion is applied
3. The radia_ngsolve integration correctly handles A field scaling

## Unit Conversion (Critical)

**Radia ALWAYS uses mm internally**, regardless of `FldUnits()` setting:

| Field | Radia Internal | With `FldUnits('m')` | NGSolve Expected |
|-------|---------------|---------------------|------------------|
| B | Tesla | Tesla | Tesla |
| H | A/m | A/m | A/m |
| **A** | **T*mm** | **T*mm** (not auto-scaled!) | **T*m** |

### Why A Needs Scaling

1. Vector potential A has dimensions [T*length]
2. Radia computes A using mm-based geometry internally
3. NGSolve differentiates in meters: `curl(A) = dA/dx [m^-1]`
4. For B = curl(A) to hold: `A_SI [T*m] = A_radia [T*mm] / 1000`

### radia_ngsolve.cpp Fix

The fix in `src/radia/radia_ngsolve.cpp` applies automatic scaling:

```cpp
// Vector potential A unit scaling:
// Radia ALWAYS uses mm internally, so A is always in T*mm
// NGSolve differentiates in meters: curl(A) = dA/dx_m
// To get correct B = curl(A), we scale A by 0.001:
double scale = (field_type == "a") ? 0.001 : 1.0;
```

## Test Results

With proper A field scaling in radia_ngsolve:

| Metric | Expected | Actual |
|--------|----------|--------|
| `\|curl(A)\| / \|B\|` ratio | ~1.0 | ~1.0 |
| Ratio variation | < 10% | < 5% |

## Running the Test

```bash
cd examples/ngsolve_integration/verify_curl_A_equals_B
python verify_curl_A_equals_B.py
```

## Output Files

- `verify_curl_A_B.vtu` - VTK file with A, curl(A), and B fields
- `verify_curl_A_B_error.vtu` - VTK file with |curl(A) - B| error field

## Workflow

1. Create hexahedral permanent magnet using ObjPolyhdr
2. Create NGSolve mesh in air region outside magnet
3. Project A onto HCurl space using RadiaField
4. Compute curl(A) using NGSolve curl() operator
5. Project B onto HDiv space using RadiaField
6. Compare |curl(A)| with |B| at test points
7. Verify ratio is consistent (~1.0)

## Key Findings

1. **Radia always uses mm internally** - `FldUnits('m')` only scales coordinate input, not A output
2. **A field requires /1000 scaling** - to convert from T*mm to T*m for correct curl(A) = B
3. **B and H fields need no scaling** - they are dimensionally correct in all unit systems
4. **The fix is in radia_ngsolve.cpp** - automatic scaling applied when field_type == "a"

---

**Last Updated**: 2025-12-27
