# Verify curl(A) = B

Verification script for the Maxwell relation B = curl(A) using Radia and NGSolve.

## Purpose

This script verifies that:
1. Vector potential A is correctly computed by Radia for ObjHexahedron permanent magnets
2. The Maxwell relation B = curl(A) holds with proper unit handling
3. The radia_ngsolve integration correctly handles coordinate unit conversion

## Unit Handling in radia_ngsolve

The `radia_ngsolve.cpp` implementation handles unit conversion as follows:

### Coordinate Conversion

| Setting | NGSolve Input | Radia Call | Notes |
|---------|---------------|------------|-------|
| `units='m'` (default) | meters | meters (coord_scale_=1.0) | Standard SI units |
| `units='mm'` | millimeters | millimeters (coord_scale_=1000.0) | For mm-based meshes |

### Field Values

All field types (B, H, A, M, phi) are returned **without additional scaling**:

| Field | Units | Notes |
|-------|-------|-------|
| B | Tesla | Magnetic flux density |
| H | A/m | Magnetic field strength |
| A | T*m | Vector potential (consistent with `FldUnits` setting) |
| M | A/m | Magnetization |
| phi | A | Magnetic scalar potential |

### Why No A-field Scaling is Needed

The implementation in `radia_ngsolve.cpp` (lines 483-490) explicitly states:

```cpp
// Vector potential A: No additional scaling needed
// Radia returns A in T*m when FldUnits('m') is set, or T*mm when FldUnits('mm')
// The numerical value is the same, but units match the FldUnits setting
// Since we use coord_scale_ to convert coords to Radia's unit system,
// the returned A is already in the correct units (T*m for NGSolve)
```

For FMM-accelerated computation (lines 719-759), the dipole formula handles units consistently:
- Dipole moment m is in A*m^2 (SI units)
- Positions are converted to match the coordinate system
- The formula `A = (mu0/4pi) * (m x r) / |r|^3` gives A in T*m when using SI units

## Running the Test

```bash
cd examples/NGSolve_Integration/verify_curl_A_equals_B
python verify_curl_A_equals_B.py
```

## Output Files

- `verify_curl_A_B.vtu` - VTK file with A, curl(A), and B fields
- `verify_curl_A_B_error.vtu` - VTK file with |curl(A) - B| error field

## Workflow

1. Create hexahedral permanent magnet using ObjHexahedron
2. Create NGSolve mesh in air region outside magnet
3. Project A onto HCurl space using RadiaField
4. Compute curl(A) using NGSolve curl() operator
5. Project B onto HDiv space using RadiaField
6. Compare |curl(A)| with |B| at test points
7. Verify ratio is approximately 1.0

## Expected Results

With the current implementation, the |curl(A)|/|B| ratio should be close to 1.0, confirming that the Maxwell relation B = curl(A) is satisfied.

| Metric | Expected | Notes |
|--------|----------|-------|
| \|curl(A)\| / \|B\| ratio | ~1.0 | Maxwell relation satisfied |
| Ratio variation | < 10-20% | Due to FE discretization error |

## Key Implementation Details

1. **Coordinate conversion** - `coord_scale_` handles m <-> mm conversion for input coordinates
2. **FldUnits setting** - Radia's `FldUnits('m')` ensures consistent field output units
3. **No explicit A scaling** - The implementation uses `scale = 1.0` for all field types
4. **FMM path** - Dipole approximation uses SI units consistently (m, A*m^2, T*m)

---

**Last Updated**: 2026-01-09
