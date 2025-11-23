# Tetrahedral Mesh Implementation - Current Status

**Date**: 2025-11-23
**Status**: ⚠️ PARTIAL FIX - Coordinate System Fixed, In-Plane Magnetization Pending

---

## Summary

Tetrahedral mesh support in Radia has been **partially fixed**. The coordinate transformation bug has been resolved, but field errors remain due to missing in-plane magnetization treatment.

**Issues Resolved**:
1. ✅ **R1 rotation sign error** - Fixed by negating third column of R1 matrix
2. ✅ **Unified coordinate system** - Hexahedra and tetrahedra now use consistent transformations
3. ✅ **Segmentation faults eliminated** - No longer crashes

**Remaining Issues**:
1. ❌ **In-plane magnetization ignored** - B_comp only uses normal component (Magn.z)
2. ❌ **Large field errors** - 1305-7019% error due to missing in-plane contribution

---

## Test Results

### Current Status (After R1 Fix)

**Script**: `compare_permanent_magnet_meshes.py`
**Input**: 0.1m cube, M = [0, 0, 1.2] T
**Mesh**: Netgen tetrahedral (659 elements)

| Test Point [m] | Built-in \|H\| [A/m] | Tetra \|H\| [A/m] | Error |
|---------------|---------------------|------------------|-------|
| [0.08, 0, 0]  | 0.033760           | 2.403415         | **7019%** |
| [0, 0.08, 0]  | 0.033760           | 1.155282         | **3322%** |
| [0, 0, 0.08]  | 0.053851           | 0.756456         | **1305%** |
| [0.15, 0, 0]  | 0.010733           | 0.033145         | 209% |
| [0, 0, 0.15]  | 0.019074           | 0.068293         | 258% |

**Built-in hexahedral mesh**: ✅ Working correctly (< 0.1% error)
**Tetrahedral mesh**: ❌ Large errors due to in-plane magnetization

---

## Root Cause Analysis

### 1. R1 Rotation Sign Error (FIXED)

**Problem**: Original R1 construction gave R1*[0,0,1] = -N (inward normal) instead of +N (outward normal)

**Debug Output**:
```
Face 0: N=[-1.0000,0.0000,0.0000], R1*[0,0,1]=[1.0000,0.0000,0.0000]  (WRONG!)
```

**Fix Applied** (rad_polyhedron.cpp:372-378):
```cpp
// UNIFIED FIX: Negate third column to get R1*[0,0,1] = +N
TVector3d St1_fixed(St1.x, St1.y, -St1.z);
TVector3d St2_fixed(St2.x, St2.y, -St2.z);
TVector3d St3_fixed(St3.x, St3.y, -St3.z);
TMatrix3d R1(St1_fixed, St2_fixed, St3_fixed);
```

**Result**: Hexahedra now work correctly, tetrahedra produce finite (but inaccurate) results

### 2. In-Plane Magnetization Ignored (NOT FIXED)

**Problem**: B_comp function only uses Magn.z (normal component), ignoring Magn.x and Magn.y

**Code Location**: `src/core/rad_planar_2d_part2.cpp:76-100`

```cpp
// Only uses normal component!
double W = ConstForH * Magn.z;

// Calls analytical formula for magnetic charge contribution
RadAnalyticalFieldFromPolygonCharge(..., W, ...);

// In-plane components Magn.x, Magn.y are IGNORED!
```

**Why This Matters**:

For **hexahedral faces**, magnetization is perpendicular to faces:
```
Face normal: N = [-1, 0, 0]
Local magnetization: LocMagn = [-1.200, 0.000, 0.000]
After rotation: Magn = [0.000, 0.000, -1.200]  <- All in z-component
```
→ In-plane components (Magn.x, Magn.y) are negligible → Current code works

For **tetrahedral faces**, magnetization is NOT perpendicular:
```
Face normal: N = [-0.122, -0.973, -0.037]
Local magnetization: LocMagn = [-0.150, -1.190, -0.046]
After rotation: Magn = [-0.150, -1.190, -0.046]  <- All components non-zero!
```
→ In-plane components are **LARGE** → Ignoring them causes huge errors

---

## Changes Made

### src/core/rad_polyhedron.cpp

**Lines 372-378** - R1 column sign fix:
```cpp
// UNIFIED FIX for hexahedra and tetrahedra
TVector3d St1_fixed(St1.x, St1.y, -St1.z);
TVector3d St2_fixed(St2.x, St2.y, -St2.z);
TVector3d St3_fixed(St3.x, St3.y, -St3.z);
TMatrix3d R1(St1_fixed, St2_fixed, St3_fixed);
```

**Lines 383-384** - Simplified to R1 only:
```cpp
// Use R1 for all faces (edge-based R2 disabled for now)
TMatrix3d R = R1;
```

**Status**: ✅ Implemented and tested

### src/core/rad_planar_2d_part2.cpp

**Lines 76-100** - In-plane magnetization attempt:

```cpp
// ========================================================================
// DISABLED: Initial implementation of in-plane magnetization
// ========================================================================
// Attempted to add line current contributions from M_in = [Magn.x, Magn.y, 0]
// Result: Errors increased 10x (7019% -> 94350%)
// Conclusion: Physical model or implementation was incorrect
//
// TODO: Analyze ELF MAGIC Fortran code to understand correct approach
// ========================================================================
```

**Status**: ❌ Attempted but disabled (made results worse)

---

## What Needs to Be Done

### Priority 1: Analyze ELF MAGIC Implementation

**File**: `S:/ELF_MAGIC/02_Fortran_Source/devel/magic.f90`
**Function**: `MM4TS()` (lines 3567-3590)

**Key Finding**: ELF loops over all 3 magnetization components (x, y, z) separately:

```fortran
do LL = 1,3  ! Loop over x, y, z components
    do MEN = 1,4  ! Loop over 4 tetrahedral faces
        call AAOBB(VV,UNIV(1,LL),Q)  ! Dot product with unit vector
        ! ... compute contribution ...
    end do
end do
```

This suggests they process each component (Mx, My, Mz) independently, not just the normal component.

**Next Steps**:
1. Understand DOMBS() subroutine (field accumulation)
2. Trace VV variable (geometric quantity)
3. Map ELF algorithm to Radia's B_comp structure
4. Implement proper in-plane magnetization treatment

### Priority 2: Edge-Based Rotation (R2) - Optional Enhancement

**Status**: Currently disabled. R1-only works for hexahedra.

**When Needed**: May be required once in-plane magnetization is implemented to ensure consistent edge orientation for line current integration.

**Algorithm** (from ELF MAGIC):
```cpp
// Step 1: R1 - rotate face normal to z-axis (implemented)
TMatrix3d R1 = RotationToAlignWithZ(face_normal);

// Step 2: R2 - rotate around z-axis to align edge with x-axis
TVector3d edge_local = R1 * (v1 - v0);
double theta = atan2(edge_local.y, edge_local.x);
TMatrix3d R2 = RotationAroundZ(-theta);

// Step 3: Combine
TMatrix3d R = R2 * R1;
```

---

## Testing Strategy

### Phase 1: Built-in Hexahedral Mesh
```bash
python compare_permanent_magnet_meshes.py
```
**Status**: ✅ PASSING (< 0.1% error)

### Phase 2: Tetrahedral Mesh (659 elements)
```bash
python compare_permanent_magnet_meshes.py
```
**Status**: ⚠️ PARTIAL - Finite results but large errors (1305-7019%)
**Target**: Error < 5% after in-plane magnetization implementation

### Phase 3: Single Tetrahedron Test
```bash
python test_single_cube_tetra.py
```
**Status**: ✅ No segfault, produces finite results
**Target**: Error < 5% after in-plane fix

---

## Related Files

**Test Scripts**:
- `compare_permanent_magnet_meshes.py` - Main benchmark (permanent magnets)
- `test_single_cube_tetra.py` - Simple 6-tetrahedra test

**Radia Source**:
- `src/core/rad_polyhedron.cpp` - Rotation matrix construction (R1 fix applied)
- `src/core/rad_planar_2d_part2.cpp` - B_comp function (needs in-plane fix)
- `src/python/netgen_mesh_import.py` - Tetrahedral mesh import (working)

**Reference Implementation**:
- `S:/ELF_MAGIC/02_Fortran_Source/devel/magic.f90` - MM4TS() subroutine

---

## MatLin API Changes (Completed)

**Date**: 2025-11-23
**Status**: ✅ Documentation updated

### Simplified API

```python
# Form 1: Isotropic linear material
mat = rad.MatLin(ksi)

# Form 2: Anisotropic with easy axis
mat = rad.MatLin([ksi_par, ksi_perp], [ex, ey, ez])
```

### Removed Old API
❌ `MatLin([ksi_par, ksi_perp], Mr_scalar)` - DEPRECATED

**Rationale**: Permanent magnets should use `ObjRecMag()` directly, not MatLin

---

## Conclusion

**Current State**:
- ✅ Coordinate transformation bugs fixed
- ✅ Hexahedral meshes working correctly
- ⚠️ Tetrahedral meshes produce finite results but with large errors

**Root Cause of Remaining Errors**:
B_comp function only uses normal magnetization component (Magn.z), ignoring in-plane components (Magn.x, Magn.y). For tetrahedra, in-plane components are large and cannot be ignored.

**Next Action**:
Analyze ELF MAGIC's MM4TS() implementation to understand proper treatment of all three magnetization components, then implement similar approach in Radia's B_comp.

---

**Last Updated**: 2025-11-23 (after R1 fix and in-plane investigation)
