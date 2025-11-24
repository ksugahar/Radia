# Radia Tetrahedral Mesh Root Cause Analysis

## Problem Statement

Radia's tetrahedral mesh support has **critical accuracy and stability issues**:

1. **Permanent magnets**: 2999% field error  
2. **Linear materials (MatLin)**: Complete solver divergence (NaN)

## Root Cause: Coordinate System Mismatch

**Location**: `radTPolygon::B_comp()` in `src/core/rad_planar_2d_part2.cpp`

**The Problem**:
```cpp
// radTPolygon stores face vertices as 2D points (x, y)
radTVect2dVect EdgePointsVector;  // TVector2d = (x, y) only!

// Basis vectors HARDCODED for XY plane
TVector3d AA(1, 0, 0);  // X-axis
TVector3d BB(0, 1, 0);  // Y-axis  
TVector3d CC(0, 0, 1);  // Z-axis (face normal)
```

**Why This Fails**:
- Tetrahedral faces have **arbitrary 3D orientations**
- Using fixed basis [1,0,0], [0,1,0], [0,0,1] is **fundamentally wrong**

**Result**: 2999% error

## Solution: Analytical Method

**Key Insight**: Construct basis from **actual face vertices** (not transformation matrix)

**Algorithm**:
1. Get 3 face vertices P1, P2, P3  
2. Construct orthonormal basis:
   - AA = P2 - P1  
   - BB = P3 - P1  
   - CC = AA × BB (cross product)
   - Orthonormalize via Gram-Schmidt
3. Transform evaluation point to local coordinates
4. Apply analytical formulas (log + atan terms)
5. Transform field back to global

**Expected Accuracy**: <10% error (comparable to hexahedral: 3.7%)

## Implementation Plan

See [TETRAHEDRAL_IMPLEMENTATION.md](../../../TETRAHEDRAL_IMPLEMENTATION.md) for full details.

**Timeline**: 5-7 days (2 days coding + 3 days testing)

## References

- **Implementation**: [rad_polyhedron.cpp](../../../src/core/rad_polyhedron.cpp#L701-L873)
- **Test results**: [COMPARISON_SUMMARY.md](COMPARISON_SUMMARY.md)

---

**Date**: 2025-11-24  
**Status**: Root cause identified, solution designed, implementation pending
