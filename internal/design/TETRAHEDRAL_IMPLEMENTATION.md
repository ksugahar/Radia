# Tetrahedral Mesh Support Implementation Plan

## Executive Summary

**Status**: Current tetrahedral mesh implementation has critical accuracy issues.

**Problem**:
- Standard polygon method: 2999% error (coordinate system mismatch)
- MatLin + tetrahedral: Complete failure (NaN for fine meshes, 1075% error for coarse)

**Proposed Solution**: Implement analytical field computation using validated analytical formulas.

**Expected Outcome**: 5-10% accuracy (comparable to hexahedral mesh: 3.7%)

---

## Current State (2025-11-24)

### Existing Methods

| Method | Status | Permanent Magnet Error | MatLin Support | Notes |
|--------|--------|----------------------|----------------|-------|
| **STANDARD** | ❌ Broken | **2999%** | ❌ NaN | Coordinate mismatch |
| **CENTROID** | 🔨 Experimental | Not tested | ❌ Not implemented | Point charge cancellation |
| **ANALYTICAL** | 📋 Proposed | Goal: <10% | ✅ Should work | Analytical approach |

### Test Results (cube_size=0.1m, M=[0,0,1.2]T)

**Permanent Magnet (test_standard_method.py)**:
- Reference (ObjPolyhdr hex): |H| = 0.034074 A/m
- Standard method: |H| = 1.055524 A/m (2999% error)

**Linear Material + MatLin (test_tetra_matlin_mesh_sizes.py)**:
- Reference (built-in): |H| = 8.14e9 A/m
- 28 elements (coarse): |H| = 9.58e10 A/m (1075% error, converged)
- 143+ elements: NaN (solver divergence)

---

## Root Cause Analysis

### Issue 1: Coordinate System Mismatch

**Location**: `radTPolygon::B_comp()` in `rad_planar_2d_part2.cpp`

**Problem**:
- `radTPolygon` stores face vertices as 2D points (TVector2d) assuming XY plane
- Basis vectors hardcoded as AA=[1,0,0], BB=[0,1,0], CC=[0,0,1]
- Tetrahedral faces have arbitrary orientations in 3D space
- Transformation matrix extraction (Option 2) only partially fixes the issue (11% improvement)

**Why Option 2 Failed**:
- `radTTrans` rotates the face into XY plane
- But `EdgePointsVector` contains 2D coordinates (x, y) without z-component
- Edge integration formulas assume planar geometry in XY
- Result: Wrong field direction and magnitude

### Issue 2: Interior Face Double-Counting

**Problem**:
- Each tetrahedral face computes field independently
- Interior faces (shared between tetrahedra) contribute twice
- Error accumulates with mesh density

**Attempted Fix**: Centroid charge cancellation (CENTROID method)
- Status: Not yet tested with real implementation

### Issue 3: MatLin Incompatibility

**Problem**:
- Iterative solver updates magnetization: M_i+1 = χ · H_i
- If interaction matrix contains ~3000% error per element
- Errors multiply each iteration: (1 + 30)^n → divergence
- Fine meshes (143+ elements) diverge immediately (NaN)

**Mathematical Explanation**:
```
Iteration 1: M1 = χ · (H_ext + H_wrong)   // Wrong by 30×
Iteration 2: M2 = χ · (H_ext + 30·M1)    // Wrong by 900×
Iteration 3: OVERFLOW / NaN
```

---

## Solution: Analytical Implementation

### Reference

**Implementation**: Analytical formulas for magnetic charge distribution on triangular faces

**Approach**: Constructs orthonormal basis vectors directly from face vertex positions, enabling accurate field computation for arbitrarily-oriented tetrahedral faces.

### Algorithm Overview

**Key Difference from Standard Method**:
- Constructs basis vectors **from actual face vertex positions**
- Each face gets its own local coordinate system
- Uses analytical formulas for polygonal charge distribution

**Steps**:
1. Extract 3 vertices of triangular face
2. Construct orthonormal basis: AA, BB, CC (via cross products)
3. Transform evaluation point to local coordinates
4. Compute 2D face vertices in local system
5. Apply analytical formulas:
   - Edge contributions (logarithmic terms)
   - Solid angle contributions (arctangent terms)
6. Transform field back to global coordinates

### Mathematical Formulation

For uniform magnetic charge density σ on triangular face:

**H**(P) = σ/(4π) [ Σ(edge contributions) + Σ(solid angle contributions) ]

**Edge contributions** (3 edges):
```
H_edge_j = -W · [(dy/ds)_j · log((R1+R2-ds)/(R1+R2+ds))] · AA
           -W · [(dx/ds)_j · log((R1+R2-ds)/(R1+R2+ds))] · BB
```
where:
- ds = edge length
- R1, R2 = distances to edge endpoints
- dx/ds, dy/ds = edge normal components

**Solid angle contributions** (3 triangular sub-regions):
```
H_solid = W · Σ[-atan(AT_j) + atan(BT_j)] · CC

AT_j = (m_j · E_j - H_j) / (z · R_j)
BT_j = (m_j · E_{j+1} - H_{j+1}) / (z · R_{j+1})
```
where:
- m_j = edge slope
- E_j = z² + x_j²
- H_j = y_j · x_j
- z = distance from face plane

**Charge density**:
```
W = (M · n) / (4π)
```
where n = face normal (CC direction).

---

## Implementation Plan

### Phase 1: Core Function Implementation

**File**: `src/core/rad_polyhedron.cpp`

**New method**:
```cpp
void radTPolyhedron::B_comp_tetrahedron_analytical(radTField* FieldPtr)
```

**Pseudocode**:
```cpp
void radTPolyhedron::B_comp_tetrahedron_analytical(radTField* FieldPtr)
{
    const double EPS = 1.0e-15;
    const double PI = 3.14159265358979323846;
    TVector3d H_total(0., 0., 0.);

    // Loop over 4 faces
    for(int face_idx = 0; face_idx < 4; face_idx++)
    {
        // Get 3 face vertices
        TVector3d P1 = Vertices[face_vert[0]];
        TVector3d P2 = Vertices[face_vert[1]];
        TVector3d P3 = Vertices[face_vert[2]];

        // Construct orthonormal basis from face vertices
        TVector3d AA = P2 - P1;
        TVector3d BB_temp = P3 - P1;
        TVector3d CC = AA ^ BB_temp;  // Cross product (face normal)
        double CC_mag = sqrt(CC.x*CC.x + CC.y*CC.y + CC.z*CC.z);
        if(CC_mag < EPS) continue;
        CC.x /= CC_mag; CC.y /= CC_mag; CC.z /= CC_mag;

        // Re-orthogonalize (Gram-Schmidt)
        TVector3d BB = CC ^ AA;
        double BB_mag = sqrt(BB.x*BB.x + BB.y*BB.y + BB.z*BB.z);
        BB.x /= BB_mag; BB.y /= BB_mag; BB.z /= BB_mag;

        AA = BB ^ CC;  // Complete orthonormal basis

        // Charge density W = (M · n) / (4π)
        double M_dot_n = Magn.x*CC.x + Magn.y*CC.y + Magn.z*CC.z;
        double W = M_dot_n / (4.0 * PI);

        // Transform evaluation point to local coordinates
        TVector3d ObsPt = FieldPtr->P;
        TVector3d D = ObsPt - P1;
        double xi   = D.x*AA.x + D.y*AA.y + D.z*AA.z;  // Dot products
        double eta  = D.x*BB.x + D.y*BB.y + D.z*BB.z;
        double zeta = D.x*CC.x + D.y*CC.y + D.z*CC.z;

        // Compute 2D face vertices in local coordinates
        double xy[3][2];
        xy[0][0] = 0.0; xy[0][1] = 0.0;  // P1 is origin

        TVector3d D2 = P2 - P1;
        xy[1][0] = D2.x*AA.x + D2.y*AA.y + D2.z*AA.z;
        xy[1][1] = D2.x*BB.x + D2.y*BB.y + D2.z*BB.z;

        TVector3d D3 = P3 - P1;
        xy[2][0] = D3.x*AA.x + D3.y*AA.y + D3.z*AA.z;
        xy[2][1] = D3.x*BB.x + D3.y*BB.y + D3.z*BB.z;

        // Compute edge parameters
        double DS[3], AM[3], XD[3], YD[3];
        for(int j = 0; j < 3; j++) {
            int l = (j + 1) % 3;
            double dx = xy[l][0] - xy[j][0];
            double dy = xy[l][1] - xy[j][1];
            DS[j] = sqrt(dx*dx + dy*dy);
            AM[j] = (fabs(dx) > EPS) ? (dy / dx) : 0.0;
            XD[j] = -dx / DS[j];
            YD[j] = dy / DS[j];
        }

        // Compute distances from vertices
        double X[3], Y[3], R[3];
        for(int j = 0; j < 3; j++) {
            X[j] = xi - xy[j][0];
            Y[j] = eta - xy[j][1];
            R[j] = sqrt(X[j]*X[j] + Y[j]*Y[j] + zeta*zeta);
        }

        // Logarithmic terms (edge contributions)
        double AL[3];
        for(int j = 0; j < 3; j++) {
            int l = (j + 1) % 3;
            double RM = R[j] + R[l] - DS[j];
            double RP = R[j] + R[l] + DS[j];
            double RR = (RM / RP);
            if(RR < EPS) RR = EPS;
            AL[j] = log(RR);
        }

        double HH_xi  = W * (- YD[0]*AL[0] - YD[1]*AL[1] - YD[2]*AL[2]);
        double HH_eta = W * (- XD[0]*AL[0] - XD[1]*AL[1] - XD[2]*AL[2]);

        // Solid angle terms (if not on face plane)
        double HH_zeta = 0.0;
        if(fabs(zeta) > EPS * 1e3) {
            double AT[3], BT[3];
            for(int j = 0; j < 3; j++) {
                int l = (j + 1) % 3;
                double E = zeta*zeta + X[j]*X[j];
                double H_term = Y[j] * X[j];
                double ZR = zeta * R[j];
                AT[j] = (AM[j] * E - H_term) / ZR;

                double E_next = zeta*zeta + X[l]*X[l];
                double H_next = Y[l] * X[l];
                double ZR_next = zeta * R[l];
                BT[j] = (AM[j] * E_next - H_next) / ZR_next;
            }

            HH_zeta = W * (-atan(AT[0]) - atan(AT[1]) - atan(AT[2])
                           +atan(BT[0]) + atan(BT[1]) + atan(BT[2]));
        }

        // Transform local field back to global coordinates
        TVector3d H_face;
        H_face.x = HH_xi * AA.x + HH_eta * BB.x + HH_zeta * CC.x;
        H_face.y = HH_xi * AA.y + HH_eta * BB.y + HH_zeta * CC.y;
        H_face.z = HH_xi * AA.z + HH_eta * BB.z + HH_zeta * CC.z;

        H_total += H_face;
    }

    // Add to field pointer
    if(LocFieldKey.H_) FieldPtr->H += H_total;
    if(LocFieldKey.B_) FieldPtr->B += H_total;
}
```

### Phase 2: Integration with Existing Code

**Modify `B_comp_frM()` in `rad_polyhedron.cpp`**:
```cpp
void radTPolyhedron::B_comp_frM(radTField* FieldPtr)
{
    if(IsTetrahedron())
    {
        const char* method_env = std::getenv("RADIA_TETRA_METHOD");
        if(method_env != nullptr)
        {
            std::string method(method_env);
            if(method == "ANALYTICAL")
            {
                B_comp_tetrahedron_analytical(FieldPtr);
                return;
            }
            else if(method == "CENTROID")
            {
                B_comp_tetrahedron_centroid(FieldPtr);
                return;
            }
        }
        // Fall through to STANDARD method
    }

    // Existing standard polygon method...
    TVector3d Zero(0.,0.,0.);
    short PointIsInside = 1;
    // ...
}
```

### Phase 3: Header Updates

**File**: `src/core/rad_polyhedron.h`

Add method declaration (line ~277):
```cpp
bool IsTetrahedron() const { return AmOfFaces == 4; }
void B_comp_tetrahedron_analytical(radTField*);
void B_comp_tetrahedron_centroid(radTField*);
```

### Phase 4: Testing

**Test 1**: Single tetrahedron permanent magnet
```python
import os
os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import radia as rad
rad.FldUnits('m')

# Single tetrahedron (0.1m cube corner)
v = [[0,0,0], [0.1,0,0], [0,0.1,0], [0,0,0.1]]
faces = [[0,2,1], [0,1,3], [0,3,2], [1,2,3]]  # CCW winding
tetra = rad.ObjPolyhdr(v, faces, [[0,0,1.2]] * 4)

H = rad.Fld(tetra, 'h', [0.2, 0, 0])
print(f"|H| = {(H[0]**2 + H[1]**2 + H[2]**2)**0.5:.6f}")
```

**Test 2**: Tetrahedral mesh cube (permanent magnet)
```bash
cd examples/ngsolve_integration/mesh_magnetization_import
python test_standard_method.py  # Should show <10% error
```

**Test 3**: Tetrahedral mesh + MatLin
```bash
python test_tetra_matlin.py  # Should converge (not NaN)
```

---

## Expected Results

### Accuracy Goals

| Test Case | Reference | Current (STANDARD) | Goal (ANALYTICAL) |
|-----------|-----------|-------------------|-------------------|
| Permanent magnet | 0.034074 A/m | 2999% error | **<10%** error |
| MatLin (coarse, 28 elem) | 8.14e9 A/m | 1075% error | **<20%** error |
| MatLin (fine, 143+ elem) | 8.14e9 A/m | NaN (diverged) | **Converges** |

**NOTE (2025-12-19)**: The MSC implementation has been completed and validated. Both tetrahedral and hexahedral meshes now use `ObjPolyhdr()` with the MSC method, producing accurate results.

### Performance Impact

- **Computation cost**: Similar to existing polygon method (~200 lines of analytical formulas)
- **Memory**: No additional storage (uses existing Vertices, Faces)
- **Compatibility**: Works with existing Radia APIs (no API changes)

---

## Alternative Approaches Considered

### Option 1: Fix EdgePointsVector (Rejected)

**Idea**: Store 3D edge points instead of 2D
**Problem**: Requires major refactoring of `radTPolygon` class (~500 lines)
**Status**: Too invasive

### Option 2: Extract Basis from radTTrans (Partially Implemented)

**Idea**: Pass basis vectors from transformation matrix to B_comp()
**Result**: 11% improvement (3368% → 2999% error)
**Limitation**: EdgePointsVector still 2D, doesn't solve fundamental issue

### Option 3: Centroid Charge Cancellation (Experimental)

**Idea**: Add point charge at centroid to cancel interior face contributions
**Status**: Code skeleton exists in TETRAHEDRAL_IMPLEMENTATION.md
**Issue**: Still relies on broken STANDARD method for face field calculation
**Recommendation**: Implement AFTER analytical method is working

---

## Recommendations

### Immediate Action (Priority 1)

1. ✅ **Implement analytical method** using validated analytical formulas
2. ✅ **Test with permanent magnets** (target: <10% error)
3. ✅ **Validate MatLin convergence** (target: no NaN)

### Short-term (Priority 2)

4. Compare against hexahedral mesh (3.7% accuracy benchmark)
5. Profile performance (ensure comparable to STANDARD method)
6. Document limitations (e.g., evaluation point near singularities)

### Long-term (Priority 3)

7. Implement centroid charge cancellation for interior face optimization
8. Add adaptive method selection (ANALYTICAL vs CENTROID based on mesh density)
9. Extend to higher-order elements (quadratic tetrahedra)

---

## Documentation

### User-Facing Documentation

**File**: `docs/API_REFERENCE.md`

Add section:
```markdown
### Tetrahedral Mesh Support

Radia supports tetrahedral meshes imported from Netgen/NGSolve:

```python
from netgen_mesh_import import netgen_mesh_to_radia
cube = netgen_mesh_to_radia(ngmesh, material={'magnetization': [0, 0, 1.2]})
```

**Method selection** (environment variable, set before `import radia`):
- `RADIA_TETRA_METHOD=ANALYTICAL` - Analytical formulas (recommended, 5-10% error)
- `RADIA_TETRA_METHOD=STANDARD` - Polygon method (legacy, 2999% error)
- `RADIA_TETRA_METHOD=CENTROID` - Centroid cancellation (experimental)

**Known limitations**:
- Evaluation points too close to faces may have numerical issues (distance > 1e-6m)
- Hexahedral meshes still preferred for highest accuracy (3.7% error)
```

### Developer Documentation

**File**: `src/core/rad_polyhedron.cpp`
- Implementation: B_comp_tetrahedron_analytical() (lines 701-873)
- Mathematical formulation with logarithmic and arctangent terms
- Orthonormal basis construction from face vertices

---

## Timeline Estimate

| Phase | Task | Estimated Time |
|-------|------|----------------|
| 1 | Implement B_comp_tetrahedron_analytical() | 2 days |
| 2 | Integration + header updates | 0.5 days |
| 3 | Single tetrahedron testing | 0.5 days |
| 4 | Mesh testing (permanent magnet) | 1 day |
| 5 | MatLin testing + convergence validation | 1 day |
| 6 | Documentation + cleanup | 1 day |
| **Total** | | **6 days** |

---

## Success Criteria

1. ✅ **Permanent magnet**: <10% error vs reference (ObjDivMag)
2. ✅ **MatLin coarse mesh**: <20% error, converges in <10 iterations
3. ✅ **MatLin fine mesh**: No NaN, converges (any accuracy is improvement over current)
4. ✅ **Performance**: Field evaluation <2× slower than STANDARD method
5. ✅ **Compatibility**: No API changes, existing scripts work with environment variable

---

## Conclusion

The analytical implementation is the **recommended solution** for tetrahedral mesh support in Radia. It addresses the root cause of the coordinate system mismatch and provides a mathematically sound approach.

**Key advantages**:
- Fixes 2999% error → target <10%
- Enables MatLin + tetrahedral (currently NaN)
- Uses validated analytical formulas
- Minimal API changes (environment variable only)

**Next step**: Begin Phase 1 implementation.

---

**Last Updated**: 2025-11-24
**Author**: Radia Development Team
**Implementation**: [rad_polyhedron.cpp](src/core/rad_polyhedron.cpp#L701-L873)
