# ELF MAGIC Fortran Code Analysis

**Date**: 2025-11-23
**Purpose**: Understand how ELF MAGIC handles tetrahedral mesh magnetization

---

## Key Subroutines

### MM4TS() - Magnetic Moment for Tetrahedra (lines 3567-3590)

```fortran
subroutine MM4TS() ! MM4T SUM
    do LL = 1,3  ! Loop over x, y, z components
        do MEN = 1,4  ! Loop over 4 tetrahedral faces
            call NPXX4T()          ! Get face vertices
            call AREAV3(AREA)      ! Compute face normal VV = edge1 × edge2
            call AAOBB(VV,UNIV(1,LL),Q)  ! Q = VV · e_LL (dot product)
            call CENT(3)           ! Compute face centroid
            AA = EE-EEC
            call AAOBB(VV,AA,A)
            if(A < 0.0d0) Q = -Q   ! Fix sign based on orientation
            if(KMOF == 2) Q = Q*PAA(NAA)  ! Apply magnetization?
            call DOMBS(Q,0)        ! Compute field contribution
        end do
    end do
end
```

**Key Observations**:
1. Loops over LL=1,2,3 (x,y,z directions)
2. For each component, loops over all 4 faces of tetrahedron
3. Computes `Q = VV · e_LL` where VV is face normal (un-normalized)
4. **Q represents the LL-th component of the face normal vector**
5. This is NOT looping over magnetization components directly!

### AREAV3() - Compute Face Area and Normal (lines 5271-5282)

```fortran
subroutine AREAV3(AREA)
    do K = 1,3
        AA(K) = XX(K,2)-XX(K,1)  ! Edge vector 1
        BB(K) = XX(K,3)-XX(K,1)  ! Edge vector 2
    end do
    call AAXBB(AA,BB,VV,AREA1)   ! VV = AA × BB (cross product)
    AREA = AREA1*HAL             ! Area = |VV|/2
end
```

**Result**: VV contains the face normal vector (un-normalized, length = 2*Area)

### DOMBS() - Distribute Over Multiple Biot-Savart (lines 3592-3663)

```fortran
subroutine DOMBS(WA,ISUM)
    ! ... setup local coordinate system AA, BB, CC ...
    ! ... transform vertices to 2D polygon in local frame ...

    W = SSS(4)*WA  ! Magnetic charge density weight

    call CQ_HT43(AA,BB,CC,YY,XY,XJJ1,FJJ2,PJJ2,W,MX,NII,M43)
    ! Compute field from polygon charge distribution
end
```

**Key**: W is the magnetic charge density weight passed to field calculation

### CQ_HT43() - Calculate Charge field for Triangles/Quads (lines 3757-3870)

```fortran
subroutine CQ_HT43(AA,BB,CC,YY,XY,XX,FGH,PJJ,W,MXX,NII,KAdo)
! CALC Q(CHARGE) T(3D) 43(QUAD,TRIA)

    ! Compute field in local coordinate system -> HH1, HH2, HH3
    ! (Lines 3845-3863: analytical formulas)

    ! Transform back to global coordinates:
    FGH(1,I) = FGH(1,I)+HH1*AA(1)+HH2*BB(1)+HH3*CC(1)
    FGH(2,I) = FGH(2,I)+HH1*AA(2)+HH2*BB(2)+HH3*CC(2)
    FGH(3,I) = FGH(3,I)+HH1*AA(3)+HH2*BB(3)+HH3*CC(3)
end
```

**Key**: AA, BB, CC are basis vectors of local coordinate system in global frame

---

## Algorithm Flow

### What ELF MAGIC Does:

```
For each magnetization component (x, y, z):
    For each tetrahedral face:
        1. Compute face normal VV = edge1 × edge2
        2. Extract component: Q = VV_x (or VV_y, or VV_z)
        3. Apply orientation correction
        4. Compute field contribution from charge density W
```

### Comparison with Radia:

**Radia's approach** (rad_polyhedron.cpp:444-458):
```cpp
// For each face:
1. Transform magnetization to local coords: Magn_local = R^T * Magn_global
2. Extract normal component: W = Magn_local.z
3. Compute field in local coords
4. Transform field back to global coords
```

**Key Difference**:
- **ELF**: Loops over global x,y,z directions first, then over faces
- **Radia**: Loops over faces, transforms magnetization to local coordinates

**Why the difference?**:
- ELF's approach may be more efficient for matrix assembly (process all x-components together)
- Radia's approach is simpler for individual face evaluation

**Physical equivalence**:

Both should give the same result because:
```
σ = M · n  (magnetic charge density)

ELF approach:
  σ = Mx*nx + My*ny + Mz*nz  (computed in 3 separate passes)

Radia approach:
  σ = M_local.z  (where z-axis of local frame = n)
```

These are mathematically equivalent.

---

## Conclusion: ELF MAGIC Analysis

**Finding**: ELF MAGIC's loop structure is an implementation detail for efficiency, not a fundamental algorithmic difference.

**Key insight**: Both ELF and Radia compute the same quantity:
```
σ = M · n
```

**For Radia tetrahedra**: The current implementation (using only Magn.z in local coordinates) is **theoretically correct**.

**Remaining question**: Why do Radia tetrahedra show 1165-6864% errors?

Possible causes:
1. R1 rotation matrix construction has subtle errors
2. Coordinate transformation (TrVectField_inv) has sign or direction issues
3. Face orientation (inward vs outward normal) is inconsistent
4. There may be additional terms needed beyond simple magnetic charge (σ = M·n)

---

## Next Steps

1. Add debug output to verify:
   - Face normal directions
   - Magnetization vectors in local coordinates
   - Magnetic charge density W values

2. Compare with hexahedral case:
   - Why do hexahedra work correctly?
   - What is different for tetrahedra?

3. Check reference implementations:
   - Verify analytical formulas in RadAnalyticalFieldFromPolygonCharge
   - Cross-check with ELF MAGIC's CQ_HT43 implementation

---

**Test Results** (2025-11-23, after R1 fix):

| Point [m] | Built-in \|H\| | Tetra \|H\| | Error % |
|-----------|---------------|------------|---------|
| [0.08, 0, 0] | 0.033760 | 2.351307 | **6864%** |
| [0, 0.08, 0] | 0.033760 | 1.172550 | **3373%** |
| [0, 0, 0.08] | 0.053851 | 0.681724 | **1165%** |
| [0.15, 0, 0] | 0.010733 | 0.469363 | 4273% |
| [0, 0, 0.15] | 0.019074 | 0.423825 | 2122% |

**Note**: Lowest error is for z-axis point [0,0,0.08], which aligns with magnetization direction [0,0,1.2].

---

**Last Updated**: 2025-11-23 11:20
