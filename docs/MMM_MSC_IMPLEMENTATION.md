# MMM with MSC (Magnetic Surface Charge) Implementation in Radia

## Summary

This document describes how Radia implements the **MMM (Magnetic Moment Method)** using the
**MSC (Magnetic Surface Charge)** approach for field computation.

**Key Finding**: Radia **already implements MSC** for both rectangular blocks and general polyhedra.
No new implementation is needed - the existing code uses the same mathematical approach as ELF_MAGIC.

## Background

### MMM (Magnetic Moment Method)
- Represents magnetic objects as distributions of magnetic moments
- Each element has uniform magnetization M
- Computes mutual interactions to determine equilibrium magnetization

### MSC (Magnetic Surface Charge Method)
- Equivalent formulation using magnetic surface charges
- Surface charge density: sigma = M . n (magnetization dot outward normal)
- Field computed from surface integral over charged polygons
- Mathematically equivalent to volume integral approach for uniform M

## Radia Implementation

### 1. Rectangular Blocks (radTRecMag)

**File**: [rad_rectangular_block.cpp](../src/core/rad_rectangular_block.cpp)

**Method**: Direct analytical formula for rectangular parallelepiped

```cpp
void radTRecMag::B_comp(radTField* FieldPtr)
{
    // Uses atan-based analytical formula
    // 8 vertices, computing atan terms for each corner combination
    T0.x = atan(TransAtans(...));  // x-component
    T0.y = atan(TransAtans(...));  // y-component
    T0.z = atan(TransAtans(...));  // z-component
    ...
}
```

**Verification**: Self-demagnetization factor for cube = 1/3 (correct)

```
Demagnetization factors for cube:
  N_xx = -Hx/Mx = 0.333333
  N_yy = -Hy/My = 0.333333
  N_zz = -Hz/Mz = 0.333333
  Sum: 1.000000 (should be 1.0)
```

### 2. General Polyhedra (radTPolyhedron)

**File**: [rad_polyhedron.cpp](../src/core/rad_polyhedron.cpp)

**Method**: MSC - Sum over polygon faces

```cpp
void radTPolyhedron::B_comp_frM(radTField* FieldPtr)
{
    // Loop over all faces
    for(int i=0; i<AmOfFaces; i++)
    {
        // Transform to local coordinates
        LocField.P = TransPtr->TrPoint_inv(SumLocField.P);
        PgnPtr->Magn = TransPtr->TrVectField_inv(Magn);

        // Compute field from this polygon face (MSC)
        PgnPtr->B_comp(&LocField);

        // Transform back to global and sum
        SumLocField.H += TransPtr->TrVectField(LocField.H);
    }
}
```

### 3. Polygon Face Field (radTPolygon)

**File**: [rad_planar_2d_part2.cpp](../src/core/rad_planar_2d_part2.cpp)

**Method**: Analytical formula for charged polygon

```cpp
void radTPolygon::B_comp(radTField* FieldPtr)
{
    // Magnetic charge density: W = (1/4pi) * M_z
    // where M_z is magnetization component normal to polygon
    double W = ConstForH * Magn.z;

    // Call analytical formula for polygon charge
    RadAnalyticalFieldFromPolygonCharge(
        AA, BB, CC, YY,           // Local coordinate basis
        EdgePointsVector,          // Polygon vertices (2D)
        obs_points,                // Observation points
        field_result,              // Output field
        W,                         // Charge density
        1,                         // Element index
        AmOfEdgePoints             // Number of vertices
    );
}
```

## Comparison: Radia vs ELF_MAGIC

| Feature | Radia | ELF_MAGIC |
|---------|-------|-----------|
| **Method** | MMM with MSC | MMM with MSC |
| **Rectangular elements** | Analytical formula (8-vertex atan) | MSC (6 quad faces -> 12 triangles) |
| **Tetrahedral elements** | MSC (4 triangular faces) | MSC (4 triangular faces) |
| **Hexahedral elements** | Analytical (radTRecMag) or MSC (radTPolyhedron) | MSC (6 quad faces) |
| **General polyhedra** | MSC (N polygon faces) | Not directly supported |
| **Polygon field formula** | log + atan (same as ELF_MAGIC) | log + atan |
| **Self-demagnetization** | Correct (N=1/3 for cube) | Correct (N=1/3 for cube) |

## Mathematical Equivalence

Both Radia and ELF_MAGIC use the same fundamental formula for field from a charged polygon:

### Tangential field (H_x, H_y in local frame):
```
H_tangential = sigma * sum_edges( -edge_dir * log((r1 + r2 - L)/(r1 + r2 + L)) )
```

### Normal field (H_z in local frame):
```
H_normal = sigma * sum_edges( atan(arg1) - atan(arg2) )
```

Where:
- sigma = M . n (surface charge density)
- r1, r2 = distances from observation point to edge endpoints
- L = edge length
- arg1, arg2 = geometric arguments involving edge parameters

## Benchmark Results (Nonlinear Hexahedral Solver)

### Radia v1.3.9 (N=10, 1000 elements, mu_r=1000)

| Solver | Time (s) | M_avg_z (A/m) | Iterations |
|--------|----------|---------------|------------|
| LU (Method 0) | 2.01 | 178,631 | 8 |
| BiCGSTAB (Method 1) | 0.54 | 178,620 | 27 |

### ELF_MAGIC (HACApK) (N=10, 1000 elements, mu_r=1000)

| Solver | Time (s) | M_avg_z (A/m) | Iterations |
|--------|----------|---------------|------------|
| LU | 4.78 | 178,664 | 8 |
| HACApK | 7.63 | 178,665 | 264 (linear) |

**Key Observations**:
1. Both codes produce nearly identical M_avg_z values (< 0.03% difference)
2. Radia LU is ~2.4x faster than ELF_MAGIC LU
3. Radia BiCGSTAB is ~14x faster than ELF_MAGIC HACApK

## Conclusion

Radia's implementation of MMM with MSC is:
1. **Mathematically equivalent** to ELF_MAGIC
2. **Already implemented** - no new development needed
3. **Faster** due to optimized analytical formulas for common shapes
4. **More general** - supports arbitrary polyhedra natively

The "MMM with MSC" approach is the standard method in both codes.
The performance difference comes from:
- Radia: C++ with OpenBLAS LAPACK, optimized for common shapes
- ELF_MAGIC: Fortran with generic MSC for all element types

---

**Last Updated**: 2025-12-11 (6-face hexahedral MSC added)
**Author**: Claude Code
**Project**: Radia Magnetic Field Computation
