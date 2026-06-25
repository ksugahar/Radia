# Image Symmetry Implementation Design

## Overview

This document describes the implementation of Image symmetry for MSC hexahedra in Radia.

**Status (updated 2026-06-23):** this is a historical design note for the pre-multipole-moment MMM
surface-charge kernel. The EIEM2 `Compute6x6BlockFast` / `Compute5x5BlockFast` kernels were removed;
current surface-charge IMA is handled by `BuildMomentSystemCore` /
`CentroidFieldGradFromFace`, while tetrahedral MMM still uses `Compute3x3BlockFast`.
Old snippets and class names are kept below only to explain the original design.

**Note**: As of 2026-01-31, `TrfMlt`, `TrfPlSym`, `TrfZerPara`, and `TrfZerPerp` have been **REMOVED** from Radia. `rad.Image()` is the only supported method for plane symmetry with MSC hexahedra.

## Why TrfMlt Was Removed

The original TrfMlt had fundamental design issues:

1. **DOF Sharing Issue**: TrfMlt shared DOFs between original and virtual elements, but MSC requires independent DOFs when the external field is perpendicular to the mirror plane.

2. **Face Permutation Required**: For x-mirror, face 1 and face 3 must be permuted. TrfMlt did not handle this.

3. **Design Philosophy**: Element-based management (independent DOFs per element) is essential for correct physics. Face-based management (shared DOFs) causes errors.

## Image Matrix Construction

The Image method constructs the mirror matrix using **image summation** during assembly:

```
N_Image[i,j] = N[i,j] + sign * N[i, mirror_j] @ P
```

Where:
- `sign = +1` for **symmetric boundary condition** (field tangent to plane)
- `sign = -1` for **antisymmetric boundary condition** (field normal to plane)

Where:
- `N[i,j]`: Interaction from element j to element i (direct)
- `N[i, mirror_j]`: Interaction from mirror image of element j to element i
- `P`: DOF permutation matrix (swaps face 1 and face 3 for x-mirror)

## Face Ordering

| DOF Index | Face Normal Direction |
|-----------|----------------------|
| 0 | Face 0 |
| 1 | Face 1 |
| 2 | Face 2 |
| 3 | Face 3 |
| 4 | Face 4 |
| 5 | Face 5 |

## Permutation Matrices

### X-Mirror (YZ plane)

Swaps DOF 1 and DOF 3:

```
P_x = [[1,0,0,0,0,0],
       [0,0,0,1,0,0],
       [0,0,1,0,0,0],
       [0,1,0,0,0,0],
       [0,0,0,0,1,0],
       [0,0,0,0,0,1]]
```

### Y-Mirror (XZ plane)

Swaps DOF 0 and DOF 2:

```
P_y = [[0,0,1,0,0,0],
       [0,1,0,0,0,0],
       [1,0,0,0,0,0],
       [0,0,0,1,0,0],
       [0,0,0,0,1,0],
       [0,0,0,0,0,1]]
```

### Z-Mirror (XY plane)

Swaps DOF 4 and DOF 5:

```
P_z = [[1,0,0,0,0,0],
       [0,1,0,0,0,0],
       [0,0,1,0,0,0],
       [0,0,0,1,0,0],
       [0,0,0,0,0,1],
       [0,0,0,0,1,0]]
```

## API Reference

### Python API

```python
import radia as rad

# Set Image symmetry with optional sign
n_ima = rad.Image(intrc, symmetry)

# Build Image interaction matrix
rad.BuildImageMatrix(intrc)
```

**Symmetry String Format**: `[+|-]axis`

| String | Description | Boundary Condition |
|--------|-------------|-------------------|
| `+x`, `x` | x-mirror (symmetric) | Field tangent to YZ plane |
| `-x` | x-mirror (antisymmetric) | Field normal to YZ plane |
| `+y`, `y` | y-mirror (symmetric) | Field tangent to XZ plane |
| `-y` | y-mirror (antisymmetric) | Field normal to XZ plane |
| `+z`, `z` | z-mirror (symmetric) | Field tangent to XY plane |
| `-z` | z-mirror (antisymmetric) | Field normal to XY plane |
| `+xy` | Quarter model (both symmetric) | |
| `-xz` | xz-mirror (antisymmetric) | |
| `+xyz` | Eighth model (all symmetric) | |

### C++ Internal API

```cpp
// Set Image symmetry configuration
// symmetry: axis flags (IMA_X, IMA_Y, IMA_Z)
// sign: +1 (symmetric) or -1 (antisymmetric)
int SetIMASymmetry(int symmetry, int sign = 1);

// Get Image-reduced element count
int GetIMAElementCount();

// Check if element is in the positive half-space for given symmetry
bool IsElementInIMARegion(int elem_idx);
```

### 2. Element Mapping

For x-mirror symmetry:
- Original elements: Elements with x_center >= 0
- Mirror elements: Same elements reflected through x=0 plane

```cpp
struct IMASymmetryConfig {
    int symmetry_flags;           // Bitfield: X=1, Y=2, Z=4
    std::vector<int> ima_to_full; // IMA element index -> full element index
    std::vector<int> mirror_map;  // IMA element index -> mirror element index in full model
    std::vector<int> perm_x;      // DOF permutation for x-mirror (size 6)
    std::vector<int> perm_y;      // DOF permutation for y-mirror (size 6)
    std::vector<int> perm_z;      // DOF permutation for z-mirror (size 6)
};
```

### 3. Modified Matrix Assembly (Dense)

```cpp
int SetupInteractMatrix_IMA()
{
    // For each IMA element pair (i, j):
    // N_IMA[i,j] = N[full_i, full_j] + N[full_i, mirror_j] @ P_x

    for(int i = 0; i < n_ima_elem; i++)
    {
        int full_i = ima_to_full[i];

        for(int j = 0; j < n_ima_elem; j++)
        {
            int full_j = ima_to_full[j];
            int mirror_j = mirror_map[full_j];

            // Compute direct interaction
            double N_direct[36];
            Compute6x6Block(full_i, full_j, N_direct);

            // Compute mirror interaction
            double N_mirror[36];
            Compute6x6Block(full_i, mirror_j, N_mirror);

            // Apply permutation: N_mirror @ P_x
            double N_mirror_perm[36];
            ApplyPermutation(N_mirror, perm_x, N_mirror_perm);

            // Sum: N_IMA = N_direct + N_mirror_perm
            double N_IMA[36];
            for(int k = 0; k < 36; k++)
                N_IMA[k] = N_direct[k] + N_mirror_perm[k];

            // Store in flattened matrix
            StoreBlock(i, j, N_IMA);
        }
    }
}
```

### 4. HACApK Adaptation

For HACApK on-demand matrix computation:

```cpp
void ComputeMatrixElement_IMA(int row, int col, double* values)
{
    // row, col are IMA element indices (reduced set)
    int full_row = ima_to_full[row];
    int full_col = ima_to_full[col];
    int mirror_col = mirror_map[full_col];

    // Compute both direct and mirror contributions
    double N_direct[36], N_mirror[36];
    Compute6x6BlockFast(full_row, full_col, N_direct);
    Compute6x6BlockFast(full_row, mirror_col, N_mirror);

    // Apply permutation and sum
    ApplyPermutation(N_mirror, perm_x, N_mirror);
    for(int k = 0; k < 36; k++)
        values[k] = N_direct[k] + N_mirror[k];
}
```

### 5. Element Position for Cluster Tree

HACApK cluster tree uses element positions. For IMA:
- Use the original element positions (not mirrored)
- The cluster tree is built on the reduced element set
- Distance calculation uses original positions

```cpp
TVector3d GetIMAElementCenter(int ima_idx)
{
    int full_idx = ima_to_full[ima_idx];
    return g3dRelaxPtrVect[full_idx]->ReturnCentrPoint();
}
```

## External Field Handling

For x-mirror IMA with coil:
- Coil also needs IMA treatment
- H_ext[i] = H_coil[full_i] (evaluated at original position)
- The mirror contribution is implicitly included via matrix IMA

## Validation

Compare with ELF_MAGIC x-mirror results:
1. Matrix should match exactly (rel diff < 0.01%)
2. Field at origin should match (-228 mT for mu=1000 case)

## Performance Considerations

### Dense Matrix (LU/BiCGSTAB)
- Memory: 4x reduction (N/2)^2 vs N^2
- Assembly: 2x more work per block (direct + mirror)
- Net: ~2x memory reduction

### HACApK
- Element count: N/2 (reduced)
- Each block computation: 2x (direct + mirror)
- ACA compression: May be slightly less effective due to image contributions
- Overall: Significant memory and time savings for large problems

## Implementation Priority

1. **Phase 1**: Dense matrix IMA (SetupInteractMatrix_IMA)
   - Implement for x-mirror only first
   - Validate against ELF_MAGIC

2. **Phase 2**: HACApK IMA (ComputeMatrixElement_IMA)
   - Modify RadHACApKMSCManager for IMA mode
   - Validate performance and accuracy

3. **Phase 3**: Multi-axis symmetry
   - Extend to xy, xz, yz, xyz symmetries
   - Quarter model (xy) and eighth model (xyz) support

## Files to Modify

1. `rad_interaction.h`: Add IMASymmetryConfig struct
2. `rad_interaction.cpp`: Add SetupInteractMatrix_IMA()
3. `rad_hacapk.h/cpp`: Add IMA support to ComputeMatrixElement
4. `radentry.cpp`: Add Python API for SetIMASymmetry
5. `radia_pybind.cpp`: Bind new functions

## Implementation Status (2026-01-30)

### Phase 1: Dense Matrix IMA - COMPLETED

**Files Modified:**
- `rad_interaction.h`: Added IMA enums, permutation arrays, and method declarations
- `rad_interaction.cpp`: Implemented SetIMASymmetry(), SetupInteractMatrix_IMA(), ApplyDOFPermutation()
- `rad_c_interface.cpp`: Added RadSetIMASymmetry() and RadBuildIMAMatrix() C interfaces
- `radentry.cpp`: Added Python API bindings
- `radia_pybind.cpp`: Added pybind11 bindings for SetIMASymmetry() and BuildIMAMatrix()

**Validation Results:**

| Model | Bz at Origin | Notes |
|-------|--------------|-------|
| Full model (52 elements) | -226.24 mT | Reference |
| IMA x-mirror (26 elements) | -226.24 mT | **Matches full model exactly** |
| ELF EIEM2 | -228.12 mT | Difference due to coil modeling |

The IMA implementation produces **identical results** to the full model (0.00% difference), confirming correct implementation.

**Python API Usage:**

```python
import radia as rad

# Build model with full geometry
hex_objects = [rad.ObjHexahedron(verts, [0,0,0]) for verts in all_vertices]
container = rad.ObjCnt(hex_objects + [coil])

# Setup Image x-mirror (symmetric BC)
intrc = rad.PreRelax(container, container)
n_ima = rad.Image(intrc, '+x')  # Returns number of Image elements
rad.BuildImageMatrix(intrc)    # Build reduced matrix

# Solve with half the DOFs
rad.Solve(container, 0.0001, 100, 0)
B = rad.Fld(container, 'b', [0, 0, 0])
```

**Example with antisymmetric BC:**

```python
# Z-mirror with antisymmetric BC (field normal to XY plane)
intrc = rad.PreRelax(container, container)
n_ima = rad.Image(intrc, '-z')  # Antisymmetric z-mirror
rad.BuildImageMatrix(intrc)
rad.Solve(container, 0.0001, 100, 0)
```

### Phase 2: HACApK Image - COMPLETED (2026-02-05)

HACApK IMA is implemented via the unified `image=` parameter in `rad.Solve()` and
`rad.BuildMatrix()`. The on-demand matrix element computation (`GetCached6x6Element`,
`GetCached3x3Element`) transparently includes IMA contributions through the
interaction matrix infrastructure.

**Key fix**: Thread-local cache invalidation for IMA transitions (2026-02-05).
Previously, stale cache values from non-IMA solves were incorrectly used for IMA
solves. Fixed by generation-based cache invalidation in `rad_hacapk.cpp`.

Verified: All three solvers (LU, BiCGSTAB, HACApK) produce identical results with
IMA quarter model (`image='+x-z'`) on the C-type electromagnet quarter-model benchmark.

### Phase 3: Multi-axis Symmetry - COMPLETED (2026-01-31)

Combined symmetries are supported via the unified `image=` string parameter:
- `image='+x'`: x-mirror only (half model)
- `image='+x-z'`: x + z mirror (quarter model)
- `image='+x+y-z'`: x + y + z mirror (eighth model)

Sign selection policy:
- `+` for field parallel to mirror plane (symmetric)
- `-` for field perpendicular to mirror plane (antisymmetric)

## References

- Yano MSC method (EIEM2 evaluation points)
