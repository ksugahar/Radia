# ELF Compatibility Guide

This document consolidates the conventions, mesh differences, matrix formulation, and verification results for achieving compatibility between Radia and ELF/MAGIC (ELF_MAGIC). It serves as the single reference for users who need to compare results or migrate between the two solvers.

**Acknowledgment**: The author thanks ELF Corporation (President: Dr. Yano) for permission to document this correspondence.

---

## 1. Unit Systems & Conventions

### 1.1 Magnetization M: A/m (Amperes per meter)

**CRITICAL**: Both ELF and Radia use **M in A/m** (same units as H), NOT Tesla.

| Quantity | Symbol | SI Unit | Notes |
|----------|--------|---------|-------|
| Magnetic field strength | H | A/m | Applied + demagnetizing field |
| Magnetization | M | A/m | M = chi * H |
| Magnetic flux density | B | Tesla (T) | B = mu_0 * (H + M) |
| Susceptibility | chi | dimensionless | chi = mu_r - 1 |
| Relative permeability | mu_r | dimensionless | mu_r = B / (mu_0 * H) |
| Permeability of vacuum | mu_0 | H/m | mu_0 = 4*pi*1e-7 |

**Key Relationships**:
```
M = chi * H           (M in A/m, H in A/m, chi dimensionless)
B = mu_0 * (H + M)    (B in Tesla)
mu_r = B / (mu_0 * H) = 1 + chi
chi = mu_r - 1
```

**Radia Python API**: `rad.MatLin()` takes **mu_r** (relative permeability, >= 1.0).
Radia internally converts to chi = mu_r - 1 before passing to the C++ solver.
```python
mat = rad.MatLin(1000)           # mu_r = 1000, internally chi = 999
mat = rad.MatLin([5000, 100], [0, 0, 1])  # anisotropic, easy axis in z
```

**Common Mistake**: Do NOT confuse M with J (magnetic polarization).
- M (magnetization): A/m
- J (magnetic polarization): J = mu_0 * M, in Tesla

### 1.2 Why M in A/m (not Tesla)?

1. **Dimensional consistency**: M = chi * H requires M and H to have same units
2. **Physical meaning**: M represents magnetic dipole moment per unit volume
3. **Standard SI**: SI system defines M in A/m
4. **ELF compatibility**: ELF uses M in A/m internally

**Conversion** (if needed):
```python
# If you have M in Tesla (actually J = magnetic polarization):
M_Am = M_Tesla / mu_0    # Convert J to M
# M_Am is in A/m

# Example: NdFeB with Br = 1.2 T
# J = Br = 1.2 T (magnetic polarization)
# M = J / mu_0 = 1.2 / (4*pi*1e-7) = 954930 A/m
```

### 1.3 BH Curve Format

Both ELF and Radia use **BH curve** (not MH curve) for nonlinear materials:

| Column | Quantity | Unit |
|--------|----------|------|
| 1 | H (field) | A/m |
| 2 | B (flux density) | Tesla |

**Example BH curve data** (Steel 1008):
```
# H [A/m]    B [T]
0.0          0.0
100.0        0.1
500.0        0.8
1000.0       1.2
5000.0       1.7
50000.0      2.0
```

**Chi computation from BH curve**:
```python
# Given BH curve point (H, B):
mu_r = B / (mu_0 * H)    # Relative permeability
chi = mu_r - 1           # Susceptibility
M = chi * H              # Magnetization in A/m
```

### 1.4 Nonlinear Solver Methods (ELF/MAGIC Compatibility)

ELF/MAGIC uses two distinct nonlinear solution methods:

| Method | Name | Description |
|--------|------|-------------|
| mucal1 | Picard Method | Fixed-point iteration with dual estimation for stability |
| mucal2 | Newton Method | Newton-Raphson using differential susceptibility |

#### mucal1: Picard Method (Relaxation)

**ELF Implementation**: `magic.f90` subroutine `mucal1(it,ccon)` (line 4181)

**Radia Implementation**: Default nonlinear iteration in `rad_relaxation_methods.cpp`

The Picard method (fixed-point iteration) updates the permeability using two estimation methods and selects the one with smaller change for stability:

**Method 1** (un1): H-value interpolation from B-H curve
```
mu_1 = B(H) / (mu_0 * H)
```
where B is interpolated from the B-H curve at the current H value.

**Method 2** (un2): (H+B)-sum based interpolation
```
# Find B-H curve point where H_curve + B_curve ≈ H_current + B_current
# Then: mu_2 = B_interp / (mu_0 * H_interp)
```

**Selection**: Apply under-relaxation alpha = 0.5 to both estimates:
```
mu_k_relaxed = (1-alpha) * mu_k + alpha * mu_old,  k=1,2
```
Select the estimate with smaller |mu_new - mu_old| to improve convergence stability.

**Convergence**: Based on relative permeability change
```
max |mu_new - mu_old| / mu_old < tolerance
```

#### mucal2: Newton-Raphson Method

**ELF Implementation**: `magic.f90` subroutine `mucal2(it,ccon)` (line 4311)

**Radia Implementation**: Newton method in `rad_relaxation_methods.cpp` (activated by `use_newton` flag)

The Newton-Raphson method uses differential susceptibility from the B-H curve tangent:
```
chi_d = (dB/dH) / mu_0 - 1
```

Modified linear system:
```
[D(1/chi_d) + G] * sigma^(k+1) = H_ext + D(1/chi_d - 1/chi) * sigma^(k)
```

**Convergence**: Based on B-field change (faster convergence but requires good initial guess)
```
max |B_new - B_old| / B_sat < tolerance
```

#### Hybrid Strategy (Radia Default)

Radia uses a hybrid approach combining both methods:

1. **Initial 10 iterations**: Picard method (mucal1-like) for stable convergence
2. **Subsequent iterations**: Newton method (mucal2-like) for faster convergence

This is controlled by:
```cpp
const int newton_start_iter = 10;
bool newton_active = ctx.use_newton && iterCount >= newton_start_iter;
```

**Implementation note**: Radia implements Method 1 and Method 2 as separate options, whereas ELF's mucal1 combines both and selects the better estimate automatically. Both approaches achieve stable convergence.

#### Solver Correspondence Table

| Aspect | ELF/MAGIC | Radia |
|--------|-----------|-------|
| Picard (Method 1) | mucal1 - un1 | Method 1: mu_1 = B(H)/(mu_0*H) |
| Picard (Method 2) | mucal1 - un2 | Method 2: (H+B)-sum interpolation |
| Newton method | mucal2 | Newton with chi_d = (dB/dH)/mu_0 - 1 |
| Switching logic | newton flag | newton_start_iter = 10 |
| Convergence (Picard) | mu change | B-field change |
| Convergence (Newton) | B change | B-field change |

---

## 2. Mesh Convention Differences (Hexahedra)

### 2.1 CHEXA (8-node Hexahedron) Element

#### Vertex Ordering (Nastran Convention)

Both solvers use the **standard CHEXA vertex ordering** (Netgen-compatible, Nastran CHEXA equivalent):

```
        v7 -------- v6
        /|          /|
       / |         / |
     v4 -------- v5  |
      |  v3 -----|-- v2
      | /        | /
      |/         |/
     v0 -------- v1

Vertex positions (for an axis-aligned cube):
  v0: (x-, y-, z-)  bottom-left-front
  v1: (x+, y-, z-)  bottom-right-front
  v2: (x+, y+, z-)  bottom-right-back
  v3: (x-, y+, z-)  bottom-left-back
  v4: (x-, y-, z+)  top-left-front
  v5: (x+, y-, z+)  top-right-front
  v6: (x+, y+, z+)  top-right-back
  v7: (x-, y+, z+)  top-left-back

Local coordinate system:
  local_x: v0 -> v1  (edge from vertex 0 to vertex 1)
  local_y: v0 -> v3  (edge from vertex 0 to vertex 3)
  local_z: v0 -> v4  (edge from vertex 0 to vertex 4)
```

**Note**: In ELF MEG files, node IDs may be listed in a different order, but the element connectivity maps back to the standard CHEXA ordering.

#### Nastran CHEXA Format (Reference Only)

**Note**: Nastran BDF import is removed from Radia. This section is kept as a reference for understanding ELF's element conventions. Use Coreform Cubit for mesh generation.

```
CHEXA  EID     PID     G1      G2      G3      G4      G5      G6      +
+      G7      G8
```

Example:
```
CHEXA  1       1       1       2       3       4       5       6       +
+      7       8
```

- EID: Element ID
- PID: Property ID
- G1-G8: Grid point (node) IDs for 8 vertices

### 2.2 Face Ordering

**ELF and Radia use different face orderings.** ELF identifies faces by **topology** (which vertices form the face), NOT by geometry (face normal direction). For sheared/rotated elements, the actual normal direction may not align with any coordinate axis.

#### ELF Face Ordering (kkh array)

**Source**: `magic.f90` lines 1122-1127 (kkh array)

| Face Index | Vertices (1-indexed) | Vertices (0-indexed) | Normal direction (axis-aligned case) |
|------------|---------------------|---------------------|--------------------------------------|
| 0          | 1, 2, 6, 5 | 0, 1, 5, 4 | y- (front) |
| 1          | 2, 3, 7, 6 | 1, 2, 6, 5 | x+ (right) |
| 2          | 3, 4, 8, 7 | 2, 3, 7, 6 | y+ (back) |
| 3          | 4, 1, 5, 8 | 3, 0, 4, 7 | x- (left) |
| 4          | 4, 3, 2, 1 | 3, 2, 1, 0 | z- (bottom) |
| 5          | 5, 6, 7, 8 | 4, 5, 6, 7 | z+ (top) |

#### Radia Face Ordering (NETGEN_FACES)

Defined in `radia_pybind.cpp`:

```cpp
static const int NETGEN_FACES[6][4] = {
    {1, 2, 3, 4},  // Face 0: z- (bottom) = {v0,v1,v2,v3}
    {2, 6, 7, 3},  // Face 1: x+ (right)  = {v1,v5,v6,v2}
    {1, 5, 6, 2},  // Face 2: y- (front)  = {v0,v4,v5,v1}
    {1, 4, 8, 5},  // Face 3: x- (left)   = {v0,v3,v7,v4}
    {3, 7, 8, 4},  // Face 4: y+ (back)   = {v2,v6,v7,v3}
    {5, 8, 7, 6}   // Face 5: z+ (top)    = {v4,v7,v6,v5}
};
// Note: NETGEN_FACES uses 1-indexed vertex IDs
```

#### Face Correspondence Table

| Radia Face | Radia Direction | Radia Vertices | ELF Face | ELF Direction | ELF Vertices (kkh) |
|------------|----------------|----------------|----------|--------------|---------------------|
| 0 | z- (bottom) | {v0,v1,v2,v3} | 4 | z- | {v3,v2,v1,v0} |
| 1 | x+ (right)  | {v1,v5,v6,v2} | 1 | x+ | {v1,v2,v6,v5} |
| 2 | y- (front)  | {v0,v4,v5,v1} | 0 | y- | {v0,v1,v5,v4} |
| 3 | x- (left)   | {v0,v3,v7,v4} | 3 | x- | {v3,v0,v4,v7} |
| 4 | y+ (back)   | {v2,v6,v7,v3} | 2 | y+ | {v2,v3,v7,v6} |
| 5 | z+ (top)    | {v4,v7,v6,v5} | 5 | z+ | {v4,v5,v6,v7} |

#### Permutation Vector and Matrix

Radia-to-ELF face mapping:
```
Permutation: Radia face -> ELF face
  0 -> 4  (z-)
  1 -> 1  (x+)
  2 -> 0  (y-)
  3 -> 3  (x-)
  4 -> 2  (y+)
  5 -> 5  (z+)

Permutation vector: radia_to_elf = [4, 1, 0, 3, 2, 5]

Reordering ELF diagonal values to Radia order:
  radia_ordered[i] = elf_diag[radia_to_elf[i]]

Permutation matrix P:
[[0 0 0 0 1 0]
 [0 1 0 0 0 0]
 [1 0 0 0 0 0]
 [0 0 0 1 0 0]
 [0 0 1 0 0 0]
 [0 0 0 0 0 1]]
```

**Important**: The permutation depends only on the face naming convention, not on the element shape (cube, cuboid, or sheared hexahedron).

### 2.3 DOF Mapping

For 6-DOF MSC (Magnetic Surface Charge) hexahedra:
- DOF 0 = sigma on face 0
- DOF 1 = sigma on face 1
- DOF 2 = sigma on face 2
- DOF 3 = sigma on face 3
- DOF 4 = sigma on face 4
- DOF 5 = sigma on face 5

Face ordering differences directly affect DOF ordering. When comparing results between solvers, apply the permutation matrix.

### 2.4 Matrix Storage Order

| Language | Storage Order | A[i,j] memory location |
|----------|--------------|------------------------|
| Fortran (ELF) | Column-major | A[j*M + i] |
| C++ (Radia) | Row-major | A[i*N + j] |

When reading ELF `.mat` files in Python/NumPy, a transpose is required:
```python
# Correct loading procedure
elf_matrix_raw = read_elf_matrix(mat_file)
elf_matrix = elf_matrix_raw.T  # Transpose for use in row-major context
```

Without the transpose:
- Diagonal elements are correct
- Off-diagonal elements are transposed
- The entire matrix ends up transposed

### 2.5 Netgen/NGSolve Compatibility

Netgen generates tetrahedral meshes by default. For hexahedral meshes:
1. Use Coreform Cubit to generate hex mesh
2. Export via `cubit_mesh_export` to Netgen format
3. Import into Radia using `ObjHexahedron` with faces in ELF order

Face ordering verification:
```python
import radia as rad
import numpy as np

# ELF kkh face ordering (0-indexed vertices)
CHEXA_FACES = [
    [0, 1, 5, 4],  # Face 0
    [1, 2, 6, 5],  # Face 1
    [2, 3, 7, 6],  # Face 2
    [3, 0, 4, 7],  # Face 3
    [3, 2, 1, 0],  # Face 4
    [4, 5, 6, 7],  # Face 5
]

# Verify diagonal block matches ELF
# If diagonal values are in same positions, face ordering is correct
```

---

## 3. Matrix Formulation & Compatibility

### 3.1 Sign Convention

Both Radia (v1.4.4 and later) and ELF use **negative diagonal elements**.

| Property | Radia | ELF |
|----------|-------|-----|
| Diagonal sign | Negative (-1.5684) | Negative (-1.5684) |
| Matrix definition | N = interaction matrix | A = N - diag(1/chi), chi = mu_r - 1 |

**Physical interpretation**:
- **Radia**: Negative diagonal (ELF-compatible; self-demagnetizing field reduces magnetization)
- **ELF**: Negative diagonal (energy formulation)

**Note**: Radia versions prior to v1.4.3 used positive diagonal elements.

### 3.2 Matrix Formulation History (Radia Refactoring)

#### Previous Radia formulation (pre-v1.4.4)

```
Matrix equation: A . M = H_ext
Matrix definition: A = -N + diag(+1/chi)

Diagonal elements:   A_ii = +1/chi - N_ii  (positive values)
Off-diagonal elements: A_ij = -N_ij
```

#### Current formulation (v1.4.4+, ELF-compatible)

```
Matrix equation: A . M = H_ext
Matrix definition: A = -N + diag(-1/chi)

Diagonal elements:   A_ii = -1/chi - N_ii  (negative values)
Off-diagonal elements: A_ij = -N_ij
```

#### Mathematical relationship between the two

```
A_ELF = A_Radia_old - 2 * diag(1/chi)
```

Or element-wise:
```
A_ELF[i,i] = A_Radia_old[i,i] - 2/chi_i
A_ELF[i,j] = A_Radia_old[i,j]  (i != j)
```

### 3.3 Interaction Matrix Storage

ELF stores the interaction matrix in **column-major** format (Fortran convention).

For conversion between ELF and Radia:
```python
# A_ELF[i,j] = N_Radia[j,i]^T for each 6x6 block
radia_block = radia_matrix[i*6:(i+1)*6, j*6:(j+1)*6]
elf_block = radia_block.T  # Transpose for column-major storage
```

### 3.4 Full Conversion Formula

ELF matrix A and Radia matrix N are related by:

```
A_ELF_raw = P @ N_Radia^T @ P^T

where:
  A_ELF_raw = Matrix read directly from ELF file (no transpose)
  N_Radia^T = Transpose of Radia's interaction matrix
  P         = Permutation matrix (see Section 2.2)
  @         = Matrix multiplication
```

**Note**: When using the transposed ELF matrix: `A_ELF_T = P @ N_Radia @ P^T`.

### 3.5 Affected Source Files

#### Core files

| File | Change |
|------|--------|
| `rad_interaction.cpp` | Diagonal sign in matrix construction |
| `rad_interaction.h` | Comment updates |
| `rad_hacapk.cpp` | Sign in UpdateDiagonal() |
| `rad_hacapk.h` | Comment updates |
| `rad_relaxation_methods.cpp` | Solver adjustments |

#### Test files

| File | Change |
|------|--------|
| `tests/test_msc_matrix.py` | Expected value sign changes |
| `examples/electromagnet/mu=1000/verify_matrix.py` | ELF comparison test |

### 3.6 Implementation Details

#### Diagonal sign change (rad_interaction.cpp)

```cpp
// Previous (pre-v1.4.4):
FlatInteractMatrix[diagIdx] = N_ii + 1.0 / chi;

// Current (ELF-compatible):
FlatInteractMatrix[diagIdx] = N_ii - 1.0 / chi;
```

#### UpdateDiagonal (rad_hacapk.cpp)

```cpp
void RadHACApKManager::UpdateDiagonal(const std::vector<double>& inv_chi) {
    // Previous: diag[i] = N_ii + inv_chi[i]
    // Current:  diag[i] = N_ii - inv_chi[i]  // ELF-compatible: -1/chi
}
```

#### Solver (rad_relaxation_methods.cpp)

RHS sign is unchanged (A . M = H_ext remains the same). Residual calculation and convergence checks verified:
- `r = H_ext - A . M`
- Jacobi preconditioner diagonal extraction
- Convergence criterion `|r| / |H_ext|`

### 3.7 Rollback Plan

A compile-time flag enables switching between formulations:

```cpp
#ifdef RADIA_ELF_COMPAT
    // ELF-compatible: negative diagonal
    diag = N_ii - inv_chi;
#else
    // Legacy: positive diagonal
    diag = N_ii + inv_chi;
#endif
```

Default is ELF-compatible mode (RADIA_ELF_COMPAT defined).

### 3.8 API Compatibility Impact

| Category | Impact |
|----------|--------|
| `GetInteractMatrix()` return values | **Breaking**: sign changed on diagonal |
| `Solve()` API | No change |
| `Fld()` computed results | No change (same physics) |
| Python user API | No change |

Versioning: minor version bump (v1.x.0 -> v1.(x+1).0) with changelog entry.

### 3.9 Python Conversion Utilities

```python
import numpy as np

def create_permutation_matrix():
    """Create Radia -> ELF permutation matrix."""
    perm = [4, 1, 0, 3, 2, 5]  # Radia face -> ELF face
    P = np.zeros((6, 6))
    for radia_i, elf_i in enumerate(perm):
        P[elf_i, radia_i] = 1.0
    return P

def radia_to_elf_matrix(radia_N):
    """Convert Radia N matrix to ELF A matrix format (raw)."""
    P = create_permutation_matrix()
    # Conversion: A_raw = P @ N^T @ P^T (v1.4.4+, sign convention is identical)
    A_elf_raw = P @ radia_N.T @ P.T
    return A_elf_raw

def elf_raw_to_radia_matrix(elf_A_raw):
    """Convert ELF raw A matrix to Radia N matrix format."""
    P = create_permutation_matrix()
    # Inverse: N^T = P^T @ A_raw @ P
    # Therefore: N = (P^T @ A_raw @ P)^T = P^T @ A_raw^T @ P
    N_radia = P.T @ elf_A_raw.T @ P
    return N_radia

def read_elf_matrix_and_convert(mat_file):
    """Read ELF file and convert to Radia format."""
    from verify_elf_radia import load_elf_matrix as read_elf_matrix
    elf_raw = read_elf_matrix(mat_file)
    # Convert each 6x6 block
    n = elf_raw.shape[0]
    n_elem = n // 6
    radia_N = np.zeros_like(elf_raw)
    for i in range(n_elem):
        for j in range(n_elem):
            i0, i1 = i*6, (i+1)*6
            j0, j1 = j*6, (j+1)*6
            radia_N[i0:i1, j0:j1] = elf_raw_to_radia_matrix(
                elf_raw[i0:i1, j0:j1])
    return radia_N
```

---

## 4. Verification Against ELF Results

### 4.1 Summary of Key Differences

| Property | ELF (Fortran) | Radia (C++) |
|----------|---------------|-------------|
| Matrix storage order | Column-major | Row-major |
| Diagonal sign | Negative | Negative (v1.4.4+, identical) |
| Face ordering | Face 0-5 (kkh array) | Face 0-5 (NETGEN_FACES) |

### 4.2 Single-Element Verification

**Single cube element** (100 mm side, mu_r = 1000):
- Maximum absolute difference: 1.7e-05
- Conclusion: agreement within numerical precision

**Single cuboid element** (100 x 90 x 120 mm, mu_r = 1000):
- Maximum absolute difference: 2.8e-05
- Conclusion: agreement within numerical precision

### 4.3 Full Model Verification

**Full model** (52 elements, cuboid):
- Maximum absolute error: 3.59e-04
- Maximum relative error: 0.0053%
- Conclusion: agreement within numerical integration precision

### 4.4 V304 Model Verification (Distorted Hexahedra)

**V304 model** (74 elements, distorted hexahedra, EIEM2 mesh):
- Diagonal comparison (reordered with radia_to_elf = [4, 1, 0, 3, 2, 5]):
  - Elements with max error < 1%: 32 of 74 (43%)
  - Elements with max error < 5%: 48 of 74 (65%)
  - Diagonal sum relative difference: 0.7-1.7%
- Primary source of residual: for face pairs with identical geometry (x-direction extrusion), Radia produces identical values while ELF produces different values
- Conclusion: topological face ordering verified successfully; residuals are caused by implementation differences (e.g., evaluation point computation)

### 4.5 Cosine Similarity (Matrix-Level)

After the diagonal sign change (Phase 1):
- **Cosine similarity: 0.918** (previously 0.35)
- Optimal reordering: `(2, 4, 3, 5, 0, 1)`

### 4.6 Remaining Sources of Discrepancy

Reasons the matrices are not perfectly identical:
1. **Nastran file vertex ordering**: Vertex ordering defined in ELF Nastran files may differ from standard CHEXA depending on mesh generation tool
2. **Evaluation point differences**: Different choices of integration evaluation points on faces
3. **Numerical precision**: Different implementations of solid angle calculation

### 4.7 Verification Test Plan

#### Unit tests

| Test | Content | Expected Result |
|------|---------|-----------------|
| Single cube | 10 cm cube, mu_r = 1000 | Diagonal elements are negative |
| Symmetry | Matrix A = A^T | Symmetric matrix |
| ELF comparison | Same mesh | Match after permutation |

#### Integration tests

| Test | Content | Expected Result |
|------|---------|-----------------|
| Linear analysis | mu_r = 1000 electromagnet | Same magnetic field as ELF |
| Nonlinear analysis | BH curve | Convergence, match with ELF |
| TrfMlt | 1/4 model symmetry | Match with full model |

#### ELF verification test code

```python
# verify_matrix.py

def compare_matrices(radia_matrix, elf_matrix):
    # After v1.4.4: permutation only (no sign flip needed)

    # Verify diagonal elements are negative
    assert np.all(np.diag(radia_matrix) < 0), "Diagonal must be negative"

    # Direct comparison with ELF (after permutation)
    P = find_permutation(radia_matrix, elf_matrix)
    radia_permuted = P @ radia_matrix @ P.T

    np.testing.assert_allclose(radia_permuted, elf_matrix, rtol=1e-4)
```

### 4.8 Recommendations

#### For Radia users

1. Radia follows the Netgen standard
2. When comparing results with ELF, apply the conversion described in Sections 2 and 3
3. Face ordering differences affect DOF (degrees of freedom) ordering

#### For future compatibility

- Radia prioritizes integration with Netgen/NGSolve
- Conversion to ELF format is provided via Python utilities
- New development should use Radia (Netgen) conventions

---

## References

- ELF/MAGIC source code: `S:\ELF_MAGIC\01_GitHub\src\legacy\magic.f90`
- Radia source code: `S:\radia\01_GitHub\src\core\rad_relaxation_methods.cpp`
- Paper reference: SA-26-010 (IEEJ Technical Meeting on Static Apparatus and Rotating Machinery, 2026)
- Analysis script: `examples/electromagnet/mu=1000/analyze_elf_face_ordering.py`
- Single element comparison: `examples/electromagnet/mu=1000/compare_single_element.py`
- ELF MEG file example: `S:\ELF_MAGIC\...\single\ELF_MAGIC.meg`

## Change History

| Date | Change |
|------|--------|
| 2025-01-30 | Initial creation (refactoring plan). |
| 2026-01-30 | Face ordering and sign convention differences identified. Diagonal sign changed to ELF-compatible (positive to negative). Sign inversion removed from conversion formula. |
| 2026-02-08 | Face ordering corrected: old [4,5,0,2,1,3] -> correct [4,1,0,3,2,5]. V304 mesh verification results added. NETGEN_FACES definition source added. |
| 2026-02-19 | Consolidated from ELF_CONVENTIONS.md, ELF_RADIA_MESH_COMPARISON.md, and REFACTOR_ELF_COMPATIBLE_MATRIX.md into single guide. |
