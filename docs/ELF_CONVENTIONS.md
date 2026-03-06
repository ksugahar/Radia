# ELF Element Conventions

This document describes the ELF (ELF_MAGIC) element conventions for compatibility with Radia.

## Unit System and Physical Quantities

### Magnetization M: A/m (Amperes per meter)

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

### BH Curve Format

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

### Why M in A/m (not Tesla)?

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

## CHEXA (8-node Hexahedron) Element

### Vertex Ordering (1-indexed, Nastran convention)

```
        8--------7
       /|       /|
      / |      / |
     5--------6  |
     |  4-----|--3
     | /      | /
     |/       |/
     1--------2

Local coordinate system:
- local_x: v0 -> v1  (edge from vertex 1 to vertex 2)
- local_y: v0 -> v3  (edge from vertex 1 to vertex 4)
- local_z: v0 -> v4  (edge from vertex 1 to vertex 5)
```

### Face Ordering (ELF Convention - kkh array)

ELF uses the following face ordering for MSC (Magnetic Surface Charge) method.
**Source**: magic.f90 lines 1122-1127 (kkh array)

| Face Index | Vertices (1-indexed) | Vertices (0-indexed) |
|------------|---------------------|---------------------|
| 0          | 1, 2, 6, 5 | 0, 1, 5, 4 |
| 1          | 2, 3, 7, 6 | 1, 2, 6, 5 |
| 2          | 3, 4, 8, 7 | 2, 3, 7, 6 |
| 3          | 4, 1, 5, 8 | 3, 0, 4, 7 |
| 4          | 4, 3, 2, 1 | 3, 2, 1, 0 |
| 5          | 5, 6, 7, 8 | 4, 5, 6, 7 |

**Key Point**: Face ordering follows kkh array definition (faces 0-5 as defined by vertex topology)

### DOF Mapping

For 6-DOF MSC hexahedra:
- DOF 0 = sigma on face 0
- DOF 1 = sigma on face 1
- DOF 2 = sigma on face 2
- DOF 3 = sigma on face 3
- DOF 4 = sigma on face 4
- DOF 5 = sigma on face 5

### Face Identification

**IMPORTANT**: ELF identifies faces by **topology** (which vertices form the face), NOT by **geometry** (face normal direction).

For sheared/rotated elements:
- Face 3 is ALWAYS the face formed by vertices {0, 4, 7, 3} (from kkh array)
- The actual normal direction may not align with any coordinate axis
- This topological definition ensures consistent DOF mapping regardless of element deformation

### Interaction Matrix Storage

ELF stores the interaction matrix in **column-major** format (Fortran convention).

For conversion between ELF and Radia:
```python
# A_ELF[i,j] = N_Radia[j,i]^T for each 6x6 block
radia_block = radia_matrix[i*6:(i+1)*6, j*6:(j+1)*6]
elf_block = radia_block.T  # Transpose for column-major storage
```

## Nastran CHEXA Format (Reference Only)

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

## Netgen/NGSolve Compatibility

Netgen generates tetrahedral meshes by default. For hexahedral meshes:
1. Use Coreform Cubit to generate hex mesh
2. Export via cubit_mesh_export to Netgen format
3. Import into Radia using ObjHexahedron with faces in ELF order

### Face Ordering Verification

To verify face ordering matches ELF kkh array:
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

## Nonlinear Solver Methods (ELF/MAGIC Compatibility)

This section describes the correspondence between Radia's nonlinear solver implementation and ELF/MAGIC's solver methods (mucal1, mucal2).

**Acknowledgment**: The author thanks ELF Corporation (President: Dr. Yano) for permission to document this correspondence.

### Method Overview

ELF/MAGIC uses two distinct nonlinear solution methods:

| Method | Name | Description |
|--------|------|-------------|
| mucal1 | Picard Method | Fixed-point iteration with dual estimation for stability |
| mucal2 | Newton Method | Newton-Raphson using differential susceptibility |

### mucal1: Picard Method (Relaxation)

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

**Selection**: Apply under-relaxation α = 0.5 to both estimates:
```
mu_k_relaxed = (1-α) * mu_k + α * mu_old,  k=1,2
```
Select the estimate with smaller |mu_new - mu_old| to improve convergence stability.

**Convergence**: Based on relative permeability change
```
max |mu_new - mu_old| / mu_old < tolerance
```

### mucal2: Newton-Raphson Method

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

### Hybrid Strategy (Radia Default)

Radia uses a hybrid approach combining both methods:

1. **Initial 10 iterations**: Picard method (mucal1-like) for stable convergence
2. **Subsequent iterations**: Newton method (mucal2-like) for faster convergence

This is controlled by:
```cpp
const int newton_start_iter = 10;
bool newton_active = ctx.use_newton && iterCount >= newton_start_iter;
```

**Implementation note**: Radia implements Method 1 and Method 2 as separate options, whereas ELF's mucal1 combines both and selects the better estimate automatically. Both approaches achieve stable convergence.

### Correspondence Table

| Aspect | ELF/MAGIC | Radia |
|--------|-----------|-------|
| Picard (Method 1) | mucal1 - un1 | Method 1: μ₁ = B(H)/(μ₀H) |
| Picard (Method 2) | mucal1 - un2 | Method 2: (H+B)-sum interpolation |
| Newton method | mucal2 | Newton with χ_d = (dB/dH)/μ₀ - 1 |
| Switching logic | newton flag | newton_start_iter = 10 |
| Convergence (Picard) | μ change | B-field change |
| Convergence (Newton) | B change | B-field change |

### References

- ELF/MAGIC source code: `S:\ELF_MAGIC\01_GitHub\src\legacy\magic.f90`
- Radia source code: `S:\radia\01_GitHub\src\core\rad_relaxation_methods.cpp`
- Paper reference: SA-26-010 (IEEJ Technical Meeting on Static Apparatus and Rotating Machinery, 2026)
