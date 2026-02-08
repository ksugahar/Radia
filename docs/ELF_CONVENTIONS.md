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
