# ELF Element Conventions

This document describes the ELF (ELF_MAGIC) element conventions for compatibility with Radia.

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

### Face Ordering (ELF Convention)

ELF uses the following face ordering for MSC (Magnetic Surface Charge) method.
**Source**: magic.f90 lines 1122-1127 (kkh array)

| Index | Face Name | Description | Vertices (1-indexed) | Vertices (0-indexed) |
|-------|-----------|-------------|---------------------|---------------------|
| 0     | y-        | y-minimum face | 1, 2, 6, 5 | 0, 1, 5, 4 |
| 1     | x+        | x-maximum face | 2, 3, 7, 6 | 1, 2, 6, 5 |
| 2     | y+        | y-maximum face | 3, 4, 8, 7 | 2, 3, 7, 6 |
| 3     | x-        | x-minimum face | 4, 1, 5, 8 | 3, 0, 4, 7 |
| 4     | z-        | z-minimum face (bottom) | 4, 3, 2, 1 | 3, 2, 1, 0 |
| 5     | z+        | z-maximum face (top) | 5, 6, 7, 8 | 4, 5, 6, 7 |

**Key Point**: Face ordering is `[y-, x+, y+, x-, z-, z+]` (NOT `[y-, x-, y+, x+, z-, z+]`)

### DOF Mapping

For 6-DOF MSC hexahedra:
- DOF 0 = sigma on y- face
- DOF 1 = sigma on x- face
- DOF 2 = sigma on y+ face
- DOF 3 = sigma on x+ face
- DOF 4 = sigma on z- face
- DOF 5 = sigma on z+ face

### Face Identification

**IMPORTANT**: ELF identifies faces by **topology** (which vertices form the face), NOT by **geometry** (face normal direction).

For sheared/rotated elements:
- Face 1 (x-) is ALWAYS the face formed by vertices {0, 4, 7, 3}
- The actual normal direction may not point exactly in -x_local
- This topological definition ensures consistent DOF mapping regardless of element deformation

### Interaction Matrix Storage

ELF stores the interaction matrix in **column-major** format (Fortran convention).

For conversion between ELF and Radia:
```python
# A_ELF[i,j] = N_Radia[j,i]^T for each 6x6 block
radia_block = radia_matrix[i*6:(i+1)*6, j*6:(j+1)*6]
elf_block = radia_block.T  # Transpose for column-major storage
```

## Nastran CHEXA Format

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

To verify face ordering matches ELF:
```python
import radia as rad
import numpy as np

# Create hexahedron with ELF face ordering
CHEXA_FACES = [
    [0, 1, 5, 4],  # y- face
    [0, 4, 7, 3],  # x- face
    [2, 3, 7, 6],  # y+ face
    [1, 2, 6, 5],  # x+ face
    [0, 3, 2, 1],  # z- face
    [4, 5, 6, 7],  # z+ face
]

# Verify diagonal block matches ELF
# If diagonal values are in same positions, face ordering is correct
```
