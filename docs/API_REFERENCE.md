# Radia Python API Reference

Complete reference for Radia Python API.

**Version**: 1.3.17
**Date**: 2025-12-30
**Original ESRF Documentation**: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html

---

## Table of Contents

- [Quick Start](#quick-start)
- [Supported Elements](#supported-elements)
- [Geometry Objects](#geometry-objects)
- [Materials](#materials)
- [Solver](#solver)
- [Field Computation](#field-computation)
- [Mesh Import](#mesh-import)
- [NGSolve Integration](#ngsolve-integration)
- [Utilities](#utilities)

---

## Quick Start

### MSC Hexahedral Example (ObjThckPgn)

```python
import radia as rad
import numpy as np

rad.FldUnits('m')
rad.UtiDelAll()

MU_0 = 4 * np.pi * 1e-7
n_div = 5
cube_size = 1.0
elem_size = cube_size / n_div

# Create 5x5x5 hexahedral mesh using ObjThckPgn
elements = []
for ix in range(n_div):
    for iy in range(n_div):
        for iz in range(n_div):
            cx = (ix + 0.5) * elem_size - cube_size / 2
            cy = (iy + 0.5) * elem_size - cube_size / 2
            cz = (iz + 0.5) * elem_size - cube_size / 2
            half = elem_size / 2

            polygon = [[cx-half, cy-half], [cx+half, cy-half],
                       [cx+half, cy+half], [cx-half, cy+half]]
            obj = rad.ObjThckPgn(cz - half, elem_size, polygon, 'z', [0, 0, 0])
            elements.append(obj)

container = rad.ObjCnt(elements)
mat = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(container, mat)

ext = rad.ObjBckg([0, 0, MU_0 * 50000])
grp = rad.ObjCnt([container, ext])
rad.Solve(grp, 0.001, 1000, 1)
```

### Tetrahedral Mesh Example (Netgen)

```python
import radia as rad
rad.FldUnits('m')

# Import NGSolve BEFORE radia modules
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

# Create tetrahedral mesh
cube = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube.mat('magnetic')
mesh = Mesh(OCCGeometry(cube).GenerateMesh(maxh=0.3))

# Import to Radia
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')
```

---

## Supported Elements

| Element Type | API | Faces | DOF | Use Case |
|--------------|-----|-------|-----|----------|
| **Extruded Polygon** | `ObjThckPgn()` | N-gon extruded | 3 | General prism shapes |
| **Hexahedron (MSC)** | `ObjHexahedron()` | 6 quad | 6 | Permanent magnets, soft iron |
| **Tetrahedron** | `ObjTetrahedron()` | 4 tri | 3 | Complex curved geometry |
| **Wedge/Prism** | `ObjPolyhdr()` + `WEDGE_FACES` | 5 | 3 | Hybrid meshes |
| **Pyramid** | `ObjPolyhdr()` + `PYRAMID_FACES` | 5 | 3 | Mesh transitions |
| **General** | `ObjPolyhdr()` | custom | 3-6 | Arbitrary polyhedra |

**DOF (Degrees of Freedom)**:
- **Hexahedra (6 faces)**: 6 DOF - Surface charge density (sigma) per face (MSC method)
- **Other elements (4-5 faces)**: 3 DOF - Magnetization vector (Mx, My, Mz)
- All meshes are expected to be generated externally (Netgen, GMSH, Cubit, etc.)

### Simplified APIs (Recommended)

```python
import radia as rad

# Tetrahedron: just provide 4 vertices (faces auto-generated)
tet_vertices = [[0,0,0], [1,0,0], [0.5,0.866,0], [0.5,0.289,0.816]]
tetra = rad.ObjTetrahedron(tet_vertices, [0, 0, 1e6])

# Hexahedron: just provide 8 vertices (faces auto-generated)
hex_vertices = [
    [-0.5,-0.5,-0.5], [0.5,-0.5,-0.5], [0.5,0.5,-0.5], [-0.5,0.5,-0.5],  # bottom
    [-0.5,-0.5,0.5], [0.5,-0.5,0.5], [0.5,0.5,0.5], [-0.5,0.5,0.5]       # top
]
hexa = rad.ObjHexahedron(hex_vertices, [0, 0, 1e6])
```

### Face Topology Constants (for advanced usage)

```python
from netgen_mesh_import import TETRA_FACES, HEX_FACES, WEDGE_FACES, PYRAMID_FACES

# TETRA_FACES (1-indexed) - used internally by ObjTetrahedron
[[1, 2, 3], [1, 4, 2], [2, 4, 3], [3, 4, 1]]

# HEX_FACES (1-indexed) - used internally by ObjHexahedron
[[1, 4, 3, 2], [5, 6, 7, 8], [1, 2, 6, 5], [3, 4, 8, 7], [1, 5, 8, 4], [2, 3, 7, 6]]
```

---

## Geometry Objects

### ObjThckPgn - Thick Polygon (Extruded 2D)

```python
obj = rad.ObjThckPgn(z_base, thickness, vertices_2d, axis, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `z_base` | float | Base position along extrusion axis |
| `thickness` | float | Extrusion length |
| `vertices_2d` | [[x,y], ...] | 2D polygon vertices (CCW) |
| `axis` | str | Extrusion axis: `'x'`, `'y'`, or `'z'` |
| `magnetization` | [Mx, My, Mz] | Initial magnetization |

```python
polygon = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]
hex_elem = rad.ObjThckPgn(-0.5, 1.0, polygon, 'z', [0, 0, 0])
```

### ObjTetrahedron - Tetrahedral Element (Recommended)

```python
obj = rad.ObjTetrahedron(vertices, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `vertices` | [[x,y,z], ...] | 4 vertex coordinates |
| `magnetization` | [Mx, My, Mz] | Initial magnetization (optional, default [0,0,0]) |

Creates a tetrahedron with face topology auto-generated internally.

```python
vertices = [[0,0,0], [1,0,0], [0.5,0.866,0], [0.5,0.289,0.816]]
tet = rad.ObjTetrahedron(vertices, [0, 0, 1e6])

# Without magnetization (for soft magnetic materials)
tet2 = rad.ObjTetrahedron(vertices)
```

**Vertex ordering**:
- v1, v2, v3: Base triangle (counter-clockwise from below)
- v4: Apex (top vertex)

### ObjHexahedron - Hexahedral Element (Recommended)

```python
obj = rad.ObjHexahedron(vertices, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `vertices` | [[x,y,z], ...] | 8 vertex coordinates |
| `magnetization` | [Mx, My, Mz] | Initial magnetization (optional, default [0,0,0]) |

Creates a hexahedron with face topology auto-generated internally.

```python
s = 0.5
vertices = [
    [-s,-s,-s], [s,-s,-s], [s,s,-s], [-s,s,-s],  # bottom face
    [-s,-s,s], [s,-s,s], [s,s,s], [-s,s,s]        # top face
]
hex_obj = rad.ObjHexahedron(vertices, [0, 0, 1e6])

# Without magnetization (for soft magnetic materials)
hex2 = rad.ObjHexahedron(vertices)
```

**Vertex ordering**:
```
       v8--------v7
      /|        /|
     / |       / |
    v5--------v6 |
    |  v4-----|--v3
    | /       | /
    |/        |/
    v1--------v2

Bottom (v1-v4): counter-clockwise from below
Top (v5-v8): directly above bottom vertices
```

### ObjPolyhdr - General Polyhedron

```python
obj = rad.ObjPolyhdr(vertices, faces, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `vertices` | [[x,y,z], ...] | 3D vertex coordinates |
| `faces` | [[v1,v2,...], ...] | Face vertex indices (**1-indexed!**) |
| `magnetization` | [Mx, My, Mz] | Initial magnetization |

Use `ObjPolyhdr` for wedge, pyramid, or custom polyhedra. For tetrahedra and hexahedra, prefer `ObjTetrahedron` and `ObjHexahedron`.

```python
from netgen_mesh_import import WEDGE_FACES
vertices = [[0,0,0], [1,0,0], [0.5,0.866,0], [0,0,1], [1,0,1], [0.5,0.866,1]]
wedge = rad.ObjPolyhdr(vertices, WEDGE_FACES, [0, 0, 1e6])
```

### ObjBckg - Uniform Background Field

```python
field_src = rad.ObjBckg([Bx, By, Bz])
```

```python
MU_0 = 4 * np.pi * 1e-7
ext = rad.ObjBckg([0, 0, MU_0 * 50000])  # 50,000 A/m in z
```

### ObjCnt - Container

```python
group = rad.ObjCnt([obj1, obj2, ...])
```

### ObjArcCur - Arc/Circular Coil

```python
coil = rad.ObjArcCur(center, radii, angles, height, n_sectors, j_azim)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `center` | [x,y,z] | Center of the arc/circle |
| `radii` | [r_min, r_max] | Inner and outer radii |
| `angles` | [phi_min, phi_max] | Start and end angles (rad) |
| `height` | float | Height of the coil cross-section |
| `n_sectors` | int | Number of azimuthal sectors |
| `j_azim` | float | Azimuthal current density (A/mm^2) |

```python
import numpy as np

# Full circular coil (R=50mm, thin cross-section)
center = [0, 0, 0]
radii = [49.5, 50.5]  # 1mm radial width
angles = [-np.pi, np.pi]  # Full circle
height = 1.0  # 1mm height
j_azim = 1000.0  # A/mm^2 (equivalent to 1000A total current)

coil = rad.ObjArcCur(center, radii, angles, height, 100, j_azim)
B = rad.Fld(coil, 'b', [0, 0, 50])  # Field on axis at z=50mm
```

**Analytical Method**: Uses elliptic integral formulas for high accuracy.
See [Elliptic Integral Formulas](#elliptic-integral-formulas-for-coils) for details.

### ObjRaceTrk - Racetrack Coil

```python
coil = rad.ObjRaceTrk(center, radii, heights, current, n_segments)
```

### ObjFlmCur - Filament Conductor (Line Current)

```python
filament = rad.ObjFlmCur([[x1,y1,z1], [x2,y2,z2], ...], current)
```

**Analytical Method**: Uses Biot-Savart law with closed-form solution.

---

## Materials

### MatLin - Linear Isotropic

```python
mat = rad.MatLin(mu_r)  # relative permeability
rad.MatApl(obj, mat)
```

```python
# Soft iron (mu_r = 1000)
mat = rad.MatLin(1000)
rad.MatApl(cube, mat)
```

### MatLin - Linear Anisotropic

```python
mat = rad.MatLin([mu_r_par, mu_r_perp], [ex, ey, ez])
```

```python
# Easy axis in z-direction
mat = rad.MatLin([5001, 101], [0, 0, 1])
```

### MatSatIsoTab - Nonlinear (B-H Table)

```python
mat = rad.MatSatIsoTab(BH_data)  # [[H, B], ...] in A/m and Tesla
```

**Input Format**: Industry-standard B-H curve (H in A/m, B in Tesla).
Radia internally converts to M-H using: M = B/mu_0 - H

```python
# B-H curve: [H (A/m), B (T)]
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
]

mat = rad.MatSatIsoTab(BH_DATA)
```

### MatSatIsoFrm - Nonlinear (Formula)

```python
mat = rad.MatSatIsoFrm([ksi1, ms1], [ksi2, ms2], [ksi3, ms3])
```

Formula: `M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3)`

```python
# Steel37 (C<0.13%)
mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])
```

### MatApl - Apply Material

```python
rad.MatApl(obj, material)
```

---

## Solver

### Solve - High-Level API (Recommended)

```python
result = rad.Solve(obj, tolerance, max_iter, method=1)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | int | Object or container |
| `tolerance` | float | Convergence threshold (0.001 = 0.1%) |
| `max_iter` | int | Maximum iterations |
| `method` | int | `0` = LU, `1` = BiCGSTAB (default) |

| Returns | Description |
|---------|-------------|
| `result[0]` | Final residual |
| `result[3]` | Number of iterations |

### Solver Selection

| Problem Size | Elements | Method | Code |
|--------------|----------|--------|------|
| Small | < 1,000 | LU | `rad.Solve(grp, 0.001, 100, 0)` |
| Medium | 1,000-10,000 | BiCGSTAB | `rad.Solve(grp, 0.001, 1000, 1)` |
| Large | > 10,000 | BiCGSTAB | `rad.Solve(grp, 0.001, 1000, 1)` |

**Iteration counts**:
- Linear materials: 1-2 iterations
- Nonlinear materials: 3-6 iterations (with B-field convergence)

### Nonlinear Convergence (v1.3.15+)

Radia uses **B-field based convergence** (mucal2) for nonlinear materials:

```
rel_change = |B_new - B_old| / B_sat
```

| Parameter | Description |
|-----------|-------------|
| `B_sat` | Saturation magnetization from BH curve |
| `tolerance` | Default 0.0001 (0.01% relative change) |

This method provides fast Newton-Raphson convergence and matches industry-standard solvers.

### Solver Tolerance Parameters

Radia provides three tolerance parameters for controlling solver behavior:

```python
# 1. Nonlinear iteration tolerance (outer loop)
#    Set via Solve() - controls when Newton-Raphson iterations stop
rad.Solve(obj, nonl_tol, max_iter, method)  # nonl_tol = 0.001 recommended

# 2. BiCGSTAB inner loop tolerance
#    Set via SetBiCGSTABTol() BEFORE Solve() - controls linear system accuracy
rad.SetBiCGSTABTol(bicg_tol)  # Default: 1e-4

# 3. H-matrix ACA tolerance (Method 2 only)
#    Set via SetHACApKParams() BEFORE Solve() - controls low-rank approximation
rad.SetHACApKParams(hmat_eps, leaf_size, eta)  # Default: 1e-4, 10, 2.0
```

| Parameter | API | Default | Description |
|-----------|-----|---------|-------------|
| `nonl_tol` | `rad.Solve(obj, nonl_tol, ...)` | 0.001 | Nonlinear convergence threshold |
| `bicg_tol` | `rad.SetBiCGSTABTol(tol)` | 1e-4 | BiCGSTAB relative residual tolerance |
| `hmat_eps` | `rad.SetHACApKParams(eps, ...)` | 1e-4 | H-matrix ACA compression tolerance |

**Example - Full solver configuration:**

```python
import radia as rad

# Configure tolerances BEFORE Solve()
rad.SetBiCGSTABTol(1e-4)           # BiCGSTAB tolerance
rad.SetHACApKParams(1e-4, 10, 2.0) # H-matrix: eps=1e-4, leaf=10, eta=2.0

# Solve with nonlinear tolerance
rad.Solve(grp, 0.001, 100, 2)      # nonl_tol=0.001, max_iter=100, method=2 (HACApK)
```

### SetBiCGSTABTol - BiCGSTAB Inner Loop Tolerance

```python
rad.SetBiCGSTABTol(tol)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `tol` | float | 1e-4 | Relative residual tolerance for BiCGSTAB |

**Notes:**
- Affects Method 1 (BiCGSTAB) and Method 2 (HACApK)
- Lower values = higher accuracy but more iterations
- Call BEFORE `rad.Solve()`

### SetHACApKParams - H-Matrix Parameters (Method 2)

```python
rad.SetHACApKParams(eps, leaf_size, eta)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `eps` | float | 1e-4 | ACA+ compression tolerance |
| `leaf_size` | int | 10 | Minimum cluster size in elements |
| `eta` | float | 2.0 | Admissibility parameter |

**Notes:**
- Only affects Method 2 (HACApK H-matrix solver)
- Lower `eps` = higher accuracy, larger ranks, more memory
- Call BEFORE `rad.Solve()`

**Parameter Rationale:**

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `eps` | 1e-4 | Balance between accuracy and compression. Lower values (1e-6, 1e-8) for higher accuracy, higher values (1e-3) for faster computation. |
| `leaf_size` | 10 | Minimum cluster size. Smaller values allow deeper tree but increase H-matrix overhead. 10 provides good balance for typical element counts. ELF-compatible default. |
| `eta` | 2.0 | Standard admissibility criterion: clusters are "well-separated" when `dist(c1,c2) >= eta * max(diam(c1), diam(c2))`. eta=2.0 is conservative, ensuring accurate low-rank approximations. Lower values (1.0) allow more aggressive compression but may reduce accuracy. |

### SetRelaxParam - Under-Relaxation Coefficient

```python
rad.SetRelaxParam(relax)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `relax` | float | 0.0 | Under-relaxation coefficient (0.0-1.0) |

**Notes:**
- Affects all solver methods (0=LU, 1=BiCGSTAB, 2=HACApK)
- `relax=0.0`: Full Newton step (default, fastest convergence when stable)
- `relax>0.0`: Damped update: `chi_new = chi_new*(1-relax) + chi_old*relax`
- Use under-relaxation (e.g., 0.2-0.5) when:
  - Convergence is slow or oscillating
  - Material has steep B-H curve
  - Problem is highly nonlinear
- Call BEFORE `rad.Solve()`

**Example:**
```python
# For difficult nonlinear problems, use under-relaxation
rad.SetRelaxParam(0.3)  # 30% damping
rad.Solve(container, 0.001, 100, 1)

# Reset to full step for normal cases
rad.SetRelaxParam(0.0)
```

### BiCGSTAB Performance

Typical solve times (nonlinear BH curve material):

| Elements | Time | Iterations |
|----------|------|------------|
| 1,000 | 0.55s | 5-6 |
| 3,375 | 7.30s | 5-6 |
| 8,000 | 51.81s | 5-6 |

---

## Field Computation

### Fld - Field at Point(s)

```python
field = rad.Fld(obj, component, point)
```

| Component | Description |
|-----------|-------------|
| `'bx'`, `'by'`, `'bz'`, `'b'` | Magnetic flux density B (T) |
| `'hx'`, `'hy'`, `'hz'`, `'h'` | Magnetic field H (A/m) |
| `'ax'`, `'ay'`, `'az'`, `'a'` | Vector potential A (T*m) |
| `'mx'`, `'my'`, `'mz'`, `'m'` | Magnetization M |

```python
B = rad.Fld(magnet, 'b', [0, 0, 0.1])  # B vector at point
Bz = rad.Fld(magnet, 'bz', [0, 0, 0.1])  # Bz component
H = rad.Fld(magnet, 'h', [0, 0, 0.1])  # H vector at point
```

### FldBatch - Batch Field Computation (v1.3.16+)

```python
result = rad.FldBatch(obj, points, method=0)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | int | Object or container |
| `points` | [[x,y,z], ...] | List of evaluation points |
| `method` | int | `0` = direct (default), `1` = FMM (future) |

| Returns | Description |
|---------|-------------|
| `result['B']` | List of [Bx, By, Bz] values (T) |
| `result['H']` | List of [Hx, Hy, Hz] values (A/m) |

```python
points = [[0, 0, 0.1], [0, 0, 0.2], [0, 0, 0.3]]
result = rad.FldBatch(magnet, points)
B_list = result['B']  # [[Bx1,By1,Bz1], [Bx2,By2,Bz2], ...]
H_list = result['H']  # [[Hx1,Hy1,Hz1], [Hx2,Hy2,Hz2], ...]
```

**Note**: More efficient than calling `Fld()` in a loop for many points.

### FldPhi - Scalar Potential Batch (v1.3.16+)

```python
phi = rad.FldPhi(obj, points)
```

Computes magnetic scalar potential phi_m at multiple points.

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | int | Object or container |
| `points` | [[x,y,z], ...] | List of evaluation points |

| Returns | Description |
|---------|-------------|
| `phi` | List of scalar values (A) |

**Limitation**: Currently returns zero for ObjPolyhdr elements (MSC method).
For single-point phi, use `rad.Fld(obj, 'phi', point)` instead.

### FldA - Vector Potential Batch (v1.3.16+)

```python
A = rad.FldA(obj, points)
```

Computes magnetic vector potential A at multiple points.

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | int | Object or container |
| `points` | [[x,y,z], ...] | List of evaluation points |

| Returns | Description |
|---------|-------------|
| `A` | List of [Ax, Ay, Az] values (T*m) |

**Limitation**: Currently returns zero for ObjPolyhdr elements (MSC method).
For single-point A, use `rad.Fld(obj, 'a', point)` instead.

### ClassifyPoints - Point Classification (v1.3.16+)

```python
result = rad.ClassifyPoints(obj, points, near_threshold=3.0)
```

Classifies evaluation points relative to mesh elements (for FMM field computation).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | int | | Object or container |
| `points` | [[x,y,z], ...] | | List of evaluation points |
| `near_threshold` | float | 3.0 | Near zone multiplier |

| Returns | Description |
|---------|-------------|
| `result['classification']` | List of int: 0=inside, 1=near, 2=far |
| `result['nearest_elem']` | List of int: index of nearest element |

```python
points = [[0, 0, 0], [0, 0, 0.1], [0, 0, 1.0]]
result = rad.ClassifyPoints(magnet, points)
# classification: [0, 1, 2] = [inside, near, far]
```

### FldLst - Field Along Line

```python
field_list = rad.FldLst(obj, component, p1, p2, n_points, 'arg')
```

### ObjM - Get Magnetization

```python
all_M = rad.ObjM(obj)  # Returns [[center, [Mx, My, Mz]], ...]
```

```python
all_M = rad.ObjM(container)
M_list = [m[1] for m in all_M]
M_avg_z = np.mean([m[2] for m in M_list])
```

---

## Mesh Import

### NGSolve Mesh Access Policy (MANDATORY)

**CRITICAL**: All NGSolve mesh access MUST use functions from `netgen_mesh_import.py`.

| Rule | Description |
|------|-------------|
| **ALWAYS** | Use `netgen_mesh_to_radia()` or `extract_elements()` |
| **NEVER** | Directly access `mesh.ngmesh.Points()`, `mesh.vertices[]`, or `el.vertices[].nr` |
| **NO EXCEPTIONS** | Applies to all scripts including examples, tests, and debugging code |

**Why?** NGSolve has TWO different indexing schemes:

| Access Method | Indexing | Valid Range |
|--------------|----------|-------------|
| `mesh.ngmesh.Points()[i]` | **1-indexed** | 1 to nv |
| `mesh.vertices[i]` | **0-indexed** | 0 to nv-1 |
| `el.vertices[i].nr` | Returns **0-indexed** | Use with `mesh.vertices[]` only |

Mixing these causes off-by-one errors that are difficult to debug.

### netgen_mesh_to_radia - Netgen Tetrahedral

```python
from netgen_mesh_import import netgen_mesh_to_radia

mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')
```

### extract_elements - Custom Processing

```python
from netgen_mesh_import import extract_elements, compute_element_centroid

elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']  # Correctly extracted coordinates
    centroid = compute_element_centroid(vertices)
```

### Available Functions in netgen_mesh_import.py

| Function | Description |
|----------|-------------|
| `netgen_mesh_to_radia()` | Convert entire mesh to Radia geometry (recommended) |
| `extract_elements()` | Extract element data for custom processing |
| `compute_element_centroid()` | Compute centroid from vertex list |
| `create_radia_tetrahedron()` | Create single Radia tetrahedron |
| `create_radia_hexahedron()` | Create single Radia hexahedron |

### create_radia_from_nastran - Nastran Import

```python
from nastran_mesh_import import create_radia_from_nastran

mag_obj = create_radia_from_nastran('model.bdf',
                                     material={'magnetization': [0, 0, 1e6]},
                                     units='m')
```

**Supported Nastran elements**: CTETRA, CHEXA, CPENTA, CPYRAM, CTRIA3

---

## NGSolve Integration

### Import Order (CRITICAL)

```python
# 1. Import radia first
import radia as rad
rad.FldUnits('m')  # REQUIRED: NGSolve uses meters

# 2. Import ngsolve BEFORE radia_ngsolve
import ngsolve
from ngsolve import *

# 3. NOW import radia_ngsolve
from radia import radia_ngsolve
```

Wrong order causes `ImportError: DLL load failed`.

### NGSolve Version Requirement

**Use NGSolve 6.2.2405 only** (6.2.2406+ has Periodic BC bug).

```bash
pip install ngsolve==6.2.2405
```

### RadiaField - CoefficientFunction

```python
cf = radia_ngsolve.RadiaField(radia_obj, field_type='b')
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `radia_obj` | int | Radia object ID |
| `field_type` | str | `'b'`, `'h'`, `'a'`, or `'m'` |

```python
# Create CoefficientFunction for B field
B_cf = radia_ngsolve.RadiaField(magnet, 'b')

# Use in NGSolve
fes = HDiv(mesh, order=2)
gf = GridFunction(fes)
gf.Set(B_cf)
```

---

## Utilities

### FldUnits - Unit System

```python
rad.FldUnits('m')   # Use meters (required for NGSolve)
rad.FldUnits('mm')  # Use millimeters (default)
rad.FldUnits()      # Get current units
```

### UtiDelAll - Clear Memory

```python
rad.UtiDelAll()
```

### UtiVer - Version

```python
version = rad.UtiVer()
```

---

## Transformations

### TrfTrsl - Translation

```python
rad.TrfTrsl(obj, [dx, dy, dz])
```

### TrfRot - Rotation

```python
rad.TrfRot(obj, [x, y, z], [nx, ny, nz], angle)
```

### TrfMlt - Multiple Copies

```python
array = rad.TrfMlt(obj, transformation, n_copies)
```

---

## Common Issues

### 1. Coordinates Off by 1000x

**Cause**: Unit mismatch (NGSolve uses meters, Radia defaults to mm)

**Solution**:
```python
rad.FldUnits('m')  # Set at start of script
```

### 2. DLL Load Failed

**Cause**: Wrong import order

**Solution**: Import ngsolve BEFORE radia_ngsolve

### 3. ObjPolyhdr Face Error

**Cause**: 0-indexed faces

**Solution**: Use **1-indexed** faces (Radia convention)

### 4. Solver Not Converging

**Solutions**:
1. Use BiCGSTAB (Method 1)
2. Increase max iterations
3. Check B-H data is monotonic
4. Verify H-M conversion: `M = B/mu_0 - H`

---

## Units

| Quantity | Unit |
|----------|------|
| Length | mm (default) or m with `FldUnits('m')` |
| B (flux density) | Tesla (T) |
| H (field) | A/m |
| M (magnetization) | A/m |
| A (vector potential) | T*m (when using `FldUnits('m')`) |
| Current | Ampere (A) |

### Internal Unit System

**IMPORTANT**: Radia ALWAYS uses millimeters (mm) internally, regardless of `FldUnits()` setting.

| Setting | Coordinate Input | B, H Output | A Output | Internal |
|---------|------------------|-------------|----------|----------|
| `FldUnits('mm')` | mm | T, A/m | T*mm | mm |
| `FldUnits('m')` | m (scaled x1000) | T, A/m | **T*mm** (needs /1000) | mm |

### Vector Potential A Unit Conversion

When using NGSolve integration with `FldUnits('m')`:

- **B, H fields**: Returned correctly in SI units (no conversion needed)
- **A field**: Returned in T*mm (requires scaling for curl(A) = B verification)

**Why A needs special handling:**

1. A is dimensionally [T*length] = [Wb/m] = [V*s/m]
2. Radia computes A using mm-based geometry: A_radia = T*mm
3. NGSolve differentiates in meters: `curl(A) = dA/dx [m^-1]`
4. For B = curl(A) to hold: `A_SI = A_radia / 1000`

**In radia_ngsolve.cpp:**

```cpp
// Vector potential A unit scaling:
// Radia ALWAYS uses mm internally, so A is always in T*mm
// NGSolve differentiates in meters: curl(A) = dA/dx_m
// To get correct B = curl(A), we scale A by 0.001:
double scale = (field_type == "a") ? 0.001 : 1.0;
```

### Maxwell Relation Verification

See `examples/ngsolve_integration/verify_curl_A_equals_B/` for a complete verification script that:

1. Creates a permanent magnet using ObjPolyhdr
2. Projects A onto HCurl space
3. Computes curl(A) using NGSolve
4. Compares with B projected onto HDiv space
5. Verifies `|curl(A)|/|B| ~= 1.0`

---

## Elliptic Integral Formulas for Coils

The magnetic field of circular current loops is computed using complete elliptic integrals of the first and second kind, K(k) and E(k). This provides analytical accuracy without numerical integration.

### Mathematical Background

For a circular current loop of radius R carrying current I, the field at cylindrical coordinates (rho, z) is:

```
k^2 = 4*R*rho / ((R+rho)^2 + z^2)

B_rho = (mu_0*I / 2*pi) * z / (rho * sqrt((R+rho)^2 + z^2)) *
        (-K(k) + (R^2 + rho^2 + z^2) / ((R-rho)^2 + z^2) * E(k))

B_z = (mu_0*I / 2*pi) * 1 / sqrt((R+rho)^2 + z^2) *
      (K(k) - (R^2 - rho^2 + z^2) / ((R-rho)^2 + z^2) * E(k))
```

The elliptic integrals are computed using the Hastings polynomial approximation, which provides accuracy to ~10^-8 relative error.

### On-Axis Field (Special Case)

For points on the axis (rho=0), the field simplifies to:

```
B_z = mu_0 * I * R^2 / (2 * (R^2 + z^2)^(3/2))
B_rho = 0
```

### Vector Potential

The azimuthal component of the vector potential A_phi is also computed analytically:

```
A_phi = (mu_0*I / pi) * sqrt(R/rho) * (1/k) * ((1 - k^2/2)*K(k) - E(k))
```

### Rectangular Cross-Section Coils

For coils with finite cross-section (radial width and height), Radia uses Gaussian quadrature to integrate the thin-loop formula over the cross-section. This maintains analytical accuracy while handling practical coil geometries.

---

## Analytical Magnet Classes (Python)

The `radia.analytical_magnet` module provides pure Python analytical field computation classes for use as background field sources. These are independent of Radia's C++ solver and can be used for:
- Background field computation with `rad.ObjBckgCF()`
- Standalone field calculations
- Verification and validation

### Available Classes

| Class | Description | B-field | H-field | A-field (vector potential) |
|-------|-------------|---------|---------|---------------------------|
| `SphericalMagnet` | Uniformly magnetized sphere | Exact dipole | Exact | Exact dipole |
| `CuboidMagnet` | Rectangular block magnet | Yang/Camacho formula | Exact | Exact (surface current) |
| `CurrentLoop` | Circular current loop | Ortner elliptic integral | Exact | Elliptic integral |
| `CylindricalMagnet` | Axially magnetized cylinder | Caciagli/Derby formula | Exact | Gaussian quadrature |
| `RingMagnet` | Hollow cylindrical magnet | Caciagli formula | Exact | Gaussian quadrature |

### Usage Examples

```python
from radia.analytical_magnet import SphericalMagnet, CuboidMagnet, CurrentLoop

# Spherical magnet (diameter 20mm, Mz = 955000 A/m)
sphere = SphericalMagnet(
    center=[0, 0, 0],      # mm
    diameter=20.0,          # mm
    magnetization=[0, 0, 955000]  # A/m
)
B = sphere.get_B([15, 0, 0])  # [Bx, By, Bz] in Tesla
H = sphere.get_H([15, 0, 0])  # [Hx, Hy, Hz] in A/m
A = sphere.get_A([15, 0, 0])  # [Ax, Ay, Az] in T*m

# Cuboid magnet (20x20x10 mm)
cuboid = CuboidMagnet(
    center=[0, 0, 0],
    dimensions=[20, 20, 10],  # mm
    magnetization=[0, 0, 955000]  # A/m
)
B = cuboid.get_B([25, 0, 0])
A = cuboid.get_A([25, 0, 0])  # Exact analytical (not dipole approximation)

# Current loop (diameter 50mm, current 100A)
loop = CurrentLoop(
    center=[0, 0, 0],
    diameter=50.0,  # mm
    current=100.0,  # A
    axis='z'
)
B = loop.get_B([0, 0, 25])
```

### Use as Background Field Source

```python
import radia as rad
from radia.analytical_magnet import CuboidMagnet

rad.FldUnits('m')

# Define permanent magnet as background field
pm = CuboidMagnet(
    center=[0, 0, 50],      # 50mm above center
    dimensions=[40, 40, 20],
    magnetization=[0, 0, 955000]
)

# Create Radia background field object
bkg = rad.ObjBckgCF(pm)  # Uses pm.__call__() which returns get_B()

# Create soft iron to solve
iron = rad.ObjHexahedron(vertices, [0, 0, 0])
mat = rad.MatLin(1000)
rad.MatApl(iron, mat)

grp = rad.ObjCnt([iron, bkg])
rad.Solve(grp, 0.001, 1000, 1)
```

### Vector Potential Verification

All classes satisfy curl(A) = B (verified numerically with < 0.01% error):

```python
# Numerical curl verification
import numpy as np
h = 0.1  # mm step
h_m = h / 1000.0  # meters

def numerical_curl(magnet, pt):
    A_px = magnet.get_A([pt[0]+h, pt[1], pt[2]])
    A_mx = magnet.get_A([pt[0]-h, pt[1], pt[2]])
    A_py = magnet.get_A([pt[0], pt[1]+h, pt[2]])
    A_my = magnet.get_A([pt[0], pt[1]-h, pt[2]])
    A_pz = magnet.get_A([pt[0], pt[1], pt[2]+h])
    A_mz = magnet.get_A([pt[0], pt[1], pt[2]-h])

    return [
        (A_py[2] - A_my[2]) / (2*h_m) - (A_pz[1] - A_mz[1]) / (2*h_m),
        (A_pz[0] - A_mz[0]) / (2*h_m) - (A_px[2] - A_mx[2]) / (2*h_m),
        (A_px[1] - A_mx[1]) / (2*h_m) - (A_py[0] - A_my[0]) / (2*h_m)
    ]

curl_A = numerical_curl(cuboid, [25, 0, 0])
B = cuboid.get_B([25, 0, 0])
# curl_A should equal B within numerical precision
```

### Key Formulas

**CuboidMagnet Vector Potential**: Uses the equivalent surface current model:
- Surface current density: K = M x n on each face
- A = (mu_0 / 4*pi) * integral_S [K / |r - r'|] dS'
- Uses Urankar (1980) / Ravaud (2009) formula for rectangular surface integration

**CurrentLoop**: Uses Ortner et al. (2023) elliptic integral formulation for both B and A fields.

---

## References

### Elliptic Integral Formulas

1. **Simpson, J.C., Lane, J.E., Immer, C.D., Youngquist, R.C.** (2001). "Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop." NASA Technical Memorandum NASA/TM-2013-217919. [NASA NTRS](https://ntrs.nasa.gov/citations/20010038494)

2. **Maxwell, J.C.** (1873). "A Treatise on Electricity and Magnetism," Vol. 2, Art. 701-706. Oxford: Clarendon Press. [Cambridge University Press Edition](https://www.cambridge.org/core/books/treatise-on-electricity-and-magnetism/130A7181ECAB0C990FBC2B88341A4141)

3. **Smythe, W.R.** (1989). "Static and Dynamic Electricity," 3rd ed., pp. 290-295. New York: Hemisphere Publishing.

### Polynomial Approximation

4. **Hastings, C., Hayward, J.T., Wong, J.P.** (1955). "Approximations for Digital Computers." Princeton University Press. [De Gruyter](https://www.degruyterbrill.com/document/doi/10.1515/9781400875597/html)

5. **Cody, W.J.** (1965). "Chebyshev Approximations for the Complete Elliptic Integrals K and E." Mathematics of Computation 19(92), pp. 105-112. [Semantic Scholar](https://www.semanticscholar.org/paper/Chebyshev-Approximations-for-the-Complete-Elliptic-Cody/e120c0220534dcee9c154478226122edf124ded5)

### Analytical Magnet Formulas

6. **Yang, Z.J., et al.** (1990). "Potential and force between a magnet and a bulk Y1Ba2Cu3O7 superconductor studied by a mechanical pendulum." Supercond. Sci. Technol. 3(12):591. - Cuboid B-field formula

7. **Camacho, J.M., Sosa, V.** (2013). "Alternative method to calculate the magnetic field of permanent magnets with azimuthal symmetry." Rev. Mex. Fis. E 59, 8-17. - Cuboid B-field validation

8. **Cichon, D.** (2019). "Stability of magnetic field computation near edges using analytical formulas." Master's thesis. - Numerical stability improvements

### Surface Current Vector Potential

9. **Urankar, L.K.** (1980). "Vector potential and magnetic field of current-carrying finite arc segment in analytical form." IEEE Trans. Magn. 16(5), 1283-1288. - Rectangular surface integral formula

10. **Ravaud, R., et al.** (2009). "Analytical calculation of the magnetic field created by permanent-magnet rings." IEEE Trans. Magn. 45(4), 1572-1576. - Surface current A-field integration

### Triangle B-field (MSC Method)

11. **Guptasarma, D.** (1999). "Computation of the time-domain response of a polarizable ground." Geophysics 64(1), 70-74. - Solid angle formula for triangle B-field

12. **van Oosterom, A., Strackee, J.** (1983). "The solid angle of a plane triangle." IEEE Trans. Biomed. Eng. 30(2), 125-126. - Efficient solid angle computation

### Potential Integrals on Triangles

13. **Carley, M.** (2013). "Potential integrals on triangles." arXiv:1201.4938. - Analytical formula for 1/r integral over triangular surfaces

### General References

14. [ESRF Radia Reference Guide](https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html)
15. [examples/cube_uniform_field/](../examples/cube_uniform_field/) - Benchmark examples

---

**Last Updated**: 2025-12-30
**License**: LGPL-2.1 (modifications), BSD-style (original RADIA from ESRF)
