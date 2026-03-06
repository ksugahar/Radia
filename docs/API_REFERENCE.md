# Radia Python API Reference

Complete reference for Radia Python API.

**Version**: 1.4.4
**Date**: 2026-01-09
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
- [VTK Export](#vtk-export)
- [ESIM (Effective Surface Impedance Method)](#esim-effective-surface-impedance-method)

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

ext = rad.ObjBckg(lambda p: [0, 0, MU_0 * 50000])  # B field in Tesla
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

### ObjBckg - Background Field (Callback)

```python
field_src = rad.ObjBckg(lambda p: [Bx, By, Bz])  # Uniform field
field_src = rad.ObjBckg(callback_function)        # Non-uniform field
```

The callback function receives a point `[x, y, z]` in current units and returns `[Bx, By, Bz]` in Tesla.

**Uniform Background Field**:
```python
MU_0 = 4 * np.pi * 1e-7
H_ext = 50000  # A/m
ext = rad.ObjBckg(lambda p: [0, 0, MU_0 * H_ext])  # B = mu_0 * H in Tesla
```

**Non-uniform Background Field** (e.g., quadrupole):
```python
def quadrupole_field(point):
    x, y, z = point
    G = 10.0  # T/m gradient
    return [G * y, G * x, 0]

ext = rad.ObjBckg(quadrupole_field)
```

**IMPORTANT**: The callback returns **B field in Tesla**, NOT H field in A/m.

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
mat = rad.MatLin([5000, 100], [0, 0, 1])
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

**Note**: `ksi` here is a **fitting parameter** for the tanh saturation formula (initial susceptibility of each term), not the same as the bulk susceptibility chi = mu_r - 1.

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
| `method` | int | `0` = LU, `1` = BiCGSTAB (default), `2` = HACApK |

| Returns | Description |
|---------|-------------|
| `result[0]` | Final residual |
| `result[3]` | Number of iterations |

### Solver Selection

| Problem Size | Elements | Method | Code |
|--------------|----------|--------|------|
| Small | < 500 | LU | `rad.Solve(grp, 0.001, 100, 0)` |
| Medium | 500-5,000 | BiCGSTAB | `rad.Solve(grp, 0.001, 1000, 1)` |
| Large | > 5,000 | HACApK | `rad.Solve(grp, 0.001, 1000, 2)` |

**HACApK Advantages** (for large problems):
- O(N log N) memory vs O(N^2) for dense
- 10-20% of dense matrix memory at 10,000+ elements
- 5-10x faster than BiCGSTAB for 10,000+ elements

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
- `relax>0.0`: Damped update: `M_new = M_computed*(1-relax) + M_old*relax`
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

### SetNewtonMethod - Newton-Raphson Method (v1.4.4+)

```python
rad.SetNewtonMethod(enabled)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enabled` | bool | False | Enable Newton-Raphson with differential susceptibility |

**Notes:**
- Uses tangent stiffness matrix based on differential susceptibility χ_d = (dB/dH)/μ₀ - 1
- Reduces linear iterations per nonlinear step compared to Picard (fixed-point) method
- Recommended for large-scale problems (>10,000 DOF) with HACApK solver (Method 2)
- Combine with `SetNewtonDamping()` for robust convergence
- Call BEFORE `rad.Solve()`

**Example:**
```python
rad.SetSolver(2)  # HACApK
rad.SetNewtonMethod(True)
rad.SetNewtonDamping(True)  # Enable line search damping
rad.Solve(container, 0.001, 100)
```

### GetNewtonMethod - Query Newton-Raphson Status (v1.4.4+)

```python
enabled = rad.GetNewtonMethod()
```

Returns `True` if Newton-Raphson is enabled, `False` otherwise.

### SetNewtonDamping - Line Search Damping (v1.4.4+)

```python
rad.SetNewtonDamping(enabled=True, max_iter=5, min_omega=0.01)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enabled` | bool | True | Enable backtracking line search |
| `max_iter` | int | 5 | Maximum line search iterations |
| `min_omega` | float | 0.01 | Minimum damping factor (1%) |

**Notes:**
- Only active when Newton-Raphson is enabled (`SetNewtonMethod(True)`)
- Uses backtracking line search: σ_new = ω × σ_trial + (1-ω) × σ_old
- Starts with ω=1.0 (full Newton step), reduces by ω *= 0.5 if residual increases
- Ensures monotonic residual reduction, improves nonlinear convergence
- Particularly effective for highly nonlinear problems
- Call BEFORE `rad.Solve()`

**Algorithm:**
1. Solve linear system: (D(1/χ_d) + G) σ_trial = RHS
2. Try full step ω=1.0
3. If residual increases: ω *= 0.5, retry with damped update
4. Repeat up to `max_iter` times
5. Force accept if ω < `min_omega`

**Example:**
```python
# Enable Newton with line search damping (recommended)
rad.SetNewtonMethod(True)
rad.SetNewtonDamping(True, max_iter=5, min_omega=0.01)
rad.Solve(container, 0.001, 100, 2)  # HACApK + Newton + Damping

# Disable damping for full Newton steps
rad.SetNewtonDamping(False)
rad.Solve(container, 0.001, 100, 2)  # HACApK + Newton (no damping)
```

### GetNewtonDampingStats - Query Damping Statistics (v1.4.4+)

```python
stats = rad.GetNewtonDampingStats()
```

Returns a dictionary with damping configuration:

| Key | Type | Description |
|-----|------|-------------|
| `'enabled'` | bool | Damping enabled status |
| `'max_iter'` | int | Maximum line search iterations |
| `'min_omega'` | float | Minimum damping factor |

**Example:**
```python
stats = rad.GetNewtonDampingStats()
print(f"Damping: {stats['enabled']}, max_iter={stats['max_iter']}, min_omega={stats['min_omega']}")
# Output: Damping: True, max_iter=5, min_omega=0.01
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
| `'p'`, `'phi'` | Scalar potential Phi (A) |
| `'mx'`, `'my'`, `'mz'`, `'m'` | Magnetization M |

```python
B = rad.Fld(magnet, 'b', [0, 0, 0.1])    # B vector at point
Bz = rad.Fld(magnet, 'bz', [0, 0, 0.1])  # Bz component
H = rad.Fld(magnet, 'h', [0, 0, 0.1])    # H vector at point
A = rad.Fld(magnet, 'a', [0, 0, 0.1])    # Vector potential A
Phi = rad.Fld(magnet, 'p', [0, 0, 0.1])  # Scalar potential Phi
```

**Potential Field Notes (v1.4.2+)**:
- **A (Vector Potential)**: Uses face-based integration `A = (1/4pi) * M x BufVect`
- **Phi (Scalar Potential)**: Uses face-based integration `Phi = (1/4pi) * M . BufVect`
- Both A and Phi are computed accurately for ObjHexahedron/ObjTetrahedron
- Maxwell relations verified: `curl(A) ∝ B`, `-grad(Phi) ∝ H`

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

**Note**: Uses face-based integration for accurate scalar potential computation (v1.4.2+).

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

**Note**: Uses face-based integration for accurate vector potential computation (v1.4.2+).

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

# 3. NOW import radia_ngsolve (separate C++ module)
import radia_ngsolve
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

## VTK/VTS Export

Export Radia magnetic field to VTS (VTK XML Structured Grid) format for visualization in ParaView.

### rad.FldVTS()

Export magnetic field on a structured 3D grid to VTS format (C++ implementation with OpenMP):

```python
import radia as rad

rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# FldVTS(obj, filename, x_range, y_range, z_range, nx, ny, nz, include_B, include_H, unit_scale)
rad.FldVTS(magnet, 'B_field.vts',
           [-0.08, 0.08], [-0.08, 0.08], [0.03, 0.12],  # x, y, z ranges
           33, 33, 19,    # grid points
           1, 0,          # include_B=True, include_H=False
           1.0)           # unit_scale (1.0 for meters)
# -> B_field.vts
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `obj` | int | required | Radia object handle |
| `filename` | str | required | Output filename (.vts) |
| `x_range` | [float, float] | required | X range [min, max] |
| `y_range` | [float, float] | required | Y range [min, max] |
| `z_range` | [float, float] | required | Z range [min, max] |
| `nx` | int | 21 | Grid points in X |
| `ny` | int | 21 | Grid points in Y |
| `nz` | int | 21 | Grid points in Z |
| `include_B` | int | 1 | Include B field (0/1) |
| `include_H` | int | 0 | Include H field (0/1) |
| `unit_scale` | float | 1.0 | Coordinate scale factor |

**Output Fields**:

| Field | Type | Unit | Description |
|-------|------|------|-------------|
| `B` | Vector (3) | T | Magnetic flux density |
| `B_magnitude` | Scalar | T | \|B\| |
| `H` | Vector (3) | A/m | Magnetic field (if include_H=1) |
| `H_magnitude` | Scalar | A/m | \|H\| (if include_H=1) |

**Example** - Complete workflow:

```python
import radia as rad

rad.FldUnits('m')

# Create permanent magnet
vertices = [
    [-0.02, -0.02, -0.01], [0.02, -0.02, -0.01],
    [0.02, 0.02, -0.01], [-0.02, 0.02, -0.01],
    [-0.02, -0.02, 0.01], [0.02, -0.02, 0.01],
    [0.02, 0.02, 0.01], [-0.02, 0.02, 0.01],
]
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])

# Export field to VTS
rad.FldVTS(magnet, 'magnet_field.vts',
           [-0.05, 0.05], [-0.05, 0.05], [0.02, 0.08],
           21, 21, 13, 1, 1, 1.0)

# View in ParaView:
# 1. File -> Open -> magnet_field.vts
# 2. Apply
# 3. Color by B_magnitude or use Glyph filter for B vectors
```

---

## ESIM (Effective Surface Impedance Method)

The ESIM module provides specialized tools for **induction heating analysis** with nonlinear magnetic materials. It implements the Effective Surface Impedance Method for computing eddy current losses and coil-workpiece coupled impedance.

**Availability**: Requires `scipy` package. Check `radia.ESIM_AVAILABLE` to verify.

### Overview

| Module | Description |
|--------|-------------|
| `ESIMCellProblemSolver` | Solve 1D cell problem for effective surface impedance |
| `BHCurveInterpolator` | B-H curve interpolation with H-dependent permeability |
| `ComplexPermeabilityInterpolator` | Complex permeability mu' - j*mu" support |
| `ESIMWorkpiece` | Workpiece geometry with surface panels |
| `InductionHeatingCoil` | Coil geometry (spiral, loop, etc.) |
| `ESIMCoupledSolver` | Coupled coil-workpiece impedance solver |
| `ESIMVTKOutput` | VTK export for visualization |

### ESIMCellProblemSolver - Cell Problem Solver

Solves the 1D boundary value problem for effective surface impedance Z(H0):

```python
from radia import ESIMCellProblemSolver

solver = ESIMCellProblemSolver(
    sigma=5e6,           # Conductivity [S/m]
    frequency=50000,     # Frequency [Hz]
    mu_r=100,            # Relative permeability (constant)
    # OR
    bh_curve=bh_data,    # H-dependent permeability [[H, B], ...]
    # OR
    complex_mu=(1000, 100)  # (mu'_r, mu"_r) for magnetic loss
)

result = solver.solve(H0=5000)  # Surface field [A/m]
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `sigma` | float | Electrical conductivity [S/m] |
| `frequency` | float | Analysis frequency [Hz] |
| `mu_r` | float | Constant real relative permeability (optional) |
| `bh_curve` | list | H-dependent B-H data [[H, B], ...] (optional) |
| `complex_mu` | tuple/list | Complex permeability - see below (optional) |

**complex_mu Parameter Formats**:

| Format | Example | Description |
|--------|---------|-------------|
| Tuple `(mu'_r, mu"_r)` | `(1000, 100)` | Constant complex mu |
| List of lists | `[[H, mu'_r, mu"_r], ...]` | H-dependent complex mu |

**Returns** (dict):

| Key | Description |
|-----|-------------|
| `Z` | Complex surface impedance [Ohm] |
| `P_prime` | Total power density [W/m^2] |
| `P_ohmic` | Ohmic (Joule) loss from sigma [W/m^2] |
| `P_magnetic` | Magnetic loss from mu" [W/m^2] |
| `delta` | Effective skin depth [m] |
| `iterations` | Number of iterations (nonlinear) |

**Example - Constant Permeability**:

```python
from radia import ESIMCellProblemSolver

solver = ESIMCellProblemSolver(
    sigma=5e6,        # 5 MS/m (hot steel)
    frequency=50000,  # 50 kHz
    mu_r=100
)

result = solver.solve(H0=5000)
print(f"Z = {result['Z'].real:.4f} + j{result['Z'].imag:.4f} Ohm")
print(f"Power density = {result['P_prime']:.0f} W/m^2")
```

**Example - Nonlinear B-H Curve**:

```python
bh_curve = [
    [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
    [2500, 1.6], [5000, 1.8], [10000, 1.95]
]

solver = ESIMCellProblemSolver(
    sigma=2e6,
    frequency=50000,
    bh_curve=bh_curve
)

result = solver.solve(H0=5000)
```

**Example - Complex Permeability (Magnetic Loss)**:

```python
# Constant complex mu
solver = ESIMCellProblemSolver(
    sigma=1e6,
    frequency=50000,
    complex_mu=(1000, 100)  # mu'_r=1000, mu"_r=100
)

# H-dependent complex mu
complex_mu_data = [
    [0, 2000, 200],      # [H, mu'_r, mu"_r]
    [1000, 1500, 150],
    [5000, 500, 50],
]
solver = ESIMCellProblemSolver(
    sigma=1e6,
    frequency=50000,
    complex_mu=complex_mu_data
)
```

### BHCurveInterpolator - B-H Curve Interpolation

Interpolates B-H curves and computes derived quantities:

```python
from radia import BHCurveInterpolator

interp = BHCurveInterpolator(bh_data)

B = interp.get_B(H=5000)           # B at given H
mu_r = interp.get_mu_r(H=5000)     # Relative permeability
dB_dH = interp.get_dB_dH(H=5000)   # Differential permeability
```

| Method | Returns |
|--------|---------|
| `get_B(H)` | B [T] at field H [A/m] |
| `get_mu_r(H)` | Relative permeability mu_r at H |
| `get_dB_dH(H)` | Differential dB/dH at H |
| `get_B_sat()` | Saturation B [T] |
| `get_H_sat()` | Saturation H [A/m] |

### ComplexPermeabilityInterpolator - Complex mu Support

For materials with magnetic hysteresis loss (ferrites, amorphous metals, laminated steel):

**Complex Permeability**: mu = mu' - j*mu"

| Component | Symbol | Physical Meaning |
|-----------|--------|------------------|
| Real part | mu' | Energy storage (reactive power) |
| Imaginary part | mu" | Energy loss (hysteresis, domain wall motion) |
| Loss tangent | tan(delta_m) = mu"/mu' | Ratio of loss to storage |

**Power loss from magnetic hysteresis**:
```
P_magnetic = (omega/2) * mu_0 * mu"_r * |H|^2  [W/m^3]
```

```python
from radia import ComplexPermeabilityInterpolator

# Constant complex mu
interp = ComplexPermeabilityInterpolator(
    complex_mu=(1000, 100)  # (mu'_r, mu"_r)
)

# H-dependent complex mu
complex_mu_data = [
    [0, 2000, 200],      # [H, mu'_r, mu"_r]
    [1000, 1000, 100],
    [10000, 200, 20],
]
interp = ComplexPermeabilityInterpolator(complex_mu=complex_mu_data)

mu_prime = interp.get_mu_prime(H=5000)   # Real part mu'
mu_double_prime = interp.get_mu_double_prime(H=5000)  # Imaginary part mu"
loss_tangent = mu_double_prime / mu_prime
```

**Typical Material Properties**:

| Material | mu'_r | mu"_r | tan(delta_m) | Application |
|----------|-------|-------|--------------|-------------|
| MnZn Ferrite (1 kHz) | 2500 | 25 | 0.01 | Power transformers |
| MnZn Ferrite (100 kHz) | 2000 | 400 | 0.2 | Switching supplies |
| NiZn Ferrite (1 MHz) | 150 | 75 | 0.5 | EMI suppression |
| Amorphous Metal | 10000 | 100 | 0.01 | High-efficiency cores |
| Laminated Steel (60 Hz) | 4000 | 40 | 0.01 | Power transformers |

### InductionHeatingCoil - Coil Geometry

Creates coil geometry for field computation:

```python
from radia import InductionHeatingCoil

# Spiral (pancake) coil
coil = InductionHeatingCoil(
    coil_type='spiral',
    center=[0, 0, 0.02],       # [m]
    inner_radius=0.03,         # 30mm
    outer_radius=0.06,         # 60mm
    pitch=0.005,               # 5mm
    num_turns=4,
    axis=[0, 0, 1],
    wire_width=0.004,          # 4mm
    wire_height=0.002,         # 2mm
    conductivity=5.8e7,        # Copper
)
coil.set_current(150)  # 150 A

# Single-turn loop coil
coil = InductionHeatingCoil(
    coil_type='loop',
    center=[0, 0, 0.015],
    radius=0.04,               # 40mm radius
    normal=[0, 0, 1],
    wire_width=0.005,
    wire_height=0.003,
)
coil.set_current(300)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `coil_type` | str | `'spiral'` or `'loop'` |
| `center` | list | Coil center [x, y, z] in meters |
| `inner_radius` | float | Inner radius (spiral only) [m] |
| `outer_radius` | float | Outer radius (spiral only) [m] |
| `radius` | float | Loop radius (loop only) [m] |
| `pitch` | float | Axial pitch per turn (spiral) [m] |
| `num_turns` | int | Number of turns (spiral) |
| `axis` / `normal` | list | Coil axis direction |
| `wire_width` | float | Wire width [m] |
| `wire_height` | float | Wire height [m] |
| `conductivity` | float | Wire conductivity [S/m] |

**Methods**:

| Method | Description |
|--------|-------------|
| `set_current(I)` | Set coil current [A] |
| `get_B_field(point)` | Compute B at point [T] |
| `get_self_inductance()` | Self-inductance [H] |
| `get_ac_resistance(freq)` | AC resistance [Ohm] |

### ESIMWorkpiece - Workpiece Geometry

Represents workpiece with surface panels for ESIM analysis:

```python
from radia import create_esim_block, create_esim_cylinder

# Rectangular block workpiece
workpiece = create_esim_block(
    center=[0, 0, -0.01],          # 10mm below origin
    dimensions=[0.12, 0.12, 0.02], # 120mm x 120mm x 20mm
    bh_curve=bh_curve,
    sigma=2e6,                      # 2 MS/m
    frequency=50000,
    panels_per_side=5
)

# Cylindrical workpiece
workpiece = create_esim_cylinder(
    center=[0, 0, 0],
    radius=0.05,                    # 50mm radius
    height=0.03,                    # 30mm height
    bh_curve=bh_curve,
    sigma=2e6,
    frequency=50000,
    panels_radial=4,
    panels_axial=3
)
```

### ESIMCoupledSolver - Coupled Impedance Solver

Solves coupled coil-workpiece system and computes impedance:

```python
from radia import ESIMCoupledSolver

solver = ESIMCoupledSolver(coil, workpiece, frequency)
result = solver.solve(tol=1e-4, max_iter=30, verbose=True)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `coil` | InductionHeatingCoil | Coil object |
| `workpiece` | ESIMWorkpiece | Workpiece object |
| `frequency` | float | Analysis frequency [Hz] |

**solve() Parameters**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tol` | 1e-4 | Convergence tolerance |
| `max_iter` | 30 | Maximum iterations |
| `verbose` | False | Print iteration info |

**Returns** (dict):

| Key | Description |
|-----|-------------|
| `P_total` | Total power to workpiece [W] |
| `Q_total` | Reactive power [var] |
| `power_factor` | Power factor |
| `impedance` | Impedance analysis dict |
| `converged` | Convergence flag |
| `iterations` | Iteration count |

**Impedance Analysis** (`result['impedance']`):

| Key | Description |
|-----|-------------|
| `L_coil_uH` | Coil inductance [uH] |
| `R_coil_mOhm` | Coil AC resistance [mOhm] |
| `R_reflected_mOhm` | Reflected resistance [mOhm] |
| `X_reflected_mOhm` | Reflected reactance [mOhm] |
| `Z_total_magnitude_mOhm` | Total |Z| [mOhm] |
| `phase_deg` | Impedance phase [deg] |
| `efficiency` | Heating efficiency |

**Example - Complete Analysis**:

```python
from radia import (
    ESIMCoupledSolver, InductionHeatingCoil, create_esim_block
)

# Steel B-H curve
bh_curve = [
    [0, 0], [100, 0.2], [500, 0.9], [1000, 1.3],
    [2500, 1.6], [5000, 1.8], [10000, 1.95],
]

# Create coil
coil = InductionHeatingCoil(
    coil_type='spiral',
    center=[0, 0, 0.02],
    inner_radius=0.03,
    outer_radius=0.06,
    pitch=0.005,
    num_turns=4,
    axis=[0, 0, 1],
    wire_width=0.004,
    wire_height=0.002,
    conductivity=5.8e7,
)
coil.set_current(150)

# Create workpiece
workpiece = create_esim_block(
    center=[0, 0, -0.01],
    dimensions=[0.12, 0.12, 0.02],
    bh_curve=bh_curve,
    sigma=2e6,
    frequency=50000,
    panels_per_side=5
)

# Solve
solver = ESIMCoupledSolver(coil, workpiece, 50000)
result = solver.solve(tol=1e-4, max_iter=30, verbose=True)

# Results
imp = result['impedance']
print(f"Power to workpiece: {result['P_total']:.0f} W")
print(f"Coil inductance: {imp['L_coil_uH']:.3f} uH")
print(f"Reflected resistance: {imp['R_reflected_mOhm']:.3f} mOhm")
print(f"Efficiency: {imp['efficiency']*100:.1f}%")

# Resonance capacitor
C_res = 1 / (solver.omega**2 * solver.L_coil)
print(f"Resonance capacitor: {C_res*1e6:.2f} uF")
```

### VTK Export Functions

Export ESIM results for visualization:

```python
from radia import (
    export_esim_workpiece_vtk,
    export_esim_coil_field_vtk,
    export_esim_combined_vtk,
)

# Export workpiece with surface fields
export_esim_workpiece_vtk(
    workpiece,
    'workpiece.vtk',
    fields=['H_tan', 'P_density', 'Z_surface']
)

# Export coil field on a plane
export_esim_coil_field_vtk(
    coil,
    'coil_field.vtk',
    plane='xy',
    z_level=0.0,
    extent=[-0.1, 0.1, -0.1, 0.1],
    resolution=50
)

# Combined export
export_esim_combined_vtk(
    coil, workpiece, result,
    'combined.vtk'
)
```

### Frequency Sweep Example

```python
from radia import ESIMCoupledSolver, InductionHeatingCoil, create_esim_block

frequencies = [10000, 25000, 50000, 100000]  # 10-100 kHz

for freq in frequencies:
    workpiece = create_esim_block(
        center=[0, 0, -0.01],
        dimensions=[0.08, 0.08, 0.015],
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=freq,
        panels_per_side=4
    )

    solver = ESIMCoupledSolver(coil, workpiece, freq)
    result = solver.solve(tol=1e-3, max_iter=20)

    imp = result['impedance']
    print(f"{freq/1000:.0f} kHz: P={result['P_total']:.0f} W, "
          f"Eff={imp['efficiency']*100:.1f}%")
```

### Physical Background

The ESIM method is based on:

1. **Cell Problem**: 1D BVP for electromagnetic field penetration into conductor
2. **Surface Impedance**: Z = E_tan / H_tan at surface
3. **Power Density**: P' = Re(Z) * |H_tan|^2

**Skin Depth**:
```
delta = sqrt(2 / (omega * mu * sigma))
```

**Surface Impedance** (linear material):
```
Z = (1 + j) / (sigma * delta) = (1 + j) * sqrt(omega * mu / (2 * sigma))
```

**Complex Permeability Loss**:
```
P_magnetic = (omega * mu_0 * mu"_r / 2) * integral |H|^2 dV
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

### TrfOrnt - Apply Transformation

```python
rad.TrfOrnt(obj, trf)  # Apply transformation to orient object
```

### Image Symmetry (Replaces TrfMlt)

**Note**: `TrfMlt` has been removed (2026-01-31). Use `rad.Image()` for plane symmetry:

```python
intrc = rad.PreRelax(container, container)
n_ima = rad.Image(intrc, '+x')  # Symmetric x-mirror (tangent field)
rad.BuildImageMatrix(intrc)
rad.Solve(container, 0.0001, 100, 0)
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
| `+xy`, `xy` | Quarter model (symmetric) | Both tangent |
| `-xz` | xz-mirror (antisymmetric) | Both normal |
| `+xyz` | Eighth model (symmetric) | All tangent |

**Matrix Construction**:
- Symmetric (+): `N_Image[i,j] = N[i,j] + N[i, mirror_j] @ P`
- Antisymmetric (-): `N_Image[i,j] = N[i,j] - N[i, mirror_j] @ P`

**Legacy API** (deprecated):
- `SetIMASymmetry()` -> Use `Image()`
- `BuildIMAMatrix()` -> Use `BuildImageMatrix()`

See [IMA_SYMMETRY_DESIGN.md](IMA_SYMMETRY_DESIGN.md) for details.

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

### 3. ObjPolyhdr Face Error (Internal API)

**Cause**: 0-indexed faces when using internal ObjPolyhdr API

**Solution**: Use **1-indexed** faces (Radia convention). For Python users, prefer `ObjHexahedron` and `ObjTetrahedron` which auto-generate faces.

### 4. Solver Not Converging

**Solutions**:
1. Use BiCGSTAB (Method 1)
2. Increase max iterations
3. Check B-H data is monotonic
4. Verify H-M conversion: `M = B/mu_0 - H`

---

## Units

### Unit System (v1.4.3+)

**IMPORTANT**: Starting from v1.4.3, Radia uses SI units (meters) internally, matching ELF.

| Quantity | Unit | Notes |
|----------|------|-------|
| Length | m (default) or mm with `FldUnits('mm')` | User-selectable |
| B (flux density) | Tesla (T) | **Fixed to SI** |
| H (field) | A/m | **Fixed to SI** |
| M (magnetization) | A/m | **Fixed to SI** |
| A (vector potential) | T*m | **Fixed to SI** |
| Current | Ampere (A) | **Fixed to SI** |
| Current density | A/m^2 (or A/mm^2) | Depends on length unit |

### Design Principles

1. **Length unit only**: `FldUnits()` controls **only the length unit** (m or mm)
2. **Field values are always SI**: B, H, and A are always in SI units (T, A/m, T*m)
3. **No field scaling**: Changing length unit does NOT change B, H, or A values

| Setting | Coordinate Input | B, H Output | A Output |
|---------|------------------|-------------|----------|
| `FldUnits('m')` (default) | meters | T, A/m | T*m |
| `FldUnits('mm')` | millimeters | T, A/m | T*m |

### Maxwell Relation: B = curl(A)

With the new SI internal units, the Maxwell relation `B = curl(A)` is satisfied without any unit conversion:

```python
import radia as rad
import numpy as np

rad.FldUnits('m')  # Default: meters

# Create magnet
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 954930])

# Get fields
point = [0.05, 0.03, 0.04]
B = rad.Fld(magnet, 'b', point)  # Tesla
A = rad.Fld(magnet, 'a', point)  # T*m

# Numerical curl
h = 1e-6
A_xp = rad.Fld(magnet, 'a', [point[0]+h, point[1], point[2]])
A_xm = rad.Fld(magnet, 'a', [point[0]-h, point[1], point[2]])
# ... (compute full curl)
# Result: |curl(A)| / |B| should be approximately 1.0
```

### Maxwell Relation Verification

See `examples/ngsolve_integration/verify_curl_A_equals_B/` for a complete verification script that:

1. Creates a permanent magnet using ObjHexahedron
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

### Potential Integrals on Triangles (PEEC)

13. **Wilton, D.R., Rao, S.M., Glisson, A.W., Schaubert, D.H., Al-Bundak, O.M., Butler, C.M.** (1984). "Potential Integrals for Uniform and Linear Source Distributions on Polygonal and Polyhedral Domains." IEEE Trans. Antennas Propag. 32(3), 276-281. - Analytical 1/R integral over triangular panels using edge-based logarithm and arctangent formulas

14. **Graglia, R.D.** (1993). "On the Numerical Integration of the Linear Shape Functions Times the 3-D Green's Function or its Gradient on a Plane Triangle." IEEE Trans. Antennas Propag. 41(10), 1448-1455. - Higher-order integration formulas

15. **Carley, M.** (2013). "Potential integrals on triangles." arXiv:1201.4938. - Alternative analytical formula for 1/r integral

### Arc Coil Analytical Formulas

16. **Kameari, A.** (1990). "Calculation of Transient 3D Eddy Current using Edge-Elements." IEEE Trans. Magn. 26(2), 466-469. - Analytical formula for arc coil B-field using incomplete elliptic integrals, and 1/r integral over rectangular cross-sections. Also provides the analytical cross-section integration formula using distances to 4 corners (R[0-3]) at each azimuthal angle.

17. **Abramowitz, M., Stegun, I.A.** (1964). "Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables." National Bureau of Standards, Chapter 17. - Incomplete elliptic integrals F(phi,m) and E(phi,m), arithmetic-geometric mean algorithm

18. **Nakata, T., Takahashi, N., Fujiwara, K.** (1990). "Summary of Results for TEAM Problem 7." COMPEL 9(2), 137-154. - Validation of arc coil formulas

19. **Piessens, R., de Doncker-Kapenga, E., Überhuber, C.W., Kahaner, D.K.** (1983). "QUADPACK: A Subroutine Package for Automatic Integration." Springer-Verlag. - Gauss-Kronrod 7-15 adaptive quadrature for azimuthal integration in near-field arc coil calculation

### General References

20. [ESRF Radia Reference Guide](https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html)
21. [examples/cube_uniform_field/](../examples/cube_uniform_field/) - Benchmark examples

---

**Last Updated**: 2026-01-15
**License**: LGPL-2.1 (modifications), BSD-style (original RADIA from ESRF)
