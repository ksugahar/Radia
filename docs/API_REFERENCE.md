# Radia Python API Reference

Complete reference for Radia Python API.

**Version**: 1.3.14
**Date**: 2025-12-15
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
| **Hexahedron (MSC)** | `ObjPolyhdr()` + `HEX_FACES` | 6 quad | 6 | Permanent magnets, soft iron |
| **Tetrahedron** | `ObjPolyhdr()` + `TETRA_FACES` | 4 tri | 3 | Complex curved geometry |
| **Wedge/Prism** | `ObjPolyhdr()` + `WEDGE_FACES` | 5 | 3 | Hybrid meshes |
| **Pyramid** | `ObjPolyhdr()` + `PYRAMID_FACES` | 5 | 3 | Mesh transitions |

**DOF (Degrees of Freedom)**:
- **Hexahedra (6 faces)**: 6 DOF - Surface charge density (sigma) per face (MSC method)
- **Other elements (4-5 faces)**: 3 DOF - Magnetization vector (Mx, My, Mz)
- All meshes are expected to be generated externally (Netgen, GMSH, Cubit, etc.)

### Face Topology Constants

```python
from netgen_mesh_import import TETRA_FACES, HEX_FACES, WEDGE_FACES, PYRAMID_FACES

# TETRA_FACES (1-indexed)
[[1, 3, 2], [1, 2, 4], [2, 3, 4], [3, 1, 4]]

# HEX_FACES (1-indexed)
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

### ObjPolyhdr - General Polyhedron

```python
obj = rad.ObjPolyhdr(vertices, faces, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `vertices` | [[x,y,z], ...] | 3D vertex coordinates |
| `faces` | [[v1,v2,...], ...] | Face vertex indices (**1-indexed!**) |
| `magnetization` | [Mx, My, Mz] | Initial magnetization |

```python
from netgen_mesh_import import TETRA_FACES
vertices = [[0,0,0], [1,0,0], [0.5,1,0], [0.5,0.5,1]]
tet = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 1e6])
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

### ObjRaceTrk - Racetrack Coil

```python
coil = rad.ObjRaceTrk(center, radii, heights, current, n_segments)
```

### ObjFlmCur - Filament Conductor

```python
filament = rad.ObjFlmCur([[x1,y1,z1], [x2,y2,z2], ...], current)
```

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

## References

1. [ESRF Radia Reference Guide](https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html)
2. [examples/cube_uniform_field/](../examples/cube_uniform_field/) - Benchmark examples

---

**Last Updated**: 2025-12-15
**License**: LGPL-2.1 (modifications), BSD-style (original RADIA from ESRF)
