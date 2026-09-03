# Radia Python API Reference

Complete reference for Radia Python API.

**Version**: 4.55.0
**Date**: 2026-05-16
**Original ESRF Documentation**: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html

> Refer to [CHANGELOG.md](../../CHANGELOG.md) for per-release additions / removals
> since v4.6.0.  The PEEC STEP-loading subsystem
> (`coil_from_cad.filaments_from_step`, classification dispatch, RMF,
> adaptive resampling, cap-centroid anchoring) was overhauled in
> v4.48.2 -> v4.55.0; see
> [docs/peec/PEEC_CONDUCTOR_MODELING_GUIDE.md](../peec/PEEC_CONDUCTOR_MODELING_GUIDE.md)
> for the architectural overview and `radia-mcp peec_inductance(topic=...)`
> for the runnable knowledge layer.

---

## Table of Contents

- [Quick Start](#quick-start)
- [Removed APIs](#removed-apis)
- [Supported Elements](#supported-elements)
- [Geometry Objects](#geometry-objects)
- [Materials](#materials)
- [Solver](#solver)
- [Field Computation](#field-computation)
- [Mesh Import](#mesh-import)
- [NGSolve Integration](#ngsolve-integration)
- [Utilities](#utilities)
- [ESIM VTK Export](#esim-vtk-export)
- [PEEC Solver](#peec-solver)
- [ESIM (Effective Surface Impedance Method)](#esim-effective-surface-impedance-method)

---

## Quick Start

### HDiv-VIM Hexahedral Example

```python
import radia as rad
import numpy as np
import ngsolve as ng
from radia.vim import soft_iron_box

rad.UtiDelAll()

MU_0 = 4 * np.pi * 1e-7
cube_size = 1.0
iron = soft_iron_box(
    center=(0.0, 0.0, 0.0),
    size=(cube_size, cube_size, cube_size),
    mu_r=1000.0,
    nsub=5,
)

ext = rad.ObjBckg(lambda p: [0, 0, MU_0 * 50000])  # B field in Tesla
grp = rad.ObjCnt([iron, ext])
with ng.TaskManager():
    result = rad.Solve(grp)
```

### Netgen Mesh Example (HDiv-VIM Workflow)

```python
import radia as rad

# Import NGSolve BEFORE radia modules
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

# Create tetrahedral mesh
cube = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube.mat('magnetic')
mesh = Mesh(OCCGeometry(cube).GenerateMesh(maxh=0.3))

# Import to Radia for fixed-magnetization field evaluation.
# Soft-iron and nonlinear magnetic-material solves should use HDiv-VIM /
# reduced-FEM workflows rather than legacy 3-DOF tetrahedral moment paths.
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')
```

---

## Removed APIs

The following APIs have been removed from Radia. Calling them will raise an error.

### Removed Functions

| Removed API | Date | Replacement | Reason |
|-------------|------|-------------|--------|
| `FldUnits()` | — | None needed | Radia always uses meters. No configuration needed. |
| `RlxPre()`, `RlxMan()`, `RlxAuto()` | — | `rad.Solve(obj, prec, maxiter, method)` | Unified solver API |
| `RlxUpdSrc()`, `SetRelaxSubInterval()` | — | `rad.Solve()` | Unified solver API |
| `TrfMlt()`, `SetIMASymmetry()`, `BuildIMAMatrix()`, `PreRelax()`, `Image()` | 2026-01-31 | `rad.Solve(image=...)` / `rad.BuildMatrix(image=...)` | Unified image symmetry parameter |
| `CndLoop`, `CndRecBlock`, `CndLoopFromHelix`, `CplMagCreate`, `CplMagSolve`, `CplMagSetFrequency`, `CndHexahedron`, `CndWire`, `CndSpiral`, `MatSIBC` | 2026-02-13 | `PEECBuilder`; HDiv-VIM / reduced FEM for magnetic cores | Legacy PEEC conductor API |
| `ObjDrwVTK()`, `exportGeometryToVTK()` | — | NGSolve WebGUI / `GmshPostExport` | Old VTK visualization removed |
| `ObjDivMag()`, `ObjDivMagPln()`, `ObjCutMag()` | — | Netgen / Cubit | Mesh operations use external tools |
| `FldVTS()` | 2026-03-22 | NGSolve + `GmshPostExport` | Field visualization removed |
| `beam_tracking` module | 2026-03-22 | `radia.xsuite_bridge` + CERN Xsuite | Old in-tree engine removed; Radia field maps feed Xsuite's spatial Boris integrator through the optional `beam` extra. |
| `radia_pyvista_viewer.py` | — | NGSolve WebGUI | Visualization removed |

### Removed Libraries

| Library | Date | Replacement | Reason |
|---------|------|-------------|--------|
| ExaFMM-t (FMM acceleration, method 3) | 2026-03-06 | HACApK (H-matrix) | FMM removed from repo |
| `GmshBuilder` | 2026-03-13 | Cubit plugin (`export gmsh`) | GMSH is visualization-only, not mesh generation |

### Removed Mesh Import Paths

| Removed Path | Status |
|-------------|--------|
| Nastran BDF import (`.bdf` / `.nas`) | **Removed**. Was never a reliable path. |
| Gmsh `.msh` import for geometry creation | **Not supported**. `.msh` is export-only (visualization / other solvers). |

The **only supported mesh input format** for Radia and NGSolve is **Netgen `.vol`**:

```
Cubit → export netgen "mesh.vol" → NGSolve Mesh("mesh.vol") → netgen_mesh_to_radia()
```

---

## Supported Elements

| Element Type | API | Faces | DOF | Use Case |
|--------------|-----|-------|-----|----------|
| **Extruded Polygon** | `ObjThckPgn()` | N-gon extruded | 3 | General prism shapes |
| **Hexahedron (face charge)** | `ObjHexahedron()` | 6 quad | 6 | Fixed-magnetization field kernels |
| **Tetrahedron** | `ObjTetrahedron()` | 4 tri | 3 | Fixed-magnetization geometry import |
| **Wedge/Prism (face charge)** | `ObjWedge()` | 5 | 5 | Fixed-magnetization hybrid meshes |
| **Pyramid (face charge)** | `ObjPyramid()` | 5 | 5 | Mesh transitions for field evaluation |
| **General** | `ObjPolyhdr()` | custom | 3/5/6 | Arbitrary polyhedra |

**DOF (Degrees of Freedom)**:
- **Face-charge elements**: 6 DOF for hexahedra, 5 DOF for wedge/pyramid (sigma per face)
- **Fixed magnetization elements**: 3-vector magnetization for permanent-field evaluation only.
- Soft-iron and nonlinear magnetic-material solves are HDiv-VIM / reduced-FEM workflows.
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

Use dedicated constructors (`ObjTetrahedron`, `ObjHexahedron`, `ObjWedge`, `ObjPyramid`) when available; reserve `ObjPolyhdr` for custom polyhedra.

```python
vertices = [[0,0,0], [1,0,0], [0.5,0.866,0], [0,0,1], [1,0,1], [0.5,0.866,1]]
wedge = rad.ObjWedge(vertices, [0, 0, 1e6])
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
result = rad.Solve(obj, tolerance, max_iter, method=0)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `obj` | int | Object or container |
| `tolerance` | float | Convergence threshold (0.001 = 0.1%) |
| `max_iter` | int | Maximum iterations |
| `method` | int | Legacy C++ relaxation method; only `0` = LU is supported |

| Returns | Description |
|---------|-------------|
| `result[0]` | Final residual |
| `result[3]` | Number of iterations |

### Solver Routing

| Model | Route | Call |
|-------|-------|------|
| Mesh-backed soft iron | FEEC HDiv-VIM with its symmetric CG policy | `rad.Solve(system)` |
| Legacy C++ relaxation object | Dense LU | `rad.Solve(obj, 0.001, 100, 0)` |

The old non-symmetric BiCGSTAB and HACApK relaxation methods were retired.
HDiv-VIM owns its own matrix compression, preconditioner, and symmetric linear
solver choices; those are not selected through the legacy `method` integer.

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
rad.Solve(obj, nonl_tol, max_iter, method=0)  # legacy C++ route

# 2. Solver parameters (all in one call)
rad.SolverConfig(bicgstab_tol=1e-4, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
```

| Parameter | Keyword | Default | Description |
|-----------|---------|---------|-------------|
| `nonl_tol` | `rad.Solve(obj, nonl_tol, ...)` | 0.001 | Nonlinear convergence threshold |
| `bicgstab_tol` | `rad.SolverConfig(bicgstab_tol=...)` | 1e-4 | BiCGSTAB relative residual tolerance |
| `hacapk_eps` | `rad.SolverConfig(hacapk_eps=...)` | 1e-4 | H-matrix ACA compression tolerance |

**Example - Full solver configuration:**

```python
import radia as rad

# Configure all solver parameters in one call BEFORE Solve()
rad.SolverConfig(
    bicgstab_tol=1e-4,
    hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0,
    newton_method=True,
    newton_damping=True, newton_damping_max_iter=5, newton_damping_min_omega=0.01,
    relax_param=0.0
)

# Solve
rad.Solve(grp, 0.001, 100, 0)  # legacy C++ LU route

# Query all settings
config = rad.GetSolverConfig()
print(config)
```

### SolverConfig - Unified Solver Parameter Configuration

```python
rad.SolverConfig(**kwargs)
```

All parameters are optional keyword arguments. Only specified parameters are changed.

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `hacapk_eps` | float | 1e-4 | ACA+ compression tolerance (Method 2 only) |
| `hacapk_leaf` | int | 10 | Minimum cluster size in elements |
| `hacapk_eta` | float | 2.0 | Admissibility parameter |
| `hmatrix_eps` | float | - | H-matrix field evaluation epsilon |
| `bicgstab_tol` | float | 1e-4 | BiCGSTAB relative residual tolerance |
| `relax_param` | float | 0.0 | Under-relaxation (0=full step, <1=damped) |
| `newton_method` | bool | False | True=Newton-Raphson, False=Picard |
| `newton_damping` | bool | True | Enable Newton line search damping |
| `newton_damping_max_iter` | int | 5 | Max line search iterations |
| `newton_damping_min_omega` | float | 0.01 | Minimum damping factor |

**HACApK Parameter Rationale:**

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `hacapk_eps` | 1e-4 | Balance between accuracy and compression. Lower values (1e-6, 1e-8) for higher accuracy. |
| `hacapk_leaf` | 10 | For face-charge hexahedra, leaf_size=10 gives roughly 60 DOF/leaf (binary tree splitting). |
| `hacapk_eta` | 2.0 | Admissibility: clusters are "well-separated" when dist >= eta * max(diam). eta=2.0 is conservative. |

**Examples:**

```python
# Configure retained low-level H-matrix parameters for APIs that use them
rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

# Under-relaxation for difficult nonlinear problems
rad.SolverConfig(relax_param=0.3)
rad.Solve(container, 0.001, 100, 0)
rad.SolverConfig(relax_param=0.0)  # Reset to full step

# Newton-Raphson with line search damping (recommended for large problems)
rad.SolverConfig(newton_method=True, newton_damping=True)
rad.Solve(container, 0.001, 100, 0)  # LU + Newton + damping
```

### GetSolverConfig - Query All Solver Parameters

```python
config = rad.GetSolverConfig()
```

Returns a dictionary with all current solver parameters:

| Key | Type | Description |
|-----|------|-------------|
| `'bicgstab_tol'` | float | BiCGSTAB tolerance |
| `'relax_param'` | float | Under-relaxation coefficient |
| `'newton_method'` | bool | Newton-Raphson enabled |
| `'newton_damping'` | bool | Newton damping enabled |
| `'newton_damping_max_iter'` | int | Max line search iterations |
| `'newton_damping_min_omega'` | float | Minimum damping factor |
| `'hacapk_stats'` | dict | H-matrix statistics (only after HACApK solve) |

**Example:**
```python
config = rad.GetSolverConfig()
print(f"BiCGSTAB tol: {config['bicgstab_tol']}")
print(f"Newton: {config['newton_method']}, damping: {config['newton_damping']}")
if 'hacapk_stats' in config:
    stats = config['hacapk_stats']
    print(f"Compression: {stats['compression']:.1%}, Memory: {stats['memory_mb']:.1f} MB")
```

---

## Parallelization (TaskManager)

Radia uses NGSolve's TaskManager for all parallelism. There is no OpenMP dependency.

### Thread Control

Thread count is determined by the TaskManager (default: all available cores).

```python
import ngsolve
ngsolve.SetNumThreads(8)  # Set thread count BEFORE solving
```

**Import order**: `import radia` before `import ngsolve`. Radia's `__init__.py`
imports ngsolve internally for DLL resolution (ngcore.dll), which initializes the
TaskManager.

### Solver-Specific Parallelization

| Solver | Method | Parallelization |
|--------|--------|-----------------|
| LU (method=0) | MKL `dgesv_` | `SuspendTaskManager` + `MKLThreadGuard` (MKL multi-threading) |
| HDiv-VIM | Symmetric CG with its selected preconditioner | Caller-owned NGSolve `TaskManager` plus native C++ kernels |

### GetSolveStats - Query Solve Statistics

```python
stats = rad.GetSolveStats()
```

Returns a dictionary after `rad.Solve()`:

| Key | Type | Description |
|-----|------|-------------|
| `'t_matrix_build'` | float | Matrix construction time [s] |
| `'t_linear_solve'` | float | Linear solver time [s] |
| `'t_lu_decomp'` | float | LU decomposition time [s] (LU only) |
| `'t_hmatrix_build'` | float | Retained low-level H-matrix timing counter |
| `'linear_iterations'` | int | Retained low-level linear iteration counter |
| `'nonl_iterations'` | int | Nonlinear iterations |
| `'num_threads'` | int | Number of TaskManager threads |
| `'taskmanager_enabled'` | bool | TaskManager was active |

### Other Parallelized Operations

| Operation | Parallelization |
|-----------|-----------------|
| Interaction matrix build | `ParallelFor` |
| Field computation (`Fld`, `FldLst`) | `ParallelFor` |
| Analytical polygon integrals | `ParallelFor` |

### NGSolve Compatibility

Radia and NGSolve share the same TaskManager instance, so they coexist naturally:

```python
import radia as rad
from ngsolve import *

with TaskManager():
    rad.Solve(grp)
    V_op = LaplaceSL(j_trial.Trace() * ds) * j_test.Trace() * ds  # Also uses TaskManager
```

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

### Fld - Unified Field Computation (Batch)

`rad.Fld()` auto-detects single-point vs batch from input shape:

```python
# Single point: shape (3,) -> returns array
B = rad.Fld(obj, 'b', np.array([0, 0, 0.1]))

# Batch: shape (N, 3) -> returns (N, 3) array
points = np.array([[0, 0, 0.1], [0, 0, 0.2], [0, 0, 0.3]])
B_batch = rad.Fld(obj, 'b', points)   # (3, 3) array of B values
H_batch = rad.Fld(obj, 'h', points)   # (3, 3) array of H values
A_batch = rad.Fld(obj, 'a', points)   # (3, 3) array of A values
phi_batch = rad.Fld(obj, 'phi', points)  # (3,) array of scalar potential
```

| Field type | Returns | Units |
|------------|---------|-------|
| `'b'` | B field | Tesla |
| `'h'` | H field | A/m |
| `'a'` | Vector potential A | T*m |
| `'phi'` | Scalar potential | A |
| `'m'` | Magnetization | A/m |

**Note**: `FldBatch`, `FldA`, `FldPhi` are removed. Use `Fld()` with appropriate field type and batch points.

### ClassifyPoints - Point Classification (v1.3.16+)

```python
result = rad.ClassifyPoints(obj, points, near_threshold=3.0)
```

Classifies evaluation points relative to mesh elements (for batch field computation).

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

### Supported Mesh Path

The **only** supported mesh import path for Radia magnetostatic analysis is:

```
Cubit → export netgen "mesh.vol" → NGSolve Mesh("mesh.vol") → netgen_mesh_to_radia()
```

Nastran BDF import and Gmsh `.msh` import are **not supported** for Radia/NGSolve input.
The only input format is Netgen `.vol`.

### NGSolve Mesh Access Policy (MANDATORY)

**CRITICAL**: All NGSolve mesh access MUST use functions from `netgen_mesh_import.py`.

| Rule | Description |
|------|-------------|
| **ALWAYS** | Use `netgen_mesh_to_radia()` or `extract_elements()` |
| **NEVER** | Directly access `mesh.ngmesh.Points()`, `mesh.vertices[]`, or `el.vertices[].nr` |

**Why?** NGSolve has TWO different indexing schemes:

| Access Method | Indexing | Valid Range |
|--------------|----------|-------------|
| `mesh.ngmesh.Points()[i]` | **1-indexed** | 1 to nv |
| `mesh.vertices[i]` | **0-indexed** | 0 to nv-1 |
| `el.vertices[i].nr` | Returns **0-indexed** | Use with `mesh.vertices[]` only |

Mixing these causes off-by-one errors that are difficult to debug.

### netgen_mesh_to_radia - Mesh to Radia Geometry

```python
from ngsolve import Mesh
from radia.netgen_mesh_import import netgen_mesh_to_radia

mesh = Mesh("mesh.vol")
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')
```

### extract_elements - Custom Processing

```python
from radia.netgen_mesh_import import extract_elements, compute_element_centroid

elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']
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
# 1. Import radia (includes RadiaField CoefficientFunction)
import radia as rad

# 2. Import ngsolve
import ngsolve
from ngsolve import *

# RadiaField is available as rad.RadiaField (no separate module needed since v2.5.0)
```

### NGSolve Version Requirement

**Use NGSolve 6.2.2603 or later**. Version 6.2.2406~6.2.2501 had a Periodic BC regression. Version 6.2.2603+ includes `curvedelements` Load, hex/prism curving, and Periodic BC fix.

```bash
pip install ngsolve>=6.2.2603
```

### RadiaField - CoefficientFunction (v2.5.0+)

`RadiaField` is an NGSolve `CoefficientFunction` subclass that evaluates Radia fields at arbitrary points. It is now part of `_radia_pybind.pyd` (no separate `radia_ngsolve.pyd` needed).

```python
cf = rad.RadiaField(radia_obj, field_type='b', origin=None, u=None, v=None, w=None,
                    precision=None, units='m')
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `radia_obj` | int | Radia object ID |
| `field_type` | str | `'b'`, `'h'`, `'a'`, `'m'`, or `'phi'` |
| `origin` | [x,y,z] | Origin for coordinate transform (optional) |
| `u`, `v`, `w` | [x,y,z] | Local axes for coordinate transform (optional) |
| `precision` | float | Field computation precision (optional) |
| `units` | str | Must be `'m'` |

**Dimensions**: `dim=3` for B/H/A/M (vector), `dim=1` for phi (scalar).

**Batch evaluation**: When NGSolve calls `gf.Set(cf)`, the batch `Evaluate(mir, result)` is used internally. This calls `rad.Fld(obj, field_type, points)` with batch points (shape N x 3), which is TaskManager-parallelized in C++.

```python
# Basic usage
B_cf = rad.RadiaField(magnet, 'b')
fes = HDiv(mesh, order=2)  # BDM2; use RT=True for Raviart--Thomas
gf = GridFunction(fes)
gf.Set(B_cf)

# Scalar potential
phi_cf = rad.RadiaField(magnet, 'phi')
fes_h1 = H1(mesh, order=2)
gf_phi = GridFunction(fes_h1)
gf_phi.Set(phi_cf)

# With coordinate transform (rotated magnet)
B_rotated = rad.RadiaField(magnet, 'b',
    origin=[0.1, 0, 0],
    u=[cos_a, sin_a, 0], v=[-sin_a, cos_a, 0])
```

#### VoxelCoefficient (`as_voxel_cf`)

Creates a trilinearly-interpolated voxel grid for fast repeated evaluation (e.g., particle trajectory tracking).

```python
B_voxel = B_cf.as_voxel_cf(mesh, resolution=61)
# B_voxel is a VoxelCoefficient - fast evaluation at any point
```

#### Point Cache

Pre-cache field values at known points for faster `gf.Set()`:

```python
B_cf.PrepareCache(integration_points)
gf.Set(B_cf)  # Uses cache
stats = B_cf.GetCacheStats()  # {'enabled': True, 'hits': ..., 'hit_rate': ...}
B_cf.ClearCache()
```

---

## Utilities

### UtiDelAll - Clear Memory

```python
rad.UtiDelAll()
```

### UtiVer - Version

```python
version = rad.UtiVer()
```

---


---

## PEEC Solver

The PEEC (Partial Element Equivalent Circuit) solver provides circuit parameter extraction from conductor geometries. The solver is implemented in C++ with MKL LAPACK/BLAS and exposed via pybind11.

### Overview

| Module | Description |
|--------|-------------|
| `peec_matrices.PyPEECBuilder` | C++ matrix builder (L, R, P, M_LS) |
| `peec_matrices.MNASolver` | C++ MNA multi-port solver (LAPACK zgesv_) |
| `peec_topology.PEECCircuitSolver` | Python API for port impedance/coupling |
| `fasthenry_parser.FastHenryParser` | FastHenry .inp file parser |
| Magnetic core coupling | HDiv-VIM / reduced FEM workflow |

### PyPEECBuilder - Matrix Construction

Build PEEC matrices from node-segment topology:

```python
from peec_matrices import PyPEECBuilder

builder = PyPEECBuilder()

# Add nodes at 3D positions
n1 = builder.add_node_at(0, 0, 0)         # Returns node ID (int)
n2 = builder.add_node_at(0.05, 0, 0)
n3 = builder.add_node_at(0.1, 0, 0)

# Add segments connecting nodes
builder.add_connected_segment(n1, n2, w=1e-3, h=1e-3, sigma=5.8e7)
builder.add_connected_segment(n2, n3, w=1e-3, h=1e-3, sigma=5.8e7)

# Multi-filament: subdivide cross-section for skin/proximity effect
builder.add_connected_segment(n1, n2, w=3e-3, h=3e-3, sigma=5.8e7, nwinc=3, nhinc=3)

# Add ports (measurement terminals)
port_id = builder.add_port(n1, n3)

# Build matrices with topology information
topo = builder.build_topology()
# Returns dict with: L, R, P, segment_nodes, n_nodes, ports, etc.
```

**add_node_at Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `x, y, z` | float | 3D position (meters) |

**add_connected_segment Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `node_from` | int | required | Start node ID |
| `node_to` | int | required | End node ID |
| `w` | float | required | Width (meters) |
| `h` | float | required | Height (meters) |
| `sigma` | float | 5.8e7 | Conductivity (S/m) |
| `seg_type` | str | 'normal' | Segment type |
| `nwinc` | int | 1 | Width subdivisions |
| `nhinc` | int | 1 | Height subdivisions |

### PEECCircuitSolver - Port Impedance

Solve for port impedance using MNA (Modified Nodal Analysis):

```python
from peec_topology import PEECCircuitSolver

solver = PEECCircuitSolver(topo)

# Single frequency
Z = solver.compute_port_impedance(freq=1e6)
R = Z.real   # Resistance (Ohm)
L = Z.imag / (2 * np.pi * 1e6)  # Inductance (H)

# Multi-port Z-matrix
Z_matrix = solver.compute_Z_matrix(freq=1e6)  # (n_ports x n_ports) complex

# Coupling coefficient
result = solver.compute_coupling_coefficient(freq=1e6)
k = result['k']      # Coupling coefficient
L1 = result['L1']    # Self-inductance port 1
L2 = result['L2']    # Self-inductance port 2
M = result['M']      # Mutual inductance

# Frequency sweep
freqs = np.logspace(2, 6, 50)
Z_sweep = solver.frequency_sweep(freqs)

# Frequency sweep with surface impedance (Zs callback)
def Zs_func(freq, n_loop):
    """Return per-filament surface impedance array."""
    omega = 2 * np.pi * freq
    delta = np.sqrt(2 / (omega * MU_0 * sigma))
    Zs_val = (1 + 1j) / (sigma * delta)
    return np.full(n_loop, Zs_val)

Z_sweep = solver.frequency_sweep(freqs, Zs_func=Zs_func)

# Multi-port frequency sweep
Z_matrix_sweep = solver.frequency_sweep_multiport(freqs)
```

**Solver Method Selection**:

```python
solver.set_solver_method(0)  # LU (LAPACK zgesv_, default)
solver.set_solver_method(1)  # BiCGSTAB (templated, MKL BLAS)

solver.set_bicgstab_params(tol=1e-10, max_iter=1000)
```

### FastHenryParser - FastHenry .inp Import

Parse FastHenry input files and solve:

```python
from fasthenry_parser import FastHenryParser

parser = FastHenryParser()

# Parse from file or string
parser.parse_file('inductor.inp')
# or:
parser.parse_string("""
.Units mm
.default sigma=5.8e7

N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1 nwinc=3 nhinc=3

.external N1 N2
.freq fmin=100 fmax=1e6 ndec=5
.end
""")

# Inspect model
print(parser.get_summary())
freqs = parser.get_frequencies()

# Convert to PEECBuilder
builder = parser.to_peec_builder()
topo = builder.build_topology()

# Or one-step solve
result = parser.solve()
# Returns: {'freqs': array, 'Z_port': array, 'R': array, 'L': array, 'topology': dict}
```

**Supported Directives**:

| Directive | Example | Description |
|-----------|---------|-------------|
| `.Units` | `.Units mm` | Length unit |
| `N<name>` | `N1 x=0 y=0 z=0` | Node definition |
| `E<name>` | `E1 N1 N2 w=1 h=1` | Segment definition |
| `.external` | `.external N1 N2` | Port definition |
| `.freq` | `.freq fmin=1e3 fmax=1e6 ndec=5` | Frequency sweep |
| `.default` | `.default w=1 h=1 sigma=5.8e7` | Default parameters |
| `.equiv` | `.equiv N1 N3` | Node merge |
| `.magnetic` | Rejected by `solve()` | Magnetic core input |

**Magnetic core inputs** are parsed only for diagnostics, and
`FastHenryParser.solve()` rejects them. Use HDiv-VIM / reduced FEM for
magnetic cores and keep PEEC for conductor/shield circuit extraction.

### Magnetic Material Coupling Policy

PEEC is conductor/shield oriented. Magnetic material coupling uses the
HDiv-VIM / reduced-FEM workflow, with field exchange at the application layer.

```python
# conductor PEEC: PEECCircuitSolver(topology_dict)
# magnetic core: radia.vim / reduced FEM
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
| `ESIMVTKOutput` | Field export for visualization |

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

### ESIM VTK Export

Export ESIM results for visualization (requires `esim_vtk_export` module):

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

### Image Symmetry (IMA)

IMA exploits mirror symmetry to reduce problem size. The `image=` parameter is
supported by the current HDiv-VIM route and by the retained LU interaction-matrix
route through `rad.Solve()` and `rad.BuildMatrix()`.

**Note**: The old API (`TrfMlt`, `SetIMASymmetry`, `BuildIMAMatrix`, `PreRelax`, `Image`) has been
removed (2026-01-31). Use the unified `image=` parameter instead.

```python
# Quarter model: x-mirror (symmetric) + z-mirror (antisymmetric)
rad.Solve(container, image='+x-z')

# Pre-build matrix with IMA (for inspection)
handle = rad.BuildMatrix(model, image='+x-z')
matrix, dof = rad.GetInteractMatrix(handle)
```

**IMA Sign Selection Policy**:

| Field vs Mirror Plane | IMA Sign | Description |
|----------------------|----------|-------------|
| Field **parallel** to mirror | **+** (symmetric) | Field tangent to mirror plane |
| Field **perpendicular** to mirror | **-** (antisymmetric) | Field normal to mirror plane |

**Symmetry Combinations**:

| `image=` | Model Reduction | Example Use |
|----------|----------------|-------------|
| `'+x'` | Half model | Single mirror |
| `'+x-z'` | Quarter model | Two mirrors (e.g., C-type electromagnet) |
| `'+x+y-z'` | Eighth model | Three mirrors |

**Matrix Construction**:
- Symmetric (+): `N_Image[i,j] = N[i,j] + N[i, mirror_j] @ P`
- Antisymmetric (-): `N_Image[i,j] = N[i,j] - N[i, mirror_j] @ P`

**IMA Boundary Element Note**:
Image symmetry is a solver contract of the HDiv-VIM / charge-Gram route.  Keep
full-model and image-model tests in sync for boundary elements lying on symmetry
planes; do not use the retired tetrahedral moment path as an exception.

---

## Common Issues

### 1. Coordinates Off by 1000x

**Cause**: Legacy scripts may have used millimeters. Radia always uses meters now.

**Solution**: Ensure all coordinates are specified in meters. `FldUnits()` has been removed — Radia always uses meters.

### 2. DLL Load Failed

**Cause**: Wrong import order or missing MKL DLLs

**Solution**: Import ngsolve before radia. Since v2.5.0, RadiaField is in the main radia module (no separate radia_ngsolve needed).

### 3. ObjPolyhdr Face Error (Internal API)

**Cause**: 0-indexed faces when using internal ObjPolyhdr API

**Solution**: Use **1-indexed** faces (Radia convention). For Python users, prefer `ObjHexahedron` and `ObjTetrahedron` which auto-generate faces.

### 4. Solver Not Converging

**Solutions**:
1. Confirm the model uses mesh-backed soft iron and therefore reaches HDiv-VIM
2. Increase the HDiv-VIM iteration limit through its named API
3. Check B-H data is monotonic
4. Verify H-M conversion: `M = B/mu_0 - H`

---

## Units

### Unit System

Radia always uses SI units (meters) internally.

| Quantity | Unit | Notes |
|----------|------|-------|
| Length | m | Always meters |
| B (flux density) | Tesla (T) | SI |
| H (field) | A/m | SI |
| M (magnetization) | A/m | SI |
| A (vector potential) | T*m | SI |
| Current | Ampere (A) | SI |
| Current density | A/m^2 | SI |

### Design Principles

1. **Always meters**: All coordinates are in meters. `FldUnits()` has been removed.
2. **Field values are always SI**: B, H, and A are always in SI units (T, A/m, T*m)
3. **No field scaling**: Field values (B, H, A) are always in SI units

### Maxwell Relation: B = curl(A)

With the new SI internal units, the Maxwell relation `B = curl(A)` is satisfied without any unit conversion:

```python
import radia as rad
import numpy as np

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

See `validation_test/ngsolve_integration/verify_curl_A_equals_B/` for a complete verification script that:

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

# Define permanent magnet as background field
pm = CuboidMagnet(
    center=[0, 0, 50],      # 50mm above center
    dimensions=[40, 40, 20],
    magnetization=[0, 0, 955000]
)

# Create Radia background field object
bkg = rad.ObjBckgCF(pm)  # Uses pm.__call__() which returns get_B()

# Create mesh-backed soft iron to solve
from radia.vim import soft_iron_box
iron = soft_iron_box(
    center=(0.0, 0.0, 0.0), size=(0.04, 0.04, 0.02), mu_r=1000.0, nsub=4)

grp = rad.ObjCnt([iron, bkg])
rad.Solve(grp)
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

### Triangle B-field / Solid Angle Kernels

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
21. [validation_test/cube_uniform_field/](../validation_test/cube_uniform_field/) - Cube uniform-field benchmark corpus

---

**Last Updated**: 2026-02-16
**License**: LGPL-2.1 (modifications), BSD-style (original RADIA from ESRF)
