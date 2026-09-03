"""
Radia library usage knowledge base for Radia MCP server.

Covers the core Radia BEM library API: geometry creation, materials,
solving, field computation, mesh import, and best practices.
"""

RADIA_OVERVIEW = """
# Radia Library Overview

Radia is a boundary element method (BEM) / integral method library for
computing magnetostatic fields in unbounded domains. Written in C++ with
Python bindings (pybind11).

## Key Advantages over FEM

- **No air mesh**: Only discretize magnetic materials, not surrounding air
- **Naturally open boundaries**: No truncation or absorbing layers
- **Exact for permanent magnets**: Analytical integrals for simple shapes
- **Efficient for thin structures**: No volume meshing of thin conductors

## Core Workflow

```
1. rad.UtiDelAll()       # Clear previous objects
2. Create geometry       # ObjHexahedron, ObjRecMag, etc.
3. Apply materials       # MatLin, MatSatIsoFrm, etc.
4. rad.Solve(...)        # Compute magnetization
5. rad.Fld(...)          # Evaluate fields
6. rad.UtiDelAll()       # Cleanup
```

Note: Radia always uses meters. All coordinates are in meters.

## Module Structure

- `import radia as rad` - Main module (re-exports from _radia_pybind)
- C++ backend: pybind11 bindings in _radia_pybind.pyd
"""

RADIA_GEOMETRY = """
# Geometry Creation

## Simple Objects

| Function | Description | DOF Type |
|----------|-------------|----------|
| `ObjRecMag(center, dims, M)` | Rectangular permanent/current magnet | fixed magnetization |
| `ObjHexahedron(vertices, M)` | 8-vertex permanent/current magnet | fixed magnetization |
| `ObjTetrahedron(vertices, M)` | 4-vertex permanent/current magnet | fixed magnetization |
| `ObjWedge(vertices, M)` | 6-vertex wedge/prism | fixed magnetization |
| `ObjPyramid(vertices, M)` | 5-vertex square-base pyramid | fixed magnetization |
| `ObjThckPgn(z, dz, polygon, axis, M)` | Extruded polygon | varies |
| `ObjPolyhdr(vertices, faces, M)` | General polyhedron | varies |

## Soft-Iron Policy

- Use `radia.vim.MeshSoftIron` / `radia.Solve(..., demag_backend="hdiv")`
  for mesh-backed soft iron.
- Raw mesh-less polyhedra remain useful for fixed magnetization and field
  evaluation, but they are not the production soft-iron solve path.
- For Cubit/Netgen meshes, keep the NGSolve mesh as the owner of material
  labels and finite-element spaces, then route the solve through HDiv-VIM.

## Examples

```python
import radia as rad
import numpy as np

rad.UtiDelAll()

# Rectangular permanent magnet: 20x20x10 mm, Br=1.2T in Z
mag = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.01], [0, 0, 955000])

# Hexahedron with 8 vertices
verts = [
    [-0.01, -0.01, -0.005], [0.01, -0.01, -0.005],
    [0.01, 0.01, -0.005], [-0.01, 0.01, -0.005],
    [-0.01, -0.01, 0.005], [0.01, -0.01, 0.005],
    [0.01, 0.01, 0.005], [-0.01, 0.01, 0.005]
]
hex_elem = rad.ObjHexahedron(verts, [0, 0, 0])
```

## Current Sources

| Function | Description |
|----------|-------------|
| `ObjArcCur(center, [r_min, r_max], [phi_min, phi_max], h, n_sec, j)` | Arc/circular coil |
| `ObjFlmCur(points, current)` | Filament (Biot-Savart) |
| `ObjRaceTrk(center, radii, heights, current, n_seg)` | Racetrack coil |

```python
# Full circular coil: R=50mm, 1mm cross-section, J=1e6 A/m^2
coil = rad.ObjArcCur([0, 0, 0], [0.0495, 0.0505],
                     [-np.pi, np.pi], 0.001, 100, 1e6)
```

## Background Field

```python
# IMPORTANT: ObjBckg requires a CALLABLE, not a list!
MU_0 = 4 * np.pi * 1e-7

# Correct: lambda or function
ext = rad.ObjBckg(lambda p: [0, 0, MU_0 * 200000])

# WRONG - will fail silently or error:
# ext = rad.ObjBckg([0, 0, MU_0 * 200000])

# Non-uniform field
def gradient_field(p):
    x, y, z = p
    G = 10.0  # T/m
    return [G * y, G * x, 0]
ext = rad.ObjBckg(gradient_field)
```

SHOWCASE NOTEBOOK: `docs/background_fields/background_fields.ipynb` -- a
spatially-varying quadrupole `ObjBckg` callback driving a nonlinear soft-iron
cube / sphere, plus a permeability sweep (executed + rendered). Durable run
metadata and script/source hashes are in `docs/background_fields/background_fields_results.json`.

## Containers

```python
container = rad.ObjCnt([obj1, obj2, obj3])
grp = rad.ObjCnt([container, background])
```
"""

RADIA_MATERIALS = """
# Material Definition

## Linear Materials

```python
# Isotropic: scalar mu_r
mat_iron = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(obj, mat_iron)

# Anisotropic: [mu_parallel, mu_perp], easy_axis
mat_aniso = rad.MatLin([5000, 100], [0, 0, 1])
rad.MatApl(obj, mat_aniso)
```

## Nonlinear Materials

### B-H Table (Industry Standard)

```python
# Format: [[H (A/m), B (T)], ...]
BH_STEEL = [
    [0.0, 0.0], [100.0, 0.1], [500.0, 0.8],
    [1000.0, 1.2], [5000.0, 1.7], [100000.0, 2.1]
]
mat = rad.MatSatIsoTab(BH_STEEL)
rad.MatApl(core, mat)
```

### Tanh Formula (Analytical)

```python
# M = ms1*tanh(ksi1*H/ms1) + ms2*tanh(ksi2*H/ms2) + ms3*tanh(ksi3*H/ms3)
mat = rad.MatSatIsoFrm(
    [1596.3, 1.1488],   # [ksi1, ms1]
    [133.11, 0.4268],   # [ksi2, ms2]
    [18.713, 0.4759]    # [ksi3, ms3]
)
rad.MatApl(core, mat)
```

### Standard Material Library

```python
# Pre-defined materials (name-based lookup)
mat = rad.MatStd('NdFeB')   # Neodymium Iron Boron
mat = rad.MatStd('Xc06')    # Low-carbon steel
rad.MatApl(obj, mat)
```

## Hysteresis Materials

### Play Model (Recommended)

```python
from radia.hysteresis_io import hys_to_play_radia

# One-step from JMAG .hys file
K, eta, f_k_tables = hys_to_play_radia('material.hys', K=20)
mat = rad.MatPlayHysteresis(K, eta, f_k_tables)
rad.MatApl(iron, mat)
```

### Energy Model (Legacy)

```python
from radia.hysteresis_io import hys_to_radia

params = hys_to_radia('material.hys', K=20, eps=1e-8)
mat = rad.MatEnergyHysteresis(**params)
rad.MatApl(iron, mat)
```

### State Management (Both Models)

```python
state = rad.MatHysSaveState(mat)       # Save play operator states
rad.MatHysRestoreState(mat, state)     # Restore to saved state
rad.MatHysCommitState(mat)             # Commit as new reference
```

## Material Application

```python
# Apply to single object
rad.MatApl(hex_elem, mat)

# Apply to container (all children)
rad.MatApl(container, mat)
```
"""

RADIA_SOLVING = """
# Solving the Magnetostatic Problem

## rad.Solve()

```python
result = rad.Solve(mesh_backed_system)
```

`rad.Solve` is the convenience dispatch. A registered mesh-backed soft-iron
body routes to FEEC HDiv-VIM. The retained C++ compatibility signature is
`rad.Solve(obj, tolerance, max_iter, method=0, image="")`; its only supported
method is dense LU (`0`).

## Solver Routing

| Model | Entry point | Controls |
|-------|-------------|----------|
| Mesh-backed soft iron | `rad.Solve(system)` or `radia.vim.Solve(mesh, ...)` | named `linear_solver`, `preconditioner`, `gram_eps`, `leaf`, `eta`, tolerances |
| Repeated HDiv geometry | `radia.vim.HDivSolver(mesh).Solve(...)` | same named controls plus operator reuse |
| Legacy C++ relaxation object | `rad.Solve(obj, tol, max_iter, 0)` | dense LU only |

Retired legacy method integers 1 and 2 must not be recommended. HACApK remains
available inside owning HDiv/PEEC/BEM operator routes; it is not a public
`rad.Solve(method=2)` selector.

## Example

```python
import ngsolve as ng
import radia as rad
from radia.vim import soft_iron_box

rad.UtiDelAll()
iron = soft_iron_box(
    center=(0.0, 0.0, 0.0), size=(0.02, 0.02, 0.02),
    mu_r=1000.0, nsub=4,
)
source = rad.ObjBckg(lambda p: [0.0, 0.0, 0.1])
system = rad.ObjCnt([iron, source])
with ng.TaskManager():
    result = rad.Solve(system)
```

## Advanced Solver Configuration

```python
from radia import vim

solver = vim.HDivSolver(mesh, order=2)
with ng.TaskManager():
    result = solver.Solve(
        mu_r=1000.0,
        H_ext=applied_field,
        tol=1e-8,
        maxit=4000,
        linear_solver="auto",
        preconditioner="auto",
    )
```

## Vector-Accurate Hysteresis

Production hysteresis uses the mesh-backed HDiv-VIM history solver. The
constitutive reference direction is B-input (`B -> H`), and the field history
is advanced on one conforming NGSolve mesh:

```python
from radia import vim

with ng.TaskManager():
    result = vim.SolveHysteresis(
        mesh,
        h_steps,
        play=(K, eta, f_k_tables),
        tol=1e-8,
        maxit=4000,
        nl_tol=1e-3,
    )
```

Use `vim.HDivSolver(mesh).SolveHysteresis(...)` when the same geometry is
continued across runs and its operator cache should be reused. Legacy
`SolverConfig(b_input_*)` flags and mesh-less material objects are not the
production soft-iron route. Numerical evidence belongs under
`validation_test/hysteresis/`.

## IMA (Image Method of Analysis) Symmetry

IMA exploits mirror symmetry to reduce problem size (half, quarter, eighth model).
The supported mesh-backed HDiv-VIM geometry routes accept `image=`. The legacy
C++ compatibility route accepts the same string only with dense LU.

```python
# Quarter model with x-mirror (symmetric) and z-mirror (antisymmetric)
# For Z-directed field: Bz is parallel to x-plane (+) and perpendicular to z-plane (-)
rad.Solve(container, image='+x-z')

# Pre-build matrix with IMA (optional, for matrix inspection)
handle = rad.BuildMatrix(model, image='+x-z')
matrix, dof = rad.GetInteractMatrix(handle)
```

### IMA Sign Selection Policy

| Field vs Mirror Plane | IMA Sign |
|----------------------|----------|
| Field **parallel** to mirror | **+** (symmetric) |
| Field **perpendicular** to mirror | **-** (antisymmetric) |

### Supported Symmetry Combinations

| `image=` | Reduction | Example |
|----------|-----------|---------|
| `'+x'` | Half model | x-mirror only |
| `'+x-z'` | Quarter model | x + z mirrors |
| `'+x+y-z'` | Eighth model | x + y + z mirrors |

### IMA Boundary Element Limitation

IMA may produce incorrect results (~0.5x magnitude) for **boundary elements** whose
faces lie ON the symmetry plane, when observation points are also on that plane.

**When IMA is Safe**:
- Elements offset from symmetry planes (not touching)
- Observation points off-plane
- Mesh-backed HDiv models with matching image signs and symmetric cuts
"""

RADIA_PARALLELIZATION = """
# Parallelization Architecture

Radia follows the NGSolve-native execution model: NGSolve TaskManager is the single
threading substrate for Radia + NGSolve workflows. There is no Radia OpenMP thread pool.
The TaskManager is a work-stealing thread pool initialized when `import radia` loads
ngcore.dll.

## Thread Control

Thread count is determined by the TaskManager. By default, it uses all available cores.
To control thread count from Python, use NGSolve's API:

```python
import ngsolve
ngsolve.SetNumThreads(8)  # Set thread count BEFORE solving
```

Or equivalently via the `with TaskManager(n)` context manager in NGSolve scripts.

**IMPORTANT**: `import radia` must come BEFORE `import ngsolve`. This is because
Radia's `__init__.py` imports ngsolve internally for DLL resolution (ngcore.dll),
which initializes the TaskManager. If ngsolve is imported first with a different
thread count, there can be conflicts.

## How Each Solver Uses Parallelism

| Solver | Method | Parallelization |
|--------|--------|-----------------|
| Legacy LU (`method=0`) | MKL `dgesv_` | `SuspendTaskManager` + `MKLThreadGuard` around the external threaded kernel |
| HDiv-VIM | symmetric named linear solver and preconditioner | caller-owned NGSolve `TaskManager` plus native charge-Gram kernels |
| PEEC/BEM H-matrix routes | formulation-owned named solver | owning API reports its backend and TaskManager behavior |

### LU Solver Threading Detail

LU is intentionally not TaskManager-parallel internally. NGSolve sets
`mkl_set_num_threads(1)` globally to prevent conflicts with TaskManager. When Radia's LU
solver calls `dgesv_`, it uses `MKLThreadGuard` RAII to temporarily re-enable
multi-threaded MKL, then `SuspendTaskManager` to avoid contention:

```cpp
// C++ internals (rad_relaxation_methods.cpp)
{
    ngcore::SuspendTaskManager stm;          // Pause TaskManager threads
    radia::MKLThreadGuard mkl_guard(nthreads); // Enable MKL multi-threading
    dgesv_(&n, &nrhs, A, &n, ipiv, b, &n, &info);
}   // Both guards restore original state on scope exit
```

### HDiv-VIM and H-matrix Threading

Call NGSolve finite-element assembly and HDiv-VIM solves under
`with ngsolve.TaskManager():`. Native charge-Gram and H-matrix kernels reuse
that execution substrate. PEEC and BEM routes own their separate solver choice;
do not infer their threading from a retired `rad.Solve` method number.

## Querying Thread Info

```python
import radia as rad

# After rad.Solve():
stats = rad.GetSolveStats()
print(stats['num_threads'])         # e.g., 8
print(stats['taskmanager_enabled']) # True
print(stats['t_matrix_build'])      # Matrix build time [s]
print(stats['t_linear_solve'])      # Linear solve time [s]
print(stats['t_lu_decomp'])         # LU decomposition time [s] (LU only)
print(stats['t_hmatrix_build'])     # retained low-level H-matrix counter
```

## Other Parallelized Operations

| Operation | File | Method |
|-----------|------|--------|
| Interaction matrix build | `rad_interaction.cpp` | `ParallelFor` |
| Field computation (Fld, FldLst) | `rad_field_unified.cpp` | `ParallelFor` |
| Analytical poly integrals | `rad_poly_analytical.cpp` | `ParallelFor` |
| Scalar potential batch | `rad_material_impl.cpp` | `ParallelFor` |
| Vector potential batch | `rad_material_impl.cpp` | `ParallelFor` |

## NGSolve Compatibility

Since Radia and NGSolve share the same TaskManager, they coexist naturally:

```python
import radia as rad  # Initializes TaskManager (loads ngcore.dll)
from ngsolve import *

with TaskManager():
    rad.Solve(grp)
    V_op = LaplaceSL(j_trial.Trace() * ds) * j_test.Trace() * ds
```

No thread conflict because both use the same underlying TaskManager instance.

Scaling and backend comparisons are validation evidence, not timeless MCP
constants. Read the current JSON artifacts under `validation_test/` before
making a performance claim.
"""

RADIA_FIELDS = """
# Field Computation

## Single Point

```python
field = rad.Fld(obj, component, [x, y, z])
```

### Component Options

| Component | Type | Unit | Description |
|-----------|------|------|-------------|
| 'b', 'bx', 'by', 'bz' | vector/scalar | T | Magnetic flux density |
| 'h', 'hx', 'hy', 'hz' | vector/scalar | A/m | Magnetic field strength |
| 'a', 'ax', 'ay', 'az' | vector/scalar | T*m | Vector potential |
| 'p', 'phi' | scalar | A | Scalar potential |
| 'm', 'mx', 'my', 'mz' | vector/scalar | A/m | Magnetization |

### Examples

```python
B = rad.Fld(mag, 'b', [0, 0, 0.05])       # [Bx, By, Bz]
Bz = rad.Fld(mag, 'bz', [0, 0, 0.05])     # Scalar Bz
H = rad.Fld(mag, 'h', [0, 0, 0.05])        # [Hx, Hy, Hz]
```

## Field Along a Line

```python
field_values = rad.FldLst(obj, 'bz', start, end, n_points, 'arg')
# Returns n_points field values along line from start to end
```

## Field Integration

```python
flux = rad.FldInt(obj, component, start, end, axis)
```

## Visualization (NGSolve-Through)

Field visualization uses NGSolve + GmshPostExport (.msh v4.1).
Geometry uses STEP files (CoilBuilder.write_step(), OCC shapes).
rad.FldVTS() is removed. Use rad.Fld() for point evaluation only.

## Querying Magnetization

```python
all_M = rad.ObjM(container)    # [[center, [Mx,My,Mz]], ...]
rad.ObjSetM(obj, [Mx, My, Mz])  # Set magnetization manually
```

## CoilBuilder: Racetrack Coil for C-Type Dipole

CoilBuilder creates racetrack coils from straight + arc segments.

### Key Parameters

- `set_start(position, orientation=None)`: start position and optional 3x3 orientation matrix.
  Default orientation: Y_local = +y (forward), X_local = +x (lateral), Z_local = +z.
  The first `add_straight()` goes along Y_local.
- `set_cross_section(width, height)`: width = radial (in loop plane), height = axial (perpendicular to loop plane).
- `add_arc(radius, arc_angle)`: radius is CENTER-LINE radius, NOT inner radius.
  Inner R = radius - width/2, Outer R = radius + width/2.
- `close()`: optimizes arc angles to close the loop.

### Racetrack Geometry (Rounded Rectangle)

For a racetrack with rounded corners (4 straights + 4 x 90-degree arcs):

```python
# Rectangle center-line half-dimensions
cl_half_a = inner_half_a + cs_width / 2  # a direction
cl_half_b = inner_half_b + cs_width / 2  # b direction

# Straight section lengths = full side minus 2 arc radii
straight_a = 2 * (cl_half_a - r_centerline)
straight_b = 2 * (cl_half_b - r_centerline)
```

### Example: C-Type Dipole Coil

```python
from radia_coil_builder import CoilBuilder
mm = 1e-3

# Coil in x-y plane (default orientation), flux along z
# Straight sections along y, arcs turn in x-y plane
r_cl = 22.5 * mm      # center-line R (inner 5mm + width/2)
cs_w = 35 * mm         # cross-section radial width
cs_h = 105 * mm        # cross-section axial height (z)

coil = (CoilBuilder(current=20000.0)
    .set_start([cl_half_x, -straight_y/2, 0])
    .set_cross_section(width=cs_w, height=cs_h)
    .add_straight(straight_y)
    .add_arc(radius=r_cl, arc_angle=90)
    .add_straight(straight_x)
    .add_arc(radius=r_cl, arc_angle=90)
    .add_straight(straight_y)
    .add_arc(radius=r_cl, arc_angle=90)
    .add_straight(straight_x)
    .add_arc(radius=r_cl, arc_angle=90)
    .close())

coil.write_step("coil.step")       # OCC STEP export
objs = coil.to_radia()             # Radia current objects
container = rad.ObjCnt(objs)
B = rad.Fld(container, 'b', [0, 0, 0])
```

### GMSH Visualization: Full Model from Quarter

```python
import gmsh
gmsh.initialize()
gmsh.merge("yoke_quarter.step")
vols = gmsh.model.getEntities(3)
# Mirror about x=0 and z=0 for full model
cx = gmsh.model.occ.copy(vols); gmsh.model.occ.mirror(cx, 1, 0, 0, 0)
cz = gmsh.model.occ.copy(vols); gmsh.model.occ.mirror(cz, 0, 0, 1, 0)
cxz = gmsh.model.occ.copy(vols)
gmsh.model.occ.mirror(cxz, 1, 0, 0, 0)
gmsh.model.occ.mirror(cxz, 0, 0, 1, 0)
gmsh.model.occ.synchronize()
gmsh.merge("coil.step")           # Add coil overlay
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Mesh.VolumeEdges", 0)
gmsh.fltk.run()
```
"""

RADIA_REMOVED_APIS = """
# Removed APIs (breaking change reference)

Radia is a research codebase with a no-backward-compat policy. When an API is
removed, it is removed outright -- no alias, no deprecation wrapper. If a user
or AI-generated script calls one of these, it will raise `AttributeError` or
similar. Use the listed replacement instead.

## Force / energy / torque (Phase C, 2026-04-16)

| Removed                 | Replacement / Status                                   |
|-------------------------|--------------------------------------------------------|
| `rad.FldEnr(dst, src, SbdPar)`     | Use `rad.FldFrc(obj, rect_shape)` (Maxwell stress tensor on rectangular surface). A new analytical pair-interaction API (`rad.AnalEn/AnalFrc/AnalTrq`) is planned as research -- see `docs/research/magnet_design/FORCE_COMPUTATION_DESIGN.ipynb`. |
| `rad.FldEnrFrc(...)`    | Same -- use `rad.FldFrc` until the new API ships.     |
| `rad.FldEnrTrq(...)`    | Same -- no direct torque API currently; derive from force * arm or wait for `rad.AnalTrq`. |

Why removed: the old implementation relied on `SubdivideItself`-based midpoint
quadrature (physical mesh splitting with user-specified `SbdPar=[kx,ky,kz,kxs,kys,kzs]`).
This is superseded by the plan to use closed-form / HDiv-compatible field
interactions combined with NGSolve's high-order `Integrate(..., order=N)`
quadrature.

Still live: `rad.FldFrc(obj, shape)` and `rad.FldFrcShpRtg(center, size)` --
the Maxwell stress tensor path is independent and unaffected.

## Mesh operations (removed 2026-01-14)

| Removed                 | Replacement                                            |
|-------------------------|--------------------------------------------------------|
| `rad.ObjDivMag(obj, k)` | Use external meshers: Cubit (`export netgen`) or Netgen OCC -> `.vol` -> `netgen_mesh_to_radia`. |
| `rad.ObjDivMagPln(obj, planes)` | Same. Plane-based splits are mesh generation, not Radia's job. |
| `rad.ObjCutMag(obj, plane)` | Same. |

Why removed: Radia is no longer a mesh generator. All meshing goes through
Cubit or Netgen; Radia consumes the resulting `.vol` via `netgen_mesh_to_radia`.

## Serialization (Phase A/B, 2026-04-15 -- 04-16)

| Removed                 | Replacement                                            |
|-------------------------|--------------------------------------------------------|
| `rad.UtiDmp(obj)` / `rad.UtiDmpPrs(s)` | No replacement. Save your Python script; rebuild the Radia object tree on next run. `.rad` save/load is not supported. |
| `rad.ObjGeoLim(obj)` | No replacement. Compute bounding boxes from the input parameters you used to construct the object (vertices, center+size, etc.). |
| `rad.ObjDrwAttr(...)` / `rad.ObjDrwVTK(...)` | Visualization goes through NGSolve + GMSH. Use `rad.RadiaField` CoefficientFunction + `GmshPostExport`. |

## Unit system (removed long ago)

| Removed                 | Replacement                                            |
|-------------------------|--------------------------------------------------------|
| `rad.FldUnits("mm")` / any `FldUnits` call | Radia always uses meters. There is no unit switch. |

## Visualization (removed)

| Removed                 | Replacement                                            |
|-------------------------|--------------------------------------------------------|
| `rad.FldVTS(...)` | Use `rad.Fld(obj, field_type, points)` for point evaluation, then export to `.msh v4.1` via `GmshPostExport` for visualization. |

## Unsafe legacy extrusions (removed 2026-07-20)

| Removed | Replacement |
|---------|-------------|
| `rad.ObjMltExtPgn(...)` | Build the extruded geometry with Netgen or Cubit, then use the NGSolve-native mesh path. |
| `rad.ObjMltExtRtg(...)` | Same; the legacy side-face construction can also become non-planar. |
| `rad.ObjMltExtTri(...)` | Same; the Triangle-backed implementation is no longer bundled. |

The C ABI declarations, definitions, and exports were deleted together with
the pybind11 and MATLAB MEX entries. No compatibility shim remains.

## Policy

When adding a new API or breaking an existing one, the committer MUST update
this section. No code-level compat shim will be added. Users and AI coders
discover the change by querying this MCP tool.
"""


RADIA_MESH_IMPORT = """
# NGSolve Mesh Import

Convert Netgen/NGSolve meshes to Radia geometry objects.

## Key Module

```python
from radia.netgen_mesh_import import netgen_mesh_to_radia
```

## Supported Element Types

| Element | Radia Object | DOF |
|---------|-------------|-----|
| Tetrahedron (ET.TET) | ObjTetrahedron | fixed magnetization export |
| Hexahedron (ET.HEX) | ObjHexahedron | fixed magnetization export |
| Wedge (ET.PRISM) | ObjWedge | 5 DOF |
| Pyramid | ObjPyramid | 5 DOF |

## Basic Usage

```python
import radia as rad
from ngsolve import *
from radia.netgen_mesh_import import netgen_mesh_to_radia

rad.UtiDelAll()

# Create/load NGSolve mesh
from netgen.occ import Box, Pnt, OCCGeometry
box = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.3))

# Convert to Radia
container = netgen_mesh_to_radia(
    mesh,
    material={'magnetization': [0, 0, 0]},
    units='m',
    material_filter='magnetic'
)

# Apply material and solve
mat = rad.MatLin(1000)
rad.MatApl(container, mat)
```

## Cubit .vol Import (Recommended)

```python
from ngsolve import Mesh
from radia.netgen_mesh_import import netgen_mesh_to_radia

# Load .vol exported from Cubit (export netgen "mesh.vol" order N)
mesh = Mesh("mesh.vol")
container = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

## Cubit Direct Import

```python
from radia.netgen_mesh_import import cubit_hex_to_radia

# After Cubit meshing
container = cubit_hex_to_radia(
    cubit_mesh_data,
    magnetization=[0, 0, 0]
)
```

## NGSolve Magnetization -> Radia Open Boundary Field Evaluation

NGSolve FEM solves M(x) inside bounded domains but cannot easily compute fields
in unbounded external regions (needs PML). Radia provides natural open boundary
field evaluation using exact analytical formulas (NOT dipole approximation).

Pipeline:
```
NGSolve FEM Solve -> M per element -> netgen_mesh_to_radia() -> Radia objects -> rad.Fld()
```

IMPORTANT: Do NOT use dipole approximation (m=M*V). Register elements as proper
Radia ObjHexahedron/ObjTetrahedron with solved magnetization. Radia's surface
charge/surface current analytical formulas are exact for constant M per element.

```python
import radia as rad
from ngsolve import *
from radia.netgen_mesh_import import netgen_mesh_to_radia

rad.UtiDelAll()

# 1. NGSolve solves nonlinear problem -> M per element

# 2. Convert mesh to Radia objects with per-element magnetization
def material_from_ngsolve(el_idx):
    M = get_element_magnetization(gf_M, mesh, el_idx)  # user function
    return {'magnetization': M.tolist()}

container = netgen_mesh_to_radia(mesh, material=material_from_ngsolve, units='m')
# No Solve() needed - M is already known from NGSolve

# 3. Evaluate field at arbitrary external points (exact analytical formulas)
B = rad.Fld(container, 'b', [0, 0, 0.1])
```

Why Radia objects, not dipoles:
- Surface charge model is exact for constant M (zero approximation error)
- No distance limitation (dipoles fail at r < 2 * element_size)
- netgen_mesh_to_radia() already supports per-element material via callable
"""

RADIA_BEST_PRACTICES = """
# Radia Best Practices

## 1. Units: Always Meters

Radia always uses meters. All coordinates must be specified in meters.
`rad.FldUnits()` has been removed. Do not call it.

## 2. Always Clean Up

```python
rad.UtiDelAll()     # At start of script AND at end
```

## 3. ObjBckg Requires Callable

```python
# CORRECT
ext = rad.ObjBckg(lambda p: [0, 0, B_value])

# WRONG - silent failure!
ext = rad.ObjBckg([0, 0, B_value])
```

## 4. Choose the Owning Solver Route

- Mesh-backed soft iron: HDiv-VIM through `vim.Solve` or `vim.HDivSolver`.
- PEEC/BEM: use that formulation's named solver and preconditioner API.
- Legacy C++ relaxation: dense LU (`method=0`) only.
- Never select a current backend by the retired legacy method values 1 or 2.

## 5. Magnetization Units

- Permanent magnets: M in A/m (e.g., 955000 A/m for NdFeB, ~1.2T)
- Conversion: M = Br / mu_0

## 6. B-H Table Format

- Always [[H in A/m, B in Tesla], ...]
- Must be monotonically increasing
- Start from [0, 0]

## 7. Hexahedra Must Be Convex

- ObjHexahedron requires convex shapes
- Non-convex hexahedra cause incorrect surface normals
- For complex shapes, use tetrahedra

## 8. NGSolve Integration

- **Import radia BEFORE ngsolve**: Radia's `__init__.py` imports ngsolve
  internally to locate ngcore.dll and initialize the TaskManager thread pool.
  If ngsolve is imported first, thread pool configuration may conflict.
- Radia and NGSolve share the same TaskManager (ngcore), so they coexist
  without thread conflicts.
- Radia always uses meters (all coordinates in meters)
- Use `netgen_mesh_to_radia()` for mesh conversion (not manual extraction)

## 9. Convergence Checks

```python
with ngsolve.TaskManager():
    result = radia.vim.Solve(mesh, mu_r=1000.0, H_ext=H_ext)
print(result["iters"], result["linear_solver"], result["preconditioner"])
```

## 10. Field Evaluation Points

- Do NOT evaluate ON element surfaces (singularity)
- Offset evaluation points slightly from boundaries
- Use rad.FldLst for systematic field profiles

## 11. Keep demonstration and validation evidence separate

`docs/**/*.ipynb` demonstrates what Radia can do. It is executed and stores its
result and WebGUI output in the notebook itself; it does not require a result
JSON sidecar and is not a benchmark.

`validation_test/` owns analytical comparisons, convergence studies,
benchmarks, and paper/conference evidence. Every such run writes a checked JSON
artifact so the result can be inspected without rerunning an expensive solve.

Why JSON:
- Human-readable + diff-able + version-control-friendly.
- Self-contained: params + results + metadata in one file, so a plot or
  comparison can be regenerated months later from the JSON alone.
- Language-agnostic (re-load in Python / MATLAB / a notebook).

What each result JSON should contain (mirrors the lab Benchmark Policy):
```python
import json, os, platform
from datetime import datetime

def save_results(path, name, problem, results):
    data = {
        "timestamp": datetime.now().isoformat(),   # when
        "hostname":  platform.node(),               # which machine
        "experiment": name,                         # what
        "problem":   problem,   # input params (mesh h, freq, material, ...)
        "results":   results,   # list of per-case dicts (the numbers)
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
```
Each per-case dict records the scalars you would otherwise have to
re-run to recover (e.g. ``{"ndof": ..., "t_solve": ..., "iterations":
..., "converged": ..., "L_total": ...}``).

Folder hygiene (the cleanup half of the rule):
- Put all artifacts of one validation campaign under one dedicated
  `validation_test/<topic>/` folder (result JSON, plots, logs, generated meshes).
- Keep the campaign folder SELF-CONTAINED so that, when the results are
  no longer needed, you can ``rm -rf <folder>`` (delete the folder
  whole) and leave NOTHING orphaned elsewhere.
- Do NOT scatter ``results_*.json`` / ``*.png`` across the repo root or
  shared dirs -- that defeats the one-shot folder delete and accumulates
  cruft.

The validation schema and its focused tests, not copied numbers in MCP
knowledge, define the required metadata.
"""


RADIA_PEEC = """
# PEEC (Partial Element Equivalent Circuit)

Radia includes a PEEC solver for conductor analysis: inductance, resistance,
and capacitance extraction from 3D conductor geometries.

## Architecture

```
FastHenry .inp file  -->  FastHenryParser  -->  PEECBuilder (C++)  -->  PEECCircuitSolver
                                                                    -->  HDiv-VIM / reduced-FEM core coupling
                                                                    -->  ShieldedPEECSolver (with BEM shield)
```

## Key Modules

| Module | Class | Purpose |
|--------|-------|---------|
| `peec_matrices.pyd` | `PEECBuilder` | C++ matrix assembly (L, R, P matrices) |
| `peec_topology.py` | `PEECCircuitSolver` | MNA nodal admittance circuit solver |
| Magnetic core coupling | HDiv-VIM / reduced FEM | Keep PEEC conductor-only |
| `peec_shielded.py` | `ShieldedPEECSolver` | Conductor + BEM shield coupling |
| `fasthenry_parser.py` | `FastHenryParser` | FastHenry .inp file parser |

## Quick Start: FastHenry Input

```python
from fasthenry_parser import FastHenryParser

parser = FastHenryParser()
parser.parse_string(\"\"\"
.Units mm
.default sigma=5.8e7

N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1 nwinc=3 nhinc=3

.external N1 N2
.freq fmin=100 fmax=1e6 ndec=5
.end
\"\"\")

result = parser.solve()
print(f"DC: R={result['R'][0]*1e3:.3f} mOhm, L={result['L'][0]*1e9:.1f} nH")
```

## Topology API (Node-Segment)

```python
from radia.peec_matrices import PEECBuilder
from radia.peec_topology import PEECCircuitSolver

builder = PEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.05, 0, 0)
n3 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_connected_segment(n2, n3, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n3)

topo = builder.build_topology()
solver = PEECCircuitSolver(topo)
Z = solver.compute_port_impedance(freq=1e6)
```

## Multi-Filament (Skin/Proximity Effect)

Use `nwinc` and `nhinc` to subdivide conductor cross-sections:

```python
# 5x5 = 25 parallel sub-filaments for skin effect
builder.add_connected_segment(n1, n2, 3e-3, 3e-3, sigma=5.8e7, nwinc=5, nhinc=5)
```

## Magnetic Cores

PEEC handles conductor / shield circuit extraction.  Magnetic cores use
the HDiv-VIM / reduced-FEM workflow and are coupled at the application
layer.  FastHenry `.magnetic` blocks are rejected by the PEEC parser.

## Accuracy Parameters

### CRITICAL: Circular Coil Segment Count

| n_seg | Circle Error | Delta_L Error (vs FEM) |
|-------|-------------|----------------------|
| 16 | ~1% near core | ~26% |
| 32 | ~0.3% | ~15% |
| **64** | **<0.1%** | **~9%** |

**Rule**: Use n_seg >= 64 for circular coils when computing coupling with nearby objects.

### Core Mesh Divisions

| Divisions | Elements | Delta_L Error |
|-----------|----------|---------------|
| 2,2,1 | 4 | ~26% |
| **3,3,2** | **18** | **~9%** |
| 5,5,3 | 75 | ~5% |

### Bessel SIBC for Circular Wire

**CRITICAL**: Use modified Bessel functions `iv` (I0, I1), NOT regular `jv` (J0, J1):

```python
from scipy.special import iv  # CORRECT: modified Bessel

k = np.sqrt(1j * omega * mu * sigma)
Z_internal = (k * length) / (2 * np.pi * radius * sigma) * (iv(0, k*radius) / iv(1, k*radius))
```

Using `jv` gives correct R_ac/R_dc but WRONG SIGN on internal inductance.

### Internal Inductance Double-Counting

When computing air inductance with Neumann formula (GMD-based):
- Neumann GMD already includes internal inductance (Li/4)
- Do NOT add SIBC internal impedance on top: `use_sibc=False`
- Only use SIBC when computing AC impedance (frequency-dependent skin effect)
"""

RADIA_NGBEM_PEEC = """
# PEEC with ngsolve.bem (NGSolve BEM)

ngsolve.bem provides Galerkin BEM assembly for PEEC impedance extraction.
The key insight is that the PEEC Loop-Star decomposition maps directly
to NGSolve function spaces: `HDivSurface` (Loop) x `SurfaceL2` (Star).

## Architecture

```
Surface Mesh (Netgen OCC)
  |
  +-> NGBEMPEECSolver (ngsbem_peec.py)
  |     - L matrix: LaplaceSL on HDivSurface (inductance)
  |     - P matrix: SingleLayerPotential on SurfaceL2 (capacitance)
  |     - M_LS:     div coupling (FEM, not BEM)
  |     - R:        sheet resistance from sigma, thickness
  |     - solve_frequency() -> Z(f) = R + jwL (MQS) or full Loop-Star
  |
  +-> ShieldBEMSIBC (ngsbem_eddy.py)
  |     - BEM + SIBC for conducting shields
  |     - compute_impedance_matrix() -> Delta_Z (shield coupling)
  |
  +-> Application-specific magnetic coupling
        - Keep PEEC/BEM impedance extraction separate from mesh-backed
          HDiv-VIM soft-iron solves.
```

## Key Modules

| Module | Class | Purpose |
|--------|-------|---------|
| `ngsbem_peec.py` | `NGBEMPEECSolver` | Galerkin PEEC: L, P, M_LS assembly + impedance sweep |
| `ngsbem_eddy.py` | `ShieldBEMSIBC` | Loop-only BEM+SIBC for conducting shields |
| `ngsbem_eddy.py` | `LoopBasisBuilder` | Div-free loop basis from face-edge topology |
| `ngsbem_interface.py` | `extract_edge_geometry` | Bridge: ngsolve.bem mesh to PEEC topology |

## Quick Start: Plate Impedance

```python
from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh
import numpy as np

# Create conductor surface mesh
mesh = create_plate_mesh(0.01, 0.01, maxh=0.003, label="conductor")

# Assemble PEEC matrices (Galerkin BEM)
solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                          sigma=5.8e7, thickness=35e-6, order=0)
solver.assemble()

# MQS frequency sweep (10 Hz - 1 MHz)
freqs = np.logspace(1, 6, 50)
Z = solver.solve_frequency(freqs, mode='mqs')
L_nH = np.imag(Z) / (2 * np.pi * freqs) * 1e9
```

## CRITICAL: .Trace() Requirement for LaplaceSL on HDivSurface

When assembling BEM operators (LaplaceSL, HelmholtzSL) on HDivSurface,
you MUST use `.Trace()` on both trial and test functions:

```python
# CORRECT - with .Trace()
j_trial, j_test = fes_J.TnT()
V_op = LaplaceSL(j_trial.Trace() * ds(label)) * j_test.Trace() * ds(label)

# WRONG - without .Trace() -> corrupted boundary-edge DOFs
V_op = LaplaceSL(j_trial * ds(label)) * j_test * ds(label)
```

**Bug symptom**: Without `.Trace()`, boundary-edge DOFs (RWG basis functions
with support on only 1 triangle) get corrupted diagonal entries (e.g., -1.1e+11
instead of ~1e-10). This causes wildly wrong inductance values.

**Diagnostic**: Extract dense matrix and check diagonal:
```python
M = extract_dense(V_op.mat, n_J)
diag = np.real(np.diag(M))
print(f"min diag = {diag.min():.3e}")  # Should be positive (~1e-10)
# If any diagonal < 0 or >> 1e-5: missing .Trace()
```

**Note**: `.Trace()` is a purely mathematical operation (tangential trace
projection for H(div) functions). It is NOT related to physical thickness.
Same `.Trace()` requirement applies to HelmholtzSL.

Fixed in NGBEMPEECSolver (2026-02-23).

## Air-Core Coil Inductance (Standalone ngsolve.bem)

Direct BEM inductance computation WITHOUT the PEEC framework.
Uses Hodge decomposition to extract the harmonic current mode on a genus-1
surface (ring/frame), then computes L via the L/R ratio formula.

### Single Loop Pattern

```python
import numpy as np
from scipy.linalg import null_space
from ngsolve import Mesh, HDivSurface, SurfaceL2, BilinearForm, ds, BND, TaskManager, InnerProduct
from ngsolve import div as ng_div
from ngsolve.bem import LaplaceSL

MU_0 = 4e-7 * np.pi

# 1. Create ring mesh (circular or rectangular frame)
mesh = create_ring_mesh(R_center, trace_width, maxh)
fes_J  = HDivSurface(mesh, order=0)   # RWG edge-based
fes_L2 = SurfaceL2(mesh, order=0)     # face-based
n_J, n_v, n_f = fes_J.ndof, mesh.nv, fes_L2.ndof

# 2. Divergence matrix D (face x edge)
u_J, q_L2 = fes_J.TrialFunction(), fes_L2.TestFunction()
bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
bf_D += ng_div(u_J.Trace()) * q_L2 * ds
bf_D.Assemble()
D = extract_rect(bf_D.mat, n_f, n_J)

# 3. Incidence matrix C (edge x vertex)
C = np.zeros((n_J, n_v))
for e_idx, edge in enumerate(mesh.edges):
    verts = list(edge.vertices)
    C[e_idx, verts[0].nr] = -1
    C[e_idx, verts[1].nr] = +1

# 4. Mass matrix M_J
u2, v2 = fes_J.TnT()
bf_M = BilinearForm(fes_J)
bf_M += InnerProduct(u2.Trace(), v2.Trace()) * ds
bf_M.Assemble()
M_J = np.real(extract_dense(bf_M.mat, n_J))

# 5. Hodge decomposition: harmonic mode = null([D; C^T @ M_J])
constraint = np.vstack([D, C.T @ M_J])
c_h = null_space(constraint, rcond=1e-10)[:, 0]  # genus-1 -> 1 mode
energy = c_h @ M_J @ c_h

# 6. BEM inductance matrix (LaplaceSL with .Trace()!)
j_trial, j_test = fes_J.TnT()
with TaskManager():
    V_op = LaplaceSL(
        j_trial.Trace() * ds(label)
    ) * j_test.Trace() * ds(label)
V_A = extract_dense(V_op.mat, n_J)
V_A_proj = np.real(c_h @ V_A @ c_h)

# 7. Inductance via L/R ratio (normalization-independent)
R_sheet   = 1.0 / (sigma * thickness)
perimeter = 2 * np.pi * R_center  # or 4*side for rectangle
R_loop    = R_sheet * perimeter / trace_width
LR_ratio  = MU_0 * V_A_proj / (R_sheet * energy)
L = LR_ratio * R_loop
```

### Key Formulas

- **L/R ratio**: `L = mu_0 * V_A_proj / (R_sheet * energy) * R_loop`
  - V_A_proj = c_h^T @ V_A @ c_h (BEM projection)
  - energy = c_h^T @ M_J @ c_h (mass matrix norm)
  - R_sheet = 1/(sigma*t), R_loop = R_sheet * perimeter / trace_width
  - This ratio is independent of c_h normalization

- **Hodge decomposition**: For genus-g surface, null([D; C^T@M_J]) has dim = g.
  Ring/frame (genus-1) has exactly 1 harmonic mode.
  The harmonic mode is both div-free (D@c_h=0) and M_J-orthogonal to gradients (C^T@M_J@c_h=0).

### Multi-Turn Solenoid

For N-turn coils, combine BEM self-inductance with Neumann mutual inductance:
```python
from scipy.special import ellipk, ellipe

def neumann_mutual(R1, R2, d):
    \"\"\"Mutual inductance of coaxial circular loops (Neumann formula).\"\"\"
    k2 = 4*R1*R2 / ((R1+R2)**2 + d**2)
    K, E = ellipk(k2), ellipe(k2)
    return MU_0 * np.sqrt(R1*R2) * ((2/np.sqrt(k2) - np.sqrt(k2))*K - 2*E/np.sqrt(k2))

L_total = N * L_self + sum(2*M[i][j] for i<j)  # Neumann sum
```

### Validated Results

- Circular ring (R=10mm, w=1mm): BEM=48.64 nH, Analytical=48.78 nH (0.3% error)
- Rectangular frame (10mm, w=1mm): BEM~24 nH, FastHenry=21.6 nH, Grover~24 nH
- Computation time: ~420 ms per mesh (no bonus_intorder needed)

See: `validation_test/peec_integration/ngsbem_peec_demo/ngbem/1_turn_coil.py`
     `validation_test/peec_integration/ngsbem_peec_demo/compute_L_final.py`

### With Conductor Shield (SIBC, Standalone)

For a coil above a conducting shield plate, use ShieldBEMSIBC:

```python
from radia.ngsbem_eddy import ShieldBEMSIBC
from 1_turn_coil import (compute_loop_inductance, create_circular_ring_mesh,
                          discretize_ring_coil, create_shield_plate_mesh)

# 1. BEM air-core inductance
mesh = create_circular_ring_mesh(R, trace_width, maxh)
L_air, info = compute_loop_inductance(mesh, R, trace_width, sigma, thickness)

# 2. Coil filament model + shield plate mesh
topo_dict = discretize_ring_coil(R, N_seg=60, z=0.0)
plate_mesh = create_shield_plate_mesh(30e-3, 2e-3, gap=0.5e-3, maxh=3e-3)

# 3. SIBC solver (also supports mu_r for magnetic conductors)
shield = ShieldBEMSIBC(plate_mesh, sigma=3.7e7, mu_r=1.0)
shield.assemble(intorder=4)

# 4. Frequency sweep: Delta_Z -> L_eff, R_eff
Delta_Z = shield.compute_impedance_matrix(freq, topo_dict)
Delta_Z_total = np.sum(Delta_Z)  # single-turn: all segments in series
L_eff = L_air + np.imag(Delta_Z_total) / omega   # L decreases (Lenz)
R_eff = R_dc + np.real(Delta_Z_total)             # R increases (eddy loss)
```

**Physics:**
- SIBC surface impedance: Zs = (1+j)*sqrt(omega*mu_r*mu_0/(2*sigma))
- Eddy currents in shield oppose incident flux (Lenz's law): L_eff < L_air
- Ohmic loss in shield adds resistance: R_eff > R_dc
- Effect depends on skin depth vs plate thickness

**Limitations:**
- SIBC captures eddy current effects only (frequency-dependent shielding + loss)
- Does NOT capture DC magnetic flux enhancement (mu_r > 1 magnetostatic effect)
- For penetrable magnetic bodies, PMCHWT formulation would be needed
  (not yet available in ngsbem; see arxiv:2408.01321 for low-freq stabilized PMCHWT)

## Shield Coupling (BEM+SIBC)

Shield eddy currents are computed using the EFIE with surface impedance:

    (Zs * M_LL + jw * mu_0 * V_LL) * I_loop = -jw * b_loop

**CRITICAL SIGN CONVENTIONS (verified 2026-02-15):**

1. **V_LL sign**: MUST be **positive** (+jw*mu_0*V_LL).
   A minus sign violates Lenz's law (A_scat would reinforce A_inc).
   Diagnostic: check dot(A_scat, wire_dir) < 0.

2. **Loop basis**: Must use FACE-EDGE divergence matrix (not edge-vertex incidence).
   The LoopBasisBuilder computes null space of div operator
   to ensure div(J) = 0 in RT0 sense.

3. **Faraday sign (Delta_Z)**:
   Delta_Z[i][j] = jw * dot(A_scat_i, dir_i) * len_i
   Physical check: Im(Delta_Z) < 0 (reduced inductance)
                   Re(Delta_Z) > 0 (added resistance)

```python
from radia.ngsbem_eddy import ShieldBEMSIBC
from radia.ngsbem_interface import extract_edge_geometry

# Create shield mesh (aluminum plate)
shield_mesh = ...  # Netgen OCC Box mesh
shield = ShieldBEMSIBC(shield_mesh, sigma=3.7e7)
shield.assemble(intorder=4)

# Extract PEEC edge geometry for coupling
edge_geom = extract_edge_geometry(conductor_mesh)
topo_dict = {
    'segment_centers': edge_geom['centers'],
    'segment_directions': edge_geom['directions'],
    'segment_lengths': edge_geom['lengths'],
}

# Compute coupling at each frequency
Delta_Z = shield.compute_impedance_matrix(freq, topo_dict)
Z_shielded = Z_air + Delta_Z  # Add to PEEC branch impedance
```

## Ferrite Core Coupling (Image Method)

For planar ferrite cores, keep the analytical image-method approximation as an
application-level model instead of reviving the removed PEEC-moment coupling
module:

```python
mu_r = 1000
Delta_L = solver.L / (mu_r + 1)
L_total = solver.L + Delta_L * (mu_r - 1)  # = L_air * 2*mu_r/(mu_r+1)
```

**Why analytical?** BEM L_air uses Galerkin RT0 surface integrals.
Mixing with filament-based Delta_L (line integrals) produces basis mismatch.
The image method scales L_air directly, preserving matrix structure.

## Slab Impedance for Thin Conductors

Standard SIBC (Zs = (1+j)/(sigma*delta)) is valid only for delta << thickness.
For finite-thickness conductors:

    Zs_slab = Zs * coth(gamma * t)

where gamma = (1+j)/delta. Limits:
- delta << t: Zs_slab -> Zs (standard SIBC)
- delta >> t: Zs_slab -> 1/(sigma*t) (DC sheet resistance)

Available in Radia PEEC pipeline: lanczos_reduction.py, veriloga_generator.py

## Stabilized EFIE Connection (Weggler) - IMPLEMENTED

PEEC Loop-Star uses the SAME product space as Weggler's stabilized EFIE:

| Stabilized EFIE | PEEC Loop-Star | ngsbem space |
|:---|:---|:---|
| A_kappa (vector SL) | Inductance L | LaplaceSL on HDivSurface |
| V_kappa (scalar SL) | Potential P (V_0) | LaplaceSL on SurfaceL2 |
| Q_kappa (coupling) | BEM Q_0 | LaplaceSL(div, SurfaceL2) |

Key: stabilized EFIE multiplies V by kappa^2 (not divides V by kappa^2):
- Stabilized: cond = O(1) for all kappa (DC to RF)
- Classical:  cond = O(kappa^{-2}) (blows up at low frequency)

### Implementation in ngsbem_peec.py (2026-02-22)

- `_build_bem_coupling()`: Assembles Q_0 via LaplaceSL product space
- `_solve_stabilized()`: Full block system with k^2*V_0
- `_solve_full_loop_star()`: Reformulated Schur with precomputed P^{-1}@M_LS
- Validated: full vs stabilized agree within 0.53% L, 0.01% R (5/5 tests PASS)

### Research Findings (Q1-Q3)

- **Q1 (k definition)**: k = omega/c_0 (free-space), NOT k^2 = -jw*mu*sigma
- **Q2 (Laplace vs Helmholtz)**: At k=0.001, ||L_Lap - L_Helm|| = 1.75e-12
  Laplace kernel works perfectly for MQS (consistent with Radia policy)
- **Q3 (Circuit extraction)**: L = mu_0*A_0, P = V_0/eps_0, R external (SIBC)
  MQS regime: k^2*V_k -> 0, system reduces to pure L+R (loop-only)

Reference: https://github.com/Weggler/docu-ngsbem/blob/main/demos/Maxwell_DtN_Stabilized.ipynb

### Product-Space vs Loop-Star (Lucy / Weggler Comparison)

Lucy (Weggler)'s product-space formulation is mathematically equivalent to
Loop-Star but automates the decomposition at the function space level:

| Aspect | Loop-Star (discrete basis) | Product-space (variational) |
|:---|:---|:---|
| Stabilization | Explicit Λ, Σ basis construction + transform | HDivSurface × SurfaceL2 separates naturally |
| κ² scaling | Manual rescaling required | Built into the formulation |
| Implementation cost | Build D, C matrices + null_space | NGSolve FES declaration only |
| Condition number | O(1) | O(1) |
| Multiply connected | Global loop detection required | Same (harmonic mode extraction needed) |

**Key insight**: The product-space approach automates Loop-Star at the function
space level. The manual D, C, M_J construction and null_space computation in
1_turn_coil.py would be replaced by a single FES declaration in Lucy's framework.
However, for genus >= 1 geometries (rings, frames), both approaches require
explicit harmonic mode handling.

### PMCHWT for Penetrable Magnetic Bodies (NOT YET IMPLEMENTED)

SIBC is a surface impedance approximation that captures eddy current effects only.
To accurately model DC permeability enhancement of magnetic bodies via BEM,
the PMCHWT (Poggio-Miller-Chang-Harrington-Wu-Tsai) formulation is required:
- Coupled EFIE + MFIE with both electric (J) and magnetic (M) surface currents
- Low-frequency stabilized version uses quasi-Helmholtz projectors for O(1)
  condition number as kappa -> 0
- Reference: arxiv:2408.01321 (2024)
  "Low-Frequency Stabilizations of the PMCHWT Equation"
- Not yet available in ngsbem

### BEM EFIE-SIBC: Known SL Eigenvalue Limitation (2026-03-28)

The BEM EFIE `Z_s*J + jw*mu0*SL(J) = -jw*A_inc` uses `A_scat = mu0*SL(J_s)`, but
the Laplace SL eigenvalue for l=1 on a sphere is R/3 (not R).  This gives denominator
`(3*Z_s + jw*mu0*R)` instead of `(Z_s + jw*mu0*R)`, producing BEM/Analytical ratios
of 0.35-1.0 depending on Z_s. Only correct for PEC (Z_s -> 0).

Fix requires MFIE `n x curl(SL)` (not available in ngsolve.bem).
**Use FEM-SIBC (fem_esim_3d.py) for finite Z_s problems.**

See: `validation_test/induction_heating/cubit_panels_legacy/efie_sibc.py`

## Loop-Star Solver Modes (FIXED 2026-02-22)

NGBEMPEECSolver supports three modes via `solve_frequency(freqs, mode=...)`:

| Mode | System Solved | Use Case |
|------|--------------|----------|
| `'mqs'` | R + jwL (loop-only) | **Standard conductor impedance** |
| `'full'` | Reformulated Schur complement | Includes displacement current |
| `'stabilized'` | Weggler's block system | Same as full, better numerics |

### mode='mqs' (Recommended for Conductors)

Loop-only: ignores capacitive (star) unknowns. Correct for good conductors
where displacement current is negligible (sigma >> omega*epsilon).

### mode='full' (Reformulated Schur)

Uses precomputed P^{-1} @ M_LS to avoid the old `P/(jw)` blowup:
```
cap_correction = jw * M_LS^T @ P^{-1} @ M_LS
Schur = (R + jwL) - cap_correction
```
P^{-1} @ M_LS is computed once during assemble(); only the jw factor varies.

### mode='stabilized' (Weggler's BEM)

Full block system with BEM Q_0 coupling (LaplaceSL):
```
[R + jwL,        jw*mu_0*Q_0^T ] [I_loop]   [V]
[jw*mu_0*Q_0,    jw*mu_0*k^2*V_0] [rho   ] = [0]
```
Well-conditioned for all frequencies (O(1) condition number).

### Full/Stabilized vs MQS: Expected Differences

For **single conductors** (no return path), full/stabilized give L values
~3x larger than MQS. This is CORRECT physics:
- Single plate has very low self-capacitance
- Capacitive impedance 1/(jwC) >> jwL adds to total impedance
- Cap/L eigenvalue ratio ~ 470,000x (ill-scaled for single conductor)

For **multi-conductor systems** (with return paths), the difference is smaller
because mutual capacitance is larger and more physically relevant.

**Recommendation**: Use `mode='mqs'` for standard PEEC impedance extraction.
Use `mode='full'` or `mode='stabilized'` when parasitic capacitance matters.

### AVOID: Classical P/(jw) Formulation

The old formulation `Z_SS = P / (1j * omega)` causes O(kappa^{-2})
condition number blow-up at low frequency. This pattern is DEPRECATED.
Use the reformulated Schur complement or stabilized mode instead.

## FEM Cross-Verification Results

PEEC+BEM results verified by independent NGSolve FEM (A-formulation):

| Quantity | FEM | PEEC | Agreement |
|----------|-----|------|-----------|
| L_air | 96.58 nH | 100.69 nH | -4.1% |
| L_core (mu_r=1000) | 102.11 nH | 105.76 nH | -3.5% |
| Delta_L_core | +5.54 nH | +5.07 nH | 9.1% |
| L_shield (varies with freq) | 72.89-89.03 nH | 83.91-96.00 nH | 7-14% |
| Analytical L_air | 97.96 nH | -- | FEM: -1.4%, PEEC: +2.8% |

Shield effect is frequency-dependent (100 Hz - 100 kHz).
See docs/solver/NGBEM_INTEGRATION_DESIGN.md for full frequency sweep results.

## Physical Checklist

1. L_air > 0 (positive inductance)
2. Delta_L_core > 0 (ferrite concentrates flux -> L increases)
3. Delta_L_shield < 0 (Lenz's law: eddy currents oppose flux -> L decreases)
4. R_shield > R_air (Ohmic loss in shield adds resistance)
5. dot(A_scat, wire_dir) < 0 (scattered potential opposes incident)
6. Classical EFIE cond -> infinity as kappa -> 0; stabilized stays O(1)

## FastHenry PEEC vs ngsolve.bem PEEC: When to Use Which (2026-02-22)

Verified on 100mm x 10mm x 1mm Cu bus bar (same geometry, Dowell skin effect):

### Performance Comparison

| Metric | FastHenry PEEC | ngsolve.bem PEEC | Ratio |
|--------|---------------|-----------|-------|
| Total time | 2.8 ms | 276 ms | **FastHenry 98x faster** |
| DOFs | 1 segment | 53 loop | 53x fewer |
| Assembly | 2.1 ms | 262 ms | 126x |
| Freq sweep (60 pts) | 0.7 ms | 13.6 ms | 19x |

### Dowell F_R Agreement (Skin Effect)

Both approaches correctly apply Dowell formula (F_R ratio within ~5-10%):
- 100 kHz: ngsolve.bem=2.14, FastHenry=2.25, analytical=2.36
- 1 MHz: ngsolve.bem=7.35, FastHenry=7.42, analytical=7.57

### Key Differences

| Aspect | FastHenry PEEC | ngsolve.bem PEEC |
|--------|---------------|-----------|
| Model type | 1D filament (w x h cross-section) | 2D surface mesh (triangles, RT0) |
| L vs frequency | Constant (single filament, no redistribution) | Decreases (captures current redistribution) |
| DC resistance | Exact: R = rho*L/A | Higher (current distributed across 2D mesh) |
| Skin effect | Dowell/Bessel formula OR nwinc/nhinc | BEM naturally + Dowell optional |
| Proximity effect | Via nwinc/nhinc multi-filament | Automatic (mesh-based) |
| Complex geometry | Straight segments only | Arbitrary 2D conductor shapes |
| Circuit extraction | Direct L,R,C output | Requires post-processing |

### Selection Guide

| Use Case | Recommended | Why |
|----------|-------------|-----|
| Bus bars, straight wires | **FastHenry** | 98x faster, same accuracy |
| Circular/helical coils | **FastHenry** | Polygon approximation sufficient |
| PCB traces (straight) | **FastHenry** | Fast, multi-filament for proximity |
| L-shaped / T-shaped conductors | **ngsolve.bem** | 2D current redistribution matters |
| Spiral inductors (2D) | **ngsolve.bem** | Complex geometry |
| Quick estimation | **FastHenry** | Always start here |
| Precision analysis | **ngsolve.bem** | After FastHenry confirms baseline |

### Practical Workflow

1. **First pass**: FastHenry for fast L,R estimation (ms-scale)
2. **If needed**: ngsolve.bem for precision (current redistribution, complex geometry)
3. **Dowell/Bessel**: Both support Zs_func callback for skin effect
4. **Shield/core coupling**: Shield coupling uses ShieldBEMSIBC; magnetic cores use HDiv-VIM / reduced FEM
"""

RADIA_EFIE_PRECONDITIONER = """
# EFIE Calderon Preconditioner (Andriulli / Schoeberl)

Calderon-type preconditioner for the Electric Field Integral Equation (EFIE),
based on preconditioning the single layer (SL) operator by the rotated SL
within the dual sequence. Implemented by Joachim Schoeberl in ngsolve.bem.

Reference: S.B. Adrian, A. Dély, D. Consoli, A. Merlini, F.P. Andriulli,
"Electromagnetic Integral Equations: Insights in Conditioning and Preconditioning",
IEEE Open Journal of Antennas and Propagation, Vol.2, pp.1143-1174, 2021.
DOI: 10.1109/OJAP.2021.3121097

## Problem: EFIE Conditioning

The standard EFIE discretization on HDivSurface produces an ill-conditioned system:

    lhs = (j*kappa) * V1 + 1/(j*kappa) * V2

where V1 = HelmholtzSL on HDivSurface, V2 = HelmholtzSL on div(HDivSurface).
Without preconditioning, GMRES converges slowly or stagnates.

## Preconditioner Formula

The preconditioner applies the rotated SL operator:

    pre = M^{-1} @ (kappa * V_rot - 1/kappa * nabla_s @ M_H1^{-1} @ V_pot @ M_H1^{-1} @ nabla_s^T) @ M^{-1}

where:
- M^{-1}: HDivSurface mass matrix inverse (sparsecholesky)
- V_rot: HelmholtzSL applied to n x u (rotated trial/test functions)
- V_pot: HelmholtzSL on scalar H1 space (potential part)
- M_H1^{-1}: H1 mass matrix inverse (sparsecholesky)
- nabla_s: Surface curl operator (Cross(grad(upot).Trace(), n))

## Spaces Required

| Space | Order | Purpose |
|-------|-------|---------|
| HDivSurface | p (e.g. 3) | EFIE unknowns (surface current) |
| H1 | p+1 (e.g. 4) | Potential part of preconditioner |

Both spaces must be complex=True for electromagnetic problems.

## Complete Working Example (from Schoeberl's EMpre.ipynb)

```python
from ngsolve import *
from netgen.occ import *
from ngsolve.bem import *
from time import time

# --- Mesh ---
face = Glue(Sphere((0,0,0),1).faces)
mesh = Mesh(OCCGeometry(face).GenerateMesh(maxh=0.2)).Curve(4)

# --- EFIE setup ---
kappa = 3*pi
d = CF((0,-kappa,0))
E_inc = exp(1j*d*CF((x,y,z))) * CF((1,0,0))

fes = HDivSurface(mesh, order=3, complex=True)
u, v = fes.TnT()

fespot = H1(mesh, order=4, complex=True)
upot, vpot = fespot.TnT()

# --- EFIE operator assembly ---
with TaskManager():
    V1 = HelmholtzSL(u.Trace()*ds(bonus_intorder=4), kappa) \\
         * v.Trace()*ds(bonus_intorder=4)
    V2 = HelmholtzSL(div(u.Trace())*ds(bonus_intorder=4), kappa) \\
         * div(v.Trace())*ds(bonus_intorder=4)

lhs = (1j*kappa) * V1.mat + 1/(1j*kappa) * V2.mat
rhs = LinearForm(E_inc*v.Trace()*ds(bonus_intorder=3)).Assemble()

# --- Preconditioner assembly ---
invMHd = BilinearForm(u.Trace()*v.Trace()*ds).Assemble() \\
         .mat.Inverse(inverse="sparsecholesky")
n = specialcf.normal(3)

with TaskManager():
    Vrot = HelmholtzSL(Cross(u.Trace(),n)*ds(bonus_intorder=2), kappa) \\
           * Cross(v.Trace(),n)*ds(bonus_intorder=2)
    Vpot = HelmholtzSL(upot*ds(bonus_intorder=2), kappa) \\
           * vpot*ds(bonus_intorder=2)

surfcurl = BilinearForm(Cross(grad(upot).Trace(), n) \\
           * v.Trace()*ds).Assemble().mat
invMH1 = BilinearForm(upot*vpot*ds).Assemble() \\
         .mat.Inverse(inverse="sparsecholesky")

pre = invMHd @ (kappa*Vrot.mat \\
      - 1/kappa*surfcurl@invMH1@Vpot.mat@invMH1@surfcurl.T) @ invMHd

# --- Solve with GMRES ---
gfj = GridFunction(fes)
with TaskManager():
    gfj.vec[:] = solvers.GMRes(A=lhs, b=rhs.vec, pre=pre,
                                maxsteps=500, tol=1e-8)
```

## Performance (Sphere, maxh=0.2, order=3, kappa=3*pi)

| Metric | Value |
|--------|-------|
| HDivSurface ndof | 10,108 |
| H1 ndof | 5,778 |
| EFIE assembly | 32.3 s |
| Preconditioner assembly | 15.5 s |
| GMRES iterations | **26** |
| GMRES time | 69.1 s |
| Final residual | 8.2e-9 |

## Key Points

1. **HDivSurface only**: Works with surface H(div) discretization (not HCurl)
2. **Small kappa regime**: Particularly effective for low-frequency / MQS
3. **.Trace() required**: All BEM operators on HDivSurface need .Trace()
4. **bonus_intorder**: Use 4 for EFIE operators, 2-3 for preconditioner
5. **Operator composition**: NGSolve's @ operator chains BEM + sparse operators
6. **Two BEM assemblies**: One for EFIE (V1, V2), one for preconditioner (Vrot, Vpot)

## Connection to PEEC / Low-Frequency Stabilization

This preconditioner addresses the same low-frequency ill-conditioning as
Weggler's stabilized EFIE (see ngsbem_peec topic). Both target the O(kappa^{-2})
condition number blow-up:

| Approach | Mechanism | Use Case |
|----------|-----------|----------|
| Weggler stabilized EFIE | Product-space reformulation | PEEC impedance extraction |
| Calderon preconditioner | Operator preconditioning via rotated SL | Full-wave EM scattering |

For pure MQS (kappa -> 0), both achieve O(1) condition number.
The Calderon preconditioner is more general (works for arbitrary kappa)
while Weggler's approach integrates directly with PEEC circuit extraction.

## Source

EMpre.ipynb by Joachim Schoeberl (2026-02). Shared for ngsolve.bem integration.
"""

RADIA_FEM_VERIFICATION = """
# NGSolve FEM Verification Reference

NGSolve FEM was used to verify PEEC+BEM results for a circular coil with
ferrite core and aluminum shield. Full results in `docs/solver/NGBEM_INTEGRATION_DESIGN.md`.

## Summary: ALL 6 CHECKS PASS

| Check | Result |
|-------|--------|
| L_air vs analytical | -1.4% (PASS, <5%) |
| L_air FEM vs PEEC | -4.1% (PASS, <5%) |
| Delta_L_core > 0 | +5.54 nH (PASS) |
| Delta_L_core FEM vs PEEC | 9.1% (PASS, <15%) |
| Shield decreases L | All frequencies (PASS) |
| Shield DeltaL trend | Matches (PASS) |

## Key FEM Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| FEM order | 2 | HCurl, nograds=True |
| air_r | 120 mm (static), 60 mm (eddy) | Sphere with Dirichlet BC |
| maxh | 8 mm (static), 4 mm (eddy) | Air mesh |
| core.faces.maxh | 2 mm | Local refinement only |
| Solver | PARDISO | NGSolve interface to Intel MKL |

## PARDISO (NGSolve Feature)

PARDISO is Intel MKL's multi-threaded direct sparse solver. It is accessed
through NGSolve's inverse interface, NOT implemented in Radia.

```python
# NGSolve usage:
gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec
```

Speedup vs UMFPACK: 3.9-7.7x (multi-threaded vs single-threaded).
"""

RADIA_SCALAR_POTENTIAL = """
# Phi-Reduced Scalar Potential (Radia + NGSolve)

## Method: Simkin-Trowbridge (1979)

The phi-reduced formulation splits the magnetic field:

    H = H_s - grad(phi)

where H_s is the source field (permanent magnets or coils) computed by Radia,
and phi is the correction potential solved by NGSolve FEM to account for iron.

This is ideal when:
- The source field (magnets, coils) is easy to compute analytically (Radia)
- Iron regions need FEM for their nonlinear/complex geometry response
- Open boundary is needed (Radia handles it naturally)

## ScalarPotentialSolver API

```python
from radia.scalar_potential_solver import ScalarPotentialSolver

solver = ScalarPotentialSolver(mesh, iron_domains='iron', mu_r=1000, order=2)

# Source from Radia object (voxel interpolation)
solver.set_source_from_radia(radia_obj, resolution=31)

# Or from callback
solver.set_source_from_callback(lambda x,y,z: (Hx, Hy, Hz), resolution=31)

# Or from CoefficientFunction directly
solver.set_source_cf(H_source_cf)

# Solve
solver.solve(method='auto')  # 'single' if mu_r<5000, 'two' otherwise

# Results
B_cf = solver.get_B()        # CoefficientFunction (Tesla)
H_cf = solver.get_H()        # CoefficientFunction (A/m)
B_hdiv = solver.project_to_hdiv()  # HDiv GridFunction (div(B)=0)
```

## Coil Modeling via Equivalent Magnetization

A solenoid with N turns, current I, height h produces the same external
field as a uniformly magnetized body with M_z = N*I/h (A/m).

For a hollow rectangular coil (bore = inner hole, outer = winding OD):
```python
M_z = N * I / h  # A/m
a_in, a_out = bore/2, outer/2

# 4 wall blocks forming the hollow rectangle
left  = rad.ObjRecMag([-(a_in+a_out)/2, 0, 0], [a_out-a_in, outer, h], [0,0,M_z])
right = rad.ObjRecMag([ (a_in+a_out)/2, 0, 0], [a_out-a_in, outer, h], [0,0,M_z])
front = rad.ObjRecMag([0,  (a_in+a_out)/2, 0], [bore, a_out-a_in, h], [0,0,M_z])
back  = rad.ObjRecMag([0, -(a_in+a_out)/2, 0], [bore, a_out-a_in, h], [0,0,M_z])
coil = rad.ObjCnt([left, right, front, back])
```

## Solver Method Selection

| mu_r Range | Method | Why |
|------------|--------|-----|
| < 5000 | `single` (reduced potential) | Simple, fast, adequate accuracy |
| >= 5000 | `two` (Simkin-Trowbridge) | Avoids cancellation in high-mu iron |

## Nonlinear Material (Newton + SymbolicEnergy)

For nonlinear B-H curves, use NGSolve's SymbolicEnergy with Newton iteration
and energy-based line search. This is MUCH better than Picard iteration
(fixed-point mu update): Newton gives quadratic convergence vs linear.

**Key difference from A-formulation**: In scalar potential, H is the primary
variable (H = H_s - grad(phi)). The energy density uses B(H), not H(B).

```python
from ngsolve import BSpline
import pandas as pd

# Load B-H data (column 0 = H in A/m, column 1 = B in Tesla)
df = pd.read_csv('BH.txt', sep='\\t')
H_data = list(df.iloc[:, 0])
B_data = list(df.iloc[:, 1])

# B as function of H (for scalar potential formulation)
bh_direct = BSpline(2, [0] + H_data, B_data)

# Energy density: w(H) = integral_0^H B(H') dH'
w_H = bh_direct.Integrate()

# Compare with A-formulation (where H(B) is used):
# bh_curve = BSpline(2, [0] + B_data, H_data)  # H(B)
# energy_dens = bh_curve.Integrate()             # w*(B) = integral H(B')dB'
```

### Formulation: Energy Minimization

The total magnetic energy functional for the Simkin formulation:

```
E(phi) = integral w(|H_s - grad(phi)|) dx  [iron, nonlinear]
       + integral (mu_0/2)|H_s - grad(phi)|^2 dx  [air, linear]
       + integral (mu_0/2)*kelvin_weight*|H_s - grad(phi)|^2 dx  [kelvin]
```

Note: H_s (source field from Radia coil) is embedded in the energy.
No separate LinearForm needed.

```python
from radia.scalar_potential_solver import ScalarPotentialSolver
from ngsolve import *

mu0 = 4e-7 * 3.14159265

solver = ScalarPotentialSolver(
    mesh, iron_domains='iron', order=2,
    kelvin_region='kelvin', kelvin_radius=AIR_R,
    kelvin_center=SPHERE_CENTER)
solver.set_source_from_radia(coil, resolution=51)

# Build SymbolicEnergy form
H_s = solver._H_source_cf  # VoxelCoefficient from Radia
fes = H1(mesh, order=2, dirichlet='outer')
phi, v = fes.TnT()

H_trial = H_s - grad(phi)
H2 = InnerProduct(H_trial, H_trial)

a = BilinearForm(fes, symmetric=True)
# Air: linear energy
a += SymbolicEnergy(mu0/2 * H2, definedon=mesh.Materials("air"))
# Iron: nonlinear energy from B-H curve
a += SymbolicEnergy(w_H(sqrt(H2 + 1e-12)),
                    definedon=mesh.Materials("iron"))
# Kelvin: linear with weight
if solver._kelvin_region:
    R = solver._kelvin_radius
    cx, cy, cz = solver._kelvin_center
    r_sq = (x-cx)**2 + (y-cy)**2 + (z-cz)**2 + 1e-30
    kelvin_w = R**2 / r_sq
    a += SymbolicEnergy(mu0/2 * kelvin_w * H2,
                        definedon=mesh.Materials("kelvin"))

# Regularization (optional, small)
a += SymbolicBFI(1e-8 * phi * v)

c = Preconditioner(a, type="bddc", inverse="sparsecholesky")
```

### Newton + Energy Line Search

```python
sol = GridFunction(fes)
sol.vec[:] = 0

au = sol.vec.CreateVector()
r = sol.vec.CreateVector()
w = sol.vec.CreateVector()
sol_new = sol.vec.CreateVector()

# No source LinearForm needed (H_s embedded in energy)
with TaskManager():
    for it in range(50):
        E0 = a.Energy(sol.vec)
        a.AssembleLinearization(sol.vec)
        a.Apply(sol.vec, au)
        r.data = -au  # residual = 0 - au (no source term)

        inv = CGSolver(mat=a.mat, pre=c.mat)
        w.data = inv * r

        err = InnerProduct(w, r)
        print(f"Newton {it}: err = {err:.2e}")
        if abs(err) < 1e-4:
            break

        # Energy line search
        sol_new.data = sol.vec + w
        E = a.Energy(sol_new)
        tau = 1
        while E > E0:
            tau *= 0.5
            sol_new.data = sol.vec + tau * w
            E = a.Energy(sol_new)

        sol.vec.data = sol_new

# Post-process
H_cf = H_s - grad(sol)
B_cf = mu0 * H_cf  # air approximation; use B(H) for iron
```

### Comparison: Newton vs Picard

| Feature | Newton + SymbolicEnergy | Picard (mu update) |
|---------|------------------------|--------------------|
| Convergence | Quadratic (~5-10 iter) | Linear (~20-50 iter) |
| Line search | Energy-based (guaranteed) | Under-relaxation (manual) |
| Jacobian | Exact (auto-differentiated) | Approximate (secant slope) |
| Remanence | Via energy functional | Needs polarization method |
| Implementation | SymbolicEnergy + BSpline | Manual mu_gf update loop |

**IMPORTANT**: For the scalar potential formulation, use `B(H)` BSpline
(not `H(B)`) because H is the primary variable. For the A-formulation,
use `H(B)` BSpline because B = curl(A) is the primary variable.

### Multi-Level Current Sweep

For sweeping multiple current levels (1000 AT to 20000 AT), reuse the
converged solution as initial guess for the next level:

```python
for AT in range(1000, 21000, 1000):
    coil = build_radia_coil(AT)
    solver.set_source_from_radia(coil, resolution=51)
    # ... rebuild a with new H_s ...
    # sol.vec keeps previous converged state -> fast convergence
    # Newton iteration (typically 3-5 iterations after first level)
```

### Reduced vs Total Potential Regions

Classical Simkin-Trowbridge uses two scalar potentials:
- **Reduced region** (contains coils): H = H_s - grad(phi_r)
- **Total region** (source-free, iron + air): H = -grad(phi_t)
  Interface matching conditions connect the two regions.

ScalarPotentialSolver uses **reduced potential everywhere**:

    H = H_s - grad(phi)   in ALL regions (iron, air, kelvin)

This works because H_s (Radia coil field) is accurately computed
via VoxelCoefficient grid interpolation. Air regions CAN be included
in either reduced or total approach without affecting accuracy.

**IMPORTANT**: H_s is a **magnetic field vector** (A/m) from Biot-Savart
(free-space field, as if no iron exists). It is NOT a scalar potential.
Compute via `rad.Fld(coil, 'h', points)`.

### Source Field H_s Computation

H_s is the coil field in FREE SPACE (no iron). It is a 3-component
vector field, NOT a scalar potential phi.

ScalarPotentialSolver computes it automatically:

```python
solver.set_source_from_radia(coil, resolution=51)
# Internally: rad.Fld(coil, 'h', grid_points) -> VoxelCoefficient
```

Resolution 31-51 points per axis is typically sufficient.

Alternative source methods:
```python
# From a callback function
solver.set_source_from_callback(
    lambda x, y, z: (Hx, Hy, Hz), resolution=31)

# From an existing NGSolve CoefficientFunction
solver.set_source_cf(H_source_coefficient_function)
```

### Nonlinear Solver API

```python
solver = ScalarPotentialSolver(
    mesh, iron_domains='iron', order=2,
    kelvin_region='kelvin', kelvin_radius=0.3,
    kelvin_center=[0, 0.07, 0])
solver.set_source_from_radia(coil, resolution=51)

# Newton (RECOMMENDED) -- quadratic convergence, ~30 iterations
phi = solver.solve_nonlinear_newton(
    bh_data=bh_data,   # [[H(A/m), B(T)], ...]
    tol=1e-4,
    maxiter=80,
    dirichlet='outer')

# Picard (fixed-point) -- linear convergence, needs under-relaxation
phi = solver.solve_nonlinear(
    bh_data=bh_data,
    tol=1e-4,
    maxiter=50,
    relax=0.3,           # 0=full step, 0.3=30% damping
    dirichlet='outer')

# Post-process
B_cf = solver.get_B()   # CoefficientFunction
H_cf = solver.get_H()
```

Newton is strongly recommended: quadratic convergence, energy-based
line search, typically ~30 iterations. Picard converges linearly and
may need significant under-relaxation (relax=0.3-0.5) for stability.

### Verified: C-type Electromagnet (NI=20000 AT)

| Pole Geometry | FEM Bz (mT) | vs BEM ~945 mT |
|---------------|-------------|----------------|
| Rectangular   | -834        | 11.8% error    |
| Y-Z taper     | -910        | 3.8% error     |
| 3D frustum    | -956        | 1.2% error     |

**Key lesson**: Nonlinear results are very sensitive to pole tip
geometry. Linear results (mu_r=1000) are robust (~1% error even
with rectangular poles). Always verify geometry against reference.

### Pole Tip Geometry: Frustum via Cubit or Netgen OCC

For tapered pole tips (truncated pyramid), use Cubit journal files
or NGSolve OCC WorkPlane with Loft operations. Do NOT use GMSH for
mesh generation (GMSH is for visualization/post-processing only).

## Examples

- `docs/ngsolve_integration/integration_basics.ipynb` - SHOWCASE: Radia+NGSolve integration (field eval, scalar/vector potential, coordinate transform)

## Reference

J. Simkin and C. W. Trowbridge, "On the use of the total scalar potential
in the numerical solution of field problems in electromagnetics,"
Int. J. Numer. Methods Eng., vol. 14, pp. 423-440, 1979.
"""

RADIA_VECTOR_POTENTIAL = """
# VectorPotentialSolver: Reduced Vector Potential (A_r) Formulation

## Formulation

B = B_s + curl(A_r)

- B_s: Biot-Savart source field from Radia coil (`rad.Fld(coil, 'b', points)`)
- A_r: Reduced vector potential (material response only)
- HCurl FE space with `nograds=True` gauge

### Linear Weak Form

    int(nu * curl(A_r) . curl(v)) dx = int((nu_0 - nu) * B_s . curl(v)) dx

- nu = 1/(mu_0 * mu_r) (reluctivity), nu_0 = 1/mu_0
- RHS non-zero only in iron domains (nu != nu_0)

### Nonlinear Newton (SymbolicEnergy + Coenergy)

Minimize reduced energy:

    E(A_r) = int w*(|B_s + curl(A_r)|) dOmega - nu_0 * int B_s . curl(A_r) dOmega

- w*(B) = integral_0^B H(B') dB' (magnetic coenergy)
- Air: w*(B) = B^2/(2*mu_0)
- Iron: w*(B) from inverted B-H curve via BSpline.Integrate()
- Source coupling term -nu_0 * B_s . curl(A_r) is CRITICAL (ensures curl(H_r)=0)

Without the source coupling term, the energy stationarity gives curl(H)=0 (wrong).
With it, stationarity gives curl(H - H_s) = 0, i.e. curl(H_r) = 0 (correct).

## Comparison with ScalarPotentialSolver (Simkin)

| Item | phi_r (Simkin) | A_r (This solver) |
|------|----------------|-------------------|
| Unknown | phi (H1, scalar) | A_r (HCurl, vector) |
| Source | H_s (rad.Fld 'h') | B_s (rad.Fld 'b') |
| Operator | grad(phi), grad(v) | curl(A), curl(v) |
| Coefficient | mu (permeability) | nu = 1/mu (reluctivity) |
| Energy | w(H) = int B(H')dH' | w*(B) = int H(B')dB' |
| Gauge | Not needed | nograds=True + eps regularization |
| DOFs | ~N_vertices | ~N_edges (3-5x more) |

## API

```python
from radia.vector_potential_solver import VectorPotentialSolver

solver = VectorPotentialSolver(mesh, iron_domains='iron', mu_r=1000.0, order=2)
solver.set_source_from_radia(coil, resolution=41)

# Linear solve
solver.solve_linear(dirichlet='outer')

# Nonlinear Newton (RECOMMENDED)
solver.solve_nonlinear_newton(bh_data=bh_data, tol=1e-4, maxiter=50, dirichlet='outer')

# Nonlinear Picard
solver.solve_nonlinear(bh_data=bh_data, tol=1e-4, maxiter=50, relax=0.3, dirichlet='outer')

# Hysteresis (Hantila polarization method, energy play model)
solver.solve_hysteresis(mat_factory, tol=1e-3, maxiter=60, alpha=500.0,
                        dirichlet='outer', relax=0.5)

# Results
B = solver.get_B()    # CoefficientFunction: B_s + curl(A_r)
H = solver.get_H()    # CoefficientFunction: nu * B
A = solver.get_A()    # GridFunction (HCurl)
B_gf = solver.project_to_hdiv()  # HDiv GridFunction
```

## Hysteresis: Hantila Polarization Method for A_r

Uses the Hantila (1975) polarization method adapted for vector potential:

    H = (B/mu_0 - R_prev) / (1 + alpha)

This derives H from B **without needing the Forward(B->H) operator**.
Only the Inverse operator `rad.MatMvsH(handle, 'm', H_vec)` is used.

- LHS (constant): nu_alpha * curl(A).curl(v) [iron] + nu_0 * curl(A).curl(v) [air]
  where nu_alpha = 1/(mu_0*(1+alpha))
- RHS (updated): (nu_0 - nu_alpha) * B_s.curl(v) + R/(1+alpha) . curl(v) [iron]
- R = M - alpha*H (polarization residual), updated each iteration
- Under-relaxation (`relax` parameter) damps oscillation near play model
  pinning thresholds where local dM/dH may exceed alpha
- Convergence: max|dB|/B_sat < tol

**A_r is the RECOMMENDED solver for B-input hysteresis models**:
- B = B_s + curl(A_r) is the primary variable -> B directly available for Play operator input
- ScalarPotentialSolver computes H = H_s - grad(phi) first, then derives B via constitutive law;
  H_s and grad(phi) cancellation degrades B accuracy (especially at high mu_r)
- Test results confirm: A_r oscillation amplitude 0.06 vs Simkin 0.21

The Hantila method is also available for ScalarPotentialSolver (less accurate for B-input models):
```python
from radia.scalar_potential_solver import ScalarPotentialSolver
solver_S = ScalarPotentialSolver(mesh, iron_domains='iron', order=2,
                                  kelvin_region='kelvin', kelvin_radius=R,
                                  kelvin_center=center)
solver_S.set_source_from_radia(coil, resolution=41)
solver_S.solve_hysteresis(mat_factory, tol=1e-3, maxiter=60,
                          alpha=500.0, dirichlet='outer', relax=0.5)
```

## Verified Results (C-type Electromagnet, Frustum Poles)

### Linear (mu_r=1000, NI=2000)
- A_r: Bz = -231.71 mT (0.1% vs BEM ~-232 mT)
- Simkin: Bz = -230.31 mT (0.7% vs BEM)
- A_r vs Simkin: 0.6%

### Nonlinear (BH curve, NI=20000)
- A_r Newton: Bz = -999.15 mT (5.7% vs BEM ~-945 mT)
- Simkin Newton: Bz = -922.40 mT (2.4% vs BEM)

### Hysteresis (Energy Play, NI=20000)
- A_r Hantila: Bz = -741 mT
- Simkin Hantila: Bz = -729 mT
- A_r vs Simkin: 1.6%

## FEM -> Radia Analytical Field Pipeline (Beam Tracking)

After solve_hysteresis(), the solved per-element magnetization M can be
exported to Radia objects for exact analytical field evaluation in the gap
(no mesh needed, no cancellation error):

```python
# 1. Solve hysteresis
solver = VectorPotentialSolver(mesh, iron_domains='iron', mu_r=1000.0, order=2)
solver.set_source_from_radia(coil)
solver.solve_hysteresis(mat_factory, alpha=500.0, relax=0.5)

# 2. Export to Radia (iron M + coil combined)
combined = solver.to_radia(coil=coil)

# 3. Exact analytical B at any point (no gap mesh needed)
B = rad.Fld(combined, 'b', [0, 0, 0])           # single point
B_batch = rad.Fld(combined, 'b', points_Nx3)     # batch (N,3)

# 4. VoxelCoefficientFunction for fast trajectory integration
from radia.radia_ngsolve import create_voxel_cf
B_voxel = create_voxel_cf(combined, 'b', mesh=mesh, resolution=61)
```

**Why this is better than FEM direct evaluation**:
- Radia surface charge formulas are **exact** for constant M per element
- No H_s - grad(phi) cancellation error
- No mesh resolution limitation in gap
- Arbitrary evaluation points without element search

**API**:
- `solver.get_M_per_element()` -> dict {el_nr: [Mx, My, Mz]} in A/m
- `solver.to_radia(coil=None)` -> Radia container handle
- Available on both VectorPotentialSolver and ScalarPotentialSolver

### Verified Results (C-type, NI=20000, Energy Play)
- FEM direct Bz at gap center: -741 mT
- Radia analytical Bz: -771 mT
- Agreement: 4.1%
- 2197 iron elements converted in 0.2s

## Test Scripts

- `validation_test/ngsolve_integration/mesh_magnetization_import/verified_ngsolve_to_radia.py` - FEM->Radia analytical pipeline
- `docs/kelvin/KELVIN_TRANSFORMATION.md` and
  `docs/kelvin/kelvin_exterior_source_and_aphi.ipynb` - maintained Kelvin
  theory and executed result view; implementation is `radia.kelvin_source`
"""

RADIA_PLAY_MODELS = """
# Play Models: Typical Radia Usage Patterns

For fixed magnetization, Radia evaluates fields from geometric objects directly.
For soft iron, use the mesh-backed HDiv-VIM route.

| Object | Use |
|--------|-----|
| ObjRecMag | optimized rectangular fixed magnet/current source |
| ObjHexahedron / ObjTetrahedron / ObjWedge / ObjPyramid | fixed magnetization on imported geometry |
| radia.vim.MeshSoftIron | mesh-backed soft iron for HDiv-VIM |

Do not build new soft-iron workflows by applying magnetic material directly to
raw imported polyhedra.  Keep the NGSolve mesh and use the HDiv-VIM path so the
same mesh/material labels can be reused by reduced FEM.

---

## Pattern A: Permanent Magnet Only (No Solve)

Fixed magnetization -- no iterative solve needed.  Fastest pattern.

```python
import radia as rad
rad.UtiDelAll()

# Rectangular PM: 20x20x10 mm, Br=1.2T in Z
# M = Br / mu_0 = 1.2 / (4*pi*1e-7) = 954930 A/m
mag = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.01], [0, 0, 954930])

# Or use ObjHexahedron with 8 vertices for arbitrary shapes
s = 0.01  # half-size in meters
verts = [[-s,-s,-s],[s,-s,-s],[s,s,-s],[-s,s,-s],
         [-s,-s,s],[s,-s,s],[s,s,s],[-s,s,s]]
mag2 = rad.ObjHexahedron(verts, [0, 0, 954930])

B = rad.Fld(mag, 'b', [0, 0, 0.03])  # No Solve() needed!
rad.UtiDelAll()
```

**Example**: the cubic_polyhedron_magnet section of `docs/simple_problems/simple_problems.ipynb`

---

## Pattern B: PM + Soft Iron (Linear)

Keep the permanent magnet and soft iron in separate conforming spaces. The PM
is an immutable magnetization source; the mesh-backed iron owns the unknowns.

```python
import numpy as np
import ngsolve as ng
from radia import vim

pm = vim.MagnetizationSource(pm_mesh, np.array([0.0, 0.0, 954930.0]))
with ng.TaskManager():
    result = vim.Solve(
        iron_mesh, mu_r=1000.0, magnetization_sources=pm)
```

Use separate PM and iron meshes when their normal magnetization is
discontinuous at a touching interface.

---

## Pattern C: Iron in Background Field

Soft iron in an externally applied field.  Classic magnetization problem.

```python
import ngsolve as ng
import radia as rad
from radia.vim import soft_iron_box

rad.UtiDelAll()
iron = soft_iron_box(
    center=(0.0, 0.0, 0.0), size=(0.02, 0.02, 0.02),
    mu_r=1000.0, nsub=3,
)
source = rad.ObjBckg(lambda p: [0.0, 0.0, 0.25])
system = rad.ObjCnt([iron, source])
with ng.TaskManager():
    result = rad.Solve(system)
B = rad.Fld(system, "b", [0.0, 0.0, 0.02])
rad.UtiDelAll()
```

The HDiv order and mesh, not an old relaxation method number, determine the
finite-element unknown space.

---

## Pattern D: Nonlinear Iron (B-H Curve)

Saturable iron with tabulated or functional B-H data.

```python
import ngsolve as ng
from radia import vim

# Tabulated B-H curve: [[H (A/m), B (T)], ...]
BH = [[0,0], [100,0.1], [500,0.8], [1000,1.2],
      [5000,1.6], [20000,1.9], [50000,2.0]]
with ng.TaskManager():
    result = vim.Solve(
        mesh, bh_table=BH, H_ext=applied_field,
        nonlinear_solver="energy-newton", nl_tol=1e-6, nl_maxit=300,
    )
```

**Example**: `docs/background_fields/sphere_in_quadrupole.py`

---

## Pattern E: Mesh-Backed Soft Iron

Use the NGSolve mesh directly for complex soft-iron geometries.

```python
import radia as rad
import radia.vim as vim

rad.UtiDelAll()

iron = vim.MeshSoftIron(mesh, mu_r=1000)
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
grp = rad.ObjCnt([iron, bkg])
with ng.TaskManager():
    result = rad.Solve(grp, demag_backend="hdiv")
rad.UtiDelAll()
```

**Rule**: imported soft iron stays mesh-backed.  Use raw imported polyhedra for
fixed magnetization export and field evaluation, not for new soft-iron solves.
**Example**: `validation_test/ngsolve_integration/mesh_magnetization_import/`
**Showcase notebook** (B=curl(A) check, RadiaField->HDiv projection, batch
evaluation, executed + rendered): `docs/ngsolve_integration/integration_basics.ipynb`

---

## Pattern F: Repeated or Large HDiv-VIM Solve

Use the persistent solver for operator reuse and configure the charge-Gram and
linear solve through named arguments.

```python
solver = vim.HDivSolver(mesh, order=2, gram_eps=1e-10, leaf=32, eta=2.0)
with ng.TaskManager():
    result = solver.Solve(
        mu_r=1000.0, H_ext=applied_field,
        linear_solver="auto", preconditioner="auto")
```

Read `result["linear_solver"]`, `result["preconditioner"]`, and the current
validation JSON before reporting which backend or threshold is fastest.
"""


RADIA_HYSTERESIS = """
# Magnetic Hysteresis: Mesh-Backed HDiv-VIM Contract

Radia's production field-history route is the B-input Play model on a
mesh-backed NGSolve HDiv space. Material identification and material-only
evaluation are separate from the spatial field solve.

## Production Entry Points

| Need | Entry point |
|------|-------------|
| One history on one mesh | `radia.vim.SolveHysteresis(mesh, h_steps, ...)` |
| Continuation and geometry reuse | `radia.vim.HDivSolver(mesh).SolveHysteresis(...)` |
| Material object | `radia.vim.PlayHysteresisMaterial(K, eta, f_k_tables)` |
| JMAG / MATLAB identification import | `radia.hysteresis_io.hys_to_play_radia` / `mat_to_play_radia` |

```python
import ngsolve as ng
from radia import vim
from radia.hysteresis_io import hys_to_play_radia

K, eta, tables = hys_to_play_radia("material.hys", K=20)
material = vim.PlayHysteresisMaterial(K, eta, tables)
h_steps = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 2.0e5],
    [0.0, 0.0, 0.0],
]

with ng.TaskManager():
    result = vim.SolveHysteresis(
        mesh,
        h_steps,
        material=material,
        tol=1e-8,
        maxit=4000,
        nl_tol=1e-3,
        nl_maxit=200,
    )
```

For repeated histories on the same geometry:

```python
solver = vim.HDivSolver(mesh, order=2)
with ng.TaskManager():
    first = solver.SolveHysteresis(h_steps_a, material=material)
    continued = solver.SolveHysteresis(h_steps_b, material=material)
```

## Contracts

- Supply exactly one of `play=(K, eta, tables)` or
  `material=PlayHysteresisMaterial(...)`.
- `h_steps` is the applied H history in A/m.
- The caller owns the surrounding `ngsolve.TaskManager`.
- Use `HDivSolver` when continuation and charge-Gram reuse matter.
- Image symmetry is not supported by the current HDiv history solver and must
  fail loudly.
- Raw `ObjRecMag` / `ObjHexahedron` plus `MatApl` is not a production
  soft-iron workflow. Mesh-less soft-iron `rad.Solve` is retired.
- The low-level `rad.MatPlayHysteresis`, `rad.MatMvsH`, and
  `MatHysSaveState/RestoreState/CommitState` functions remain useful for
  constitutive material checks. They do not replace the mesh-backed field
  solver.
- The energy-stop model is a separate intentional model with stricter shape
  constraints; do not silently substitute it for B-input Play.

## Evidence

Do not copy benchmark times, iteration counts, or accuracy percentages into MCP
knowledge. Read the checked JSON produced by the relevant validation campaign:

- `validation_test/hysteresis/test_binput_hdiv.py`
- `validation_test/hysteresis/test_real_material_hysteresis.py`
- `validation_test/hysteresis/test_loop_pollution_binput.py`
- `validation_test/hysteresis/test_energy_stop_irreversible_pm.py`
- `validation_test/hysteresis/bench_hysteresis_step.py`

Fast regression coverage belongs in `tests/test_hdiv_vim_hysteresis_rt2.py`
and the focused hysteresis API tests. Numerical and publication claims belong
in validation JSON, not in this executable manual.
"""

RADIA_ESIM = """
# ESIM (Effective Surface Impedance Method)

ESIM computes frequency- and H-dependent surface impedance Zs(H, omega)
for nonlinear magnetic conductors by solving a 1D cell problem.

## Physics

## Impedance Model Selection Rules

| Model | Type | Material input | Use case |
|-------|------|---------------|----------|
| **SIBC** | Linear | mu_r, sigma | Known constant permeability |
| **Dowell** | Linear | mu_r, sigma | Slab geometry, analytical |
| **ESIM** | Nonlinear | BH curve, sigma | Steel, saturation, Karl iteration |

- SIBC/Dowell: specify mu_r (scalar). BH curve NOT used.
- ESIM: specify BH curve (H[A/m], B[T] table). mu_r NOT used (derived from BH).
- Do NOT pass both mu_r and BH curve to the same model.

Standard SIBC uses constant Zs = (1+j)/(sigma*delta) (linear mu).
ESIM solves the 1D penetration problem with nonlinear mu(H):

```
rho * d^2H/dz^2 + j*omega*mu(|H|)*H = 0
H(0) = H0 (surface), H(inf) = 0 (bulk)
```

This yields Zs(H0) = E(0)/H0 which depends on the surface field amplitude.

## Solvers

| Solver | Domain | Use Case |
|--------|--------|----------|
| `ESIMCellProblemSolver` | Semi-infinite | Thick conductors (delta << thickness) |
| `ESIMFiniteSlabSolver` | Finite slab [0, a] | Transition region (delta ~ thickness) |

## Quick Start

```python
from radia import generate_esi_table_from_bh_curve, ESITable

# BH curve: [[H(A/m), B(T)], ...]
bh_curve = [[0, 0], [100, 0.5], [500, 1.2], [2000, 1.5], [50000, 2.0]]

# Generate ESI table for steel workpiece
esi_table = generate_esi_table_from_bh_curve(
    bh_curve, sigma=5e6, frequency=10000, n_points=50
)

# Query at specific surface H
Zs = esi_table.get_impedance(H0=1000)           # complex [Ohm]
P_loss, Q_react = esi_table.get_power_loss(1000) # W/m^2, var/m^2

# Save/load
esi_table.save('steel_10kHz.esi')
esi_table = ESITable.load('steel_10kHz.esi')
```

## Finite Slab (Thin Conductors)

For conductors where delta ~ thickness, use `ESIMFiniteSlabSolver`:

```python
from radia.esim_cell_problem import ESIMFiniteSlabSolver

solver = ESIMFiniteSlabSolver(
    half_thickness=0.005,   # 5mm half-thickness [m]
    bh_curve=bh_curve,
    sigma=5e6,              # [S/m]
    frequency=10000         # [Hz]
)
result = solver.solve(H0=1000)
# result['Z']: surface impedance with coth(gamma*t) correction
# result['R_ac_over_R_dc']: AC/DC resistance ratio
# result['H_profile']: H(z) depth profile
```

## Complex Permeability

ESIM supports complex mu = mu' - j*mu" for magnetic losses
(hysteresis loss, grain-boundary eddy currents):

```python
# Constant complex permeability
solver = ESIMCellProblemSolver(
    complex_mu=(500, 50),  # mu'_r=500, mu"_r=50
    sigma=5e6, frequency=10000
)

# H-dependent complex permeability
complex_mu_data = [
    [100, 800, 30],   # [H(A/m), mu'_r, mu"_r]
    [1000, 500, 50],
    [5000, 200, 20],
]
solver = ESIMCellProblemSolver(
    complex_mu=complex_mu_data,
    sigma=5e6, frequency=10000
)
```

## Usage Contexts

### 1. PEEC Conductor Impedance (Radia)

ESIMCoupledSolver uses ESI tables for PEEC conductor impedance:
```python
from radia.esim_coupled_solver import ESIMCoupledSolver
solver = ESIMCoupledSolver(topology, bh_curve, sigma, freq)
Z_port = solver.compute_port_impedance()
```

### 2. NGSolve FEM Robin BC

ESIM's Zs(H) can be used as a nonlinear Robin BC in NGSolve A-Phi
formulation (see ngsolve_usage topic "darwin" for details):
```python
a += (1/Zs_cf) * A.Trace() * N.Trace() * ds("workpiece")  # Robin BC
# Picard iteration updates Zs from ESI table at each step
```

### 3. BEM + SIBC Shield (ShieldBEMSIBC)

For conducting shields with ngsolve.bem BEM:
```python
from radia.ngsbem_eddy import ShieldBEMSIBC
shield = ShieldBEMSIBC(mesh, sigma=3.7e7, mu_r=1.0)
```

## Linear vs ESIM Comparison

| Property | Linear SIBC | ESIM |
|----------|-------------|------|
| mu_r | Constant | H-dependent (B-H curve) |
| Zs | `(1+j)/(sigma*delta)` | Numerical 1D solve |
| Complex mu | No | Yes (mu' - j*mu") |
| Solver | One-shot | Picard iteration |
| Best for | Cu, Al | Steel, iron, ferrite |

## 4. BEM Coil + ESIM Workpiece (Induction Heating)

Coupled analysis: BEM computes coil current, ESIM computes workpiece heating.

Pipeline:
1. BEM (source/sink saddle point EFIE) -> coil surface current J
2. Biot-Savart from J -> H at workpiece surface panels
3. ESIM cell problem (or Dowell analytical) -> Z_s(H), P', Q' per panel
4. Integration over surface -> total R, P, Q, effective impedance

```python
from impedance_esim import run

# Torus coil (R=30mm) + steel cylinder workpiece at center
result = run(
    material='steel',       # Nonlinear BH curve
    frequency=50000,        # 50 kHz
    R=0.030,                # Coil major radius [m]
    a=0.003,                # Coil wire radius [m]
    wp_radius=0.010,        # Workpiece radius [m]
    wp_height=0.030,        # Workpiece height [m]
)

# Scale to actual current (results are per unit current I=1A)
I = 100.0  # Amperes
P_heating = result['P_total'] * I**2       # Total heating power [W]
R_wp = result['R_effective']               # Workpiece resistance [Ohm]
L_coil = result['L_coil']                  # Coil inductance [H]
```

Typical results (1-turn coil, steel, 50 kHz):
- P = 11.5 W at 100A (scales as N^2 * I^2 for N-turn coil)
- Skin depth = 0.040 mm (steel, high mu)
- Power factor cos(phi) ~ 0.08 (resonant compensation needed)

Model selection:
- **ESIM**: Nonlinear BH, Picard iteration, accurate for steel/iron/ferrite
- **Dowell**: Analytical Z_s = (rho/a)*gamma*a*tanh(gamma*a), linear materials only

### Cubit Panel Workpiece Block

The Cubit inductance panel auto-detects a `workpiece` block:

```
block N add volume <workpiece_volume_id>
block N name "workpiece"
```

When present, the panel shows ESIM/Dowell settings:
- Model: ESIM / Dowell
- Material: Steel / Copper / Aluminum (sigma auto-set)
- Frequency, sigma, half-thickness
- Results: L (coil), R (workpiece), P, Q, skin depth, |Z_total|

### Verification

ESIM verified against NGSolve H1 FEM (p=4) as independent method:

| Test | Max error |
|------|-----------|
| Linear Z_s vs analytical | 0.25% |
| Linear Z_s vs NGSolve FEM | 0.04% |
| Nonlinear Z_s vs NGSolve FEM + Picard | 0.78% |

Run: `python validation_test/induction_heating/cubit_panels_legacy/verify_esim.py`

## Key Files

| File | Class/Function | Purpose |
|------|---------------|---------|
| `esim_cell_problem.py` | `ESIMCellProblemSolver` | Semi-infinite 1D solver |
| `esim_cell_problem.py` | `ESIMFiniteSlabSolver` | Finite thickness solver |
| `esim_cell_problem.py` | `ESITable` | Zs(H) lookup with interpolation |
| `esim_cell_problem.py` | `generate_esi_table_from_bh_curve()` | Main entry point |
| `esim_coupled_solver.py` | `ESIMCoupledSolver` | PEEC + ESIM coupled solver |
| `esim_workpiece.py` | `ESIMWorkpiece` | 3D workpiece with ESI tables |
| `panels/calc_inductance.py` | `_compute_workpiece_impedance()` | Panel ESIM/Dowell backend |
| `examples/.../impedance_esim.py` | `run()` | Standalone BEM+ESIM analysis |
| `examples/.../verify_esim.py` | `test1..4()` | ESIM vs NGSolve FEM verification |
"""

RADIA_BUILD_AND_RELEASE = """
# Build, Release, and PyPI Publishing

## Build Pipeline

```
Build.ps1 (MSVC + MKL + NGSolve)
  |-> _radia_pybind.pyd (main C++ extension, required)
  |-> peec_matrices.pyd  (PEEC matrix assembly, optional)
  |-> cln_core.pyd       (Lanczos MOR, optional)
  +-> All copied to src/radia/
```

### Prerequisites

- **Visual Studio 2022** (MSVC compiler)
- **Intel MKL** from `python -m pip install mkl-devel`
- **NGSolve / Netgen** at the exact versions pinned in `pyproject.toml`
- **Python 3.12** with pybind11

### Build Commands

```powershell
# Standard build
pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1

# Clean rebuild
pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -Rebuild

# Build + run tests
pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -Test

# Verbose output
pwsh -NoProfile -ExecutionPolicy Bypass -File Build.ps1 -Verbose
```

### NGSolve for CI Runner

The self-hosted CI runner (NETWORK SERVICE account) cannot access the mapped lab share.
NGSolve must be copied to C:\\NGSolve locally:

```powershell
# Run as Administrator when NGSolve is rebuilt:
robocopy "<configured-ngsolve-source>" C:\\NGSolve /MIR
```

## CI/CD Pipeline (GitHub Actions)

```
git push (main or v* tag)
  |
  v
CI (build-test.yml)
  |-> Verify NGSolve at C:\\NGSolve
  |-> Build.ps1 -Verbose (MSVC + MKL)
  |-> Verify _radia_pybind.pyd exists
  |-> pytest -m basic (quick tests)
  |-> Build_Wheel.ps1 -DryRun (build wheel, verify, no upload)
  |-> Upload artifacts: radia-pyd, radia-wheel, test-results
  |-> Upload exact ref context: ref type/name + SHA + run ID + tag snapshot
  |
  v
Release (release.yml) -- triggered by CI success
  |
  +-> [main branch] upload-binaries
  |     Upload .pyd to GitHub Releases (tag: binaries)
  |
  +-> [exact tag-ref CI only] verify per-run ref context
        |-> [v* tag on the same SHA] publish-pypi
        Download wheel artifact
        Publish to PyPI via OIDC Trusted Publishers
        (pypa/gh-action-pypi-publish, no token needed)
```

### PyPI Publishing is AUTOMATIC

When you push a version tag (e.g., `v2.5.0`), the pipeline:
1. CI builds and tests
2. If CI passes, Release workflow publishes wheel to PyPI
3. Uses OIDC Trusted Publishers (no API token stored)

## Release Checklist

### 1. Prepare Release

```python
# 1. Bump version in TWO files (must match):
#    - pyproject.toml: version = "X.Y.Z"
#    - src/radia/__init__.py: __version__ = "X.Y.Z"

# 2. Update CHANGELOG.md

# 3. Build locally and verify
powershell -ExecutionPolicy Bypass -File Build.ps1 -Rebuild -Test
```

### 2. Verify Wheel

```python
import zipfile
whl = zipfile.ZipFile('dist/radia-X.Y.Z-cp312-cp312-win_amd64.whl')
for info in whl.infolist():
    if info.filename.endswith('.pyd'):
        print(f'{info.filename}: {info.file_size} bytes')
# Must contain radia/_radia_pybind.pyd (> 2 MB)
# Must NOT contain any .dll files (MKL policy)
```

### 3. Tag and Push

```bash
git add pyproject.toml src/radia/__init__.py CHANGELOG.md
git commit -m "Release vX.Y.Z: description"
git tag vX.Y.Z
git push origin main vX.Y.Z
```

### 4. Monitor

```bash
# gh-free (No GitHub CLI policy): use the REST helper, not `gh run`
python tools/check_ci.py --branch main        # Check CI status
python tools/check_ci.py --sha <sha> --watch  # Watch a commit to completion
pip install radia==X.Y.Z                       # Verify after publish
```

## Wheel Build Details (Build_Wheel.ps1)

1. Clean dist/ and build/ directories
2. `python -m build --wheel`
3. **Remove any bundled DLLs** (MKL policy enforcement)
4. Repack wheel with correct platform tag: `cp312-cp312-win_amd64`
5. `twine check` for metadata validation

**MKL DLL Policy**: Wheel MUST NOT bundle Intel MKL DLLs.
Users install MKL via pip dependency: `mkl>=2024.2.0`.
DLLs go to `{sys.prefix}/Library/bin/`, loaded by `__init__.py`.

## Policy Lint (policy-lint.yml)

Runs on every push. Checks:
1. No `FldUnits()` calls (removed API)
2. No binary files tracked in git
3. No Helmholtz kernel in C++ core (Laplace-only)
4. No `CblasColMajor` (row-major policy)
5. No generated files at repo root
6. No legacy import paths
7. Every example directory has README.md

## Package Structure

```
src/radia/
  __init__.py           # DLL path setup + re-export
  _radia_pybind.pyd     # Main C++ extension (required)
  peec_matrices.pyd     # PEEC matrix assembly (optional)
  cln_core.pyd          # CLN transient solver (optional)
  *.py                  # Python utility modules
  # NO .dll files (MKL loaded from pip install location)
```

## Troubleshooting

- **CI fails "NGSolve not found"**: Run `robocopy` to sync C:\\NGSolve
- **Wheel too large**: Check for accidentally bundled .dll files
- **Import fails on user machine**: Ensure `pip install mkl` was installed
- **PyPI publish fails**: Check OIDC Trusted Publishers config on PyPI
- **Binary upload fails**: Ensure `binaries` release exists on GitHub
"""


RADIA_MULTILEVEL_SIMULATOR = """
# Multi-Level Simulator

Radia provides 3 levels of EM analysis in one repository.
All share the same geometry and coordinate system for cross-validation.

## The Three Levels

| Level | Method | Speed | Accuracy | Use Case |
|-------|--------|-------|----------|----------|
| **1. PEEC** | Filament + Neumann | seconds | ~5% | Design exploration, parametric |
| **2. NGSBEM** | Surface BEM (Laplace SL) | minutes | ~1% | Detailed analysis, validation |
| **3. FEM** | Volume (NGSolve) | hours | reference | Final verification, nonlinear |

## Level 1: PEEC (fast screening)

```python
from fasthenry_parser import FastHenryParser
result = FastHenryParser().parse_string(inp).solve()  # L, R, Z(f)
```
- Sub-second for 100 frequency points
- Magnetic core via Delta_L from the Radia HDiv-VIM soft-iron route
- Limitation: filament approximation (no skin/proximity in conductor)

## Level 2: NGSBEM (detailed)

```python
from radia.ngsbem_peec import NGBEMPEECSolver
solver = NGBEMPEECSolver(mesh, order=0, sigma=5.8e7)
solver.assemble(intorder=6)
Z = solver.solve_frequency(1e6)
```
- Surface current captures skin/proximity effects
- Eddy current: scalar FEM-BEM (mu_r=1) or vector FEM-BEM (any mu_r)
- Limitation: dense O(N^2), ~10k DOF max for direct

## Level 3: FEM (final verification)

- Any geometry, nonlinear, adaptive refinement
- Kelvin transformation for open boundaries
- Via NGSolve / esim_coupled_solver.py

## Cross-Validation

Any two levels validate each other on the same geometry:
- PEEC vs NGSBEM: ~5% (filament vs surface current)
- PEEC vs Analytical: ~0.1% (same approximation)

## Magnetic Core Across Levels

| Level | Core Method | Eddy | Nonlinear |
|-------|------------|------|-----------|
| 1 PEEC | Radia HDiv-VIM (Delta_L) | No | Yes |
| 2 NGSBEM | Scalar/Vector FEM-BEM | Yes | No |
| 3 FEM | Volume FEM | Yes | Yes |

## Design Workflow

1. **Explore** (PEEC): Sweep 100 core positions in 10s → select top 5
2. **Validate** (NGSBEM): Z(f) at 20 freqs in 5min → confirm trends
3. **Sign-off** (FEM): Nonlinear solve → final report
"""

RADIA_PEEC_CORE_PITFALLS = """
# PEEC + Magnetic Core: Common Pitfalls

## Critical Issues (will produce wrong results silently)

1. **Coordinates in meters**: Radia always uses meters. `60mm` = `0.06`, not `60`.

2. **Call `rad.UtiDelAll()` first**: Radia keeps global state. Previous objects persist.

3. **NGSBEM: Use `Glue(wire.faces)` for surface mesh**:
   Volume mesh causes BEM cond=1e17 (singular). Surface-only gives cond=1e4.
   ```python
   geo = OCCGeometry(Glue(wire.faces))   # CORRECT
   # geo = OCCGeometry(wire)              # WRONG (volume mesh)
   ```

4. **NGSBEM: Set maxh <= min_cross_section / 2**:
   For 1mm wire: `maxh=0.5e-3`. Larger creates elongated triangles → bad SL entries.

5. **Soft iron is mesh-backed HDiv-VIM**:
   keep material labels on the NGSolve mesh and inspect the HDiv charge map,
   charge-Gram H-matrix stats, and nonlinear iteration metadata.

6. **No mesh-less soft-iron shortcut**:
   raw Radia polyhedra are fine for fixed magnetization and field evaluation,
   but new soft-iron workflows should use `radia.vim.MeshSoftIron`.

7. **Image symmetry must be verified with a full-model check**:
   on a truly symmetric mesh, reduced-image `rad.Fld` and an explicitly mirrored
   full model should agree to near roundoff at field probes.

8. **Loop port**: Don't use `add_port(n1, n1)`. Split the loop:
   ```python
   n1 = builder.add_node_at(x, y, z)    # port A
   n1b = builder.add_node_at(x, y, z)   # port B (same position!)
   builder.add_port(n1, n1b)
   ```

9. **Hex vertex order still matters for fixed-magnet export**:
   bottom CCW (v0-v3), top CCW (v4-v7).  For soft iron, prefer the NGSolve mesh
   route so element orientation is handled by the FEM layer.

10. **Nonlinear metadata is part of the result**:
    store tolerance, iteration count, convergence flag, and material-state
    update in JSON artifacts.

## Solver Selection

```
Core conducting?
 No  → Radia HDiv-VIM ('radia') [ferrite, laminated steel, nonlinear]
 Yes → mu_r > 1?
        No  → Scalar FEM-BEM ('fembem')     [Al/Cu shield]
        Yes → Vector FEM-BEM ('vector_fembem') [solid steel]
```
"""

RADIA_MAGNETIC_CORE_SOLVER_GUIDE = """
# Magnetic Core Solver Selection Guide

When coupling PEEC conductors with magnetic cores, choose the solver based on
material properties and frequency range.

## Decision Matrix

| Core Type | Solver | core_model | Key Limitation |
|-----------|--------|------------|----------------|
| Ferrite (high mu_r, sigma~0) | Radia HDiv-VIM | 'radia' | No eddy currents |
| Laminated steel (low eff. sigma) | Radia HDiv-VIM + effective mu | 'radia' | Static Delta_L |
| Solid steel (high mu_r + sigma) | Vector FEM-BEM | 'vector_fembem' | Linear only |
| Al/Cu shield (mu_r=1, high sigma) | Scalar FEM-BEM | 'fembem' | mu_r=1 only |
| Nonlinear (B-H curve) | Radia HDiv-VIM | 'radia' | No eddy currents |

## Radia HDiv-VIM vs NGSBEM (vector FEM-BEM)

| Aspect | Radia HDiv-VIM | NGSBEM vector FEM-BEM |
|--------|-----------|----------------------|
| Eddy currents | No | Yes |
| Nonlinear | Yes (B-H) | No |
| Domain | Unbounded (no air mesh) | Unbounded (BEM) |
| Unknowns | H(div) mesh DOFs | FEM volume DOFs |
| Acceleration | HACApK charge Gram | Dense/FMM |
| Best regime | DC / low freq / nonlinear | AC / eddy current |

## API

```python
# Keep the soft-iron state mesh-backed and exchange fields through
# NGSolve GridFunction / CoefficientFunction data.
# PEEC stays conductor/shield oriented; magnetic cores use HDiv-VIM /
# reduced FEM at the application layer.
```
"""

RADIA_HDIV_SOFT_IRON = """
# HDiv Soft-Iron Implementation Notes

The current Radia soft-iron route is mesh-backed HDiv-VIM:

- keep the NGSolve mesh as the owner of material and boundary labels;
- create `radia.vim.MeshSoftIron(mesh, mu_r=... | bh_table=...)`;
- solve with `rad.Solve(..., demag_backend="hdiv")`;
- evaluate the resulting model with `rad.Fld` or pass fields back to NGSolve.

For diagnostics, inspect the charge map B, the charge Gram H-matrix stats,
nonlinear iteration metadata, and image-symmetry materialization.  For truly
symmetric meshes, a reduced image model and an explicit full model should agree
to near roundoff at field probes.
"""


RADIA_FEM_KELVIN_CUBIT = """
# Kelvin Transform: Periodic Kelvin for All Formulations

## Fundamental Result (Sugahara, IEEE Trans. Magn. 2022)

The Kelvin transformation conserves the conformal symmetry of Maxwell's equations.
ALL material properties (mu, epsilon, sigma) transform with the SAME factor:

| Dimension | Material scaling | Formulation-independent |
|-----------|-----------------|------------------------|
| **2D** (circle inversion) | In-plane: **1** (no change), Out-of-plane: **(a/r')^4** | Yes |
| **3D** (sphere inversion) | ALL components: **(a/r')^2** | Yes |

This is a property of electromagnetism, not of any particular formulation.
All formulations (A, H, phi, E) use the same Kelvin weight.

### 3D Material Properties in Exterior Domain

```
mu'_r / mu_r = mu'_theta / mu_theta = mu'_phi / mu_phi = (a/r')^2
sigma'_r / sigma_r = sigma'_theta / sigma_theta = sigma'_phi / sigma_phi = (a/r')^2
```

where a = Kelvin radius, r' = distance from exterior domain center.

## Periodic Kelvin Implementation (Verified)

Two separate mesh domains connected by Periodic identification:
- Interior sphere (origin): physical domain, standard material
- Exterior sphere (offset): Kelvin-mapped exterior, modified material

```python
# 3D HCurl (A-formulation): verified +1.2% L error
# ALL material constants scale by (a/r')^2 (conformal symmetry)
# mu_kelvin = mu0 * (a/r')^2, nu_kelvin = nu0 * (a/r')^2
# sigma_kelvin = sigma0 * (a/r')^2
r_prime_sq = (x - offset)**2 + y**2 + z**2
kelvin_fac = R_K**2 / (r_prime_sq + 1e-20)  # (a/r')^2

nu_cf = mesh.MaterialCF({
    "kelvin": NU_0 * kelvin_fac,
    "air": NU_0, "coil": NU_0,
}, default=NU_0)

# Periodic identification between sphere surfaces
kelvin_int_face.Identify(kelvin_ext_face, "periodic",
                          IdentificationType.PERIODIC)
fes = Periodic(HCurl(mesh, order=1, dirichlet_bbnd="GND"))

# GND at exterior center (maps to physical infinity)
# bonus_intorder=4 for spatially varying nu
a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
```

## Performance: Periodic vs Dirichlet vs Shell

| Method | L error | Extra DOFs | Truncation |
|--------|---------|------------|------------|
| Dirichlet at R=120mm | **+7%** | 0 | PEC wall error |
| Kelvin exterior domain (DEPRECATED) | +1.2% | +20-30% | Small (R_outer) |
| **Periodic Kelvin** | **+1.5%** | **0** | **None (exact)** |

Periodic Kelvin is the recommended method. Shell approach is deprecated.

## Geometry Pattern (OCC)

```python
# Interior: physical domain at origin
inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
inner_air = inner_sphere - coil_shape - workpiece_shape
inner_air.name = "air"

# Exterior: Kelvin-mapped domain at offset
offset = 2.5 * R_K  # separation distance
outer_sphere = Sphere(Pnt(offset, 0, 0), R_K)
outer_sphere.name = "kelvin"

# GND at exterior center
gnd = Vertex(Pnt(offset, 0, 0)); gnd.name = "GND"

# Periodic identification BEFORE Glue
int_face = [f for f in inner_air.faces if f.name == "kelvin_int"][0]
ext_face = [f for f in outer_sphere.faces if f.name == "kelvin_ext"][0]
int_face.Identify(ext_face, "periodic", IdentificationType.PERIODIC)

shape = Glue([inner_air, coil, outer_sphere, gnd])
```

## SIBC + Kelvin for Eddy Current (Induction Heating)

```python
# Robin BC on workpiece surface:
robin = -1j * omega / Z_s
a_bf += robin * u.Trace() * v.Trace() * ds("wp_surface")

# ESIM Karl iteration for nonlinear Z_s:
sol = esim_solver.solve(H_t_rms)
Z_s = relax * sol['Z'] + (1 - relax) * Z_s_old
```

## Cubit Workflow for Quarter-Sphere Periodic Kelvin

CRITICAL: The operation order matters. Mesh air FIRST, then create kelvin and copy mesh.

### Correct Order (verified 2026-03-30)

```
1. Build yoke hex mesh in Cubit (e.g., from journal)
2. Create air sphere at model center, webcut for symmetry (z=0, x=0)
3. Imprint/merge air with yoke
4. Tet mesh air volumes (creates triangle mesh on hemisphere boundary)
5. Create kelvin sphere at offset (e.g., offset_x = R * 10)
   - webcut with z=0 and x=offset_x for quarter
   - delete pieces: keep x > offset_x AND z > 0
6. copy mesh surface <int_hemis> onto surface <ext_hemis>
   source curve <src_c> source vertex <src_v>
   target curve <dst_c> target vertex <dst_v>
7. Tet mesh kelvin volumes (copied surface constrains boundary)
8. Create blocks: yoke, air, kelvin, kelvin_int, kelvin_ext
9. export netgen -> Scale(mm_to_m) -> add_periodic_kelvin
```

### Why This Order

- `copy mesh surface` requires the SOURCE surface to already have triangles
- If you try to create both spheres first and mesh surfaces independently,
  `imprint` with yoke changes interior hemisphere topology and copy fails:
  "Source surface and target surface must be topologically identical"

### Hemisphere Surface Detection

Distinguish curved hemisphere from flat symmetry faces:
- Check vertex distance to sphere center (on sphere? distance ~ R)
- ALSO check bounding box span: hemisphere spans all 3 axes,
  flat faces (z=0, x=0) have zero span in one axis

### Kelvin Offset for IdentifyPeriodicBoundaries

The offset is the TRANSLATION vector from int to ext sphere, NOT the ext center:
```python
# Int sphere at (0, sy, 0), ext at (offset_x, sy, 0)
# Translation = (offset_x, 0, 0), NOT (offset_x, sy, 0)
kelvin_offset = (offset_x * scale, 0.0, 0.0)
add_periodic_kelvin(mesh, kelvin_offset)
```

detect_kelvin_offset() returns the centroid of kelvin region (= ext center),
which is WRONG for asymmetric placement. Compute offset directly.

### Symmetry BC for Omega-Reduced (Quarter C-type Dipole)

```python
# Omega is ODD about z=0 (source ~ H0*z is odd)
#   -> Omega=0 on z=0 (sym_normal, Dirichlet)
# Omega is EVEN about x=0 (symmetric coil + yoke)
#   -> dOmega/dn=0 on x=0 (natural, no constraint)
# For A-formulation, swap: sym_tangential gets Dirichlet
dir_parts = []
if "sym_normal" in boundaries:
    dir_parts.append("sym_normal")  # Omega: Dirichlet
# sym_tangential: natural (no constraint for Omega)
```

### Verify Periodic BC

```python
freedof_before = sum(1 for d in fes_base.FreeDofs() if d)
fes = Periodic(fes_base)
freedof_after = sum(1 for d in fes.FreeDofs() if d)
assert freedof_after < freedof_before, "Periodic BC not working!"
```

## Air Gap Mesh (CRITICAL for Accuracy)

For C-type magnets with narrow gaps, the air in the gap MUST have its own
fine mesh, separate from the far-field air mesh. Without this, the FEM solution
can be 30-50% too low.

### Pattern: air_gap box (from Keiko Sugahara's CEFC 2020 work)

```
# In Cubit: create a small box around the gap region
create brick x 42 y 73 z 13      # slightly larger than gap
move volume {id} x 21 y 0 z 6.5  # center on gap

# Imprint with air sphere and yoke
imprint volume all
merge volume all

# Mesh air_gap with fine size (2mm for 10mm gap)
volume {air_gap_vid} size 2.0
mesh volume {air_gap_vid}

# Mesh outer air with coarse size
volume {air_vid} size 37.0
mesh volume {air_vid}
```

### Naming Convention

| Block | Mesh size | Purpose |
|-------|-----------|---------|
| `air_gap` | 1-3mm (gap/5) | Gap region, fine mesh for field resolution |
| `air` | R_kelvin/5 | Far-field air, coarse |
| `kelvin` | R_kelvin/3 | Kelvin exterior domain |

### Why This Matters

The magnetic field in the gap changes rapidly over the gap width (5-10mm).
Without resolving this gradient:
- Linear mu=1000: Bz = -709 mT (without air_gap) vs -879 mT (with air_gap)
- Reference target: Bz near -1 T for the detailed C-yoke benchmark
- Factor 1.24x improvement from gap mesh alone

### Pyramid Element Warning

Cubit `imprint` between hex volumes and tet volumes creates pyramid transition
elements. These DO NOT support order >= 2 in NGSolve (DOF count mismatch).
For order=2, either:
- Use pure tet mesh for all air regions (no hex-tet transition)
- Or avoid imprint and use conformal tet-only near yoke interface

## Known Issues

1. **Cubit scipy conflict**: Import scipy BEFORE cubit on Windows
2. **NGSolve DLL conflict**: Import ngsolve BEFORE adding Cubit to sys.path
3. **bonus_intorder=4**: Required for Kelvin domain integration accuracy
4. **Gauge regularization**: Small mass term (1e-8 * nu0 * u*v dx)
5. **Verify Periodic**: Check FreeDofs reduced after Periodic() wrapping
6. **Cubit copy mesh**: Fails if source surface has no mesh or topology differs
7. **Pyramid order>=2**: Cubit imprint creates pyramids; NGSolve order=2 fails on them

## Reference

Sugahara, K., "Electromagnetic Analysis of Eddy Current Testing With
Kelvin Transformation," IEEE Trans. Magn., 2022.
"""

RADIA_NGSBEM_MQS_LIMITS = """
# ngsolve.bem MQS Limitations for Eddy Current SIBC (2026-03-28)

## Critical: BEM Formulations Have Severe Limitations in MQS-SIBC

When using ngsolve.bem for eddy current analysis with Surface Impedance Boundary
Condition (SIBC) in the MQS regime, most BEM formulations fail or are limited.

### Formulation Summary

| Formulation | Works in MQS? | Limitation |
|-------------|---------------|------------|
| EFIE-SIBC (LaplaceSL) | Only Z_s/(jw*mu0*R) < 0.1 | SL eigenvalue R/3 != R, factor-of-3 error |
| MFIE tangential (LaplaceDL) | PEC only | No Z_s dependence (always gives J_pec) |
| MFIE normal (Sugahara) | All Z_s | Not in ngsolve.bem (needs Biot-Savart matrix) |
| PMCHWT-SIBC | Impossible | jw*eps*SL(M) ~ 10^{-14} in MQS |
| FEM-SIBC | All Z_s | Not BEM (uses NGSolve FEM + Robin penalty) |

### When to Use What

- **PEC conductors**: BEM EFIE or MFIE tangential (ngsolve.bem)
- **Copper at high frequency** (Z_s/(jw*mu0*R) < 0.1): BEM EFIE-SIBC
- **Steel, ferrite, any finite Z_s**: FEM-SIBC (`fem_esim_3d.py`)
- **Circuit extraction (L, R)**: Radia PEEC pipeline (not affected by this limitation)

### Impact on Radia Multi-Level Simulator

Level 1 (PEEC) and Level 2 (ngsolve.bem inductance extraction for PEC/low-Z_s)
are unaffected. The limitation only affects eddy current SIBC problems at Level 2
with high-Z_s materials. Use FEM-SIBC (Level 3) for these cases.

### Verification

- `validation_test/induction_heating/cubit_panels_legacy/verify_sphere_sibc.py` -- Sphere benchmark
- `validation_test/induction_heating/cubit_panels_legacy/fem_esim_3d.py` -- FEM-SIBC reference
"""

RADIA_MATLAB_MEX = """
# MATLAB / Python / NGSolve bridge

The production MATLAB boundary is `radia_mex`, backed by the same C++ Radia
and NGSolve-facing implementation used by the Python stack.  The live command
inventory is exposed by the `matlab_radia_mex_contract` MCP tool and by
`radia.quickCheck()` in MATLAB.

The bridge is intentionally contract-based:

* Python keeps pybind11 objects and NumPy arrays.
* MATLAB uses MEX handles plus numeric matrices, vectors, and diagnostic
  structs.
* NGSolve remains the owner of meshes, finite-element spaces, Piola maps,
  curved-element transformations, and orientation.  MATLAB may dump assembled
  matrices and metadata, but it does not reimplement those transformations.
* `radia.ngsolve.CoefficientFunction` keeps native NGSolve coefficient
  expression trees behind checked MEX handles.  Constants, arithmetic,
  scaling, metadata, and mapped physical-point evaluation are available
  without a Python process.
* `radia.ngsolve.GridFunction` keeps native NGSolve GridFunctions behind
  checked MEX handles.  MATLAB can create real or complex H1/HCurl/HDiv
  spaces, exchange DoF vectors, interpolate a native coefficient function,
  and obtain a GridFunction-backed CoefficientFunction view.
* `radia.ngsolve.Vector` keeps a GridFunction component or an independent
  native `BaseVector` work copy behind a checked MEX handle.  `setZero`,
  `scale`, `axpy`, `dot`, and `norm` stay in C++; `values()` and `setValues()`
  are explicit MATLAB observation/control boundaries.
* `radia.ngsolve.Mesh`, `FESpace`, `BilinearForm`, and `Matrix` keep the
  NGSolve object graph alive behind checked MEX handles.  A form assembles a
  built-in real integrator once; the resulting matrix can export 1-based sparse
  triplets, create native vectors, apply matvecs, and construct a free-DoF
  inverse without rebuilding the mesh.
* The persistent matrix boundary is intentionally explicit rather than a
  claim of full NGSolve Python parity: arbitrary Python callbacks, tensor-valued
  bilinear forms, general preconditioners, and solver objects
  remain outside the current MEX slice. Scalar CoefficientFunction weighting
  is available for the built-in real/complex volume integrators, and real or
  complex volume and boundary CoefficientFunction right-hand sides are
  assembled natively for H1, HCurl, and HDiv spaces.
* `radia.ngsolve.GridFunction.fromFESpace` shares an existing native space
  without reloading the mesh.  `radia.ngsolve.LinearForm` exposes real or
  complex constant-source and native CoefficientFunction volume or boundary
  RHS assembly for H1, HCurl, and HDiv spaces as native vector views, including
  lifetime retention after the form wrapper is released.
* The native MEX gateway links the NGSolve C++ libraries directly and does not
  link or start Python.  The legacy `radentry.cpp` compatibility layer is
  compiled with Python callback support disabled for MEX; numeric Radia calls
  remain available and callback objects fail explicitly at the boundary.
* Fixed reduced IH and HCurl Eddy Bubble/CLN models use
  `simulink.state_space.create/info/step/reset/destroy`.  Matrices and the
  initial state cross MATLAB/MEX once at block start, while each Simulink step
  transfers only the input and output.  Moving height-family interpolation
  remains a MATLAB S-function because its operator changes with position.
* `radia.optuna.Study` stores trial, parameter, intermediate-value, and user
  attribute tables in a MAT-file.  `bestValue`, `bestParams`, and
  `bestSolution` expose the persisted single-objective best after restart;
  `paretoFront` is the multi-objective route.  `radia.optuna.SimulinkRunner`
  uses `SimulationInput -> sim -> score -> Study.tell`.

This is a shared numerical interface, not shared Python/MATLAB object identity.
The MATLAB MEX path does not launch `python.exe`, but the current pip-provided
Windows `libngsolve.dll` has a transitive `python312.dll` dependency. A truly
Python-DLL-free deployment requires rebuilding NGSolve/Netgen without Python
support.
The MCP contract audits retired unsafe constructors and fails if
`ObjMltExtPgn`, `ObjMltExtRtg`, or `ObjMltExtTri` reappears in pybind11, MEX, or
the legacy C ABI source and export table.
"""


def get_radia_documentation(topic: str = "all") -> str:
    """Return Radia usage documentation by topic."""
    topics = {
        "overview": RADIA_OVERVIEW,
        "geometry": RADIA_GEOMETRY,
        "materials": RADIA_MATERIALS,
        "solving": RADIA_SOLVING,
        "parallelization": RADIA_PARALLELIZATION,
        "fields": RADIA_FIELDS,
        "removed_apis": RADIA_REMOVED_APIS,
        "mesh_import": RADIA_MESH_IMPORT,
        "best_practices": RADIA_BEST_PRACTICES,
        "peec": RADIA_PEEC,
        "ngsbem_peec": RADIA_NGBEM_PEEC,
        "efie_preconditioner": RADIA_EFIE_PRECONDITIONER,
        "fem_verification": RADIA_FEM_VERIFICATION,
        "scalar_potential": RADIA_SCALAR_POTENTIAL,
        "vector_potential": RADIA_VECTOR_POTENTIAL,
        "play_models": RADIA_PLAY_MODELS,
        "hysteresis": RADIA_HYSTERESIS,
        "esim": RADIA_ESIM,
        "build_and_release": RADIA_BUILD_AND_RELEASE,
        "hdiv_soft_iron": RADIA_HDIV_SOFT_IRON,
        "magnetic_core_guide": RADIA_MAGNETIC_CORE_SOLVER_GUIDE,
        "peec_core_pitfalls": RADIA_PEEC_CORE_PITFALLS,
        "multilevel": RADIA_MULTILEVEL_SIMULATOR,
        "ngsbem_mqs_limits": RADIA_NGSBEM_MQS_LIMITS,
        "fem_kelvin_cubit": RADIA_FEM_KELVIN_CUBIT,
        "matlab_mex": RADIA_MATLAB_MEX,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(topics.keys())}"
        )
