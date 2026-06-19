"""
Netgen/NGSolve high-order curving workflow knowledge base.

This is the most important knowledge module. It covers:
- Export path: APREPRO command (C++) `export netgen`
- Workflow: Cubit geometry -> mesh -> export netgen -> NGSolve Mesh()
- CallbackGeometry and ACIS curving (compact_netgen)
- .vol file as sole interface between Cubit and NGSolve
- Troubleshooting guide
"""

WORKFLOW_OVERVIEW = """
# Netgen High-Order Curving: Workflow Overview

## Two Export Paths (produce identical results)

### Path A: APREPRO Command (recommended, fast)
```
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
# -> mesh.vol (with curvedelements section) + mesh.vol.json (CAD reference)
```
Uses NetgenCurver (compact_netgen C++ static link). No Python, no DLL dependency.

### Path B: Python (reference, deprecated)
Path B (`extract_curved_mesh`) has been removed. Use Path A (`export netgen`) for all workflows.

Both paths use **CallbackGeometry** to delegate surface/edge projection
to Cubit's ACIS kernel via `closest_point_trimmed`. No STEP files,
no OCC geometry, no SetGeomInfo needed.

## .vol as Sole Interface

```
Cubit (ACIS geometry) -> export netgen -> .vol (self-contained)
                                           |
                              NGSolve: Mesh("mesh.vol")
                              (no Cubit, no STEP needed)
```

The .vol file contains: mesh points, volume elements, surface elements,
material labels, boundary labels, and curvedelements section (high-order
curving coefficients). NGSolve reads it without any geometry file.

## Choose Your Workflow

1. **Is your geometry planar (no curved surfaces)?**
   -> Use any export format. Curving is not needed.
   -> Simplest: `export netgen "mesh.vol" order 1` + `Mesh("mesh.vol")`

2. **Do you only need 2nd order (not 3rd+)?**
   -> Use `export netgen`:
   ```
   Cubit -> mesh -> export netgen "mesh.vol" order 2 -> Mesh("mesh.vol")
   ```

3. **Do you need 3rd order or higher?**
   -> Use `export netgen` with higher order:
   ```
   Cubit -> mesh -> export netgen "mesh.vol" order 3 -> Mesh("mesh.vol")
   ```
   Works for ANY geometry shape — cylinder, sphere, torus, cone,
   Boolean operations, freeform surfaces, etc. Supports order 1-5.

## Accuracy: p-Convergence Results (Verified 2026-04-02)

All shapes tested with ACIS CallbackGeometry + edge snapping:

| Shape | Surfaces | Curves | p=2 V err | p=3 V err | p=5 V err | p=5 A err |
|-------|----------|--------|-----------|-----------|-----------|-----------|
| Sphere | 1 | 0 | -0.023% | +0.002% | -0.000003% | -0.000003% |
| Cylinder | 3 | 2 | -0.003% | +0.001% | -0.000001% | -0.000001% |
| Frustum (cone) | 3 | 2 | -0.009% | +0.002% | -0.000006% | -0.000009% |
| Torus | 1 | 0 | -0.026% | +0.004% | -0.000018% | -0.000011% |
| Box with hole | 7 | 14 | +0.004% | -0.001% | +0.000003% | +0.000003% |

Key: p=5 achieves 10^-5 to 10^-6 % error for ALL shapes, matching OCC native accuracy.

| Method | p=2 Error | p=5 Error | Max Order | Complexity |
|--------|-----------|-----------|-----------|------------|
| export netgen (ACIS) | ~0.003-0.03% | ~1e-6% | 5 (tet/hex/wedge) | Low |
| OCC mesh.Curve() | ~0.003-0.03% | ~1e-6% | 5+ (tet) | Low |
| 1st order (no curving) | ~1.4% | N/A | 1 | None |

## Key Principle

`export netgen` uses Cubit's ACIS kernel for surface projection via
CallbackGeometry (compact_netgen C++ static link). The mesh curving is done
entirely in the C++ plugin — no STEP files, no OCC geometry, no Python
dependency. The .vol file contains the curved mesh ready for NGSolve.
"""

WORKFLOW_EXPORT_CURVED = """
# export netgen APREPRO Command Reference

## Signature

```python
import tempfile
from ngsolve import Mesh
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `order` | int (1-5) | Polynomial order for mesh curving (1=linear, 2=quadratic, ..., 5=quintic) |
| `overwrite` | flag | Overwrite existing file |

## Returns

`.vol` file with curvedelements section + companion `.vol.json` with CAD reference values.

## How It Works Internally

1. Reads mesh topology (nodes, elements, surface tris/quads) from Cubit
2. Extracts **1D segment elements** on geometry curves (edges between surfaces)
   - Each segment has surfnr1/surfnr2 (0-based FD indices of adjacent surfaces)
   - Arc-length normalized dist parameter for each endpoint
3. Creates a `netgen.meshing.Mesh` with volume + surface + segment elements
4. Uses **CallbackGeometry** with ACIS callbacks (C++ static link):
   - `project_func`: Project point onto Cubit surface (ACIS closest_point_trimmed)
   - `normal_func`: Surface normal at a point (ACIS normal_at)
   - `edge_project_func`: Project point onto Cubit curve (ACIS curve.closest_point_trimmed)
5. Calls `BuildCurvedElements(order)`:
   - Surface nodes: projected via `PointBetween` (surface callback)
   - Edge nodes (on curves between surfaces): projected via `PointBetweenEdge` (curve callback)
   - Requires segments with correct surfnr1/surfnr2 for `use_edge` flag
6. Writes the curved mesh to `.vol` file with curvedelements section

### Critical Implementation Details (Edge Snapping)

- **1D segments are required**: Without them, `BuildCurvedElements` does not set
  `use_edge[edgenr]=1` and edge curving is skipped entirely.
- **surfnr1/surfnr2 must be 0-based** FD indices (not -1). BuildCurvedElements passes
  these directly to `PointBetweenEdge` as surfi1/surfi2.
- **CallbackGeometry receives 0-based surfnr**: The `PointBetweenEdge` implementation
  must check `surfi >= 0` (not `> 0`) and convert to 1-based before calling
  `edge_project_func` (which uses 1-based Cubit surface/curve indices).
- **Without edge snapping**: Sphere (1 surface, no edges) works perfectly, but
  cylinder (3 surfaces, 2 edges) area does not p-converge (stuck at -0.4%).
- **With edge snapping**: Both sphere and cylinder converge to 1e-6% at p=5.

## Why ACIS (Not OCC)?

- **No seam lines**: ACIS represents a cylinder as 1 surface; OCC splits
  it at the seam into 2 faces. ACIS has no parametric discontinuities.
- **No STEP needed**: Direct access to Cubit's geometry kernel, no file exchange.
- **Universal**: Works for ANY surface type — analytic, BSpline, freeform,
  Boolean results. Not limited to cylinder/sphere/torus/cone.
- **Exact**: ACIS surface projection is exact for the underlying geometry.

## Example: Basic Usage

```python
import tempfile
from ngsolve import Mesh

# Create and mesh geometry in Cubit
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Export with curving
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Verify
from ngsolve import Integrate, CF
import math
R, H = 0.5, 2.0
expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Example: Any Geometry Shape

```python
import tempfile
from ngsolve import Mesh

# Complex Boolean geometry — no special handling needed
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create cylinder height 4 radius 0.3")
cubit.cmd("subtract volume 2 from volume 1")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Works for any geometry — cylinder holes, fillets, chamfers, BSplines...
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

## Requirements

- Coreform Cubit 2025.12+ with Radia plugin installed (`cubit-plugin-install`)
- NGSolve 6.2.2603+ (curvedelements Load, hex/prism curving)
"""

WORKFLOW_CYLINDER = """
# Example: Cylinder

## Step-by-Step

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

R = 0.5   # Radius
H = 2.0   # Height

# Step 1: Create geometry in Cubit
cubit.cmd(f"create cylinder height {H} radius {R}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "domain"')

# Step 3: Export with curving (no STEP, no OCC, no SetGeomInfo!)
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Note

No STEP export/reimport needed. No OCCGeometry. No set_cylinder_geominfo().
`export netgen` handles everything via Cubit's ACIS kernel.
"""

WORKFLOW_SPHERE = """
# Example: Sphere

## Step-by-Step

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

R = 0.5

# Step 1: Create geometry in Cubit
cubit.cmd(f"create sphere radius {R}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Step 3: Export with curving
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = 4/3 * math.pi * R**3
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```
"""

WORKFLOW_TORUS = """
# Example: Torus

## Step-by-Step

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

R_MAJOR = 1.0
R_MINOR = 0.3

# Step 1: Create geometry in Cubit
cubit.cmd(f"create torus major {R_MAJOR} minor {R_MINOR}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.08")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Step 3: Export with curving
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = 2 * math.pi**2 * R_MAJOR * R_MINOR**2
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Half-Torus Coil (IH Sample)

For induction heating, a half-torus coil with source/sink terminals:

```python
# Native torus + webcut = clean curving (1 toroidal surface)
cubit.cmd("create torus major radius 0.11 minor radius 0.01")
cubit.cmd("webcut volume 1 with plane xplane noimprint nomerge")
cubit.cmd("delete volume 2")  # Keep half-torus

cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "coil"')
cubit.cmd('sideset 1 add surface 2')
cubit.cmd('sideset 1 name "source"')
cubit.cmd('sideset 2 add surface 3')
cubit.cmd('sideset 2 name "sink"')
```

NOTE: Do NOT use `sweep` to create torus geometry for high-order meshing.
`create torus` (native ACIS) produces a single toroidal surface that
curves correctly. `sweep` splits the surface at z=0, which can cause
cross-projection issues in ACIS closest_point_trimmed at high order.
Use `create torus` + `webcut` instead.
"""

WORKFLOW_COMPLEX = """
# Complex Geometry (Boolean Operations)

With `export netgen`, complex geometries require NO special workflow.
Boolean operations, multiple curved surfaces, freeform surfaces — all
are handled automatically by the ACIS kernel.

## Example: Brick with Cylindrical Hole

```python
import tempfile
from ngsolve import Mesh, Integrate, CF
import math

BRICK_SIZE = 2.0
R_HOLE = 0.3

# Step 1: Create geometry with Boolean operations in Cubit
cubit.cmd(f"create brick x {BRICK_SIZE} y {BRICK_SIZE} z {BRICK_SIZE}")
cubit.cmd(f"create cylinder height {BRICK_SIZE*2} radius {R_HOLE}")
cubit.cmd("subtract volume 2 from volume 1")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Step 3: Export with curving — works for any geometry
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)

# Step 4: Verify
expected_vol = BRICK_SIZE**3 - math.pi * R_HOLE**2 * BRICK_SIZE
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Why No Special Workflow?

The old workflow required:
- Creating geometry in OCC (not Cubit)
- Calling name_occ_faces() for face name mapping
- STEP export from OCC
- STEP reimport in Cubit with 'noheal'
- export_netgen_with_names() for name-based mapping
- set_*_geominfo() for each curved surface type

With `export netgen`, ALL of this is eliminated. The ACIS kernel handles
surface projection for any geometry directly, regardless of complexity.
"""

WORKFLOW_GMSH_2ND_ORDER = """
# Alternative: Netgen .vol 2nd Order Workflow (Simplest)

If you only need 2nd order elements and don't need 3rd order or higher,
the APREPRO export netgen command provides a simple workflow.

## Step-by-Step

```python
from ngsolve import Mesh, Integrate, CF

# Step 1: Create geometry and mesh in Cubit
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Step 2: Register blocks
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')

# Step 3: Export to Netgen .vol with order 2
cubit.cmd('export netgen "mesh.vol" order 2 overwrite')

# Step 4: Read into NGSolve
mesh = Mesh("mesh.vol")
# Done! No geometry reference needed at compute time.
```

## Advantages
- No geometry reference needed at compute time
- Very simple workflow (APREPRO command)
- Good accuracy (~0.003%)

## Limitations
- Supports order 1-5 (arbitrary order via ACIS CallbackGeometry)
"""

WORKFLOW_ACCURACY = """
# Accuracy Guide: Choosing the Right Order

## Volume Error by Method and Order

| Method | Order 1 | Order 2 | Order 3 | Order 4 | Order 5 |
|--------|---------|---------|---------|---------|---------|
| No curving | ~1.4% | - | - | - | - |
| export netgen | ~1.4% | ~0.003% | ~0.0004% | ~0.00005% | ~0.000006% |

## When Higher Order Matters

- **Structural/thermal FEM**: Order 2 usually sufficient
- **Electromagnetics (curl-curl)**: Order 2-3 recommended
- **BEM inductance extraction**: Order 3 recommended (surface accuracy critical)
- **Acoustic/wave propagation**: Order 3-5 for dispersion control
- **Geometry verification only**: Order 2 is fine

## Mesh Size vs Order Trade-off

For a target accuracy, you can either:
- **h-refinement**: More elements, keep order low
- **p-refinement**: Fewer elements, increase order

High-order curving (order 3+) is most beneficial when:
- Geometry has high curvature
- Coarse meshes are needed (computational cost)
- High accuracy is required on curved boundaries

## Verification Pattern

Always verify accuracy after curving:

```python
import math
from ngsolve import Integrate, CF, BND

# Volume check
expected_vol = math.pi * R**2 * H  # Exact volume
computed_vol = Integrate(CF(1), mesh)
vol_error = abs(computed_vol - expected_vol) / expected_vol * 100

# Surface area check (optional)
expected_area = 2 * math.pi * R * H + 2 * math.pi * R**2
computed_area = Integrate(CF(1), mesh, VOL_or_BND=BND)
area_error = abs(computed_area - expected_area) / expected_area * 100

print(f"Volume error: {vol_error:.4f}%")
print(f"Area error: {area_error:.4f}%")
```
"""

WORKFLOW_FIELD_SOLVE = """
# High-Order Hex in an Actual FIELD Solve (not just volume)

Curving improves not only the VOLUME integral but a solved FIELD quantity. Verified on the
coaxial-capacitor Laplace problem: a hex-meshed annular TUBE (inner radius a, outer b,
length L), exported `order N`, carried into an NGSolve Laplace solve, gives the per-length
capacitance against the exact closed form

    C/L = 2 pi eps / ln(b/a)        (eps = 1 here).

Because the exact solution is purely radial V(r)=ln(b/r)/ln(b/a), the flat ends are exact
Neumann (dV/dz = 0) -- NO end effects -- so the only error is the hex geometry's chord error
on the CURVED cylindrical electrodes, which curving removes.

## Result (a=1, b=2, L=2; 707 hexes; FIELD order fixed, GEOMETRY order varied)

| geom order | volume err | C/L err  |
|------------|------------|----------|
| 1 (faceted)| +0.049%    | +0.854%  |
| 2 (curved) | +0.001%    | +0.002%  |
| 3 (curved) | +0.000%    | +0.000%  |

Faceted hex mis-states the capacitance by ~0.85% (chord error on the round electrodes);
curved hex (order 2-3) nails it. Curving matters for FIELDS, not just geometry checks.

## Boundary conditions without exported sideset names

`export netgen` may not carry your Cubit sidesets as named NGSolve boundaries. Robust trick:
select boundaries by GEOMETRY in NGSolve and impose Dirichlet by a penalty. For the tube,
pick the cylindrical electrodes by their RADIAL normal (|n_z| < 0.5) -- this excludes the
flat ends (|n_z| ~ 1), which must stay Neumann. Do NOT select the electrodes by radius alone:
the flat-end corner rings also have r ~ a or b and pinning them to V=1/0 over a finite width
wrongly distorts the field (gave a stuck +7.5% error until switched to the normal selector).

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction, grad, dx, ds,
                     Integrate, IfPos, sqrt, x, y, specialcf)
mesh = Mesh("tube_o3.vol")                 # high-order hex from Cubit
r = sqrt(x*x + y*y); n = specialcf.normal(3)
cyl    = IfPos(0.25 - n[2]*n[2], 1.0, 0.0) # |n_z|<0.5 -> cylindrical face only
on_in  = cyl*IfPos(0.5*(a+b) - r, 1.0, 0.0)
mask   = cyl                               # penalise both electrodes
fes = H1(mesh, order=3); u, v = fes.TnT(); alpha = 1e7
A = BilinearForm(fes, symmetric=True); A += grad(u)*grad(v)*dx + alpha*mask*u*v*ds
f = LinearForm(fes); f += alpha*on_in*v*ds  # target V=1 inner, 0 outer
A.Assemble(); f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = A.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*f.vec
C_per_L = Integrate(grad(gfu)*grad(gfu), mesh, order=10) / L   # = 2W/V^2 / L, V=1
```

Same pattern (swap eps->mu/sigma/k, or curl-curl for magnetics) carries any field
solve onto curved high-order hex. Loader rule from `troubleshooting` still applies: read the
curving with `Mesh()` + high `order=` quadrature, never `m.Curve()` (no CAD ref in the .vol).

## EIGENVALUE solves on the high-order hex mesh (not just source-driven fields)

The loaded curved-hex `.vol` also supports EIGENVALUE problems, where curving matters most
(eigenvalues are sensitive to the boundary shape). The Laplace-Dirichlet spectrum
`-nabla^2 u = lambda u, u=0 on the surface` of a meshed BALL has the exact lowest eigenvalue
`(pi/R)^2` (radial s-mode `sin(pi r/R)/r`), then `(4.493409/R)^2` (l=1, triply degenerate).

Build a hex SPHERE (`volume 1 scheme sphere` O-grid), `block 1 add volume all`,
`export netgen ... order 3`; in a separate ngsolve process load the `.vol` AS-IS (no
`mesh.Curve`) and call
`radia_mcp.radia_ngsolve.waveguide.laplace_dirichlet_eigenvalues(mesh, n, order=3)`. On a coarse
56-hex order-3 sphere (R=0.5) this gives lambda_1 = 39.450 vs exact `4 pi^2 = 39.478` (rel err
7e-4) and lambda_2 = 80.71 vs the l=1 mode 80.76 -- the curved hex recovers the round-domain
spectrum on very few elements. These eigenvalues are also the modal DECAY rates of the transient
heat equation (`multiphysics.solve_heat_transient`, `T_n ~ exp(-alpha lambda_n t)`).

The same works on a hex CYLINDER (`create cylinder height Lz radius R`): the Dirichlet spectrum
is `lambda = (j_mn/R)^2 + (p pi/Lz)^2` (Bessel-zero radial + axial); on a 544-hex order-3 cylinder
(R=0.5, Lz=1) lambda_1 = 33.0022 vs the exact `(j01/R)^2+(pi/Lz)^2 = 33.0024` (rel err 4e-6), with
lambda_2/lambda_3 matching the (0,1,2) and (1,1,1) modes -- the curved lateral surface is carried
by the order-3 export.
"""

TROUBLESHOOTING = """
# Troubleshooting High-Order Curving

## export netgen() Fails or Produces Wrong Results

### Symptom: RuntimeError during Curve()

**Cause 1**: Missing boundary blocks
```
Fix: Ensure both volume and surface blocks are registered:
     cubit.cmd("block 1 add tet all")
     cubit.cmd("block 2 add tri all")
```

**Cause 2**: Mesh quality too poor for high-order curving
```
Fix: Reduce element size or improve mesh quality:
     cubit.cmd("volume all size 0.05")  # Smaller elements
     cubit.cmd("smooth volume all")     # Improve quality
```

## Empty Mesh (0 Elements)

**Cause**: No blocks registered
```
Fix: cubit.cmd("block 1 add tet all")
     cubit.cmd("block 2 add tri all")
```

## Missing Boundary Elements

**Cause**: Only volume element block, no surface element block
```
Fix: cubit.cmd("block 2 add tri all")  # Add boundary elements
```

## Volume Error > 1%

**Cause**: 1st order mesh without curving
```
Fix: Use export netgen(cubit, order=2) or higher
     Or use Gmsh 2nd order alternative
```

## Netgen Import Error: "No module named 'netgen'"

```
Fix: Use system Python with CUBIT_PATH environment variable.
     Cubit's bundled Python cannot import ngsolve.
     System Python with CUBIT_PATH can access BOTH Cubit API and NGSolve.

     # Step 1: Set CUBIT_PATH
     set CUBIT_PATH="C:/Program Files/Coreform Cubit 2025.12/bin"

     # Step 2: Run with system Python (which has NGSolve installed)
     python my_script.py
```

In the script:
```python
import sys, os

# CRITICAL: Import NGSolve BEFORE Cubit to avoid DLL conflicts
import ngsolve
from ngsolve import Mesh

cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])
```

**Key insight**: By using system Python with `CUBIT_PATH`, scripts can access both
the Cubit API and NGSolve/Netgen simultaneously. This is essential for the
export netgen() workflow.

## NGSolve/Cubit DLL Conflict: Import Order Matters

When using system Python with both NGSolve and Cubit, **NGSolve MUST be imported
BEFORE Cubit**. If Cubit is imported first, its bundled DLLs (VTK, etc.) conflict
with NGSolve's Netgen library, causing `ImportError: initialization failed` on
`from netgen import libngpy`.

**Correct import order** (NGSolve first):
```python
import ngsolve                    # MUST be first - loads Netgen DLLs cleanly
from ngsolve import Mesh

import sys, os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit                      # Safe: Netgen DLLs already loaded
cubit.init(['cubit', '-nojournal', '-batch'])
```

**Wrong import order** (causes DLL conflict):
```python
import sys
sys.path.append("C:/Program Files/Coreform Cubit 2025.12/bin")
import cubit                      # Loads Cubit's bundled VTK DLLs
import ngsolve                    # FAILS - Netgen can't initialize
```

**Root cause**: Cubit bundles its own versions of VTK and other shared libraries.
When `import cubit` executes, these DLLs are loaded into the process. When NGSolve
subsequently tries to load Netgen's `libngpy`, the already-loaded Cubit DLLs
conflict with Netgen's expected library versions, causing initialization failure.

**Rule of thumb**: Always `import ngsolve` (and any `from netgen...` imports) at
the very top of the script, before adding Cubit to `sys.path` or importing `cubit`.

## Standalone cubit.init() Segfaults in the Radia Panel (verified 2026-06)

**Symptom**: From system Python, `import cubit; cubit.init([...])` prints the banner,
auto-plays `site-packages/radia/panels/startup.py` ("[Radia] Panel debug log: ..."),
then crashes with exit code -1073741819 (0xC0000005 access violation) BEFORE your first
`cubit.cmd()` runs. Happens with or without `-nographics` and regardless of the ngsolve
import order above. Cause: when the Radia Cubit *panel* plugin is installed, its
`startup.py` is auto-played on init and segfaults under headless embedded Python (it
expects a GUI main window). The single-process recipe above only works when that panel
plugin is absent.

**Fix -- robust two-process pattern** (use this whenever the Radia panel is installed):

  1. EXPORT in a child process via the real Cubit executable -- it degrades the panel
     gracefully ("Cubit main window not found -- Radia Export menu not installed", then
     continues) instead of crashing:

         coreform_cubit.exe -nographics -batch -nojournal mesh_export.py

     In a script launched this way, `cubit` is PRE-INJECTED -- do NOT `import cubit` or
     call `cubit.init()` (that re-inits and can crash). Just use `cubit.cmd(...)`.

  2. LOAD the resulting `.vol` in a SEPARATE ngsolve-only process (no `import cubit`).

  This split also sidesteps the NGSolve/Cubit DLL conflict entirely -- the two libraries
  never share a process, so import order stops mattering.

**`-batch` .py playback is LINE-ORIENTED**: multi-line Python compound statements
(for/if/try blocks, multi-line parenthesized tuples) break with
`SyntaxError: '(' was never closed`. Keep every statement on ONE physical line.
`print()` output may be swallowed -- write results to a file and read it back.

**Exit code 2 / "ERROR: Errors found during session." is benign teardown noise** from the
panel's "main window not found"; verify success from the OUTPUT FILE, not the exit code.

Verified 2026-06-14: hex cylinder (R=0.5, H=2), order-3 `export netgen` (572 hexes,
770 -> 25368 curved nodes) -> NGSolve volume 1.57082 vs pi*R^2*H = 1.57080, error
0.0017% (vs ~0.4% for a 1st-order hex cylinder: curving cuts the error ~250x).
"""

DELETED_APIS = """
# Deleted APIs (Replaced by export netgen APREPRO command)

The following APIs have been **completely removed**.
Do NOT use them — they no longer exist.

## Removed Functions

| Function | Replacement |
|----------|-------------|
| `export_NetgenMesh()` | `export netgen` APREPRO command |
| `export_netgen()` (alias) | `export netgen` APREPRO command |
| `export_netgen_with_names()` | `export netgen` APREPRO command |
| `extract_curved_mesh()` | `export netgen` APREPRO command |
| `name_occ_faces()` | Not needed (no OCC geometry) |
| `set_cylinder_geominfo()` | Not needed (ACIS handles curving) |
| `set_sphere_geominfo()` | Not needed |
| `set_torus_geominfo()` | Not needed |
| `set_cone_geominfo()` | Not needed |
| `compute_cylinder_uv()` | Not needed |
| `compute_sphere_uv()` | Not needed |
| `compute_torus_uv()` | Not needed |
| `compute_cone_uv()` | Not needed |

## Why They Were Removed

The old workflow required multiple steps:
1. STEP export from Cubit (ACIS -> STEP)
2. STEP reimport into Cubit (to match OCC face topology)
3. OCC geometry loading (OCCGeometry(step_file))
4. Mesh export with geometry reference (export_netgen(geometry=geo))
5. Per-shape UV computation (set_*_geominfo())
6. Manual mesh.Curve(order)

`export netgen` replaces ALL of these steps with a single APREPRO command.
It uses CallbackGeometry to delegate surface projection to Cubit's ACIS
kernel directly, without any STEP file exchange or OCC geometry.

## Migration Guide

### Old Code (REMOVED)
```python
# This code NO LONGER WORKS — all these functions are deleted
geo = OCCGeometry("cylinder.step")
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)    # DELETED
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=R, height=H)  # DELETED
mesh = Mesh(ngmesh)
mesh.Curve(3)
```

### New Code
```python
import tempfile
from ngsolve import Mesh
# Single APREPRO command replaces everything
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

### Old Complex Workflow (REMOVED)
```python
# This code NO LONGER WORKS
cubit_mesh_export.name_occ_faces(shape)
shape.WriteStep("geometry.step")
geo = OCCGeometry("geometry.step")
cubit.cmd('import step "geometry.step" noheal')
ngmesh = cubit_mesh_export.export_netgen_with_names(cubit, geo)
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=R_HOLE, height=H)
mesh = Mesh(ngmesh)
mesh.Curve(3)
```

### New Code
```python
import tempfile
from ngsolve import Mesh
# Create geometry directly in Cubit, no OCC needed
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create cylinder height 4 radius 0.3")
cubit.cmd("subtract volume 2 from volume 1")
# ... mesh and blocks ...
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```
"""


def get_netgen_documentation(workflow: str = "overview") -> str:
	"""Return Netgen workflow documentation by topic."""

	KELVIN_AUTO = """
# Kelvin Auto-Add in Cubit Workflow (2026-04-14)

## Overview

Kelvin open-boundary transformation is automatically added when the user
clicks "Radia-NGSolve" -> OK in the Cubit GUI. Kelvin is added automatically
-- there is no separate "Kelvin Transform" menu item.

## How It Works (register_toolbar.py)

1. User creates physical geometry (coil + air + optional workpiece hole)
2. User clicks "Radia-NGSolve" -> selects analysis mode -> OK
3. register_toolbar.py checks: is "kelvin" block already present?
   - YES -> skip, proceed to export
   - NO  -> auto-detect R, symmetry, and call add_kelvin_cubit()
4. export netgen -> .vol (with Kelvin + periodic identification)
5. Launch analysis window (calc_fem.py reads .vol with Kelvin)

## Auto-Detection Logic

### Sphere Radius (R)
- Find the "air" block volumes
- Find the largest-area surface on those volumes -> outer boundary
- R = max vertex distance from origin on that surface
- Works for full sphere, hemisphere (1/2), quarter sphere (1/4), octant (1/8)

### Symmetry Planes
- Check all vertices of air volumes
- If ALL vertices have x >= 0 AND some are at x = 0 -> "x" symmetry
- Same for y, z
- Result: [] (full), ["z"] (1/2), ["x","z"] (1/4), ["x","y","z"] (1/8)

### Offset Direction
- auto_offset_direction(symmetry) picks a free axis
- ["z"] -> offset in x; ["x","z"] -> offset in y

## What add_kelvin_cubit() Does

1. Creates exterior sphere (same R, at offset position)
2. Webcuts for symmetry planes
3. Copies mesh from interior sphere surface -> exterior sphere (1:1 nodes)
4. Tet-meshes the exterior sphere volumes
5. Assigns blocks ("kelvin"), sidesets ("kelvin_int", "kelvin_ext")
6. Creates GND vertex + nodeset at Kelvin center

## CRITICAL: Kelvin Domain Must Be Tet

Spherical geometry is best approximated by triangular high-order elements.
Hex elements on sphere surfaces introduce systematic geometry error.
Always use `scheme tetmesh` for Kelvin volumes.

## User Requirements

The user only needs to provide:
- Coil geometry (block "coil" with source/sink sidesets)
- Air sphere (block "air") containing the coil
- Optional: workpiece hole (sideset "sibc") -- subtracted from air, NOT meshed

The user does NOT need to:
- Create Kelvin geometry manually
- Know the Kelvin sphere radius or offset
- Specify symmetry planes
- Add Kelvin manually (it is auto-added on Solve)

## .jou Example (Minimal)

```python
# User only creates physical geometry:
reset
create sphere radius 0.06     # air sphere
sweep surface 1 axis ...      # coil (inside air sphere)
subtract volume <coil> from volume <air> keep_tool
imprint all; merge all
mesh volume all
block 1 add volume <coil>; block 1 name "coil"
block 2 add volume <air>;  block 2 name "air"
sideset 1 add surface <gap1>; sideset 1 name "source"
sideset 2 add surface <gap2>; sideset 2 name "sink"
# That's it. Kelvin is added automatically on "Solve".
```
"""

	topics = {
		"overview": WORKFLOW_OVERVIEW,
		"export netgen": WORKFLOW_EXPORT_CURVED,
		"simple_cylinder": WORKFLOW_CYLINDER,
		"simple_sphere": WORKFLOW_SPHERE,
		"simple_torus": WORKFLOW_TORUS,
		"complex": WORKFLOW_COMPLEX,
		"complex_named": WORKFLOW_COMPLEX,  # Alias for backward compat
		"accuracy": WORKFLOW_ACCURACY,
		"field_solve": WORKFLOW_FIELD_SOLVE,
		"capacitance": WORKFLOW_FIELD_SOLVE,
		"kelvin_auto": KELVIN_AUTO,
		"gmsh_2nd_order": WORKFLOW_GMSH_2ND_ORDER,
		"troubleshooting": TROUBLESHOOTING,
		"deleted_apis": DELETED_APIS,
		# Legacy aliases that redirect to new content
		"setgeominfo_api": DELETED_APIS,
		"seam_problem": DELETED_APIS,
		"tolerance_tuning": DELETED_APIS,
		"uv_math": DELETED_APIS,
		"multi_surface": WORKFLOW_COMPLEX,
		"freeform": WORKFLOW_COMPLEX,
		"simple_cone": WORKFLOW_EXPORT_CURVED,
	}

	workflow = workflow.lower().strip()
	if workflow == "all":
		main_topics = [
			WORKFLOW_OVERVIEW, WORKFLOW_EXPORT_CURVED,
			WORKFLOW_CYLINDER, WORKFLOW_SPHERE, WORKFLOW_TORUS,
			WORKFLOW_COMPLEX, WORKFLOW_GMSH_2ND_ORDER,
			WORKFLOW_ACCURACY, WORKFLOW_FIELD_SOLVE, TROUBLESHOOTING, DELETED_APIS,
		]
		return "\n\n".join(main_topics)
	elif workflow in topics:
		return topics[workflow]
	else:
		return (
			f"Unknown workflow: '{workflow}'. "
			f"Available: all, {', '.join(k for k in topics if k not in ('complex_named', 'setgeominfo_api', 'seam_problem', 'tolerance_tuning', 'uv_math', 'multi_surface', 'freeform', 'simple_cone'))}"
		)
