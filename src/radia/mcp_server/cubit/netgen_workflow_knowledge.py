"""
Netgen/NGSolve high-order curving workflow knowledge base.

This is the most important knowledge module. It covers:
- The export_NGSolveCurvedMesh() API (the ONLY mesh export function for curved meshes)
- Workflow: Cubit geometry -> mesh -> export_NGSolveCurvedMesh(order=N) -> done
- CallbackGeometry and ACIS UV curving
- Gmsh 2nd order alternative workflow
- Troubleshooting guide
"""

WORKFLOW_OVERVIEW = """
# Netgen High-Order Curving: Workflow Overview

## Background: How export_NGSolveCurvedMesh() Works

`export_NGSolveCurvedMesh()` is the single unified function for exporting Cubit meshes
to NGSolve with high-order curving. It uses **CallbackGeometry** (from the
ksugahar/netgen fork) to delegate surface projection directly to Cubit's
**ACIS kernel** via Python callbacks. This eliminates the need for:

- STEP files (no STEP export/reimport)
- OCC geometry (no OCCGeometry needed)
- SetGeomInfo / UV parameter computation
- Manual per-shape geominfo functions (set_cylinder_geominfo, etc.)
- Seam line workarounds (ACIS has no seam lines)

The overall workflow is simply:

```
1. Create geometry in Cubit
2. Mesh with Cubit
3. mesh = Mesh(extract_curved_mesh(cubit, order=3))
4. Done — mesh is already curved, ready for NGSolve
```

## Choose Your Workflow

1. **Is your geometry planar (no curved surfaces)?**
   -> Use any export format. Curving is not needed.
   -> Simplest: `export_Gmesh(version="2.2")` + `ReadGmsh()`
   -> Or: `export_NGSolveCurvedMesh(cubit, order=1)` for Netgen mesh

2. **Do you only need 2nd order (not 3rd+)?**
   -> You have two options:
   a. **export_NGSolveCurvedMesh()** (recommended):
   ```
   Cubit -> mesh -> export_NGSolveCurvedMesh(cubit, order=2) -> done
   ```
   b. **Gmsh 2nd order** workflow (no geometry reference needed):
   ```
   Cubit -> block element type tetra10 -> export_Gmesh(version="2.2") -> ReadGmsh()
   ```

3. **Do you need 3rd order or higher?**
   -> Use `export_NGSolveCurvedMesh()`:
   ```
   Cubit -> mesh -> export_NGSolveCurvedMesh(cubit, order=3) -> done
   ```
   Works for ANY geometry shape — cylinder, sphere, torus, cone,
   Boolean operations, freeform surfaces, etc.

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
| extract_curved_mesh (ACIS) | ~0.003-0.03% | ~1e-6% | 5+ (tet) | Low |
| OCC mesh.Curve() | ~0.003-0.03% | ~1e-6% | 5+ (tet) | Low |
| Gmsh 2nd order (ReadGmsh) | ~0.001% | N/A | 2 | Low |
| 1st order (no curving) | ~1.4% | N/A | 1 | None |

## Key Principle

`export_NGSolveCurvedMesh()` uses Cubit's ACIS kernel for surface projection via
CallbackGeometry. The mesh curving is done entirely in-memory — no STEP
files, no OCC geometry, no intermediate file I/O. The function returns
an `ngsolve.Mesh` object that is already curved to the requested order.
"""

WORKFLOW_EXPORT_CURVED = """
# export_NGSolveCurvedMesh() API Reference

## Signature

```python
mesh = Mesh(extract_curved_mesh(cubit, order=3, surface_only=False, split_quads=False))
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `order` | int | 3 | Polynomial order for mesh curving (1=linear, 2=quadratic, 3=cubic, ...) |
| `surface_only` | bool | False | If True, export only surface elements (for BEM) |
| `split_quads` | bool | False | If True, split quad elements into triangles |

## Returns

`ngsolve.Mesh` — already curved to the requested order. Ready for FEM/BEM use.

## How It Works Internally

1. Reads mesh topology (nodes, elements, surface tris/quads) from Cubit
2. Extracts **1D segment elements** on geometry curves (edges between surfaces)
   - Each segment has surfnr1/surfnr2 (0-based FD indices of adjacent surfaces)
   - Arc-length normalized dist parameter for each endpoint
3. Creates a `netgen.meshing.Mesh` with volume + surface + segment elements
4. Uses **CallbackGeometry** to register Python callbacks:
   - `project_func`: Project point onto Cubit surface (ACIS closest_point_trimmed)
   - `normal_func`: Surface normal at a point (ACIS normal_at)
   - `edge_project_func`: Project point onto Cubit curve (ACIS curve.closest_point_trimmed)
5. Calls `BuildCurvedElements(order)`:
   - Surface nodes: projected via `PointBetween` (surface callback)
   - Edge nodes (on curves between surfaces): projected via `PointBetweenEdge` (curve callback)
   - Requires segments with correct surfnr1/surfnr2 for `use_edge` flag
6. Returns the curved `netgen.meshing.Mesh`

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
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh

# Create and mesh geometry in Cubit
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Export with curving — that's it!
mesh = Mesh(extract_curved_mesh(cubit, order=3))

# Verify
from ngsolve import Integrate, CF
import math
R, H = 0.5, 2.0
expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Example: Surface-Only (for BEM)

```python
# Surface mesh for BEM inductance extraction
mesh = Mesh(extract_curved_mesh(cubit, order=3, surface_only=True))
```

## Example: Any Geometry Shape

```python
# Complex Boolean geometry — no special handling needed
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create cylinder height 4 radius 0.3")
cubit.cmd("subtract volume 2 from volume 1")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Works for any geometry — cylinder holes, fillets, chamfers, BSplines...
mesh = Mesh(extract_curved_mesh(cubit, order=3))
```

## Requirements

- Coreform Cubit 2025.3+
- ksugahar/netgen fork with CallbackGeometry support
- System Python with both NGSolve and Cubit access (CUBIT_PATH)
"""

WORKFLOW_CYLINDER = """
# Example: Cylinder

## Step-by-Step

```python
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh
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
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

# Step 3: Export with curving (no STEP, no OCC, no SetGeomInfo!)
mesh = Mesh(extract_curved_mesh(cubit, order=3))

# Step 4: Verify
expected_vol = math.pi * R**2 * H
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```

## Note

No STEP export/reimport needed. No OCCGeometry. No set_cylinder_geominfo().
export_NGSolveCurvedMesh() handles everything via Cubit's ACIS kernel.
"""

WORKFLOW_SPHERE = """
# Example: Sphere

## Step-by-Step

```python
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh
from ngsolve import Mesh, Integrate, CF
import math

R = 0.5

# Step 1: Create geometry in Cubit
cubit.cmd(f"create sphere radius {R}")

# Step 2: Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Step 3: Export with curving
mesh = Mesh(extract_curved_mesh(cubit, order=3))

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
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh
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
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Step 3: Export with curving
mesh = Mesh(extract_curved_mesh(cubit, order=3))

# Step 4: Verify
expected_vol = 2 * math.pi**2 * R_MAJOR * R_MINOR**2
vol = Integrate(CF(1), mesh)
print(f"Volume error: {abs(vol-expected_vol)/expected_vol*100:.4f}%")
```
"""

WORKFLOW_COMPLEX = """
# Complex Geometry (Boolean Operations)

With export_NGSolveCurvedMesh(), complex geometries require NO special workflow.
Boolean operations, multiple curved surfaces, freeform surfaces — all
are handled automatically by the ACIS kernel.

## Example: Brick with Cylindrical Hole

```python
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh
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
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Step 3: Export with curving — works for any geometry
mesh = Mesh(extract_curved_mesh(cubit, order=3))

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

With export_NGSolveCurvedMesh(), ALL of this is eliminated. The ACIS kernel handles
surface projection for any geometry directly, regardless of complexity.
"""

WORKFLOW_GMSH_2ND_ORDER = """
# Alternative: Gmsh 2nd Order Workflow (Simplest)

If you only need 2nd order elements and don't need 3rd order or higher,
this is an alternative approach with excellent accuracy (~0.001%).

## Step-by-Step

```python
from cubit_netgen_bridge import extract_curved_mesh
from ngsolve import Mesh
from netgen.read_gmsh import ReadGmsh
from ngsolve import Mesh

# Step 1: Create geometry and mesh in Cubit
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Step 2: Register blocks with 2nd order elements
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 element type tetra10")    # Convert to 2nd order
cubit.cmd("block 2 add tri all")
cubit.cmd("block 2 element type tri6")       # 2nd order boundary

# Step 3: Export to Gmsh v2.2
cubit.cmd('export gmsh "mesh.msh" overwrite')

# Step 4: Read into NGSolve
mesh = Mesh(ReadGmsh("mesh.msh"))
# Done! No geometry reference needed.
```

## Advantages
- No geometry reference needed
- Very simple workflow
- Good accuracy (~0.001%)

## Limitations
- Limited to 2nd order (no 3rd, 4th, 5th...)
- Cubit places mid-edge nodes, not exact CAD placement
- For 3rd+ order, use `export_NGSolveCurvedMesh()` instead

## Comparison

| Feature | export_NGSolveCurvedMesh() | Gmsh 2nd Order |
|---------|----------------|----------------|
| Max order | Unlimited | 2 |
| Geometry reference | ACIS (automatic) | None |
| Complexity | Low | Low |
| Best for | Any order, high accuracy | Quick 2nd order |
"""

WORKFLOW_ACCURACY = """
# Accuracy Guide: Choosing the Right Order

## Volume Error by Method and Order

| Method | Order 1 | Order 2 | Order 3 | Order 4 | Order 5 |
|--------|---------|---------|---------|---------|---------|
| No curving | ~1.4% | - | - | - | - |
| Gmsh 2nd order | - | ~0.001% | - | - | - |
| export_NGSolveCurvedMesh | ~1.4% | ~0.003% | ~0.0004% | ~0.00005% | ~0.000006% |

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

TROUBLESHOOTING = """
# Troubleshooting High-Order Curving

## export_NGSolveCurvedMesh() Fails or Produces Wrong Results

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
Fix: Use export_NGSolveCurvedMesh(cubit, order=2) or higher
     Or use Gmsh 2nd order alternative
```

## Netgen Import Error: "No module named 'netgen'"

```
Fix: Use system Python with CUBIT_PATH environment variable.
     Cubit's bundled Python cannot import ngsolve.
     System Python with CUBIT_PATH can access BOTH Cubit API and NGSolve.

     # Step 1: Set CUBIT_PATH
     set CUBIT_PATH="C:/Program Files/Coreform Cubit 2025.3/bin"

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
export_NGSolveCurvedMesh() workflow.

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
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")
import cubit                      # Loads Cubit's bundled VTK DLLs
import ngsolve                    # FAILS - Netgen can't initialize
```

**Root cause**: Cubit bundles its own versions of VTK and other shared libraries.
When `import cubit` executes, these DLLs are loaded into the process. When NGSolve
subsequently tries to load Netgen's `libngpy`, the already-loaded Cubit DLLs
conflict with Netgen's expected library versions, causing initialization failure.

**Rule of thumb**: Always `import ngsolve` (and any `from netgen...` imports) at
the very top of the script, before adding Cubit to `sys.path` or importing `cubit`.
"""

DELETED_APIS = """
# Deleted APIs (Replaced by export_NGSolveCurvedMesh)

The following APIs have been **completely removed** along with cubit_mesh_export.py.
Do NOT use them — they no longer exist.

## Removed Functions

| Function | Replacement |
|----------|-------------|
| `export_NetgenMesh()` | `export_NGSolveCurvedMesh()` |
| `export_netgen()` (alias) | `export_NGSolveCurvedMesh()` |
| `export_netgen_with_names()` | `export_NGSolveCurvedMesh()` |
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

export_NGSolveCurvedMesh() replaces ALL of these steps with a single function call.
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
# Single function call replaces everything
mesh = Mesh(extract_curved_mesh(cubit, order=3))
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
# Create geometry directly in Cubit, no OCC needed
cubit.cmd("create brick x 2 y 2 z 2")
cubit.cmd("create cylinder height 4 radius 0.3")
cubit.cmd("subtract volume 2 from volume 1")
# ... mesh and blocks ...
mesh = Mesh(extract_curved_mesh(cubit, order=3))
```
"""


def get_netgen_documentation(workflow: str = "overview") -> str:
	"""Return Netgen workflow documentation by topic."""
	topics = {
		"overview": WORKFLOW_OVERVIEW,
		"export_NGSolveCurvedMesh": WORKFLOW_EXPORT_CURVED,
		"simple_cylinder": WORKFLOW_CYLINDER,
		"simple_sphere": WORKFLOW_SPHERE,
		"simple_torus": WORKFLOW_TORUS,
		"complex": WORKFLOW_COMPLEX,
		"complex_named": WORKFLOW_COMPLEX,  # Alias for backward compat
		"accuracy": WORKFLOW_ACCURACY,
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
			WORKFLOW_ACCURACY, TROUBLESHOOTING, DELETED_APIS,
		]
		return "\n\n".join(main_topics)
	elif workflow in topics:
		return topics[workflow]
	else:
		return (
			f"Unknown workflow: '{workflow}'. "
			f"Available: all, {', '.join(k for k in topics if k not in ('complex_named', 'setgeominfo_api', 'seam_problem', 'tolerance_tuning', 'uv_math', 'multi_surface', 'freeform', 'simple_cone'))}"
		)
