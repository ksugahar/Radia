# export_NetgenMesh - Cubit to Netgen Mesh Export

## Overview

`extract_curved_mesh()` creates a `netgen.meshing.Mesh` object directly from Cubit mesh data. When combined with a geometry file, it enables high-order curved mesh generation using `mesh.Curve(order)` method.

This function bridges Cubit's powerful mesh generation capabilities with Netgen/NGSolve's high-order finite element analysis.

**Note**: `netgen.meshing.Mesh` is the core mesh data structure. `ngsolve.Mesh` is a wrapper/view for FEM analysis.

## Key Features

- Direct conversion from Cubit mesh to `netgen.meshing.Mesh` object (no intermediate files)
- Geometry attachment via STEP/BREP/IGES files for `mesh.Curve()` support
- Support for all common 3D element types: Tet, Hex, Wedge, Pyramid
- Block names preserved as Materials and Boundary names in NGSolve
- High-order curving (order 2, 3, 4, 5, ...) via Netgen's geometry-based algorithm

## Design Philosophy

### Why 1st-Order Mesh + Geometry-Based Curving?

This approach separates mesh topology from geometry approximation:

1. **Cubit's role**: Generate high-quality 1st-order mesh (element shapes, transitions, boundaries)
2. **Netgen's role**: Add high-order nodes that conform to the exact geometry

**Advantages over 2nd-order Cubit export:**
- Arbitrary curve orders (2, 3, 4, 5, ...) - not limited to quadratic
- Nodes placed exactly on CAD geometry surfaces
- No node ordering conversion needed for high-order nodes
- Works with all element types (Tet, Hex, Wedge, Pyramid)

### Comparison with Gmsh Export

| Feature | extract_curved_mesh() | radia_export gmsh |
|---------|-----------------|-------------------|
| Max order | Unlimited (via Curve) | 2nd order |
| Intermediate file | None | .msh file |
| Geometry reference | STEP/BREP/IGES | None |
| High-order accuracy | Exact (CAD-based) | Mesh-based |
| Wedge elements | Supported | Supported |
| Pyramid elements | Supported | Supported |

## Recommended: Netgen .vol Workflow

For computation, use the `radia_export netgen` command to produce `.vol` files directly. This supports arbitrary-order curving (order 1-5) and preserves material/boundary labels.

### Basic Workflow

```python
import cubit
from ngsolve import Mesh

# 1. Create and mesh geometry
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")

# 2. Register blocks
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("sideset 1 add surface all")
cubit.cmd('sideset 1 name "boundary"')

# 3. Export and load
cubit.cmd('radia_export netgen "mesh.vol" order 3 overwrite')
mesh = Mesh("mesh.vol")
```

### Volume Accuracy by Curve Order

| Order | Volume Error (sphere) |
|-------|-----------------------|
| 1     | ~1.4%                 |
| 2     | ~0.001%               |
| 3     | ~1e-5%                |
| 5     | ~1e-8%                |

## Usage

### Basic Usage

```python
import cubit
import radia_cubit_mesh
import ngsolve

# Create mesh in Cubit
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'sphere'")
cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'boundary'")

# Export geometry for Curve() support
cubit.cmd('export step "geometry.step" overwrite')

# Convert to Netgen mesh
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, geometry_file="geometry.step")

# Wrap with NGSolve for FEM analysis
mesh = ngsolve.Mesh(ngmesh)

# Apply high-order curving
mesh.Curve(3)  # 3rd order curved elements

# Use in NGSolve
print(mesh.GetMaterials())   # ('sphere',)
print(mesh.GetBoundaries())  # ('boundary',)
```

### Without Geometry (No Curve Support)

```python
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit)
mesh = ngsolve.Mesh(ngmesh)
# mesh.Curve() will not work without geometry
```

### Mixed Element Types

```python
# Hex mesh
cubit.cmd("volume 1 scheme map")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'hex_region'")

# Tet mesh (with pyramid transition)
cubit.cmd("volume 2 scheme tetmesh")
cubit.cmd("mesh volume 2")
cubit.cmd("block 2 add tet all")
cubit.cmd("block 2 name 'tet_region'")
cubit.cmd("block 3 add pyramid all")
cubit.cmd("block 3 name 'pyramid_region'")

# Wedge mesh (sweep)
cubit.cmd("volume 3 scheme sweep")
cubit.cmd("mesh volume 3")
cubit.cmd("block 4 add wedge all")
cubit.cmd("block 4 name 'wedge_region'")

# Export all
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, "geometry.step")
```

## Function Signature

```python
def extract_curved_mesh(order: int = 2, geometry_file: str = None) -> netgen.meshing.Mesh:
    """Export Cubit mesh to Netgen mesh format.

    Args:
        cubit: Cubit Python interface object
        geometry_file: Path to geometry file (.step, .stp, .brep, .iges) for
                       mesh.Curve() support. If None, mesh is created without
                       geometry reference.

    Returns:
        netgen.meshing.Mesh: Netgen mesh object ready for use with NGSolve
    """
```

## Supported Elements

### 3D Volume Elements

| Element Type | Cubit | Netgen | Notes |
|--------------|-------|--------|-------|
| Tetrahedron | TET4 | TET | 4-node |
| Hexahedron | HEX8 | HEX | 8-node |
| Wedge/Prism | WEDGE6 | PRISM | 6-node |
| Pyramid | PYRAMID5 | PYRAMID | 5-node |

### 2D Boundary Elements

| Element Type | Cubit | Netgen | Notes |
|--------------|-------|--------|-------|
| Triangle | TRI3 | TRIG | 3-node |
| Quadrilateral | QUAD4 | QUAD | 4-node |

### 1D Boundary Elements

| Element Type | Cubit | Netgen | Notes |
|--------------|-------|--------|-------|
| Edge | EDGE2 | SEGMENT | 2-node |

## High-Order Curving with mesh.Curve()

The key advantage of `extract_curved_mesh()` is the ability to use Netgen's `mesh.Curve(order)` method for high-order geometry approximation.

### How It Works

1. **Cubit generates a 1st-order mesh** - Linear elements with nodes only at vertices
2. **Geometry is attached** - STEP file provides exact surface definitions
3. **mesh.Curve(order) is called** - Netgen adds high-order nodes on boundaries that conform to the geometry

**Important**: Only 1st-order elements are transferred from Cubit. High-order nodes are generated by Netgen based on the geometry, not by Cubit's 2nd-order meshing.

### Curve Orders

| Order | Element Nodes (Tet) | Description |
|-------|---------------------|-------------|
| 1 | 4 | Linear (original mesh) |
| 2 | 10 | Quadratic |
| 3 | 20 | Cubic |
| 4 | 35 | Quartic |
| 5 | 56 | Quintic |

### Example: Comparing Orders

```python
# Create mesh and attach geometry
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, "geometry.step")
mesh = ngsolve.Mesh(ngmesh)

# Apply different curve orders
for order in [2, 3, 4, 5]:
    mesh.Curve(order)
    print(f"Order {order}: mesh curved successfully")
```

## Node Ordering Conversion

### Vertex Ordering: GMSH/Cubit/Nastran/VTK vs Netgen

All four external formats (GMSH, Cubit, Nastran, VTK) use the same vertex ordering.
Netgen uses a different ordering for Hex, Prism, and Pyramid. The reordering
arrays are derived from `netgen/read_gmsh.py` and are self-inverse (the same
permutation converts in both directions).

| Element | GMSH/Cubit vertex order | Netgen vertex order | Reorder array |
|---------|------------------------|---------------------|---------------|
| Tet | [0,1,2,3] | [0,1,2,3] | identity |
| Hex | [0,1,2,3,4,5,6,7] | [0,1,5,4,3,2,6,7] | `{0,1,5,4,3,2,6,7}` |
| Prism | [0,1,2,3,4,5] | [0,2,1,3,5,4] | `{0,2,1,3,5,4}` |
| Pyramid | [0,1,2,3,4] | [3,2,1,0,4] | `{3,2,1,0,4}` |
| Tri | [0,1,2] | [0,1,2] | identity |
| Quad | [0,1,2,3] | [0,1,2,3] | identity |

Reference coordinates for Hex vertices:

```
GMSH/Cubit:  0=(0,0,0) 1=(1,0,0) 2=(1,1,0) 3=(0,1,0) 4=(0,0,1) 5=(1,0,1) 6=(1,1,1) 7=(0,1,1)
Netgen:      0=(0,0,0) 1=(1,0,0) 2=(1,0,1) 3=(0,0,1) 4=(0,1,0) 5=(1,1,0) 6=(1,1,1) 7=(0,1,1)
```

### High-Order Mid-Edge Node Ordering (GMSH = Nastran = VTK)

All three output formats use identical mid-edge node ordering.
`EdgeTables` in `MeshData.hpp` encodes this shared convention.

**TET10** (6 edges):

| HO node | Edge (vertex pair) | EdgeTables index |
|---------|-------------------|-----------------|
| 4 | {0,1} | 0 |
| 5 | {1,2} | 1 |
| 6 | {0,2} | 2 |
| 7 | {0,3} | 3 |
| 8 | {1,3} | 4 |
| 9 | {2,3} | 5 |

**HEX20** (12 edges):

| HO node | Edge (vertex pair) | EdgeTables index |
|---------|-------------------|-----------------|
| 8 | {0,1} | 0 |
| 9 | {1,2} | 1 |
| 10 | {2,3} | 2 |
| 11 | {3,0} | 3 |
| 12 | {4,5} | 4 |
| 13 | {5,6} | 5 |
| 14 | {6,7} | 6 |
| 15 | {7,4} | 7 |
| 16 | {0,4} | 8 |
| 17 | {1,5} | 9 |
| 18 | {2,6} | 10 |
| 19 | {3,7} | 11 |

**PRISM15/WEDGE15** (9 edges):

| HO node | Edge (vertex pair) | EdgeTables index |
|---------|-------------------|-----------------|
| 6 | {0,1} | 0 |
| 7 | {1,2} | 1 |
| 8 | {2,0} | 2 |
| 9 | {3,4} | 3 |
| 10 | {4,5} | 4 |
| 11 | {5,3} | 5 |
| 12 | {0,3} | 6 |
| 13 | {1,4} | 7 |
| 14 | {2,5} | 8 |

**PYRAMID13** (8 edges):

| HO node | Edge (vertex pair) | EdgeTables index |
|---------|-------------------|-----------------|
| 5 | {0,1} | 0 |
| 6 | {1,2} | 1 |
| 7 | {2,3} | 2 |
| 8 | {3,0} | 3 |
| 9 | {0,4} | 4 |
| 10 | {1,4} | 5 |
| 11 | {2,4} | 6 |
| 12 | {3,4} | 7 |

**TRI6** (3 edges): {0,1}, {1,2}, {2,0}

**QUAD8** (4 edges): {0,1}, {1,2}, {2,3}, {3,0}

### Netgen Internal Edge Ordering (different from GMSH)

Netgen's internal edge ordering differs from the GMSH/Nastran/VTK convention.
This is derived from the HO node reorder arrays in `netgen/read_gmsh.py`:

| Element | GMSH HO reorder (import) | Netgen edge order |
|---------|--------------------------|-------------------|
| TET10 | `[4,6,7,5,9,8]` | {0,1},{0,2},{0,3},{1,2},{2,3},{1,3} |
| HEX20 | `[8,16,10,12,13,19,15,14,9,11,18,17]` | (complex, see read_gmsh.py) |
| PRISM15 | `[7,6,9,8,11,10,13,12,14]` | (differs from GMSH) |
| QUAD8 | `[4,6,7,5]` | {0,1},{2,3},{3,0},{1,2} |
| TRI6 | `[4,5,3]` | {0,1},{1,2},{0,2} |

The Netgen internal ordering does NOT affect our export because `edge_ho_nodes_`
uses Cubit node IDs (physical identifiers) as keys, not local vertex indices.
`build_ho_conn_nc` looks up edges by Cubit node ID pairs, so the lookup is
format-independent and always finds the correct edge regardless of which
vertex ordering was used during HO node generation.

### Implementation: Vertex Reorder in NetgenCurver

`build_netgen_mesh()` applies GMSH-to-Netgen vertex reordering when adding
volume elements to the Netgen mesh (required for correct `BuildCurvedElements`
and `.vol` export):

```cpp
static const int hex_gmsh_to_ng[8]    = {0, 1, 5, 4, 3, 2, 6, 7};
static const int prism_gmsh_to_ng[6]  = {0, 2, 1, 3, 5, 4};
static const int pyr_gmsh_to_ng[5]    = {3, 2, 1, 0, 4};

// Hex example:
for (int ng_k = 0; ng_k < 8; ng_k++)
    el[ng_k] = cubit_nid_to_ng_pi[elem.conn[hex_gmsh_to_ng[ng_k]]];
```

No reordering is needed for Tet (identity), Tri (identity), or Quad (identity).

## Block to Material/Boundary Mapping

- **3D element blocks** → NGSolve Materials
- **2D element blocks** → NGSolve Boundaries
- **Block names** are preserved

```python
# In Cubit
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'steel'")
cubit.cmd("block 2 add tri all in surface 1")
cubit.cmd("block 2 name 'dirichlet_bc'")

# In NGSolve
mesh.GetMaterials()   # ('steel',)
mesh.GetBoundaries()  # ('dirichlet_bc',)
```

## Examples

Example scripts are available in the `examples/radia_cubit_mesh/netgen/` folder:

| File | Description |
|------|-------------|
| `occ_cubit_workflow.py` | **Recommended**: OCC → Cubit → Netgen workflow for curved surfaces |
| `netgen_sphere_example.py` | Basic sphere mesh with high-order curving |
| `netgen_mixed_elements_example.py` | Mixed element types (Hex, Tet, Wedge, Pyramid) |
| `netgen_poisson_example.py` | Complete FEM workflow: Poisson equation on sphere |
| `netgen_high_order_convergence.py` | Convergence study with different curve orders |
| `netgen_heat_conduction.py` | Heat conduction analysis |
| `netgen_linear_elasticity.py` | Linear elasticity analysis |
| `netgen_eigenvalue.py` | Eigenvalue problem |

### Running Examples

```bash
python examples/radia_cubit_mesh/netgen/occ_cubit_workflow.py
python examples/radia_cubit_mesh/netgen/netgen_sphere_example.py
python examples/radia_cubit_mesh/netgen/netgen_poisson_example.py
```

## Complete Example: FEM Analysis

```python
import os, sys
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
import radia_cubit_mesh
from ngsolve import *

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Create geometry and mesh
cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Define blocks
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'domain'")
cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'dirichlet'")

# Export geometry
cubit.cmd('export step "sphere.step" overwrite')

# Convert to Netgen mesh, then wrap with NGSolve
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, "sphere.step")
mesh = Mesh(ngmesh)
mesh.Curve(3)

# Solve Poisson equation: -Laplacian(u) = 1
fes = H1(mesh, order=3, dirichlet="dirichlet")
u, v = fes.TnT()

a = BilinearForm(fes)
a += grad(u)*grad(v)*dx
a.Assemble()

f = LinearForm(fes)
f += 1*v*dx
f.Assemble()

gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

# Visualize
from ngsolve.webgui import Draw
Draw(gfu)
```

## Recommended Workflow: OCC → Cubit → Netgen

### The Seam Line Problem

When Cubit creates cylindrical or curved geometry and exports to STEP, OCC (OpenCASCADE) splits these surfaces at seam lines (typically y=0). This causes `mesh.Curve(2+)` to fail or produce artifacts on elements crossing the seam.

| Workflow | Seam Split | Curve(2+) |
|----------|------------|-----------|
| Cubit → STEP → OCC | ❌ Split | ❌ Fails |
| **OCC → STEP → Cubit → OCC** | ✅ No split | ✅ Works |

### Solution: Create Geometry in OCC First

The recommended workflow for curved surfaces:

1. **Create geometry in OCC** (no seam splitting)
2. **Export to STEP**
3. **Import into Cubit** for meshing (ACIS preserves OCC topology)
4. **Export mesh with original OCC geometry reference**

```python
import os, sys
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

from netgen import occ
from netgen.occ import OCCGeometry
import cubit
import radia_cubit_mesh
from ngsolve import Mesh

# Step 1: Create geometry in OCC
cyl = occ.Cylinder(occ.Pnt(0, 0, -1), occ.Vec(0, 0, 1), r=0.5, h=2)
geo = OCCGeometry(cyl)
print(f"OCC geometry: {len(geo.shape.faces)} faces")  # 3 faces (no seam)

# Step 2: Export to STEP
cyl.WriteStep("geometry.step")

# Step 3: Import to Cubit and mesh
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd('import step "geometry.step" heal')

# IMPORTANT: Mesh surfaces first for accurate node positions
cubit.cmd("surface all scheme trimesh")
cubit.cmd("surface all size 0.15")
cubit.cmd("mesh surface all")

cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Step 4: Export with OCC geometry reference
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, geometry=geo)

# Now Curve() works!
mesh = Mesh(ngmesh)
mesh.Curve(3)  # OK!
```

### When to Use This Workflow

Use the OCC → Cubit → Netgen workflow when:
- Your geometry contains **cylindrical, conical, or toroidal surfaces**
- You need **high-order curving** (`Curve(2)` or higher)
- Standard STEP export from Cubit causes `Curve()` failures

For **planar or simple geometries**, the standard workflow (Cubit → STEP → Netgen) works fine.

### Example Script

See `examples/radia_cubit_mesh/netgen/occ_cubit_workflow.py` for a complete working example.

## Alternative: SetDeformation Approach

When `mesh.Curve()` doesn't work correctly (due to missing `geominfo` UV parameters), you can use NGSolve's `SetDeformation` to manually project boundary elements onto the exact geometry surface.

### Why SetDeformation?

The `mesh.Curve()` method requires `Element2D.geominfo` (UV parameters) to project high-order nodes onto geometry. When importing meshes from Cubit, these parameters are not set, causing `Curve(2+)` to produce incorrect results.

`SetDeformation` bypasses this limitation by:
1. Defining a deformation field that projects points to the exact surface
2. Applying the deformation to boundary DOFs only
3. Using NGSolve's standard H1 interpolation for smooth blending

### Cylinder Example

```python
from ngsolve import Mesh, GridFunction, VectorH1, Integrate, BND, CF
from ngsolve import x, y, z, sqrt, IfPos

# After creating mesh from Cubit...
mesh = Mesh(ngmesh)
mesh.Curve(1)  # Linear mesh only

# Define deformation for cylinder (radius R, z-axis aligned)
R = 0.5
fes = VectorH1(mesh, order=4)  # Order 3-4 recommended
deform = GridFunction(fes)

# Project (x,y,z) to (R*x/r, R*y/r, z) where r = sqrt(x²+y²)
r_xy = sqrt(x*x + y*y)
scale = IfPos(r_xy - 0.01, R/r_xy - 1, 0)  # Avoid division by zero
deform_cf = CF((scale * x, scale * y, 0))

# Apply to cylinder surface only
deform.Set(deform_cf, definedon=mesh.Boundaries('face_0'))
mesh.SetDeformation(deform)

# Now integration is accurate!
area = Integrate(CF(1), mesh, VOL_or_BND=BND)  # Exact!
```

### Accuracy by Deformation Order

| Order | Area Error | Volume Error |
|-------|------------|--------------|
| 1 | 0.07% | 0.004% |
| 2 | 0.001% | 0.002% |
| 3 | 0.0001% | 0.00001% |
| 4+ | ~0% | ~0% |

### When to Use SetDeformation

- When `mesh.Curve()` produces incorrect results (inflated areas/volumes)
- When you know the exact mathematical definition of your geometry
- For cylindrical, spherical, or other analytically defined surfaces
- When you need high accuracy without modifying Netgen source code

### Example Script

See `examples/radia_cubit_mesh/netgen/setdeformation_curving.py` for complete examples including:
- `apply_cylinder_deformation()` - For cylindrical surfaces
- `apply_sphere_deformation()` - For spherical surfaces

## Troubleshooting

### mesh.Curve() Fails

**Cause 1**: Geometry file not provided or not matching mesh.

**Solution**: Ensure the STEP file is exported from the same geometry:
```python
cubit.cmd('export step "geometry.step" overwrite')
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, "geometry.step")
```

**Cause 2**: Seam line problem on curved surfaces (cylinder, cone, etc.).

**Solution**: Use the OCC → Cubit → Netgen workflow described above.


### Empty Mesh (ne=0)

**Cause**: Elements not added to blocks.

**Solution**: Use `block X add <element_type> all in volume Y`:
```python
# Wrong
cubit.cmd("block 1 volume 1")

# Correct
cubit.cmd("block 1 add tet all in volume 1")
```

### Missing Boundaries

**Cause**: Surface elements not in separate blocks.

**Solution**: Create blocks for boundary elements:
```python
cubit.cmd("block 2 add tri all in surface 1")
cubit.cmd("block 2 name 'inlet'")
cubit.cmd("block 3 add tri all in surface 2")
cubit.cmd("block 3 name 'outlet'")
```

### Import Errors

**Cause**: NGSolve/Netgen not installed or not in path.

**Solution**: Install NGSolve and ensure it's accessible:
```bash
pip install ngsolve
```

### Geometry Not Found

**Cause**: STEP file path incorrect or file doesn't exist.

**Solution**: Use absolute path or ensure file exists:
```python
import os
step_file = os.path.abspath("geometry.step")
cubit.cmd(f'export step "{step_file}" overwrite')
ngmesh = radia_cubit_mesh.extract_curved_mesh(cubit, step_file)
```

## Netgen vs NGSolve Mesh

| | netgen.meshing.Mesh | ngsolve.Mesh |
|---|---|---|
| Role | Core mesh data structure | Wrapper/view for FEM |
| Usage | Mesh manipulation, generation | Finite element analysis |
| Returned by | `extract_curved_mesh()` | `ngsolve.Mesh(ngmesh)` |

## Requirements

- Coreform Cubit 2025.3 or later
- NGSolve/Netgen (with OCC support for geometry)
- Python 3.8+

## See Also

- [NGSolve Documentation](https://docu.ngsolve.org/)
- [Netgen Mesh Generation](https://docu.ngsolve.org/latest/netgen_tutorials/)
- [Cubit Documentation](https://cubit.sandia.gov/)
