# export_NetgenMesh - Cubit to Netgen Mesh Export

## Overview

`export_NetgenMesh()` creates a `netgen.meshing.Mesh` object directly from Cubit mesh data. When combined with a geometry file, it enables high-order curved mesh generation using `mesh.Curve(order)` method.

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

| Feature | export_NetgenMesh() | export_Gmsh_ver2() |
|---------|-----------------|-------------------|
| Max order | Unlimited (via Curve) | 2nd order |
| Intermediate file | None | .msh file |
| Geometry reference | STEP/BREP/IGES | None |
| High-order accuracy | Exact (CAD-based) | Mesh-based |
| Wedge elements | Supported | Supported |
| Pyramid elements | Supported | Supported |

## Recommended: Gmsh 2nd Order Workflow

For most use cases, **Gmsh 2nd order elements with `ReadGmsh()`** provides the best accuracy with minimal complexity.

### Why Gmsh 2nd Order?

| Element Type | 1st Order Error | 2nd Order Error |
|--------------|-----------------|-----------------|
| Tet10        | ~1.4%           | **0.001%**      |
| Hex20        | ~1.6%           | **0.002%**      |
| Wedge15      | ~1.6%           | ~1.1% (see note)|

**Note**: Wedge15 mid-edge nodes are not projected to curved surfaces by Cubit. Use SetDeformation for Wedge elements on curved surfaces.

### Basic Workflow

```python
import cubit
import cubit_mesh_export
from netgen.read_gmsh import ReadGmsh
from ngsolve import Mesh

# 1. Create and mesh geometry
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")

# 2. Register blocks with 2nd order elements
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")

# 3. Export and load
cubit_mesh_export.export_Gmsh_ver2(cubit, "mesh.msh")
mesh = Mesh(ReadGmsh("mesh.msh"))
```

### Supported 2nd Order Element Types

| Cubit Type | Gmsh Type | Nodes |
|------------|-----------|-------|
| tetra10    | Tet10     | 10    |
| hex20      | Hex20     | 20    |
| wedge15    | Wedge15   | 15    |
| pyramid13  | Pyramid13 | 13    |
| tri6       | Tri6      | 6     |
| quad8      | Quad8     | 8     |

## Usage

### Basic Usage

```python
import cubit
import cubit_mesh_export
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
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, geometry_file="geometry.step")

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
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit)
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
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, "geometry.step")
```

## Function Signature

```python
def export_NetgenMesh(cubit, geometry_file: str = None) -> netgen.meshing.Mesh:
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

The key advantage of `export_NetgenMesh()` is the ability to use Netgen's `mesh.Curve(order)` method for high-order geometry approximation.

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
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, "geometry.step")
mesh = ngsolve.Mesh(ngmesh)

# Apply different curve orders
for order in [2, 3, 4, 5]:
    mesh.Curve(order)
    print(f"Order {order}: mesh curved successfully")
```

## Node Ordering Conversion

`export_NetgenMesh()` automatically handles node ordering differences between Cubit and Netgen:

| Element | Cubit Order | Netgen Order |
|---------|-------------|--------------|
| Tet | [0,1,2,3] | [0,1,2,3] |
| Hex | [0,1,2,3,4,5,6,7] | [0,1,5,4,3,2,6,7] |
| Wedge | [0,1,2,3,4,5] | [0,2,1,3,5,4] |
| Pyramid | [0,1,2,3,4] | [3,2,1,0,4] |

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

Example scripts are available in the `examples/cubit/netgen/` folder:

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
python examples/cubit/netgen/occ_cubit_workflow.py
python examples/cubit/netgen/netgen_sphere_example.py
python examples/cubit/netgen/netgen_poisson_example.py
```

## Complete Example: FEM Analysis

```python
import os, sys
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
import cubit_mesh_export
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
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, "sphere.step")
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
import cubit_mesh_export
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
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)

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

See `examples/cubit/netgen/occ_cubit_workflow.py` for a complete working example.

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

See `examples/cubit/netgen/setdeformation_curving.py` for complete examples including:
- `apply_cylinder_deformation()` - For cylindrical surfaces
- `apply_sphere_deformation()` - For spherical surfaces

## Troubleshooting

### mesh.Curve() Fails

**Cause 1**: Geometry file not provided or not matching mesh.

**Solution**: Ensure the STEP file is exported from the same geometry:
```python
cubit.cmd('export step "geometry.step" overwrite')
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, "geometry.step")
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
ngmesh = cubit_mesh_export.export_NetgenMesh(cubit, step_file)
```

## Netgen vs NGSolve Mesh

| | netgen.meshing.Mesh | ngsolve.Mesh |
|---|---|---|
| Role | Core mesh data structure | Wrapper/view for FEM |
| Usage | Mesh manipulation, generation | Finite element analysis |
| Returned by | `export_NetgenMesh()` | `ngsolve.Mesh(ngmesh)` |

## Requirements

- Coreform Cubit 2025.3 or later
- NGSolve/Netgen (with OCC support for geometry)
- Python 3.8+

## See Also

- [NGSolve Documentation](https://docu.ngsolve.org/)
- [Netgen Mesh Generation](https://docu.ngsolve.org/latest/netgen_tutorials/)
- [Cubit Documentation](https://cubit.sandia.gov/)
