# Gmsh Export Examples

Export Cubit mesh to Gmsh format (.msh).

## Why GMSH v4.1

The Radia Cubit plugin emits GMSH v4.1 only. v4.1 is the modern, entity-block
format used by current solvers:

- Includes the `$Entities` section (geometry topology: points, curves,
  surfaces, volumes)
- Groups nodes and elements by entity block, which preserves the
  geometry-to-mesh association
- Supports explicit DIM control for 2D meshes (normal orientation)
- Is the default format for FEniCS and other modern FEM/BEM tools

**For NGSolve computation**, use `export netgen "mesh.vol"` instead of
.msh. The .vol format supports arbitrary-order curving (1-5) and preserves
material/boundary labels natively. The `.msh` export is maintained for GMSH
visualization only.

## Usage

```python
import cubit

cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Gmsh v4.1 export
cubit.cmd('export gmsh "mesh.msh" overwrite')
```

## NGSolve Computation

**For NGSolve FEM computation, use `export netgen` (.vol) instead of .msh.**

```python
import cubit
from ngsolve import Mesh

# 1. Create and mesh in Cubit
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# 2. Register blocks
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')

# 3. Export .vol for computation, .msh for GMSH visualization
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
cubit.cmd('export gmsh "mesh.msh" order 2 overwrite')  # visualization only
mesh = Mesh("mesh.vol")
```

See [the Netgen/NGSolve workflow](../netgen/README.md) for details.

## DIM Parameter

The v4.1 format supports dimension control:

```python
# Auto-detect dimension (default)
cubit.cmd('export gmsh "mesh.msh" overwrite')

# Force 2D mode (normals to +z, z=0)
cubit.cmd('export gmsh "plate.msh" dimension 2 overwrite')

# Force 3D mode
cubit.cmd('export gmsh "solid.msh" dimension 3 overwrite')
```

## Sample Files

| File | Description |
|------|-------------|
| `cube.msh` | Gmsh v4.1 format (1st order) |
| `cube_2nd_order.msh` | v4.1 with 2nd order elements (TET10, TRI6) |

## Regenerate Samples

```bash
"${CUBIT_PATH:-<Coreform Cubit 2025.8+>/bin}/python3/python.exe" gmsh_export_example.py
```

## See Also

- [docs/cubit/export_Gmsh.md](../../../docs/cubit/export_Gmsh.md) - Gmsh export reference
- [Netgen/NGSolve workflow](../netgen/README.md) - NGSolve integration
