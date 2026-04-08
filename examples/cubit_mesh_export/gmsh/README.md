# Gmsh Export Examples

Export Cubit mesh to Gmsh format (.msh).

## v2.2 vs v4.1 Format

| Feature | v2.2 (`radia_export gmsh`) | v4.1 (`radia_export gmsh version 4`) |
|---------|-------------------------|-------------------------|
| Format version | 2.2 | 4.1 |
| Structure | Flat lists | Entity blocks |
| $Entities section | No | Yes |
| DIM parameter | No | Yes |
| 2D normal control | No | Yes |
| File size | Smaller | Larger |
| **NGSolve/Netgen** | **Supported** | Not recommended |

### When to Use Which

**Use v2.2 when:**
- GMSH visualization of exported meshes
- Maximum compatibility with older software
- Smaller file size needed

**Use v4.1 when:**
- Need geometry topology ($Entities section)
- 2D meshes requiring normal orientation control
- FEniCS or other modern solvers

**For NGSolve computation**, use `radia_export netgen "mesh.vol"` instead of .msh. The .vol format supports arbitrary-order curving (1-5) and preserves material/boundary labels.

## Usage

```python
import cubit

cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Gmsh v2.2 (recommended for NGSolve)
cubit.cmd('radia_export gmsh "mesh.msh" overwrite')

# Gmsh v4.1 (with $Entities section)
cubit.cmd('radia_export gmsh "mesh.msh" version 4 overwrite')
```

## NGSolve Computation

**For NGSolve FEM computation, use `radia_export netgen` (.vol) instead of .msh.**

The `.vol` format supports arbitrary-order curving (order 1-5) and preserves material/boundary labels natively. The `.msh` export is maintained for GMSH visualization only.

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
cubit.cmd('radia_export netgen "mesh.vol" order 3 overwrite')
cubit.cmd('radia_export gmsh "mesh.msh" order 2 version 2 overwrite')  # visualization only
mesh = Mesh("mesh.vol")
```

See [examples/netgen/README.md](../netgen/README.md) for detailed NGSolve workflow.

## v4.1 DIM Parameter

The v4.1 format supports dimension control:

```python
# Auto-detect dimension (default)
cubit.cmd('radia_export gmsh "mesh.msh" version 4 overwrite')

# Force 2D mode (normals to +z, z=0)
cubit.cmd('radia_export gmsh "plate.msh" version 4 dimension 2 overwrite')

# Force 3D mode
cubit.cmd('radia_export gmsh "solid.msh" version 4 dimension 3 overwrite')
```

## Sample Files

| File | Description |
|------|-------------|
| `cube_v2.msh` | Gmsh v2.2 format (1st order) |
| `cube_v4.msh` | Gmsh v4.1 format (1st order) |
| `cube_2nd_order.msh` | v2.2 with 2nd order elements (TET10, TRI6) |

## Regenerate Samples

```bash
"${CUBIT_PATH:-C:/Program Files/Coreform Cubit 2025.3/bin}/python3/python.exe" gmsh_export_example.py
```

## See Also

- [docs/cubit/export_Gmsh_ver2.md](../../../docs/cubit/export_Gmsh_ver2.md) - v2.2 documentation
- [docs/cubit/export_Gmsh_ver4.md](../../../docs/cubit/export_Gmsh_ver4.md) - v4.1 documentation
- [examples/netgen/README.md](../netgen/README.md) - NGSolve integration
