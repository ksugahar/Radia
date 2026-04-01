# Gmsh Export Examples

Export Cubit mesh to Gmsh format (.msh).

## v2.2 vs v4.1 Format

| Feature | v2.2 (`radia export gmsh`) | v4.1 (`radia export gmsh version 4`) |
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
- **NGSolve/Netgen integration** (ReadGmsh supports v2.2)
- Maximum compatibility with older software
- Simple mesh transfer without geometry information
- Smaller file size needed

**Use v4.1 when:**
- Need geometry topology ($Entities section)
- 2D meshes requiring normal orientation control
- FEniCS or other modern solvers

## Usage

```python
import cubit

cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# Gmsh v2.2 (recommended for NGSolve)
cubit.cmd('radia export gmsh "mesh.msh" overwrite')

# Gmsh v4.1 (with $Entities section)
cubit.cmd('radia export gmsh "mesh.msh" version 4 overwrite')
```

## NGSolve Integration

**NGSolve's `ReadGmsh()` supports Gmsh v2.2 format.**

This is the recommended workflow for Cubit to NGSolve:

```python
import cubit
from netgen.read_gmsh import ReadGmsh
from ngsolve import Mesh

# 1. Create and mesh in Cubit
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# 2. Register blocks with 2nd order elements
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")

# 3. Export to Gmsh v2.2 and load into NGSolve
cubit.cmd('radia export gmsh "mesh.msh" overwrite')
mesh = Mesh(ReadGmsh("mesh.msh"))
```

### Accuracy with 2nd Order Elements

| Element Type | 1st Order | 2nd Order |
|--------------|-----------|-----------|
| Tet | ~1.4% | **0.001%** |
| Hex | ~1.6% | **0.002%** |

See [examples/netgen/README.md](../netgen/README.md) for detailed NGSolve workflow.

## v4.1 DIM Parameter

The v4.1 format supports dimension control:

```python
# Auto-detect dimension (default)
cubit.cmd('radia export gmsh "mesh.msh" version 4 overwrite')

# Force 2D mode (normals to +z, z=0)
cubit.cmd('radia export gmsh "plate.msh" version 4 dimension 2 overwrite')

# Force 3D mode
cubit.cmd('radia export gmsh "solid.msh" version 4 dimension 3 overwrite')
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
