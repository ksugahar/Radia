# MEG Export Examples

Export Cubit mesh to ELF/MAGIC MEG format (.meg).

## DIM Parameter

Controls the analysis dimension:

| DIM | Description | Coordinate System |
|-----|-------------|-------------------|
| `'T'` | 3D (Three-dimensional) | X, Y, Z |
| `'K'` | 2D Planar (Kartesian) | X, Y (Z=0) |
| `'R'` | Axisymmetric (Rotational) | R (X), Z (Y=0 plane) |

```python
# 3D mesh
cubit.cmd('radia export meg "mesh.meg" overwrite')

# 2D planar mesh (XY plane)
cubit.cmd('radia export meg "plate.meg" overwrite')

# Axisymmetric mesh (R-Z plane, rotated around Z axis)
cubit.cmd('radia export meg "axisym.meg" overwrite')
```

## Block Names = ELF Element Names

**Important**: Block names must match ELF magnetic material element naming convention.

### ELF Element Naming Convention

Format: `MMB<nodes><dim>` where:
- `MMB` = Magnetic Material Block
- `<nodes>` = Number of nodes (3, 4, 6, 8)
- `<dim>` = Dimension (T, K, R)

| DIM | Triangle | Quad | Tet | Wedge | Hex |
|-----|----------|------|-----|-------|-----|
| `'T'` (3D) | - | - | **MMB4T** | **MMB6T** | **MMB8T** |
| `'K'` (2D) | **MMB3K** | **MMB4K** | - | - | - |
| `'R'` (Axisym) | **MMB3R** | **MMB4R** | - | - | - |

### Usage Example

```python
# 3D tetrahedral mesh
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'MMB4T'")  # Must match ELF naming
cubit.cmd('radia export meg "mesh.meg" overwrite')

# 2D planar triangular mesh
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3K'")  # Must match ELF naming
cubit.cmd('radia export meg "plate.meg" overwrite')

# Axisymmetric triangular mesh
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3R'")  # Must match ELF naming
cubit.cmd('radia export meg "axisym.meg" overwrite')
```

## MGR2 Parameter (Spatial Nodes)

Optional boundary condition nodes outside the mesh domain:

```python
spatial_nodes = [
    [0.0, 0.0, 10.0],   # Far-field point 1
    [0.0, 0.0, -10.0],  # Far-field point 2
]
cubit.cmd('radia export meg "mesh.meg" overwrite')
```

## Limitations

- **1st order elements only**: 2nd order elements are not supported.

## Sample Files

| File | Description |
|------|-------------|
| `cube_tet.meg` | 3D tetrahedral mesh (DIM='T', MMB4T) |
| `sphere_hex.meg` | 3D hexahedral mesh (DIM='T', MMB8T) |
| `plate_2d.meg` | 2D planar mesh (DIM='K', MMB3K) |
| `axisym.meg` | Axisymmetric mesh (DIM='R', MMB3R) |
| `with_spatial_nodes.meg` | 3D mesh with MGR2 spatial nodes |

## Regenerate Samples

```bash
"${CUBIT_PATH:-C:/Program Files/Coreform Cubit 2025.3/bin}/python3/python.exe" meg_export_example.py
```

## See Also

- [docs/cubit/export_meg.md](../../../docs/cubit/export_meg.md) - Full documentation
