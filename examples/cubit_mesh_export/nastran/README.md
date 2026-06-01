# Nastran BDF Export Examples

Export Cubit mesh to Nastran Bulk Data Format (.bdf).

## DIM Parameter

Controls the mesh dimension and element types:

| DIM | Description | Elements |
|-----|-------------|----------|
| `"3D"` | Volume elements | CTETRA, CHEXA, CPENTA, CPYRAM |
| `"2D"` | Surface elements | CTRIA3, CQUAD4 (normals oriented to +z) |

```python
# 3D solid mesh
cubit.cmd('export jmag_nastran "solid.bdf" dimension 3 overwrite')

# 2D shell mesh (normals reoriented to +z, z-coordinates set to 0)
cubit.cmd('export jmag_nastran "plate.bdf" dimension 2 overwrite')
```

## PYRAM Parameter

Controls pyramid element handling in hybrid hex/tet meshes:

| PYRAM | Output | Use Case |
|-------|--------|----------|
| `True` | CPYRAM (5-node pyramid) | Standard Nastran solvers |
| `False` | Degenerate CHEXA (8-node hex with duplicate nodes) | JMAG compatibility |

**Background**: When hex and tet regions meet, pyramid elements bridge the interface. Some solvers (e.g., JMAG) cannot import CPYRAM and interpret them as degenerate CHEXA.

```python
# Standard export with CPYRAM elements
cubit.cmd('export jmag_nastran "mesh.bdf" overwrite')

# For JMAG: convert pyramids to degenerate hex
cubit.cmd('export jmag_nastran "mesh_jmag.bdf" nopyramid overwrite')
```

## Element Mapping

| Cubit | Nastran (3D) | Nastran (2D) |
|-------|--------------|--------------|
| Tet | CTETRA | - |
| Hex | CHEXA | - |
| Wedge | CPENTA | - |
| Pyramid | CPYRAM / CHEXA* | - |
| Tri | - | CTRIA3 |
| Quad | - | CQUAD4 |
| Edge | CROD | CROD |
| Node | CMASS | CMASS |

*Pyramid exported as CPYRAM when `PYRAM=True`, degenerate CHEXA when `PYRAM=False`.

## Limitations

- **1st order elements only**: 2nd order elements (TETRA10, HEX20) are not supported.

## Sample Files

| File | Description |
|------|-------------|
| `cube_3d.bdf` | 3D tetrahedral mesh |
| `plate_2d.bdf` | 2D triangular mesh |
| `mixed_with_pyramid.bdf` | Mixed hex/tet with CPYRAM |
| `mixed_for_jmag.bdf` | Mixed hex/tet with degenerate CHEXA |

## Regenerate Samples

```bash
"${CUBIT_PATH:-<Coreform Cubit 2025.8+>/bin}/python3/python.exe" nastran_export_example.py
```

## See Also

- [docs/cubit/export_Nastran.md](../../../docs/cubit/export_Nastran.md) - Full documentation
