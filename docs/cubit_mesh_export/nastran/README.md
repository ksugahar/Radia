# Nastran BDF Export Examples

Export Cubit mesh to Nastran Bulk Data Format (.bdf).

## DIM Parameter

Controls the mesh dimension and element types:

| DIM | Description | Elements |
|-----|-------------|----------|
| `"3D"` | Volume elements | CTETRA, CHEXA, CPENTA, CPYRAM |
| `"2D"` | Surface elements | CTRIA3/CTRIA6, CQUAD4/CQUAD8; original coordinates preserved |

```python
# 3D solid mesh
cubit.cmd('export nastran_bdf "solid.bdf" dimension 3 overwrite')

# 2D shell mesh; all three GRID coordinates are preserved
cubit.cmd('export nastran_bdf "plate.bdf" dimension 2 overwrite')
```

## PYRAM Parameter

Controls pyramid element handling in hybrid hex/tet meshes:

| PYRAM | Output | Use Case |
|-------|--------|----------|
| `True` | CPYRAM (5-node pyramid) | Standard Nastran solvers |
| `False` | Degenerate CHEXA (8-node hex with duplicate nodes) | importers without CPYRAM support |

**Background**: When hex and tet regions meet, pyramid elements bridge the interface. Some solvers cannot import CPYRAM and interpret them as degenerate CHEXA.

```python
# Standard export with CPYRAM elements
cubit.cmd('export nastran_bdf "mesh.bdf" overwrite')

# For importers without CPYRAM support: convert pyramids to degenerate hex
cubit.cmd('export nastran_bdf "mesh_no_pyramid.bdf" nopyramid overwrite')
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

*Pyramid exported as CPYRAM when `PYRAM=True`, degenerate CHEXA when `PYRAM=False`.

## Contract boundary

- Orders 1 and 2 are supported for all element families in the table.
- Blocks produce PSOLID/PSHELL, sidesets produce collision-free PSHELL, and
  nodesets produce SET1 cards.
- MAT cards are intentionally omitted because the exporter cannot infer real
  constitutive data. Assign materials in the receiving application.
- `export jmag_nastran` was retired. Journals still using it must be updated to `export nastran_bdf`; the wheel gate fails if the retired token reappears in the compiled `.ccm`.

## Sample Files

| File | Description |
|------|-------------|
| `cube_3d.bdf` | 3D tetrahedral mesh |
| `plate_2d.bdf` | 2D triangular mesh |
| `mixed_with_pyramid.bdf` | Mixed hex/tet with CPYRAM |
| `mixed_no_pyramid.bdf` | Mixed hex/tet with degenerate CHEXA |

## Regenerate Samples

```bash
"${CUBIT_PATH:-<Coreform Cubit 2025.8+>/bin}/python3/python.exe" nastran_export_example.py
```

## See Also

- [docs/cubit/export_Nastran.md](../../../docs/cubit/export_Nastran.md) - Full documentation
