# export_Nastran

Export a blocked Coreform Cubit mesh as a Nastran Bulk Data Format (`.bdf`)
mesh-interchange deck.

## Command

```text
export nastran_bdf "mesh.bdf" [order <1|2>] [dimension <2|3>] [nopyramid] [overwrite]
```

`export nastran_bdf` is the primary, solver-neutral command.
`export jmag_nastran` is a deprecated compatibility alias for existing
journals. Do not use Cubit's built-in `export nastran` when this plugin format
and its high-order curving are required.

Mesh elements must belong to Cubit blocks. Sidesets and nodesets are exported
in addition to the blocked mesh.

## Options

| Option | Default | Description |
|---|---:|---|
| `order 1` | yes | Linear Nastran elements |
| `order 2` |  | Quadratic edge nodes generated through NetgenCurver |
| `dimension 3` | yes | Export compatible 3D and 2D blocked elements |
| `dimension 2` |  | Export only 2D blocked elements; reject a volume-only selection |
| `nopyramid` | off | Write each pyramid as a degenerate CHEXA for consumers that do not accept CPYRAM |
| `overwrite` | off | Replace an existing output file |

The exporter preserves all three GRID coordinates in both dimension modes. It
does not flatten an off-plane shell to `z=0` and does not silently reinterpret
volume elements as zero-thickness shells.

## Element cards

| Cubit element | Order 1 | Order 2 |
|---|---|---|
| Tetrahedron | CTETRA4 | CTETRA10 |
| Hexahedron | CHEXA8 | CHEXA20 |
| Wedge/prism | CPENTA6 | CPENTA15 |
| Pyramid | CPYRAM5 | CPYRAM13 |
| Triangle | CTRIA3 | CTRIA6 |
| Quadrilateral | CQUAD4 | CQUAD8 |

One-dimensional and zero-dimensional entities are not part of this export
contract.

## Blocks, sidesets, nodesets, and materials

- A 3D block receives a `PSOLID` property whose PID is the Cubit block ID.
- A 2D block receives a `PSHELL` property whose PID is the Cubit block ID.
- A sideset is written as boundary shell elements with a `PSHELL`. If its
  Cubit ID collides with a block PID, the exporter assigns a distinct PID and
  records it in the BDF comment.
- A nodeset is written as a machine-readable fixed-field `SET1` card as well as
  a descriptive comment.
- Property cards reference the corresponding Cubit group ID as a material ID,
  but the exporter deliberately writes no `MAT1`, `MAT10`, or other physical
  material card. Add the real constitutive data in JMAG, COMSOL, or another
  receiving application. The result is an interchange mesh, not a complete
  Nastran analysis deck.

## Examples

```python
import cubit

cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "solid"')
cubit.cmd('export nastran_bdf "solid.bdf" order 2 dimension 3 overwrite')
```

```python
cubit.cmd("block 10 add tri all")
cubit.cmd('block 10 name "plate"')
cubit.cmd('export nastran_bdf "plate.bdf" order 2 dimension 2 overwrite')
```

For a hybrid HEX/TET mesh whose receiver cannot import CPYRAM:

```python
cubit.cmd('export nastran_bdf "hybrid.bdf" dimension 3 nopyramid overwrite')
```

## Validation

The real-Cubit regression covers TET, HEX, WEDGE, PYRAMID, TRI, QUAD,
off-plane 2D coordinates, property-ID collisions, SET1 nodesets, the historical
alias, and rejection of volume-only `dimension 2` output:

```powershell
python validation_test\cubit\test_nastran_export.py
python validation_test\cubit\test_export_combinations.py
```

The second script creates an order/dimension corpus. Parse that corpus with an
independent reader and save a JSON report with:

```powershell
python -m pip install "pyNastran==1.4.1"
python validation_test\cubit\validate_nastran_with_pynastran.py `
  --output-json C:\temp\nastran-interchange.json
```

The independent gate checks GRID coordinates, element-to-node references,
element-to-property references, property kinds, and SET1 cards. Missing MAT
cards are reported as the expected boundary between a mesh-interchange deck and
a complete structural analysis deck.

## See also

- [export_Gmsh.md](export_Gmsh.md)
- [export_NetgenMesh.md](export_NetgenMesh.md)
- [Function_Reference.md](Function_Reference.md)
