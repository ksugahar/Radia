# export_meg

Export mesh to ELF/MAGIC MEG format for electromagnetic field analysis.

## Syntax

```
export meg "filename.meg" [threed|twod|axisymmetric] [labels "1:MMB,2:MWL,..."] [overwrite]
```

> **Note (2026-04-24)**: Earlier releases registered this command as
> bare `export meg` under the assumption that Cubit owned the `.meg`
> format. In practice Cubit has no built-in `meg` keyword, so the bare
> `export meg` was rejected with `Unrecognized Identifier: 'meg'`.
> The command is now `export meg`, consistent with the other plugin
> exporters (`export gmsh / netgen / vtk / femeem / meg`, plus
> `export jmag_nastran`). ELF block-name → prefix auto-detection (see "ELF Element
> Type Labels" below) is unchanged.

Block names define ELF physics element types (3-character prefix).
No explicit block assignment is required — all meshed elements are exported.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `threed` | yes | 3D analysis (DIM suffix = `T`) |
| `twod` | | 2D planar analysis (DIM = `K`, z forced to 0) |
| `axisymmetric` | | Axisymmetric analysis (DIM = `R`, y forced to 0) |
| `labels "..."` | | Per-block label override (format: `blockID:PREFIX,...`) |
| `overwrite` | off | Overwrite existing file |

## ELF Element Type Labels

Each element gets a 5-character type string: `PREFIX` + `NODECOUNT` + `DIM`.

Example: `MMB8T` = Magnetic Body, 8-node hex, 3D.

### Block Naming Convention

The plugin extracts the ELF prefix from the first 3 characters of each Cubit block name (case-insensitive):

| Prefix | Description |
|--------|-------------|
| `MMB` | Magnetic body (iron core, ferromagnetic) |
| `MMS` | Magnetic shell |
| `MMT` | Magnetic thin (line/surface elements) |
| `MMP` | Magnetic permanent magnet |
| `MWL` | Nonlinear magnet with fixed local axis |
| `MWV` | Nonlinear magnet with direction vectors |
| `MCO` | Current conductor |

If the block name does not match any valid prefix, `MMB` is used as default.
Blocks named `AIR` also map to `MMB` (air is still a magnetic region in ELF).

### Label Override

The `labels` option overrides block-name-based prefix detection:

```
export meg "mesh.meg" labels "1:MMB,2:MWL,3:MCO" overwrite
```

## Supported Elements

1st order only. Node counts are auto-detected:

| Cubit Element | ELF Nodes | Notes |
|---------------|-----------|-------|
| HEX8 | 8 | Standard hex |
| TET4 | 4 | Standard tet |
| WEDGE6 | 6 | Prism |
| PYRAMID5 | 8 | Exported as degenerate hex (apex repeated 4x) |
| QUAD4 | 4 | Surface element |
| TRI3 | 3 | Surface element |

### Pyramid Handling

Pyramids are converted to degenerate 8-node hexahedra by repeating the apex node:
`n0 n1 n2 n3 n4 n4 n4 n4`. This is standard ELF convention.

## File Format

```
BOOK  MEP  3.50
* ELF/MESH VERSION 7.3.0
* SOLVER = ELF/MAGIC
MGSC 0.001
* NODE
MGR1 1 0 0.0 0.0 0.0
MGR1 2 0 1.0 0.0 0.0
...
* ELEMENT K
MMB8T 1 0 1 1 2 3 4 5 6 7 8
MWL4T 2 0 2 9 10 11 12
...
* NODE
MGR2 1 0 0.5 0.5 0.5
...
BOOK  END
```

### MGR2 Spatial Nodes

Nodesets or sidesets named `SPACE` (case-insensitive) are exported as `MGR2` spatial nodes. These are used by ELF for field evaluation points.

## Usage Examples

### 3D Iron Core + Coil

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])

# Create geometry and mesh
cubit.cmd("create brick x 10 y 10 z 10")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

# Name blocks with ELF prefixes
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd('block 1 name "MMB_core"')

cubit.cmd('export meg "model.meg" overwrite')
```

### 2D Planar Analysis

```python
cubit.cmd('export meg "model_2d.meg" twod overwrite')
```

### Axisymmetric Analysis

```python
cubit.cmd('export meg "model_axi.meg" axisymmetric overwrite')
```

## Compatibility

- ELF/MAGIC
- JMAG (via MEG import)

## See Also

- [export_Nastran.md](export_Nastran.md) — Alternative for JMAG (BDF format)
- [Function_Reference.md](Function_Reference.md) — All plugin commands
