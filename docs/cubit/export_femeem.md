# export_femeem

Export mesh to FEMEEM format (Gifu University 3D FEM solver).

## Syntax

```
export femeem "dirname" [scale <value>] [overwrite]
```

Creates a directory containing the four files required by FEMEEM.
Only 1st-order tetrahedral elements are exported; non-tet elements (hex, wedge, pyramid) are skipped with a warning.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `scale <value>` | 1.0 | Coordinate scale factor (coordinates stored as `coord * scale`) |
| `overwrite` | off | Overwrite existing directory |

## Output Files

| File | Description |
|------|-------------|
| `in.dat` | Header + tet connectivity (2 per line) + scaled node coordinates |
| `sin.dat.B` | Boundary conditions (type 8 edges from sidesets) + element block assignments + material constants |
| `sina.dat` | Analysis control parameters template (AMGTYP=1 for Cubit mesh) |
| `d3` | File reference index |

## Supported Elements

**Tetrahedra only** (TET4). All other element types are skipped.

If the mesh contains no tetrahedra, the export fails with an error.

## Boundary Conditions

Sidesets are exported as fixed boundary edges (type 8) in `sin.dat.B`.
Each triangle edge in a sideset face becomes a boundary condition entry.

Node IDs are renumbered to contiguous 1-based indices.

## Material Assignment

Cubit block IDs are mapped to FEMEEM material numbers (1-based, in order of appearance). Material constants are initialized to placeholder values (`1.0 1.0 1.0 0.0`) — edit `sin.dat.B` to set actual material properties.

## Usage Example

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Blocks for material assignment
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd('block 1 name "iron"')

# Sidesets for boundary conditions
cubit.cmd("sideset 1 add surface all")
cubit.cmd('sideset 1 name "outer_boundary"')

# Export with mm scale (FEMEEM expects mm)
cubit.cmd('export femeem "femeem_output" scale 1000 overwrite')
```

This creates:
```
femeem_output/
  in.dat
  sin.dat.B
  sina.dat
  d3
```

## in.dat Format

```
  <npoint>  <nelem>       0       0       <scale>
  <n1><n2><n3><n4><n5><n6><n7><n8>    (2 tets per line, 8d format)
  ...
  <x><y><z>                            (scaled coordinates, 24.12E format)
  ...
```

## Limitations

- Tetrahedra only (1st order)
- No high-order element support
- Material properties are placeholder values — manual editing required
- Boundary conditions are type 8 (fixed) only

## See Also

- [export_NetgenMesh.md](export_NetgenMesh.md) — For NGSolve FEM workflows
- [Function_Reference.md](Function_Reference.md) — All plugin commands
