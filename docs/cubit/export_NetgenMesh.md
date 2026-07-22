# export_NetgenMesh — Cubit to Netgen .vol Export

Export mesh to Netgen `.vol` format with high-order curved elements (order 1-5).

## Syntax

```
export netgen "filename.vol" [order <1-5>] [overwrite]
```

No block assignment required — all meshed elements are exported automatically.
Sidesets become boundary labels. Block names become material labels.

> **NOTE**: `export netgen` is provided by the Radia Cubit plugin.
> Cubit has no built-in Netgen exporter, so there is no command conflict.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `order 1` | | Linear elements |
| `order 2` | yes | Quadratic — edge mid-nodes via NetgenCurver + ACIS projection |
| `order 3` | | Cubic — edge + face interior nodes |
| `order 4` | | Quartic |
| `order 5` | | Quintic |
| `overwrite` | off | Overwrite existing file |

### Output Files

| File | Description |
|------|-------------|
| `file.vol` | Netgen mesh with `curvedelements` section (text format) |
| `file.vol.json` | Companion JSON with CAD reference values for consistency checks |

## How It Works

The C++ plugin (`cubit_mesh_export.ccm`) performs the entire export without Python:

1. **Extract** — `MeshData::extract(order)` reads linear mesh from Cubit's `MeshExportInterface`
2. **Curve** — `NetgenCurver` builds a `netgen.meshing.Mesh` and calls `BuildCurvedElements(order)` using an ACIS `CallbackGeometry` that projects nodes onto the exact CAD surfaces
3. **Label** — Block names → material labels, sideset names → boundary labels
4. **Save** — `ng_mesh->Save(filename)` writes `.vol` text format (including `curvedelements` section)
5. **JSON** — CAD reference volumes, areas, and edge lengths from Cubit geometry kernel

No STEP export, no Python subprocess, no intermediate files.

## Basic Workflow

```python
import cubit
from ngsolve import Mesh

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create cylinder height 2 radius 0.5")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")

# Export with order 3 curving
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')

# Load in NGSolve — no STEP file needed
mesh = Mesh("mesh.vol")
print(mesh.GetMaterials())   # Material labels from block names
print(mesh.GetBoundaries())  # Boundary labels from sideset names
```

## Volume Accuracy by Curve Order

Measured on a sphere mesh (tet):

| Order | Volume Error |
|-------|-------------|
| 1     | ~1.4%       |
| 2     | ~0.001%     |
| 3     | ~1e-5%      |
| 5     | ~1e-8%      |

## Geometry Caveats for High Accuracy

### Avoid `unite volume all` on coil-like swept assemblies

For geometries built from many swept/lofted segments that share wire-shaped
cross-sections (e.g. multi-turn coils, helical bundles), **do not run
`unite volume all` before export**. The unite step merges end-cap disks at
shared terminals and introduces a polar-NURBS face whose parameter center
is singular (all UVs at one radial coordinate map to a single 3D point).
This degenerate parameterization causes the projection callback to produce
large displacement spikes (~wire diameter) that pollute neighboring
curved-tet mid-edge nodes.

Observed on a 3-turn coil (382 loft segments, wire ø6.3 mm, mesh size 6):

| Workflow | p=2 volume error | Projection rejects | max node disp |
|---|---|---|---|
| `unite volume all` before export | −0.67% | 10,115 | 61.4 mm (disk) |
| No `unite`, 382 independent volumes | **−0.31%** | 751 | 3.64 mm |

**Recommended**: mesh each lofted body independently and let the `.vol`
export preserve multi-volume topology. Downstream solvers (NGSolve, Kelvin)
handle multi-domain meshes natively; there is no accuracy penalty from
keeping volumes split, and the curving quality is nearly 2× better.

### Keep mesh size smaller than local feature size

For a wire-like cross-section with diameter `d`, pick tet size `≲ d / 2`
so that mesh triangles do not straddle the full cross-section. When
`mesh_size ≈ d` (e.g. `size 6` on ø6.3 mm wire), interior tet midpoints
land inside the volume and the surface projection becomes ill-conditioned.

## Supported Elements

### Volume Elements

| Cubit Element | Netgen | Order 2 | Order 3 | Order 5 |
|---------------|--------|---------|---------|---------|
| TET4          | TET    | TET10   | TET20   | TET56   |
| HEX8          | HEX    | HEX20   | HEX64   | —       |
| WEDGE6        | PRISM  | PRISM15 | —       | —       |
| PYRAMID5      | PYRAMID| PYR13   | —       | —       |

### Surface Elements

| Cubit Element | Netgen |
|---------------|--------|
| TRI3          | TRIG   |
| QUAD4         | QUAD   |

### Edge Elements

| Cubit Element | Netgen  |
|---------------|---------|
| EDGE2         | SEGMENT |

## Block and Sideset Mapping

### Block Names → Materials

```python
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd('block 1 name "iron"')
cubit.cmd("block 2 add tet all in volume 2")
cubit.cmd('block 2 name "air"')

cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
# mesh.GetMaterials() -> ('iron', 'air')
```

### Sideset Names → Boundaries

```python
cubit.cmd("sideset 1 add surface 1 2 3")
cubit.cmd('sideset 1 name "dirichlet"')
cubit.cmd("sideset 2 add surface 4 5")
cubit.cmd('sideset 2 name "neumann"')

# mesh.GetBoundaries() -> ('dirichlet', 'surface_6', 'neumann', ...)
```

Surfaces without an explicit sideset get auto-generated names (`surface_N`).
These fallback names are useful diagnostics only. Production meshes must name
solver-facing sidesets explicitly; `check-vol --strict-labels` rejects
`surface_N` / `Surface_N` and similar generated material names.

### Kelvin Auto-Detection

When blocks named `air` and `kelvin` are present, the plugin automatically:
- Labels the shared air|kelvin interface as `kelvin_int`
- Labels the outer kelvin boundary (DomainOut=0) as `kelvin_ext`
- Creates periodic identification pairs between inner and outer vertices

This enables Kelvin infinite element analysis without manual sideset setup.

## Companion JSON

The `.vol.json` file contains CAD reference values for mesh quality validation:

```json
{
  "materials": {"iron": 5.236e-04, "air": 1.047e-03},
  "boundaries": {"dirichlet": 3.142e-02, "neumann": 6.283e-03},
  "edges": {"curve_1": 3.142e-01},
  "n_elements": 10359,
  "n_points": 2071,
  "order": 3,
  "export_time_s": 1.234
}
```

Run `check-vol` after export and before handing the mesh to a solver or
Simulink block. The sibling sidecar is auto-discovered:

```bash
check-vol model.vol --strict-labels --report-json run/vol_check.json
```

Use a versioned `radia.vol-label-contract.v1` file with `--contract` for an
application/mode's required and allowed materials, boundaries, BBND labels,
and BBBND labels. `--json model.vol.json` is available when the exact sidecar
must be mandatory rather than auto-discovered.

## DomainIn/DomainOut

The plugin resolves `FaceDescriptor.DomainIn`/`DomainOut` from Cubit's surface-volume topology. This is required for:
- Multi-material problems (correct material assignment on shared faces)
- Periodic identification
- Kelvin boundary conditions

## In-Memory Curving

The old top-level `cubit_mesh_export.extract_curved_mesh` helper is retired
from the public API.  Use `export netgen` for production workflows; it handles
curving, labels, Kelvin detection, and companion JSON in one Cubit command.

The low-level `cubit_mesh_curver.build_curved_mesh` pybind module remains an
internal implementation detail of the `.ccm` plugin.  It expects already
extracted mesh arrays plus geometry callbacks and is not the user-facing
entry point.

## Troubleshooting

### "NetgenCurver not available"

NetgenCurver is compiled into the plugin via `compact_netgen` (static link). Ensure the plugin `.ccm` was built with `-DNETGEN_SRC_DIR` pointing to a valid Netgen source tree.

### "Interrupt Detected" on AddPoint

ABI mismatch between the plugin's Netgen and a separately installed Netgen DLL. Rebuild with `compact_netgen` to eliminate the dependency.

### Incorrect Boundary Labels

Check that sidesets are defined before export. Surfaces not in any sideset get auto-generated names (`surface_N`).

### Hex Boundary Layer Failure

If `export netgen` fails on hex boundary layers, export as Gmsh order 1 (`export gmsh "mesh.msh"`) as a fallback.

## See Also

- [export_Gmsh.md](export_Gmsh.md) — Gmsh v4.1 export (order 1-3)
- [export_Nastran.md](export_Nastran.md) — Nastran BDF export (order 1-2)
- [Function_Reference.md](Function_Reference.md) — All plugin commands
