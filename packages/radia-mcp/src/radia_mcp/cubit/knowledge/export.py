"""
Export documentation for the radia-mcp Cubit server.

Curved Netgen ``.vol`` export runs through the MCP tool ``cubit_export_vol``:
an isolated worker loads the fork NGSolve, Cubit and the native
``export netgen`` command plugin.
"""

EXPORT_OVERVIEW = """
# Cubit .vol Export - Overview

## One export path

| Route | Call | Runs in |
|-------|------|---------|
| Native plugin in an isolated worker | `cubit_export_vol(journal, "mesh.vol", order=N)` | `python -I -B` worker + `cubit_mesh_export.dll` |

`cubit_export_vol` plays the journal in an isolated external Python
(`-I -B`) that loads the fork NGSolve, then Cubit with the native
`export netgen` command plugin (`cubit_mesh_export.dll`, built locally from
`packages/radia-mcp/native/cubit_vol_export`; `plugin_dir` or
`RADIA_CUBIT_PLUGIN_DIR`).  The worker runs `import ngsolve` first so the plugin binds to
the same nglib, plays the journal in batch/nographics Cubit (`-noinitfile`),
optionally adds the Kelvin sphere, then runs
`export netgen "<out>" order N [periodic args] overwrite`.  The plugin writes
the `.vol` and `.vol.json`; the tool returns the companion record and the
sha256.  A missing Cubit, worker runtime or plugin DLL fails loudly.

Named periodic export is opt-in and requires all nine `periodic_*` fields: name
`kt`, exact master/slave boundary labels, coordinate unit `m`, BC patch-pair
list, translation x/y/z and positive finite matching distance. BC numbers are
not FaceDescriptor or CAD surface IDs. Every point must match uniquely across
all incident patches; ambiguity or incomplete coverage fails the export.

Kelvin open boundary: `add_kelvin=True` (with `kelvin_air`, `kelvin_block`,
`kelvin_mesh`, `kelvin_sym_x/y/z` = off|bn|ht) adds the exterior sphere before
export; skipped when a Kelvin block exists; a missing air block or an invalid
symmetry value raises.

## Supported Formats

| Command | Format | Max Order | Notes |
|---------|--------|-----------|-------|
| `cubit_export_vol(journal, "f.vol", order=N)` | Netgen .vol (+ .vol.json) | 1-5 | Recommended for NGSolve |

Order 2+ uses Netgen `Mesh.Curve` with Cubit surface/curve projection callbacks.

## API

```python
# Netgen .vol (recommended for NGSolve FEM computation)
# Curved .vol export runs outside Cubit: save this model as a journal,
# then call the MCP tool cubit_export_vol(journal_path, "mesh.vol", order=3).

# Netgen .vol supports order 1-5
# Curved .vol export runs outside Cubit: save this model as a journal,
# then call the MCP tool cubit_export_vol(journal_path, "mesh.vol", order=5).
```

## Companion JSON (.vol.json)

`cubit_export_vol` writes a companion JSON with CAD reference values:
```json
{
  "materials": {"iron": 5.24e-04, "air": 3.38e-03},
  "boundaries": {"surface_1": 3.14e-02, "coil": 1.26e-02},
  "edges": {"curve_1": 3.14e-01},
  "n_elements": 10359, "n_points": 2071, "order": 3
}
```

`calc_verify_vol.py` reads this JSON for per-label consistency checks.

## Common Usage Pattern

```python
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

cubit.cmd("create sphere radius 0.05")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size auto factor 5")
cubit.cmd("mesh volume 1")

# Blocks define material labels
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "sphere"')

# Sidesets define boundary labels (optional)
cubit.cmd("sideset 1 add surface 1")
cubit.cmd('sideset 1 name "outer"')

# Export (surface elements extracted automatically from volume faces)
# Curved .vol export runs outside Cubit: save this model as a journal,
# then call the MCP tool cubit_export_vol(journal_path, "sphere.vol", order=3).
```

**Note**: `block 2 add tri all` is NOT needed for cubit_export_vol.
Surface elements are extracted from volume element faces automatically.

## IH (BEM) Inductance: Sideset Setup

IH (BEM) uses surface elements only (TRI). Source/sink terminal faces
must be defined as **sidesets** (not blocks) so they become boundary
labels in the .vol file.

```python
# Coil geometry
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "coil"')

# Terminal faces: sidesets -> boundary labels in .vol
cubit.cmd("sideset 1 add surface 3")   # source terminal
cubit.cmd('sideset 1 name "source"')
cubit.cmd("sideset 2 add surface 5")   # sink terminal
cubit.cmd('sideset 2 name "sink"')

# Export (sidesets become boundary labels)
# Curved .vol export runs outside Cubit: save this model as a journal,
# then call the MCP tool cubit_export_vol(journal_path, "coil.vol", order=2).
# python calc_inductance.py --vol coil.vol --source source --sink sink
```

**Important**: Both volume .vol and surface-only .vol work for IH (BEM).
The BEM solver uses BND elements only. Volume elements are ignored.
"""

EXPORT_CURVED = """
# Curved Mesh Export (cubit_export_vol)

```python
import tempfile
from ngsolve import Mesh
vol_path = tempfile.mktemp(suffix='.vol')
# Curved .vol export runs outside Cubit: save this model as a journal,
# then call the MCP tool cubit_export_vol(journal_path, "{vol_path}", order=3).
mesh = Mesh(vol_path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `order` | int (1-5) | Polynomial order for mesh curving |
| `overwrite` | flag | Overwrite existing file |

**Returns**: `.vol` file with embedded curvedelements section. Load with `Mesh(vol_path)`.

## How It Works

The plugin (`ExportNetgenCommand`) extracts blocks/sidesets/nodesets from
Cubit (`MeshData`), builds the first-order Netgen mesh, curves it with Netgen
callback geometry projected onto Cubit surfaces and curves (`NetgenCurver`,
`callbackgeom`, `netgen_curving/geometric_refit` for order >= 3), applies the optional
named periodic identification and writes the `.vol` plus `.vol.json`.

## Compatibility and warnings

- The `.vol` carries `edgesegmentsgi3` and needs NGSolve 6.2.2605 or newer.
- Pyramids need the patched Netgen 6.2.2608 runtime (apex clamp and
  base-edge blend); `check-vol` rejects them on any other runtime.  See
  `netgen_workflow_guide` ("Pyramids require the patched Netgen").
- The worker interpreter is the `worker_python` argument, else
  `$RADIA_CUBIT_WORKER_PYTHON`, else the server's interpreter; it is started
  with `-I -B`, so `radia_mcp` and the Netgen build must be installed there.

## Key Design Decision

The command exports **1st order elements** internally, then curves
them to the requested order using ACIS surface projection. Even if
blocks contain TET10, only corner nodes are used initially. High-order
nodes are placed exactly on ACIS CAD surfaces by BuildCurvedElements(order).

## Supported Elements

| Element | Nodes | Netgen Type |
|---------|-------|-------------|
| TET4 | 4 | Tet |
| HEX8 | 8 | Hex |
| WEDGE6 | 6 | Prism |
| PYRAMID5 | 5 | Pyramid |
| TRI3 | 3 | Trig |
| QUAD4 | 4 | Quad |
| EDGE2 | 2 | Segment |

## Node Ordering Conversion (Cubit -> Netgen)

| Element | Cubit Order | Netgen Order |
|---------|-------------|-------------|
| TET | [0,1,2,3] | [0,1,2,3] (same) |
| HEX | [0,1,2,3,4,5,6,7] | [0,1,5,4,3,2,6,7] |
| WEDGE | [0,1,2,3,4,5] | [0,2,1,3,5,4] |
| PYRAMID | [0,1,2,3,4] | [3,2,1,0,4] |

## Workflows

For high-order curving, see `netgen_workflow_guide()` tool.
For simple 2nd order without geometry, use `cubit_export_vol(journal, "mesh.vol", order=2)`.

## Deleted Predecessors

`export_NetgenMesh()`, `export_netgen()`, `export_netgen_with_names()`,
`set_*_geominfo()`, `name_occ_faces()` are all removed.
`cubit_export_vol` replaces all of them.
"""

EXPORT_EXODUS = """
# Exodus II Export

```python
export_exodus(cubit, FileName, overwrite=True, large_model=False)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .exo file path |
| `overwrite` | bool | True | Overwrite existing file |
| `large_model` | bool | False | Use 64-bit integers for large meshes |

## How It Works (Exception to Design Philosophy)

Unlike other export functions, `export_exodus()` uses Cubit's **built-in**
`export mesh` command internally:
```python
cubit.cmd('export mesh "filename.exo" overwrite')
```

This is because Exodus II is Cubit's native format, and the built-in command
provides full fidelity including all element types, nodesets, sidesets, etc.

## large_model Parameter

| Value | Integer Size | Max Nodes/Elements | Use When |
|-------|-------------|-------------------|----------|
| `False` | 32-bit | ~2 billion | Most meshes |
| `True` | 64-bit | Unlimited | Very large meshes (>2B nodes/elements) |

```python
# Standard export
cubit.cmd('export mesh "mesh.exo" overwrite')

# Large model (64-bit integers)
cubit.cmd('export mesh "mesh.exo" overwrite large')
```

The `large` flag is passed to Cubit's `export mesh` command, which writes
the Exodus file with 64-bit integer IDs for node and element numbering.

## Supported Elements

All Cubit element types, including:
- 0D: NODE
- 1D: BAR, BAR2, BAR3
- 2D: TRI3, TRI6, TRI7, QUAD4, QUAD8, QUAD9
- 3D: TET4, TET10, TET11, HEX8, HEX20, HEX27, WEDGE6, WEDGE15, PYRAMID5, PYRAMID13

## Use Cases

- Archival of Cubit meshes with full fidelity
- Downstream solvers: Sierra, MOOSE, Albany, ARIA
- Multi-physics applications requiring block/nodeset/sideset structure
- Meshes with any element order (1st, 2nd, or higher)
"""

EXPORT_COMPARISON = """
# Export Format Comparison

Netgen `.vol` is the only mesh output of the Cubit export route.

| Use Case | Command | Why |
|----------|---------|-----|
| NGSolve FEM (any order) | `cubit_export_vol(journal, "f.vol", order=N)` | Order 1-5 via the native plugin |
| Cubit-native archival | `export mesh "f.exo"` (Cubit built-in) | Exodus II, Cubit's own command |

The checked `.vol` is the only solver input; run `check-vol` before analysis.
"""


EXPORT_DECISION_GUIDE = """
# Export Format Decision Guide

## "I want to use NGSolve / Netgen for FEM"

-> Use `cubit_export_vol` (recommended for any order):
  ```python
  import tempfile
  from ngsolve import Mesh
  cubit.cmd("block 1 add volume all")
  vol_path = tempfile.mktemp(suffix='.vol')
  # Curved .vol export runs outside Cubit: save this model as a journal,
  # then call the MCP tool cubit_export_vol(journal_path, "{vol_path}", order=3).
  mesh = Mesh(vol_path)
  ```
  Works for ANY geometry shape. No STEP files, no OCC, no SetGeomInfo.

## "I need multi-physics (MOOSE, Sierra)"

-> Use `export_exodus()`
- Supports all element types and orders
- Full block definitions preserved
- Large model support: `export_exodus(cubit, "mesh.exo", large_model=True)`
- Cubit-native format with highest fidelity

## Performance & Feature Summary

| Format | Max Order | 2D | 3D | HO Method |
|--------|-----------|----|----|-----------|
| netgen (.vol) | 5 | No | Yes | Native plugin: Netgen callback geometry + Cubit projection |

## Troubleshooting: High-Order Export

- `kind="environment"` from `cubit_export_vol`: the worker interpreter cannot
  import ngsolve/netgen, the plugin DLL is missing from `plugin_dir`, or Cubit
  was not found.  Run `cubit_doctor` (check `export_worker`).
- `status="error"` with the worker's message: the journal or `export netgen`
  failed in Cubit, or a label/periodic/Kelvin contract was rejected.  Fix the
  model, then export to a new path; existing `.vol` files are never overwritten.
"""


P_CONVERGENCE_KNOWLEDGE = """
# p-Convergence Testing and Known Issues

## What is p-Convergence?

When exporting a curved mesh at increasing polynomial orders (p=1..5),
the volume and area errors vs CAD geometry should decrease monotonically.
Each order should gain ~2-3 digits of accuracy.

Example (Cylinder Tet, p=1..5):
```
Order  Nodes  Volume error [%]   Area error [%]
  1       56   -2.15e+00          -3.41e+00
  2      300   -4.82e-02          -1.03e-01
  3      938   -2.11e-04          -7.28e-04
  4     2051   -8.24e-07          -3.15e-06
  5     3720   -2.67e-09          -1.42e-08
```

Only models with CURVED geometry show p-convergence. Planar models
(brick, webcut without curves) are exact at order 1.

## Verified Test Models (p-convergence passes)

| Model | Elements | Geometry | Element Types |
|-------|----------|----------|---------------|
| Cylinder Tet | 142 | Cylindrical surface | Tet |
| Cylinder Tet+BL | 199 | Cylindrical + BL | Tet + Wedge |
| Acorn (cyl+sphere) | 1935 | Sphere + cylinder | Tet + Wedge |
| Wedge Only | 144 | Cylindrical | Wedge |
| Cylinder Hex webcut | 96 | Cylindrical with cuts | Hex |
| Hex Only (planar) | 56 | Planar (exact at p=1) | Hex |
| Tet+Hex+Pyramid | 12 | Planar (exact at p=1) | Mixed |

## Known Issue: ACIS Loft Surface (05_loft)

**Problem**: Loft surfaces (rectangle to circle) in ACIS have
`closest_point_trimmed()` that overshoots the actual surface.

**Symptom**: Area error ~0.6% that does NOT converge with p-refinement.
Volume converges slowly but does converge.

**Workaround**: None. This is an ACIS/Cubit kernel limitation.
Volume convergence is acceptable. Use this geometry for volume-only tests.

## Known Behavior: Wedge/BL Volume = 0 in .vol.json

When boundary layer (BL) creates wedge elements, the .vol.json companion
file may report `0.000000e+00` volume for the BL material. This is because
Cubit's `RefVolume::measure()` returns the ACIS volume of the parent body,
and BL-generated thin layers may not have an independent ACIS volume.

This is NOT a bug. The mesh itself has correct wedge geometry.
NGSolve `Integrate(CF(1), mesh)` will return the correct volume.

## Known Issue: Hex BL on Cylinder (03) Export Failure

Hexahedral boundary layer on a cylinder with `scheme map` may fail
to export via `cubit_export_vol`. The combination of mapped hex meshing
with boundary layers on curved surfaces can produce elements that
NetgenCurver cannot process.

**Workaround**: Use tet meshing with BL (creates wedge prisms instead).

## BND Integrate Area Mismatch

`Integrate(CF(1), mesh, BND)` integrates ALL boundary surface elements,
including internal block-to-block interfaces. The .vol.json reports only
EXTERNAL surface areas from Cubit's ACIS geometry.

To get per-boundary area matching .vol.json, use per-label integration:
```python
for bnd_name, cad_area in cad["boundaries"].items():
    ng_area = Integrate(CF(1), mesh, BND,
                        definedon=mesh.Boundaries(bnd_name))
    error = (ng_area - cad_area) / cad_area * 100
```

Do NOT compare `Integrate(CF(1), mesh, BND)` (total) against
`sum(cad["boundaries"].values())` -- the total BND includes internal faces.

## .vol.json Companion File

Every `cubit_export_vol` produces a companion .vol.json:
```json
{
  "materials": {"sphere": 5.236e-04},
  "boundaries": {"surface_1": 3.142e-02},
  "edges": {"curve_1": 3.142e-01},
  "n_elements": 10359, "n_points": 2071, "order": 3
}
```

Use for:
- Automated p-convergence regression tests
- Standalone consistency checks without Cubit (`check-vol model.vol`; the
  sibling sidecar is auto-discovered when present)
- Strict application label checks (`check-vol model.vol --contract labels.json
  --strict-labels --report-json run/vol_check.json`)

Run the checker after export and before solver/Simulink initialization.
Material constants are validated by the application DesignSpec/configuration,
not inferred from `.vol` label names.
"""


def get_export_documentation(format: str = "all") -> str:
	"""Return export documentation by format name."""
	topics = {
		"overview": EXPORT_OVERVIEW,
		"curved": EXPORT_CURVED,
		"netgen": EXPORT_CURVED,  # Alias: old name redirects to cubit_export_vol
		"exodus": EXPORT_EXODUS,
		"comparison": EXPORT_COMPARISON,
		"decision_guide": EXPORT_DECISION_GUIDE,
		"p_convergence": P_CONVERGENCE_KNOWLEDGE,
	}

	format = format.lower().strip()
	if format == "all":
		return EXPORT_OVERVIEW + "\n\n" + EXPORT_COMPARISON
	elif format in topics:
		return topics[format]
	else:
		return (
			f"Unknown format: '{format}'. "
			f"Available: all, {', '.join(topics.keys())}"
		)
