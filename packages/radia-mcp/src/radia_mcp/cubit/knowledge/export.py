"""
Export format documentation for Radia mesh export MCP server.

Provides API reference, parameter tables, supported element types,
and usage guidance for each mesh export format.
"""

EXPORT_OVERVIEW = """
# Radia Mesh Export - Export Functions Overview

## Two Export Paths (must produce identical results for tet meshes)

| Path | Method | Speed | Use Case |
|------|--------|-------|----------|
| **Path A** (C++) | `cubit.cmd('export netgen/gmsh/...')` | Fast | APREPRO journal, GUI menu |

Path A is the recommended export path. Use `export netgen "mesh.vol" order N` for all workflows.

## Supported Formats

| Command | Format | Max Order | Notes |
|---------|--------|-----------|-------|
| `export netgen "f.vol" order N` | Netgen .vol (+ .vol.json) | 1-5 | Recommended for NGSolve |
| `export gmsh "f.msh" order N` | Gmsh v4.1 (.msh + .geo/.opt) | 1-3 | Raw data plus launch companion; open `.geo` |
| `export jmag_nastran "f.bdf" order N` | Nastran BDF (.bdf) | 1-2 | nopyramid for JMAG |
| `export vtk "f.vtk" order N` | VTK Legacy (.vtk) | 1-2 | ParaView visualization |
| `export meg "f.meg"` | ELF/MAGIC MEG | 1 | Block names define ELF prefixes |
| `export femeem "dir"` | FEMEEM (Gifu Univ.) | 1 (tet only) | Creates directory with 4 files |

Order 2+ uses NetgenCurver (compact_netgen BuildCurvedElements + ACIS projection).
No fallback to HighOrderMesh (removed).

## API

```python
# Netgen .vol (recommended for NGSolve FEM computation)
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
# -> produces mesh.vol + mesh.vol.json (CAD reference values)

# Gmsh v4.1 raw display data (Radia post launches the .geo companion)
cubit.cmd('export gmsh "mesh.msh" order 3 overwrite')
# -> mesh.msh + mesh.geo + mesh.geo.opt + mesh.msh.opt

# Netgen .vol supports order 1-5
cubit.cmd('export netgen "mesh.vol" order 5 overwrite')
```

## Companion JSON (.vol.json)

`export netgen` writes a companion JSON with CAD reference values:
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
cubit.cmd('export netgen "sphere.vol" order 3 overwrite')
```

**Note**: `block 2 add tri all` is NOT needed for export netgen.
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
cubit.cmd('export netgen "coil.vol" order 2 overwrite')

# Compute (no Cubit needed):
# python calc_inductance.py --vol coil.vol --source source --sink sink
```

**Important**: Both volume .vol and surface-only .vol work for IH (BEM).
The BEM solver uses BND elements only. Volume elements are ignored.
"""

EXPORT_GMSH_V2 = """
# Gmsh v4.1 Export (order 1-3)

```
export gmsh "filename.msh" [order <1-3>] [dimension <2|3>] [overwrite]
```

v4.1 is the only supported GMSH format. No block assignment required.
For Radia post-processing, treat this `.msh` as raw mesh/data and open the
generated `.geo` launch recipe. The command should attach `case.geo`,
`case.geo.opt`, and `case.msh.opt`; a plain `case.opt` is not an Explorer
auto-load contract. Do not make `.msh` double-click the primary user-facing
post-processing contract.

## Supported Elements

| Order 1 | Order 2 | Order 3 | Gmsh Type |
|---------|---------|---------|-----------|
| TET4 | TET10 (11) | TET20 (29) | 4/11/29 |
| HEX8 | HEX20 (17) | HEX32 (99) | 5/17/99 |
| WEDGE6 | WEDGE15 (18) | not supported | 6/18 |
| PYRAMID5 | PYR13 (19) | PYR21 (125) | 7/19/125 |
| TRI3 | TRI6 (9) | TRI10 (21) | 2/9/21 |
| QUAD4 | QUAD8 (16) | QUAD12 (39) | 3/16/39 |

Order 4-5 not supported in Gmsh export (use `export netgen`).
Wedge/Prism limited to order 2 (Gmsh limitation).

## Use Cases

- **GMSH raw display data**: View mesh in GMSH GUI, preferably through `.geo`
- **NGSolve**: Use `export netgen` instead (order 1-5)

```python
# Gmsh visualization (order 1-3)
cubit.cmd('export gmsh "mesh.msh" order 2 overwrite')

# NGSolve FEM (recommended: .vol)
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
```
"""

EXPORT_GMSH_V4 = """
# Gmsh v4.1 Format Policy

v4.1 is the only supported GMSH format (v2.2 removed, lab-wide standard).

## Data Flow

| Direction | Format | Purpose | Tool |
|-----------|--------|---------|------|
| **Input** (-> NGSolve) | **.vol** | Mesh import into NGSolve | `export netgen "mesh.vol" order N` -> `Mesh("mesh.vol")` |
| **Output** (Cubit ->) | **.msh v4.1 + .geo/.opt** | Mesh visualization/post launch | `export gmsh "mesh.msh" order N` + `.geo` recipe |
| **Output** (NGSolve ->) | **.msh v4.1 + .geo/.opt** | Field visualization/post launch | `GmshPostExport.write()` + `.geo` recipe |

The ONLY input format for NGSolve is `.vol`. `.msh` is output-only (Cubit mesh export or NGSolve field visualization).
The standard file to open for Radia post-processing is `.geo`; `.msh` association is optional raw mesh/data inspection.

### GmshPostExport Methods

| Method | Format | High-order | Use case |
|--------|--------|-----------|----------|
| `write(filename)` | v4.1 | Yes (any order) | Field visualization in GMSH GUI |
| `write_mesh(filename)` | v4.1 | Yes (any order) | Mesh only (no field data) |

## When to Use Which

- **`export netgen "mesh.vol"`**: For NGSolve FEM computation (any order, recommended)
- **`export gmsh "mesh.msh"`**: For Cubit mesh visualization data in GMSH (order 1-3)
- **`GmshPostExport.write()`**: For NGSolve field result data in GMSH
- **`.geo` companion**: Standard Radia post-processing launch artifact
  with exact `case.geo.opt` sidecar
- **`case.msh.opt` companion**: Raw mesh/data inspection sidecar when a
  user intentionally opens the `.msh` directly
- A plain `case.opt` is not auto-loaded when either `case.geo` or
  `case.msh` is opened
"""

EXPORT_CURVED = """
# Curved Mesh Export (export netgen)

```python
import tempfile
from ngsolve import Mesh
vol_path = tempfile.mktemp(suffix='.vol')
cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
mesh = Mesh(vol_path)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `order` | int (1-5) | Polynomial order for mesh curving |
| `overwrite` | flag | Overwrite existing file |

**Returns**: `.vol` file with embedded curvedelements section. Load with `Mesh(vol_path)`.

## How It Works

Uses NetgenCurver (compact_netgen C++ static link) with CallbackGeometry
to delegate surface projection to Cubit's ACIS kernel:

1. Reads mesh topology (nodes, elements) from Cubit blocks
2. Creates 1st order netgen.meshing.Mesh
3. Registers CallbackGeometry with ACIS surface projection
4. Calls BuildCurvedElements(order) using the CallbackGeometry
5. Writes .vol with curvedelements section + companion .vol.json

No STEP files, no OCC geometry, no SetGeomInfo needed.

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
For simple 2nd order without geometry, use `export netgen "mesh.vol" order 2`.

## Deleted Predecessors

`export_NetgenMesh()`, `export_netgen()`, `export_netgen_with_names()`,
`set_*_geominfo()`, `name_occ_faces()` are all removed.
`export netgen` replaces all of them.
"""

EXPORT_NASTRAN = """
# Nastran BDF Export (order 1-2)

```
export jmag_nastran "filename.bdf" [order <1|2>] [dimension <2|3>] [nopyramid] [overwrite]
```

No block assignment required.

| Parameter | Default | Description |
|-----------|---------|-------------|
| order | 1 | 1=CTETRA/CHEXA, 2=CTETRA10/CHEXA20 (via NetgenCurver) |
| dimension | 3 | 2D (CTRIA3/CQUAD4) or 3D |
| nopyramid | off | Convert pyramids to degenerate hex (JMAG compatible) |

## Supported Elements

| Order 1 | Order 2 | Card |
|---------|---------|------|
| CTETRA (4) | CTETRA (10) | Tetrahedron |
| CHEXA (8) | CHEXA (20) | Hexahedron |
| CPENTA (6) | CPENTA (15) | Wedge/Prism |
| CPYRAM (5) | CPYRAM (13) | Pyramid |
| CTRIA3 (3) | CTRIA6 (6) | Triangle |
| CQUAD4 (4) | CQUAD8 (8) | Quadrilateral |

**IMPORTANT**: Use `export jmag_nastran`, NOT `export nastran` (Cubit built-in has different format).

Nastran BDF is export-only. It is NOT a supported input path for Radia/NGSolve.
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

## Decision Matrix: Which Format to Use

| Use Case | Recommended Format | Why |
|----------|-------------------|-----|
| NGSolve FEM (any order) | `export netgen "f.vol" order N` | Order 1-5 via ACIS CallbackGeometry |
| GMSH visualization | `export gmsh "f.msh" order N` | v4.1, order 1-3 |
| JMAG / Nastran solver | `export jmag_nastran "f.bdf" order N` | Order 1-2, nopyramid for JMAG |
| ParaView visualization | `export vtk "f.vtk" order N` | Order 1-2 |
| ELF/MAGIC solver | `export meg "f.meg"` | Order 1, ELF element type labels |
| FEMEEM solver | `export femeem "dir"` | Order 1, tet only |

## Feature Comparison

| Feature | netgen | gmsh | nastran | vtk | meg | femeem |
|---------|--------|------|---------|-----|-----|--------|
| Max order | 5 | 3 | 2 | 2 | 1 | 1 |
| Tet | Yes | Yes | Yes | Yes | Yes | Yes |
| Hex | Yes | Yes | Yes | Yes | Yes | No |
| Wedge | Yes | Yes (o2) | Yes | Yes | Yes | No |
| Pyramid | Yes | Yes | Yes | Yes | degen hex | No |
| Labels | block+sideset | block+sideset | block | block | block+sideset | block |
| 2D mode | No | Yes | Yes | Yes | Yes (twod/axi) | No |
| Companion | .vol.json | .geo + .geo.opt launch recipe | - | - | - | d3 |

## export netgen vs Gmsh for NGSolve

| Aspect | export netgen | export gmsh |
|--------|---------------------|-------------------|
| Max order | 5 | 3 (wedge: 2) |
| NGSolve input | Yes (.vol is the ONLY input) | No (.msh not supported) |
| Best for | NGSolve FEM/BEM | GMSH visualization raw data; open `.geo` for Radia post |
"""


EXPORT_DECISION_GUIDE = """
# Export Format Decision Guide

## "I want to use NGSolve / Netgen for FEM"

-> Use `export netgen` (recommended for any order):
  ```python
  import tempfile
  from ngsolve import Mesh
  cubit.cmd("block 1 add volume all")
  vol_path = tempfile.mktemp(suffix='.vol')
  cubit.cmd(f'export netgen "{vol_path}" order 3 overwrite')
  mesh = Mesh(vol_path)
  ```
  Works for ANY geometry shape. No STEP files, no OCC, no SetGeomInfo.

- **Alternative (order 1-3)** -> Use `export gmsh` for GMSH visualization:
  ```python
  cubit.cmd('export gmsh "mesh.msh" order 3 overwrite')
  # Supports order 1-3 (wedge limited to order 2). Order 4-5: error.
  # For NGSolve use export netgen instead
  ```

## "I need structural FEA (Nastran / JMAG)"

-> Use `export jmag_nastran` (order 1-2):
  ```python
  cubit.cmd('export jmag_nastran "mesh.bdf" order 2 overwrite')
  cubit.cmd('export jmag_nastran "mesh.bdf" order 2 nopyramid overwrite')  # JMAG
  ```

### JMAG-specific: Pyramid Element Problem

JMAG **cannot read standard CPYRAM** (5-node pyramid) elements. When Cubit
generates mixed hex-tet meshes, pyramid transition elements appear at
the interface. Two solutions:

1. **Use nopyramid** (recommended): Writes pyramids as degenerate CHEXA
   ```python
   cubit.cmd('export jmag_nastran "mesh.bdf" nopyramid overwrite')
   ```

2. **Use pure tet mesh**: Avoids pyramids entirely
   ```python
   cubit.cmd("volume all scheme tetmesh")
   # ... mesh and export ...
   ```

**When to give up on hex**: If Cubit cannot webcut/decompose the geometry
into sweepable sub-volumes, switch to `tetmesh`. Forcing hex on complex
geometry leads to poor quality or failed meshing.

## "I need multi-physics (MOOSE, Sierra)"

-> Use `export_exodus()`
- Supports all element types and orders
- Full block definitions preserved
- Large model support: `export_exodus(cubit, "mesh.exo", large_model=True)`
- Cubit-native format with highest fidelity

## "I want Gmsh visualization"

-> Use `export gmsh` for GMSH GUI viewing (v4.1, order 1-3):
  ```python
  cubit.cmd('export gmsh "mesh.msh" order 2 overwrite')
  ```
- **Note**: For NGSolve FEM, use `export netgen "mesh.vol"` instead

## Performance & Feature Summary

| Format | Max Order | 2D | 3D | HO Method |
|--------|-----------|----|----|-----------|
| netgen (.vol) | 5 | No | Yes | NetgenCurver + ACIS |
| gmsh (.msh) | 3 | Yes | Yes | NetgenCurver + ACIS |
| nastran (.bdf) | 2 | Yes | Yes | NetgenCurver + ACIS |
| vtk (.vtk) | 2 | Yes | Yes | NetgenCurver + ACIS |
| meg (.meg) | 1 | Yes | Yes | - |
| femeem (dir) | 1 (tet) | No | Yes | - |

## IMPORTANT: `export jmag_nastran` (NOT `export nastran`)

Cubit has a **built-in** `export nastran` command (e.g., `export nastran "f.bdf" overwrite everything`).
Radia's Nastran export uses a DIFFERENT command name to avoid conflict:

```python
# CORRECT: Radia's export (supports order 2, nopyramid, block labels)
cubit.cmd('export jmag_nastran "mesh.bdf" order 2 dimension 3 overwrite')

# WRONG: Cubit's built-in (different format, no order 2, no nopyramid)
cubit.cmd('export nastran "mesh.bdf" overwrite everything')
```

## Coil APREPRO Command

Generate coil STEP from a Python script and import into Cubit:

```python
# In Cubit command line or .jou file:
cubit.cmd('coil "my_coil.py"')                          # generate + import
cubit.cmd('coil "my_coil.py" output "C:/out/coil.step"') # custom output path
cubit.cmd('coil "my_coil.py" noimport')                  # STEP only, no import

# my_coil.py must define build_coil() -> CoilBuilder:
# def build_coil():
#     from radia.coil_builder import CoilBuilder
#     cb = CoilBuilder(current=1000)
#     cb.set_start([0, 0, 0])
#     cb.set_cross_section(width=0.02, height=0.02)
#     cb.add_straight(length=0.1, tilt=0)
#     cb.add_arc(radius=0.05, arc_angle=180, tilt=0)
#     ...
#     return cb
```

Requires external Python 3.12 with NGSolve/OCC (not Cubit's embedded 3.10).
Set `RADIA_PYTHON` env var to override Python path.

## Troubleshooting: High-Order Export

### "Interrupt Detected" on AddPoint (NetgenCurver crash)

**Cause**: ABI mismatch between ccm plugin and nglib.dll at runtime.
The ccm was built against one version of Netgen, but a different nglib.dll
is loaded from the plugins/ directory.

**Fix**: Rebuild ccm with **compact_netgen** (static link, no nglib.dll dependency).
```bash
cmake ...  # netgen sources in-repo (compact_netgen/netgen_src/), no external path needed
```

### HEX20/TET10 has wrong node count in .msh/.bdf (order 2)

**Symptom**: GMSH crashes opening .msh, or HEX20 element has 8 nodes instead of 20.

**Cause** (fixed 2026-04-05): Volume element internal edges were not registered
in `edge_ho_nodes_`. Only surface element edges were registered, so internal
edges (not shared with any surface) had no HO nodes in the .msh/.bdf connectivity.
The .vol export was unaffected (uses Netgen internal mesh, not edge_ho_nodes_).

**Verification**: Check HEX20 node count in .msh file:
```python
# In exported .msh, HEX20 (type 17) must have 20 nodes per element
# Fields: elem_id type_id n_tags tag1 tag2 node1..node20
# If only 8 nodes -> bug (edge_ho_nodes_ not registered for volume edges)
```

## Cubit GUI Menu Structure

```
Menu bar: File Edit View Display Tools Export_Mesh Help Solve
                                       (PySide6)         (PySide6)
Export Mesh:                           Solve:
  Netgen Vol (.vol)...                   Radia-NGSolve...
  GMSH...                               Generate Coil...
  Nastran BDF...                         --------
  VTK...                                 Reload Panels
  MEG...
  FEMEEM...
```

- **Export Mesh**: PySide6 menu/dialogs calling C++ APREPRO export commands
- **Solve**: Python register_toolbar.py (subprocess to external Python 3.12)
- **Mesh p-convergence demo**: `docs/cubit_mesh_export/netgen/p_convergence_demo.ipynb`
  uses explicit `cubit.cmd(...)` export commands and is not a Cubit menu action
- Settings saved to `AppData/Roaming/Radia/export_settings.json`
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
to export via `export netgen`. The combination of mapped hex meshing
with boundary layers on curved surfaces can produce elements that
NetgenCurver cannot process.

**Workaround**: Use tet meshing with BL (creates wedge prisms instead).
Or use `export gmsh` which does not require NetgenCurver for order 1.

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

## Known Limitation: VTK Quadratic Pyramid (GMSH read)

GMSH API cannot read VTK cell type 27 (VTK_QUADRATIC_PYRAMID, 13 nodes).
This causes "Unknown type of cell 27" error in docs mesh-evaluation round-trip
when the model contains pyramids and VTK order 2 is tested.

This is a GMSH limitation, not a VTK export bug. The VTK file itself is
valid and can be opened in ParaView or other VTK readers.

VTK cell types used by Radia export:
| Element | Order 1 | Order 2 |
|---------|---------|---------|
| Tet     | 10      | 24 (VTK_QUADRATIC_TETRA) |
| Hex     | 12      | 25 (VTK_QUADRATIC_HEXAHEDRON) |
| Wedge   | 13      | 26 (VTK_QUADRATIC_WEDGE) |
| Pyramid | 14      | 27 (VTK_QUADRATIC_PYRAMID) -- GMSH cannot read |
| Tri     | 5       | 22 (VTK_QUADRATIC_TRIANGLE) |
| Quad    | 9       | 23 (VTK_QUADRATIC_QUAD) |

## .vol.json Companion File

Every `export netgen` produces a companion .vol.json:
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
- Consistency checks without Cubit (`check-vol model.vol --json model.vol.json`)
- Cross-format verification (export -> GMSH reload -> compare vs .vol.json)
"""


def get_export_documentation(format: str = "all") -> str:
	"""Return export documentation by format name."""
	topics = {
		"overview": EXPORT_OVERVIEW,
		"gmsh_v2": EXPORT_GMSH_V2,
		"gmsh_v4": EXPORT_GMSH_V4,
		"curved": EXPORT_CURVED,
		"netgen": EXPORT_CURVED,  # Alias: old name redirects to export netgen
		"nastran": EXPORT_NASTRAN,
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
