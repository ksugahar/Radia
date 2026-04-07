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
| **Path B** (Python) | `extract_curved_mesh(cubit, order=N)` | Slower | Python scripting, reference |

Path A and Path B MUST produce bit-identical curving for tet meshes.
Run `test_vol_multi_geometry.py` after any NetgenCurver or bridge.py change.

## Supported Formats

| Command | Format | Order 1 | Order 2 | Order 3-5 |
|---------|--------|---------|---------|-----------|
| `export netgen "f.vol" order N` | Netgen .vol (+ .vol.json) | Yes | Yes | Yes |
| `export gmsh "f.msh" version 2` | Gmsh v2.2 (.msh) | Yes | Yes | Yes |
| `export gmsh "f.msh" version 4` | Gmsh v4.1 (.msh) | Yes | Yes | Yes |
| `export radia_nastran "f.bdf"` | Nastran BDF (.bdf) | Yes | Yes | Yes |
| `export vtk "f.vtk"` | VTK Legacy (.vtk) | Yes | Yes | Yes |
| `extract_curved_mesh(cubit, order=N)` | ngsolve.Mesh (in-memory) | No (>=2) | Yes | Yes |

All formats use NetgenCurver (compact_netgen BuildCurvedElements) for curving.
No fallback to HighOrderMesh (removed).

## API

```python
# Netgen .vol (recommended for NGSolve FEM computation)
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
# -> produces mesh.vol + mesh.vol.json (CAD reference values)

# Gmsh v2.2 (for GMSH visualization)
cubit.cmd('export gmsh "mesh.msh" order 2 version 2 overwrite')

# Python path (reference, slower)
from cubit_mesh_export import extract_curved_mesh
ng_mesh = extract_curved_mesh(cubit, order=3)
ng_mesh.Save("mesh.vol")
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
# Gmsh v2.2 Export

```python
export_Gmesh(cubit, FileName)  # version="2.2" is default
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .msh file path |

## Supported Elements

| 1st Order | 2nd Order | Gmsh Type |
|-----------|-----------|-----------|
| TET4 | TET10 | 4 / 11 |
| HEX8 | HEX20 | 5 / 17 |
| WEDGE6 | WEDGE15 | 6 / 18 |
| PYRAMID5 | PYRAMID13 | 7 / 19 |
| TRI3 | TRI6 | 2 / 9 |
| QUAD4 | QUAD8 | 3 / 16 |
| EDGE2 | EDGE3 | 1 / 8 |

## Use Cases

- **NGSolve/Netgen integration**: Export .vol via `export netgen` (recommended)
- **GMSH visualization**: View mesh in GMSH GUI
- **2nd order elements**: Good accuracy for simple curving workflows

```python
# Recommended: export netgen .vol (any order, best accuracy)
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
from ngsolve import Mesh
mesh = Mesh("mesh.vol")

# Alternative: Gmsh v2.2 (2nd order only, for GMSH visualization)
cubit.cmd("block 1 element type tetra10")
cubit.cmd('export gmsh "mesh.msh" overwrite')
```
"""

EXPORT_GMSH_V4 = """
# Gmsh v4.1 Export

```python
export_Gmesh(cubit, FileName, version="4.1", DIM="auto")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .msh file path |
| `DIM` | str | "auto" | "auto", "2D", or "3D" |

## DIM Options

| Value | Description |
|-------|-------------|
| `"auto"` | Auto-detect (3D if volume elements exist) |
| `"2D"` | Orient normals to +z, z-coordinates set to 0 |
| `"3D"` | No normal orientation |

## v2.2 vs v4.1: Format Version Policy

| Direction | Format | Purpose | Tool |
|-----------|--------|---------|------|
| **Input** (-> NGSolve) | **.vol** | Mesh import into NGSolve | `export netgen "mesh.vol" order N` -> `Mesh("mesh.vol")` |
| **Output** (NGSolve ->) | **v4.1** | Field visualization in GMSH | `GmshPostExport.write()` -> GMSH GUI |
| **Output** (NGSolve ->) | **v2.2** | High-order mesh exchange | `GmshPostExport.write_v22()` |

| Feature | v2.2 | v4.1 |
|---------|------|------|
| $Entities section | No | Yes |
| DIM parameter | No | Yes |
| Post-processing (NodeData) | Basic | **Recommended** |
| Physical Groups | Basic | Structured |
| High-order elements (Tri6, Tet10, Tri10, ...) | **Yes (any order)** | **Yes (any order)** |
| Element type codes | Same as v4.1 | Same as v2.2 |

**Key rule**: Element type codes are identical in both versions.
High-order elements (Tri10=21 for order 3, Tri15=23 for order 4, etc.) work in both.

### GmshPostExport Methods

| Method | Format | High-order | Use case |
|--------|--------|-----------|----------|
| `write(filename)` | v4.1 | Yes (any order) | Field visualization in GMSH GUI |
| `write_v22(filename)` | v2.2 | Yes (any order) | High-order mesh exchange |
| `write_mesh(filename)` | v4.1 | Yes (any order) | Mesh only (no field data) |

## When to Use Which

- **`export netgen "mesh.vol"`**: For NGSolve FEM computation (recommended, any order)
- **`GmshPostExport.write()`**: For field visualization in GMSH GUI (v4.1)
- **`export_Gmesh(version="2.2")`**: For direct Cubit mesh output to GMSH (2nd order max)
"""

EXPORT_CURVED = """
# Curved Mesh Export (export_NGSolveCurvedMesh)

```python
from cubit_mesh_export import extract_curved_mesh
from ngsolve import Mesh
mesh = Mesh(extract_curved_mesh(cubit, order=3, surface_only=False, split_quads=False))
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `order` | int | 3 | Polynomial order for mesh curving |
| `surface_only` | bool | False | Export surface elements only (for BEM) |
| `split_quads` | bool | False | Split quad elements into triangles |

**Returns**: `netgen.meshing.Mesh` object (**already curved**). Wrap with `Mesh()` for NGSolve use.

## How It Works

Uses CallbackGeometry (upstream Netgen) to delegate surface projection to
Cubit's ACIS kernel:

1. Reads mesh topology (nodes, elements) from Cubit blocks
2. Creates 1st order netgen.meshing.Mesh
3. Registers CallbackGeometry with Cubit ACIS surface projection
4. Calls mesh.Curve(order) using the CallbackGeometry
5. Returns the curved ngsolve.Mesh

No STEP files, no OCC geometry, no SetGeomInfo needed.

## Key Design Decision

This function exports **1st order elements** internally, then curves
them to the requested order using ACIS surface projection. Even if
blocks contain TET10, only corner nodes are used initially. High-order
nodes are placed exactly on ACIS CAD surfaces by mesh.Curve(order).

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
`export_NGSolveCurvedMesh()` replaces all of them.
"""

EXPORT_NASTRAN = """
# Nastran BDF Export

```python
export_nastran(cubit, FileName, DIM="3D", PYRAM=True)
# Aliases: export_Nastran
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .bdf file path |
| `DIM` | str | "3D" | "3D" or "2D" |
| `PYRAM` | bool | True | Pyramid element handling |

## DIM Options

| Value | Elements |
|-------|----------|
| `"3D"` | CTETRA, CHEXA, CPENTA, CPYRAM |
| `"2D"` | CTRIA3, CQUAD4 (normals oriented to +z) |

## PYRAM Options

| Value | Output | Use Case |
|-------|--------|----------|
| `True` | CPYRAM (5-node) | Standard Nastran |
| `False` | Degenerate CHEXA (8-node with repeated nodes) | JMAG compatibility |

## Limitation

**1st order elements only.** Uses `get_connectivity()`.
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
| NGSolve FEM (any order) | `export netgen "f.vol" order N` | Arbitrary order via ACIS CallbackGeometry |
| NGSolve FEM (any order, Python) | `export_NGSolveCurvedMesh()` | In-memory, same curving as Path A |
| GMSH visualization | `export_Gmesh(version="2.2")` | GMSH GUI viewing |
| JMAG solver | `export_nastran()` | PYRAM=False for degenerate hex |
| Cubit-native archival | `export_exodus()` | Full fidelity, all features |

## Feature Comparison

| Feature | curved | gmsh_v2 | gmsh_v4 | nastran | exodus |
|---------|--------|---------|---------|---------|--------|
| 1st order | Yes | Yes | Yes | Yes | Yes |
| 2nd order | Yes | Yes | Yes | No | Yes |
| 3rd+ order | Yes | No | No | No | Yes |
| In-memory | Yes | No | No | No | No |
| BlockID metadata | N/A | Yes | Yes | Yes | Yes |
| 2D support | No | No | Yes | Yes | No |

## export_NGSolveCurvedMesh vs Gmsh for NGSolve

| Aspect | export_NGSolveCurvedMesh() | Gmsh v2.2 (export_Gmesh) |
|--------|----------------|---------------------------|
| Max order | Unlimited | 2nd order |
| Accuracy at order 2 | ~0.003% | ~0.001% |
| Accuracy at order 3+ | ~0.0004% | N/A |
| Complexity | Low (single function call) | Low (block element type) |
| Geometry needed | ACIS (automatic) | No |
| Best for | Any order FEM/BEM | Standard 2nd order FEM |
"""


EXPORT_DECISION_GUIDE = """
# Export Format Decision Guide

## "I want to use NGSolve / Netgen for FEM"

-> Use `export_NGSolveCurvedMesh()` (recommended for any order):
  ```python
  from cubit_mesh_export import extract_curved_mesh
  from ngsolve import Mesh
  cubit.cmd("block 1 add tet all")
  cubit.cmd("block 2 add tri all")
  mesh = Mesh(extract_curved_mesh(cubit, order=3))
  ```
  Works for ANY geometry shape. No STEP files, no OCC, no SetGeomInfo.

- **Alternative for 2nd order only** -> Use `export_Gmesh(version="2.2")` for GMSH visualization:
  ```python
  cubit.cmd("block 1 add tet all")
  cubit.cmd("block 1 element type tetra10")
  cubit.cmd("block 2 add tri all")
  cubit.cmd("block 2 element type tri6")
  cubit.cmd('export gmsh "mesh.msh" overwrite')
  # For GMSH visualization only; for NGSolve use export netgen instead
  ```

## "I need structural FEA (Nastran / JMAG)"

-> Use `export_nastran()`
- 3D: `export_nastran(cubit, "mesh.bdf", DIM="3D")`
- 2D: `export_nastran(cubit, "mesh.bdf", DIM="2D")`
- **Note**: 1st order elements only

### JMAG-specific: Pyramid Element Problem

JMAG **cannot read standard CPYRAM** (5-node pyramid) elements. When Cubit
generates mixed hex-tet meshes, pyramid transition elements appear at
the interface. Two solutions:

1. **Use PYRAM=False** (recommended): Writes pyramids as degenerate CHEXA
   ```python
   cubit.cmd('export radia_nastran "mesh.bdf" nopyramid overwrite')
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

-> Use `export_Gmesh(version="2.2")` for GMSH GUI viewing
- Or `export_Gmesh(version="4.1")` for full v4.1 with $Entities
- **Note**: For NGSolve FEM, use `export netgen "mesh.vol"` instead

## Performance & Feature Summary

| Format | File Size | Max Order | 2D | 3D | In-Memory |
|--------|-----------|-----------|----|----|-----------|
| curved | N/A | Unlimited | No | Yes | Yes |
| gmsh_v2 | Medium | 2nd | No | Yes | No |
| gmsh_v4 | Medium | 2nd | Yes | Yes | No |
| nastran | Medium | 2nd | Yes | Yes | No |
| exodus | Medium | All | No | Yes | No |

## IMPORTANT: `export radia_nastran` (NOT `export nastran`)

Cubit has a **built-in** `export nastran` command (e.g., `export nastran "f.bdf" overwrite everything`).
Radia's Nastran export uses a DIFFERENT command name to avoid conflict:

```python
# CORRECT: Radia's export (supports order 2, nopyramid, block labels)
cubit.cmd('export radia_nastran "mesh.bdf" order 2 dimension 3 overwrite')

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
#     from radia.radia_coil_builder import CoilBuilder
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
cmake -DNETGEN_SRC_DIR=/path/to/netgen/source ...  # upstream or fork
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
                                       (C++ .ccl)        (Python)
Export Mesh:                           Solve:
  Netgen Vol (.vol)...                   Radia-NGSolve...
  GMSH...                               Generate Coil...
  Nastran BDF...                         --------
  VTK...                                 Reload Panels
  --------
  Mesh Evaluation...
```

- **Export Mesh**: C++ .ccl component (Qt5 dialogs, APREPRO commands)
- **Solve**: Python register_toolbar.py (subprocess to external Python 3.12)
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
		"netgen": EXPORT_CURVED,  # Alias: old name redirects to export_NGSolveCurvedMesh
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
