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
| `export nastran "f.bdf"` | Nastran BDF (.bdf) | Yes | Yes | Yes |
| `export vtk "f.vtk"` | VTK Legacy (.vtk) | Yes | Yes | Yes |
| `export meg "f.meg"` | MEG/ELF (.meg) | Yes | No | No |
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

- **NGSolve/Netgen integration**: `ReadGmsh()` supports v2.2 format
- **2nd order elements for NGSolve**: Best approach for simple curving

```python
# NGSolve integration
from netgen.read_gmsh import ReadGmsh
from ngsolve import Mesh

cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order
cubit.cmd('export gmsh "mesh.msh" overwrite')
mesh = Mesh(ReadGmsh("mesh.msh"))  # ~0.001% volume error
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
| **Input** (-> NGSolve) | **v2.2** | Mesh import into NGSolve | `GmshPostExport.write_v22()` or `export_Gmesh(version="2.2")` -> `ReadGmsh()` |
| **Output** (NGSolve ->) | **v4.1** | Field visualization in GMSH | `GmshPostExport.write()` -> GMSH GUI |
| **Output** (NGSolve ->) | **v2.2** | High-order mesh exchange | `GmshPostExport.write_v22()` |

| Feature | v2.2 | v4.1 |
|---------|------|------|
| $Entities section | No | Yes |
| DIM parameter | No | Yes |
| NGSolve ReadGmsh() | **Supported** | Not supported |
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
| `write_v22(filename)` | v2.2 | Yes (any order) | NGSolve ReadGmsh() import, mesh exchange |
| `write_mesh(filename)` | v4.1 | Yes (any order) | Mesh only (no field data) |

## When to Use Which

- **`GmshPostExport.write_v22()`**: For high-order mesh export to v2.2 (ReadGmsh() compatible, any order)
- **`GmshPostExport.write()`**: For field visualization in GMSH GUI (v4.1)
- **`export_Gmesh(version="2.2")`**: For direct Cubit mesh output (2nd order max, no curving)
- Do NOT use v4.1 for NGSolve import (ReadGmsh cannot parse it)
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

Uses CallbackGeometry (netgen fork) to delegate surface projection to
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
For simple 2nd order without geometry, consider Gmsh export with `ReadGmsh()`.

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

EXPORT_MEG = """
# MEG Export (ELF/MAGIC)

```python
export_meg(cubit, FileName, DIM='T', MGR2=None)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .meg file path |
| `DIM` | str | 'T' | 'T' (3D), 'K' (2D planar), 'R' (axisymmetric) |
| `MGR2` | list | None | Spatial nodes [[x,y,z], ...] |

## DIM Options

| Value | Description | Coordinate System |
|-------|-------------|-------------------|
| `'T'` | 3D | X, Y, Z |
| `'K'` | 2D Planar | X, Y (Z=0) |
| `'R'` | Axisymmetric | R (X), Z |

## Block Name → ELF Element Name (Critical!)

In MEG format, **block names directly become ELF element type identifiers**.
The module takes the first 4 characters of the block name and appends the DIM
letter to form the ELF element name.

### How It Works Internally

```python
name = cubit.get_exodus_entity_name("block", block_id)
# Element line: "{name[0:4]}{DIM} {elem_id} 0 {block_id} {node1} {node2} ..."
```

For example, if block name is `"MMB8"` and DIM is `'T'`:
→ Element name in MEG file is `MMB8T`

### Standard ELF Element Names

| DIM | Tri (3-node) | Quad (4-node) | Tet (4-node) | Wedge (6-node) | Hex (8-node) |
|-----|-------------|--------------|-------------|---------------|-------------|
| `'T'` | - | - | `MMB4T` | `MMB6T` | `MMB8T` |
| `'K'` | `MMB3K` | `MMB4K` | - | - | - |
| `'R'` | `MMB3R` | `MMB4R` | - | - | - |

### Setting Block Names for MEG

**The block name MUST start with the correct ELF prefix** (first 4 characters):

```python
# 3D tet mesh for ELF
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'MMB4'")    # → MMB4T in MEG file

# 3D hex mesh for ELF
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'MMB8'")    # → MMB8T in MEG file

# 2D planar tri mesh
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3'")    # → MMB3K in MEG file (DIM='K')

# Axisymmetric tri mesh
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3'")    # → MMB3R in MEG file (DIM='R')
```

### Common Mistake

```python
# WRONG: Generic block name → ELF cannot identify element type
cubit.cmd("block 1 name 'domain'")   # → "doma" + "T" = "domaT" (not valid!)

# RIGHT: Use ELF naming convention
cubit.cmd("block 1 name 'MMB4'")     # → "MMB4" + "T" = "MMB4T" (valid!)
```

## Limitation

**1st order elements only.** Uses `get_connectivity()` (corner nodes only).
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
| NGSolve FEM (any order) | `export_NGSolveCurvedMesh()` | Arbitrary order via ACIS CallbackGeometry |
| NGSolve FEM (2nd order, simple) | `export_Gmesh(version="2.2")` | ReadGmsh() supports v2.2, ~0.001% accuracy |
| Gmsh solver | `export_Gmesh(version="4.1")` | Full v4.1 with $Entities |
| JMAG solver | `export_nastran()` | PYRAM=False for degenerate hex |
| ELF/MAGIC solver | `export_meg()` | Native ELF element names |
| Cubit-native archival | `export_exodus()` | Full fidelity, all features |

## Feature Comparison

| Feature | curved | gmsh_v2 | gmsh_v4 | nastran | meg | exodus |
|---------|--------|---------|---------|---------|-----|--------|
| 1st order | Yes | Yes | Yes | Yes | Yes | Yes |
| 2nd order | Yes | Yes | Yes | No | No | Yes |
| 3rd+ order | Yes | No | No | No | No | Yes |
| In-memory | Yes | No | No | No | No | No |
| BlockID metadata | N/A | Yes | Yes | Yes | Yes | Yes |
| 2D support | No | No | Yes | Yes | Yes | No |

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

- **Alternative for 2nd order only** -> Use `export_Gmesh(version="2.2")`:
  ```python
  cubit.cmd("block 1 add tet all")
  cubit.cmd("block 1 element type tetra10")
  cubit.cmd("block 2 add tri all")
  cubit.cmd("block 2 element type tri6")
  cubit.cmd('export gmsh "mesh.msh" overwrite')
  mesh = Mesh(ReadGmsh("mesh.msh"))
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
   cubit.cmd('export nastran "mesh.bdf" nopyramid overwrite')
   ```

2. **Use pure tet mesh**: Avoids pyramids entirely
   ```python
   cubit.cmd("volume all scheme tetmesh")
   # ... mesh and export ...
   ```

**When to give up on hex**: If Cubit cannot webcut/decompose the geometry
into sweepable sub-volumes, switch to `tetmesh`. Forcing hex on complex
geometry leads to poor quality or failed meshing.

## "I need ELF/MAGIC electromagnetic solver"

-> Use `export_meg()`
- 3D: `export_meg(cubit, "mesh.meg", DIM='T')`
- 2D planar: `export_meg(cubit, "mesh.meg", DIM='K')`
- Axisymmetric: `export_meg(cubit, "mesh.meg", DIM='R')`
- Block names become ELF element type names (MMB4T, MMB8T, etc.)
- **Note**: 1st order elements only

## "I need multi-physics (MOOSE, Sierra)"

-> Use `export_exodus()`
- Supports all element types and orders
- Full block definitions preserved
- Large model support: `export_exodus(cubit, "mesh.exo", large_model=True)`
- Cubit-native format with highest fidelity

## "I want Gmsh solver integration"

-> Use `export_Gmesh(version="4.1")` for full v4.1 format
- Includes $Entities section
- DIM parameter: `"2D"` or `"3D"`
- **Note**: For NGSolve, use `export_Gmesh(version="2.2")` instead (ReadGmsh uses v2.2)

## Performance & Feature Summary

| Format | File Size | Max Order | 2D | 3D | In-Memory |
|--------|-----------|-----------|----|----|-----------|
| curved | N/A | Unlimited | No | Yes | Yes |
| gmsh_v2 | Medium | 2nd | No | Yes | No |
| gmsh_v4 | Medium | 2nd | Yes | Yes | No |
| nastran | Medium | 1st | Yes | Yes | No |
| meg | Small | 1st | Yes | Yes | No |
| exodus | Medium | All | No | Yes | No |
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
		"meg": EXPORT_MEG,
		"exodus": EXPORT_EXODUS,
		"comparison": EXPORT_COMPARISON,
		"decision_guide": EXPORT_DECISION_GUIDE,
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
