"""
Cubit Python scripting knowledge base for mesh export.

Covers block registration, element order control, mesh schemes,
STEP exchange, initialization, and common mistakes.

Includes content migrated from Radia project's cubit_meshing_knowledge.py.
"""

CUBIT_OVERVIEW = """
# Cubit Mesh Generation Overview

Coreform Cubit provides structured and unstructured meshing with export to
multiple formats via the `cubit_mesh_export` module.

## Typical Workflow

```
1. Create geometry (Cubit commands or STEP import)
2. Set mesh scheme (tetmesh, map, sweep, etc.)
3. Set element size or interval count
4. Generate mesh
5. Register blocks
6. Export to desired format
```

## Element Types

| Element | Cubit Type | Nodes (1st) | Nodes (2nd) | Best For |
|---------|-----------|-------------|-------------|----------|
| Tet | TET | 4 | 10 | Complex geometry, auto-mesh |
| Hex | HEX | 8 | 20/27 | Structured grids, high accuracy |
| Wedge | WEDGE | 6 | 15 | Transition elements |
| Pyramid | PYRAMID | 5 | 13 | Hex-tet transition |
| Tri | TRI | 3 | 6 | Surface mesh |
| Quad | QUAD | 4 | 8/9 | Structured surface mesh |

## Why Cubit for FEM?

- **Structured hex meshing**: Higher accuracy per element
- **Interval control**: Precise element count along curves
- **Sweep meshing**: Efficient for extruded geometries
- **Multi-block**: Different mesh densities per region
- **Quality control**: Built-in mesh quality metrics
"""

CUBIT_BLOCKS = """
# Block Registration

## Why Blocks are Required

All `cubit_mesh_export` functions read mesh data from blocks.
**No blocks = no export.** This is the most common mistake.

## Mesh Element Blocks (Recommended)

```python
# 3D domain elements
cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')

# 2D boundary elements
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')
```

## Geometry Blocks (Alternative)

Blocks can contain geometry instead of mesh elements:

```python
cubit.cmd("block 1 add volume 1")    # All 3D elements in volume 1
cubit.cmd("block 2 add surface 1")   # All 2D elements on surface 1
```

| Block Contains | Elements Returned |
|----------------|-------------------|
| Volume | tet, hex, wedge, pyramid (3D only) |
| Surface | tri, quad (2D only) |
| Curve | edge (1D only) |
| Vertex | node (0D only) |

## Important: Element Order and Geometry Blocks

For 2nd order conversion, you MUST use mesh element blocks:

```python
# This works:
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")   # Converts to 2nd order

# This does NOT work for 2nd order:
cubit.cmd("block 1 add volume 1")
cubit.cmd("block 1 element type tetra10")   # No effect on geometry blocks!
```

## Multiple Blocks

Use separate blocks for different materials or boundary conditions:

```python
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd('block 1 name "iron"')
cubit.cmd("block 2 add tet all in volume 2")
cubit.cmd('block 2 name "air"')
cubit.cmd("block 3 add tri all in surface 1")
cubit.cmd('block 3 name "dirichlet"')
```

## Mixed Element Warning

If a block contains multiple 3D element types (e.g., tet + hex), a warning
is displayed with 2nd order conversion instructions. For mixed blocks,
separate `element type` commands are needed:

```python
cubit.cmd("block 1 element type hex20")
cubit.cmd("block 1 element type tetra10")
```
"""

CUBIT_ELEMENT_ORDER = """
# Element Order Control

## Creating 2nd Order Elements

1. Create mesh (generates 1st order by default)
2. Add elements to block
3. Set element type to convert

```python
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")   # Now 2nd order
```

## Element Type Commands

| 1st Order | 2nd Order | Command |
|-----------|-----------|---------|
| TET4 | TET10 | `block X element type tetra10` |
| HEX8 | HEX20 | `block X element type hex20` |
| HEX8 | HEX27 | `block X element type hex27` |
| WEDGE6 | WEDGE15 | `block X element type wedge15` |
| PYRAMID5 | PYRAMID13 | `block X element type pyramid13` |
| TRI3 | TRI6 | `block X element type tri6` |
| QUAD4 | QUAD8 | `block X element type quad8` |
| QUAD4 | QUAD9 | `block X element type quad9` |
| EDGE2 | EDGE3 | `block X element type bar3` |

## get_connectivity vs get_expanded_connectivity

```python
cubit.cmd("block 1 element type tetra10")
tet_id = cubit.get_block_tets(1)[0]

# get_connectivity: ALWAYS returns corner nodes only
nodes_1st = cubit.get_connectivity("tet", tet_id)
print(len(nodes_1st))   # Always 4

# get_expanded_connectivity: Returns ALL nodes
nodes_all = cubit.get_expanded_connectivity("tet", tet_id)
print(len(nodes_all))   # 10 for TET10, 4 for TET4
```

## Which Export Functions Use Which API

| Export Function | API | Respects Element Type |
|----------------|-----|----------------------|
| `export_gmsh_v2()` | `get_expanded_connectivity()` | Yes |
| `export_gmsh_v4()` | `get_expanded_connectivity()` | Yes |
| `export_exodus()` | `get_expanded_connectivity()` | Yes |
| `export_NGSolveCurvedMesh()` | `get_connectivity()` | No (curves via ACIS) |
| `export_nastran()` | `get_connectivity()` | No (1st order only) |
| `export_meg()` | `get_connectivity()` | No (1st order only) |

## Design: Why export_NGSolveCurvedMesh Uses 1st Order Internally

`export_NGSolveCurvedMesh()` exports 1st order elements then curves them via ACIS
CallbackGeometry. This enables arbitrary orders (2, 3, 4, 5, ...) with
exact placement on ACIS CAD surfaces. Cubit's role is topology (1st order
mesh) AND surface projection; Netgen's role is p-refinement via Curve().
"""

CUBIT_MESH_SCHEMES = """
# Mesh Schemes

## Tetrahedral (Auto)

```python
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")       # Target edge length
cubit.cmd("mesh volume all")
```

Best for: Complex geometry, automatic meshing.

## Map (Structured Quad/Hex)

```python
cubit.cmd("surface 1 scheme map")
cubit.cmd("curve 1 interval 10")
cubit.cmd("mesh surface 1")
```

Best for: Simple geometries with mappable topology.

## Sweep (Extruded)

```python
cubit.cmd("volume 1 scheme sweep source surface 1 target surface 2")
cubit.cmd("mesh volume 1")
```

Best for: Extruded or revolution geometries.

## Interval Control

```python
# Set number of elements along a curve
cubit.cmd("curve 1 interval 10")

# Set element size
cubit.cmd("volume 1 size 0.05")

# Gradation (size varies)
cubit.cmd("volume 1 sizing function type skeleton")
```

## Mesh Quality

```python
cubit.cmd("quality volume all shape")
cubit.cmd("quality volume all aspect ratio")
cubit.cmd("quality volume all jacobian")
```

## Hex Meshing Resolution Guidelines

| Resolution | Interval | Typical Elements | Use Case |
|-----------|----------|-----------------|----------|
| Coarse | 1-2 | 50-200 | Quick validation |
| Medium | 3-5 | 500-2000 | Standard analysis |
| Fine | 6-10 | 2000-10000 | High accuracy |
| Very Fine | 10+ | 10000+ | Convergence study |
"""

CUBIT_STEP_EXCHANGE = """
# STEP Import/Export

## Role of STEP in the Workflow

STEP files are used for geometry exchange between tools. Note that
for high-order curving with `export_NGSolveCurvedMesh()`, STEP files are NOT needed
(export_NGSolveCurvedMesh uses Cubit's ACIS kernel directly via CallbackGeometry).

STEP is still useful for:
- Importing geometry from external CAD tools
- Sharing geometry between Cubit and other applications

## Exporting STEP from Cubit

```python
cubit.cmd('export step "geometry.step" overwrite')
```

## Importing STEP into Cubit

```python
# With healing (recommended for simple shapes)
cubit.cmd('import step "geometry.step" heal')

# Without healing (recommended for name-based workflow)
cubit.cmd('import step "geometry.step" noheal')
```

## heal vs noheal

| Option | Behavior | Use When |
|--------|----------|----------|
| `heal` | Repairs topology, merges surfaces | Simple shapes, STEP reimport workflow |
| `noheal` | Preserves original topology and names | Name-based workflow, preserving OCC face names |

## STEP Reimport Pattern (Legacy, No Longer Needed)

**NOTE**: The STEP reimport pattern is NO LONGER NEEDED for high-order curving.
`export_NGSolveCurvedMesh()` uses Cubit's ACIS kernel directly, so no STEP exchange
or OCC topology matching is required.

STEP reimport is only relevant if you need to match Cubit topology to
an external OCC-based tool (not for export_NGSolveCurvedMesh).
"""

CUBIT_INITIALIZATION = """
# Cubit Python Initialization

## Standard Boilerplate

```python
import sys, os

# Add Cubit to Python path (uses CUBIT_PATH env var if set)
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Add cubit_mesh_export to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cubit_mesh_export
```

## With NGSolve (System Python + CUBIT_PATH)

By using **system Python** with the `CUBIT_PATH` environment variable, scripts can
access **both** the Cubit API and NGSolve/Netgen simultaneously. This is essential
for the `export_NGSolveCurvedMesh()` workflow, because Cubit's bundled Python
cannot import ngsolve.

**CRITICAL: Import NGSolve BEFORE Cubit.** Cubit bundles its own VTK DLLs which
conflict with NGSolve's Netgen library. If Cubit is imported first, NGSolve fails
with `ImportError: initialization failed` on `from netgen import libngpy`.

```bash
# Set CUBIT_PATH once (e.g., in your environment or before running)
set CUBIT_PATH="C:/Program Files/Coreform Cubit 2025.3/bin"
python my_script.py
```

```python
import sys, os

# Step 1: Import NGSolve FIRST (before Cubit!) to avoid DLL conflicts
import ngsolve
from netgen.occ import OCCGeometry
from ngsolve import Mesh
import cubit_mesh_export

# Step 2: Then import Cubit via CUBIT_PATH
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])
```

**Wrong order** (causes DLL conflict):
```python
# DON'T DO THIS — Cubit's VTK DLLs will break Netgen
sys.path.append(cubit_path)
import cubit                  # Loads Cubit's bundled DLLs
import ngsolve                # FAILS: Netgen can't initialize
```

## Reset Between Operations

```python
cubit.cmd("reset")   # Clear all geometry and mesh
```

## Batch Mode Flags

| Flag | Description |
|------|-------------|
| `-nojournal` | Don't create journal file |
| `-batch` | No GUI, batch mode |
| `-nographics` | Disable graphics (faster) |

## CUBIT_PATH Environment Variable

Set `CUBIT_PATH` to avoid hardcoding the Cubit installation path in scripts:

```bash
# Windows
set CUBIT_PATH="C:/Program Files/Coreform Cubit 2025.3/bin"

# Linux/Mac
export CUBIT_PATH="/opt/Coreform-Cubit-2025.3/bin"
```

**Key benefit**: System Python with `CUBIT_PATH` can access both the Cubit API
and NGSolve/Netgen, enabling the `export_NGSolveCurvedMesh()` high-order curving workflow.
Cubit's bundled Python cannot import ngsolve.
"""

CUBIT_MIXED_ELEMENTS = """
# Mixed Element Types

## Warning System

The module automatically warns when blocks contain multiple 3D element types:

```
WARNING: Block 1 ('mixed') contains multiple 3D element types: hex, tet
  To convert all elements to 2nd order, you need to issue separate commands:
    block 1 element type hex20
    block 1 element type tetra10
```

## Why It Matters

A single `element type` command only converts one type:
```python
cubit.cmd("block 1 element type tetra10")  # Only converts tets, not hexes
```

## Solution

Issue separate commands for each element type:
```python
cubit.cmd("block 1 element type hex20")
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 1 element type wedge15")
```
"""

CUBIT_COMMON_MISTAKES = """
# Common Cubit Scripting Mistakes

## 1. CRITICAL: No Block Registration Before Export

```python
# WRONG: Export without blocks -> empty file!
cubit.cmd("mesh volume 1")
cubit_mesh_export.export_gmsh_v2(cubit, "mesh.msh")

# RIGHT: Register blocks first
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")
cubit_mesh_export.export_gmsh_v2(cubit, "mesh.msh")
```

## 2. CRITICAL: Geometry Block with 2nd Order Conversion

```python
# WRONG: Geometry block, element type has no effect
cubit.cmd("block 1 add volume 1")
cubit.cmd("block 1 element type tetra10")   # Does nothing!

# RIGHT: Mesh element block
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 element type tetra10")   # Works!
```

## 3. HIGH: Using get_connectivity for 2nd Order

```python
# WRONG: Always returns 4 nodes for tet
nodes = cubit.get_connectivity("tet", tet_id)        # 4 nodes

# RIGHT: Returns 10 nodes for TET10
nodes = cubit.get_expanded_connectivity("tet", tet_id)  # 10 nodes
```

## 4. HIGH: Element Type Before Adding to Block

```python
# WRONG: Set type before add
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 1 add tet all")

# RIGHT: Add first, then set type
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")
```

## 5. HIGH: Using Deleted APIs (export_netgen, SetGeomInfo)

```python
# WRONG: These functions no longer exist
geo = OCCGeometry("cyl.step")
ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=0.5, height=2.0)

# RIGHT: Use export_NGSolveCurvedMesh() — single function, no STEP, no OCC
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3)
```

## 6. MODERATE: Missing Boundary Element Block

```python
# INCOMPLETE: Volume block only
cubit.cmd("block 1 add tet all")

# COMPLETE: Both volume and surface blocks
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")   # Needed for boundary conditions
```

## 7. CRITICAL: Block Registration Order for Source/Sink (BEM Inductance)

When using `set duplicate block elements on` with source/sink blocks for BEM
inductance extraction, register source/sink BEFORE the boundary block.
If boundary (all tris) is registered first, subsequent source/sink blocks
that are subsets of boundary may fail to register elements.

```python
# WRONG: boundary first -> source/sink subset registration may fail
cubit.cmd("set duplicate block elements on")
cubit.cmd("block 1 add tri all")
cubit.cmd('block 1 name "boundary"')
cubit.cmd("block 2 add tri in surface 3")   # May get 0 elements!
cubit.cmd('block 2 name "source"')

# RIGHT: source/sink first, then boundary (superset)
cubit.cmd("set duplicate block elements on")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "conductor"')
cubit.cmd("block 2 add tri in surface 3")
cubit.cmd('block 2 name "source"')
cubit.cmd("block 3 add tri in surface 5")
cubit.cmd('block 3 name "sink"')
cubit.cmd("block 4 add tri all")            # Superset last
cubit.cmd('block 4 name "boundary"')
```

In export_NGSolveCurvedMesh, the last block wins for bcname priority.
So even with `set duplicate block elements on`, the block registration
order determines which label an element gets:
- source/sink blocks registered AFTER boundary -> source/sink wins (correct)
- source/sink blocks registered BEFORE boundary -> boundary wins (wrong)

Best practice: register source/sink LAST (highest block ID), or
register boundary block first then source/sink blocks after.

## 8. MODERATE: Hardcoded Absolute Paths

```python
# WRONG: Breaks on other machines
sys.path.insert(0, "path/to/Radia/src/radia")

# RIGHT: Use relative paths
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
```

## 9. HIGH: Mixed Hex-Tet Mesh for JMAG Without PYRAM=False

```python
# WRONG: JMAG cannot read CPYRAM elements
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 add pyramid all")
cubit_mesh_export.export_nastran(cubit, "mesh.bdf")  # Pyramids as CPYRAM!

# RIGHT: Use PYRAM=False for JMAG
cubit_mesh_export.export_nastran(cubit, "mesh.bdf", PYRAM=False)
# Or use pure tet mesh to avoid pyramids entirely
```

## 10. LOW: Using 'modify mesh volume X order 2'

```python
# WRONG: May not work reliably
cubit.cmd("modify mesh volume 1 order 2")

# RIGHT: Use block element type
cubit.cmd("block 1 element type tetra10")
```
"""

CUBIT_BLOCKS_ONLY_POLICY = """
# Blocks-Only Policy (No Nodesets or Sidesets)

## Core Rule

This module uses **blocks only** for mesh export. Nodesets and sidesets are **never used**.

## Why?

In Cubit, **only blocks** support element order specification:

```python
cubit.cmd("block 1 element type tetra10")   # Works: converts to 2nd order
# nodesets and sidesets have NO equivalent command
```

To maintain consistent element order control across all export formats
(Gmsh, VTK, Nastran, MEG, Netgen, Exodus), we standardize on blocks exclusively.

## Boundary Conditions with Blocks

Boundary conditions are defined using **separate blocks for boundary elements**,
not using nodesets or sidesets:

```python
# Domain elements (3D)
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 1 element type tetra10")    # 2nd order

# Boundary elements (2D)
cubit.cmd("block 2 add tri all in surface 1")
cubit.cmd('block 2 name "boundary"')
cubit.cmd("block 2 element type tri6")       # 2nd order

# Multiple boundary conditions
cubit.cmd("block 3 add tri all in surface 2")
cubit.cmd('block 3 name "dirichlet"')
cubit.cmd("block 4 add tri all in surface 3")
cubit.cmd('block 4 name "neumann"')
```

## Common Mistakes

```python
# WRONG: Using nodesets (this module ignores them)
cubit.cmd("nodeset 1 add surface 1")     # NOT supported

# WRONG: Using sidesets (this module ignores them)
cubit.cmd("sideset 1 add surface 1")     # NOT supported

# RIGHT: Use blocks for everything
cubit.cmd("block 2 add tri all in surface 1")
cubit.cmd('block 2 name "boundary"')
```

## Design Rationale

| Feature | Blocks | Nodesets | Sidesets |
|---------|--------|----------|----------|
| Element order control | Yes | No | No |
| Named groups | Yes | Yes | Yes |
| 2nd order conversion | Yes | N/A | N/A |
| Export support | All formats | None | None |
"""

# Migrated from Radia project: hex and tet meshing workflows
CUBIT_HEX_WORKFLOW = """
# Hexahedral Mesh Workflow

## Example: Simple Brick

```python
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Create geometry
cubit.cmd('create brick x 1 y 1 z 1')

# Set mesh intervals for structured hex mesh
cubit.cmd('curve all interval 4')

# Mesh
cubit.cmd('volume all scheme auto')
cubit.cmd('mesh volume all')

# Register blocks
cubit.cmd('block 1 add hex all')
cubit.cmd('block 2 add quad all')
```

## Example: Multi-body Union

```python
# C-type geometry (multiple bricks united)
cubit.cmd('create brick x 0.1 y 0.3 z 0.1')  # Vertical leg 1
cubit.cmd('create brick x 0.3 y 0.1 z 0.1')  # Top bar
cubit.cmd('create brick x 0.1 y 0.3 z 0.1')  # Vertical leg 2

# Position parts
cubit.cmd('move volume 1 x -0.1 y 0')
cubit.cmd('move volume 2 x 0 y 0.1')
cubit.cmd('move volume 3 x 0.1 y 0')

# Union
cubit.cmd('unite volume all')

# Mesh
cubit.cmd('curve all interval 2')
cubit.cmd('volume all scheme auto')
cubit.cmd('mesh volume all')
```

## When to Give Up on Hex and Switch to Tet

Not all geometries can be hex-meshed in Cubit. If webcut and decomposition
fail to produce mappable/sweepable sub-volumes, switch to tetmesh:

```python
# Try hex first
cubit.cmd("volume all scheme auto")
cubit.cmd("mesh volume all")

# If meshing fails or produces poor quality, switch to tet
cubit.cmd("delete mesh")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.1")
cubit.cmd("mesh volume all")
```

**Warning**: Switching from hex to tet may introduce pyramid transition
elements at hex-tet interfaces. Some solvers (notably JMAG) cannot read
standard pyramid elements — use `PYRAM=False` for Nastran export.

## JMAG Compatibility

JMAG requires degenerate hex pyramids, not standard CPYRAM elements:

```python
# For JMAG: use PYRAM=False
cubit_mesh_export.export_nastran(cubit, "mesh.bdf", PYRAM=False)
```

With `PYRAM=False`, pyramid elements are written as degenerate CHEXA
(8-node hex with repeated nodes) instead of CPYRAM (5-node), which JMAG
can read.

## Tips

- Use `curve X interval N` for precise element counts
- `volume all scheme auto` selects best scheme automatically
- For sweepable volumes, Cubit auto-detects and uses sweep
- Check quality: `cubit.cmd('quality volume all shape')`
"""

CUBIT_TET_WORKFLOW = """
# Tetrahedral Mesh Workflow

## Example: Sphere

```python
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Create sphere
cubit.cmd('create sphere radius 0.5')

# Set tetrahedral meshing
cubit.cmd('volume 1 scheme tetmesh')
cubit.cmd('volume 1 size 0.05')    # Target edge length
cubit.cmd('mesh volume 1')

# Register blocks
cubit.cmd('block 1 add tet all')
cubit.cmd('block 2 add tri all')

# Check
n_tets = cubit.get_tet_count()
print(f"Tets: {n_tets}")
```

## Tips for Tetrahedral Meshing

1. **Webcut complex shapes**: Split into convex subdomains for better quality
   ```python
   cubit.cmd('webcut volume 1 with plane xplane')
   ```
2. **Control element size**: `volume N size S` where S is target edge length
3. **Quality check**: `quality volume all shape`
4. **Gradation**: `volume N sizing function type skeleton` for size grading
"""


CUBIT_2D_MESH_WORKFLOW = """
# 2D Mesh Workflow

## Overview

Several export formats support 2D mesh export:
- `export_nastran(cubit, "file.bdf", DIM="2D")`
- `export_meg(cubit, "file.meg", DIM='K')` (planar)
- `export_meg(cubit, "file.meg", DIM='R')` (axisymmetric)
- `export_gmsh_v4(cubit, "file.msh", DIM="2D")`

## Surface Creation

```python
# Planar rectangle
cubit.cmd("create surface rectangle width 2 height 2 zplane")

# Circle
cubit.cmd("create surface circle radius 1 zplane")

# From imported geometry
cubit.cmd('import step "plate.step"')
```

## Mesh Schemes for Surfaces

| Scheme | Command | Element Type | Best For |
|--------|---------|-------------|----------|
| trimesh | `surface N scheme trimesh` | TRI3 | Complex 2D shapes |
| paving | `surface N scheme paving` | QUAD4 | Structured-like quad |
| map | `surface N scheme map` | QUAD4 | Simple rectangular domains |

```python
# Tri mesh
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")

# Quad mesh (paving)
cubit.cmd("surface 1 scheme paving")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")
```

## Block Registration for 2D

```python
# Triangular elements
cubit.cmd("block 1 add tri all")
cubit.cmd('block 1 name "plate"')

# Quadrilateral elements
cubit.cmd("block 1 add quad all")
cubit.cmd('block 1 name "plate"')

# Edge elements (boundary)
cubit.cmd("block 2 add edge all in curve all")
cubit.cmd('block 2 name "boundary"')
```

## 2D Export Examples

### Nastran 2D
```python
cubit.cmd("create surface rectangle width 2 height 2 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")
cubit.cmd('block 1 name "plate"')
cubit_mesh_export.export_nastran(cubit, "plate.bdf", DIM="2D")
```

### MEG 2D Planar (DIM='K')
```python
cubit.cmd("create surface rectangle width 2 height 2 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3K'")
cubit_mesh_export.export_meg(cubit, "plate.meg", DIM='K')
```

### MEG Axisymmetric (DIM='R')
```python
cubit.cmd("create surface rectangle width 2 height 2 yplane")
cubit.cmd("surface 1 move 1 0 1")   # Positive X quadrant (r > 0)
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")
cubit.cmd("block 1 name 'MMB3R'")
cubit_mesh_export.export_meg(cubit, "axisym.meg", DIM='R')
```

### Gmsh v4 2D
```python
cubit.cmd("create surface rectangle width 2 height 2 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("surface 1 size 0.3")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")
cubit_mesh_export.export_gmsh_v4(cubit, "plate.msh", DIM="2D")
```

## 2D Normal Orientation

For 2D exports (Nastran DIM="2D", Gmsh v4 DIM="2D", MEG DIM='K'):
- Normals are oriented in +Z direction
- Z coordinates may be set to 0

For MEG axisymmetric (DIM='R'):
- Coordinates are R (mapped from X), Z
- Geometry must be in XZ plane (Y=0) with X > 0
"""

CUBIT_TROUBLESHOOTING = """
# Troubleshooting & Debugging Guide

## Problem: Export Produces Empty File (0 Elements)

**Cause**: No blocks registered before calling export function.
**Fix**:
```python
# MUST register blocks before any export
cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")
```
**Diagnostic**: Check block contents:
```python
print(f"Tets in block 1: {len(cubit.get_block_tets(1))}")
```

## Problem: 2nd Order Elements Not Appearing

**Cause 1**: Geometry block used instead of mesh element block.
```python
# WRONG: Geometry block - element type command has no effect
cubit.cmd("block 1 add volume 1")
cubit.cmd("block 1 element type tetra10")  # Silently does nothing!

# RIGHT: Mesh element block
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")  # Converts to 2nd order
```

**Cause 2**: Element type set before adding elements.
```python
# WRONG: Order matters
cubit.cmd("block 1 element type tetra10")  # Set type first
cubit.cmd("block 1 add tet all")           # Add after - type is reset

# RIGHT
cubit.cmd("block 1 add tet all")           # Add first
cubit.cmd("block 1 element type tetra10")  # Then set type
```

**Cause 3**: Using get_connectivity instead of get_expanded_connectivity.
```python
# get_connectivity ALWAYS returns 1st order nodes
nodes = cubit.get_connectivity("tet", tet_id)  # Always 4 nodes

# get_expanded_connectivity returns actual nodes
nodes = cubit.get_expanded_connectivity("tet", tet_id)  # 10 for TET10
```

**Diagnostic**: Check element order:
```python
tet_id = cubit.get_block_tets(1)[0]
n1 = len(cubit.get_connectivity("tet", tet_id))
n2 = len(cubit.get_expanded_connectivity("tet", tet_id))
print(f"1st order nodes: {n1}, Actual nodes: {n2}")
# If n1 == n2 == 4, still 1st order. If n2 == 10, it's 2nd order.
```

## Problem: Using Old SetGeomInfo / export_netgen APIs

**These APIs have been deleted.** Use `export_NGSolveCurvedMesh()` instead.

```python
# OLD (DELETED - does not work):
# geo = OCCGeometry("geometry.step")
# ngmesh = cubit_mesh_export.export_netgen(cubit, geometry=geo)
# cubit_mesh_export.set_cylinder_geominfo(ngmesh, ...)

# NEW: Single function call
mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=3)
```

## Problem: Mesh Quality Too Low

**Diagnostic**:
```python
cubit.cmd("quality volume all shape")
cubit.cmd("quality volume all aspect ratio")
```

**Fixes**:
```python
# Reduce element size
cubit.cmd("volume 1 size 0.05")

# Smooth after meshing
cubit.cmd("smooth volume 1")

# Webcut complex regions
cubit.cmd("webcut volume 1 with plane xplane")
```

## Problem: "No module named 'cubit'" / "No module named 'netgen'"

**cubit not found**: Add Cubit's bin directory to sys.path using CUBIT_PATH:
```python
import os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
```

Or set the environment variable before running:
```bash
set CUBIT_PATH="C:/Program Files/Coreform Cubit 2025.3/bin"
python my_script.py
```

**netgen not found**: NGSolve must be installed in the Python environment.
Use **system Python** (which has NGSolve installed) with `CUBIT_PATH` to access
both Cubit API and NGSolve simultaneously. Cubit's bundled Python cannot import
ngsolve, so system Python + CUBIT_PATH is the recommended approach for the
`export_NGSolveCurvedMesh()` workflow.

## Problem: NGSolve ImportError After Importing Cubit (DLL Conflict)

**Symptom**: `ImportError: initialization failed` on `from netgen import libngpy`,
or similar DLL-related errors when importing NGSolve after Cubit.

**Cause**: Cubit bundles its own VTK and other shared libraries. When `import cubit`
runs, these DLLs are loaded into the process. When NGSolve subsequently tries to
load Netgen's `libngpy`, the already-loaded Cubit DLLs conflict with Netgen's
expected library versions, causing initialization failure.

**Fix**: Always import NGSolve BEFORE Cubit:
```python
# CORRECT: NGSolve first, then Cubit
import ngsolve                    # Loads Netgen DLLs cleanly
from netgen.meshing import Mesh
from netgen.occ import OCCGeometry

import sys, os
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit                      # Safe: Netgen DLLs already loaded
cubit.init(['cubit', '-nojournal', '-batch'])
```

```python
# WRONG: Cubit first causes DLL conflict
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")
import cubit                      # Loads Cubit's VTK DLLs
import ngsolve                    # FAILS — Netgen can't initialize
```

**Rule**: In any script that uses both NGSolve and Cubit, put all `import ngsolve`
and `from netgen...` statements at the top, before any Cubit-related path
manipulation or `import cubit`.

## Problem: JMAG Cannot Read Nastran File (Pyramid Elements)

**Cause**: JMAG does not support standard CPYRAM (5-node pyramid) elements.
When Cubit produces mixed hex-tet meshes, pyramid transition elements are generated.

**Fix**: Use `PYRAM=False` to write pyramids as degenerate CHEXA:
```python
cubit_mesh_export.export_nastran(cubit, "mesh.bdf", PYRAM=False)
```

**Alternative**: Avoid pyramids entirely by using pure tet mesh:
```python
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add tet all")
cubit_mesh_export.export_nastran(cubit, "mesh.bdf")
```

## Problem: Hex Meshing Fails on Complex Geometry

**Cause**: Not all geometries are decomposable into sweepable/mappable sub-volumes.

**Fixes** (try in order):
1. Webcut to simplify:
   ```python
   cubit.cmd("webcut volume 1 with plane xplane")
   cubit.cmd("merge all")
   cubit.cmd("volume all scheme auto")
   cubit.cmd("mesh volume all")
   ```
2. Use `scheme auto` (Cubit may choose a mix):
   ```python
   cubit.cmd("volume all scheme auto")
   ```
3. Give up on hex and use tet:
   ```python
   cubit.cmd("volume all scheme tetmesh")
   ```

## Problem: export_NGSolveCurvedMesh Returns Mesh with 0 Boundary Elements

**Cause**: No surface element block (block with tri/quad).
```python
# Add boundary element block
cubit.cmd("block 2 add tri all")   # For tet meshes
cubit.cmd("block 2 add quad all")  # For hex meshes
```
"""


CUBIT_DESIGN_PHILOSOPHY = """
# Design Philosophy: Why This Module Exists

## Core Idea

This module does **NOT** use Cubit's built-in export commands
(`export genesis`, `export mesh`, `cubit.cmd("export ...")`, etc.) for most formats.
Instead, it reads mesh data via the **Cubit Python API** and constructs output
files directly in Python.

## How It Works

```
Cubit Python API                    Output file
  get_block_id_list()      ──┐
  get_block_tets(id)       ──┤
  get_connectivity()       ──┼──>  Python constructs  ──>  .msh / .vtk / .bdf / .meg
  get_expanded_connectivity()─┤
  get_nodal_coordinates()  ──┘
```

Each export function:
1. Iterates through all blocks
2. Gets element IDs and connectivity via Cubit API
3. Gets node coordinates via `get_nodal_coordinates()`
4. Writes the output file format directly

## Why Not Use Built-in Export?

Cubit's built-in export covers a limited set of formats (Exodus, STEP, SAT, etc.).
For formats like Gmsh, VTK, Nastran, MEG, and Netgen, there is no built-in export
command. This module fills that gap by providing format-specific writers that read
from Cubit's mesh data.

## Exception: export_exodus()

`export_exodus()` is the **only** function that uses a Cubit command internally:
```python
cubit.cmd('export mesh "filename.exo" overwrite')
```
This is because Exodus II is Cubit's native format and the built-in command
provides full fidelity (all element types, nodesets, sidesets, etc.).

## API Summary

| API Function | Returns | Used By |
|-------------|---------|---------|
| `get_block_id_list()` | All block IDs | All export functions |
| `get_block_tets(id)` | Tet IDs in block | All formats with tet support |
| `get_block_hexes(id)` | Hex IDs in block | All formats with hex support |
| `get_connectivity(type, id)` | Corner nodes only | nastran, meg, netgen |
| `get_expanded_connectivity(type, id)` | All nodes (incl. mid-edge) | gmsh, vtk, vtu, exodus |
| `get_nodal_coordinates(node_id)` | (x, y, z) | All export functions |

## For AI Assistants: Script Generation Rule

When generating Cubit export scripts, **never** use:
```python
# WRONG: No such built-in command for these formats
cubit.cmd('export gmsh "mesh.msh"')     # Does not exist
cubit.cmd('export vtk "mesh.vtk"')      # Does not exist
cubit.cmd('export nastran "mesh.bdf"')  # Does not exist
```

**Always** use the module functions:
```python
# RIGHT: Use cubit_mesh_export functions
cubit_mesh_export.export_gmsh_v2(cubit, "mesh.msh")
cubit_mesh_export.export_nastran(cubit, "mesh.bdf", DIM="3D")
```
"""

CUBIT_APREPRO_JOURNAL = """
# Running Python Scripts from Cubit (APREPRO Journal)

## The `play` Command

Cubit can execute Python scripts directly via the journal system:

```
Cubit> play "export_mesh.py"
```

This is the standard way to run `cubit_mesh_export` scripts from within Cubit's
GUI or command line. The `play` command:

1. Executes the Python file using Cubit's embedded Python interpreter
2. The `cubit` module is already available (no `import cubit` needed in some contexts)
3. All `cubit.cmd()` calls execute directly in the running Cubit session

## Script Template for `play`

```python
# export_mesh.py - Run with: play "export_mesh.py"
import cubit_mesh_export

# Mesh is already generated in Cubit GUI
# Just register blocks and export

cubit.cmd("block 1 add tet all")
cubit.cmd('block 1 name "domain"')
cubit.cmd("block 2 add tri all")
cubit.cmd('block 2 name "boundary"')

# Export to desired format
cubit_mesh_export.export_gmsh_v2(cubit, "mesh.msh")
```

## Prerequisites

The `cubit_mesh_export` module is part of the Radia package (`src/radia/cubit_mesh_export.py`):

```bash
# Install Radia (includes cubit_mesh_export)
pip install radia

# Or add Radia's src/radia to sys.path for direct import
sys.path.insert(0, "path/to/Radia/src/radia")
import cubit_mesh_export
```

## APREPRO Variables

Cubit's journal system supports APREPRO variable substitution:

```
# In Cubit journal file (.jou)
{radius = 0.5}
{height = 2.0}
create cylinder height {height} radius {radius}
volume 1 scheme tetmesh
volume 1 size 0.1
mesh volume 1
play "export_mesh.py"
```

APREPRO variables are expanded before the command is executed. However,
APREPRO variables do NOT propagate into Python scripts run via `play`.
To pass parameters, use environment variables or file-based communication.

## Journal + Python Hybrid Workflow

A common pattern uses a `.jou` file for geometry/meshing and a `.py` file for export:

```
# workflow.jou - Run with: cubit -batch -nographics -nojournal workflow.jou
reset
import step "geometry.step" heal
volume all scheme tetmesh
volume all size 0.1
mesh volume all
play "export_mesh.py"
```

```bash
# Execute from command line (use CUBIT_PATH or full path)
"%CUBIT_PATH%\\coreform_cubit.exe" -batch -nographics -nojournal workflow.jou
# Or with full path:
"C:\\Program Files\\Coreform Cubit 2025.3\\bin\\coreform_cubit.exe" -batch -nographics -nojournal workflow.jou
```

## Batch Mode vs GUI Mode

| Mode | Command | `cubit` Available | Use Case |
|------|---------|-------------------|----------|
| GUI | `play "script.py"` in command panel | Yes (pre-imported) | Interactive |
| Batch | `coreform_cubit -batch ... journal.jou` | Yes (pre-imported) | Automated |
| Standalone | `python script.py` (with sys.path) | Need `import cubit` | Development |
"""


CUBIT_MESH_TO_GEOMETRY = """\
# Mesh-to-Geometry Conversion (Mesh2Acis)

Convert mesh elements back to ACIS geometry volumes. Useful for creating
geometry from imported meshes (e.g., for Kelvin transformation or re-meshing).

## Core Pattern: Element Nodes -> Vertices -> Surfaces -> Volume

```python
import cubit

def create_quad(v1, v2, v3, v4):
    cubit.cmd(f"create surface vertex {v1} {v2} {v3} {v4}")

def create_tri(v1, v2, v3):
    cubit.cmd(f"create surface vertex {v1} {v2} {v3}")

# Iterate over mesh elements and reconstruct as geometry
volume_list = cubit.get_entities("volume")
for vol_id in volume_list:
    # Hexahedra: 8 nodes -> 6 quad faces -> 1 volume
    for hid in cubit.get_volume_hexes(vol_id):
        nodes = cubit.get_connectivity("hex", hid)
        verts = []
        for nid in nodes:
            x = cubit.get_nodal_coordinates(nid)
            cubit.create_vertex(x[0], x[1], x[2])
            verts.append(cubit.get_last_id("vertex"))
        # Create 6 quad faces
        for face_nodes in [(0,1,2,3), (4,5,6,7), (3,2,6,7),
                           (0,1,5,4), (2,1,5,6), (3,0,4,7)]:
            create_quad(*[verts[i] for i in face_nodes])
        sids = [cubit.get_last_id("surface") - 5 + i for i in range(6)]
        cubit.cmd(f"create volume surface {' '.join(map(str, sids))}")

    # Tetrahedra: 4 nodes -> 4 tri faces -> 1 volume
    for tid in cubit.get_volume_tets(vol_id):
        nodes = cubit.get_connectivity("tet", tid)
        verts = []
        for nid in nodes:
            x = cubit.get_nodal_coordinates(nid)
            cubit.create_vertex(x[0], x[1], x[2])
            verts.append(cubit.get_last_id("vertex"))
        for face_nodes in [(0,1,2), (3,0,1), (3,1,2), (3,2,0)]:
            create_tri(*[verts[i] for i in face_nodes])
        sids = [cubit.get_last_id("surface") - 3 + i for i in range(4)]
        cubit.cmd(f"create volume surface {' '.join(map(str, sids))}")

    # Also available: cubit.get_volume_pyramids(), cubit.get_volume_wedges()
```

## Important: Renumber nodes before processing

```python
cubit.cmd('renumber node all in Volume all start_id 1 uniqueids')
```

## Key API Functions

| Function | Returns |
|----------|---------|
| `cubit.get_entities("volume")` | List of all volume IDs |
| `cubit.get_volume_hexes(vol_id)` | Hex element IDs in volume |
| `cubit.get_volume_tets(vol_id)` | Tet element IDs in volume |
| `cubit.get_volume_pyramids(vol_id)` | Pyramid element IDs |
| `cubit.get_volume_wedges(vol_id)` | Wedge element IDs |
| `cubit.get_connectivity("hex", eid)` | Node IDs (1st order) |
| `cubit.get_nodal_coordinates(nid)` | (x, y, z) tuple |
| `cubit.create_vertex(x, y, z)` | Creates vertex, use get_last_id |
| `cubit.get_last_id("vertex")` | Most recently created entity ID |
"""

CUBIT_PARAMETRIC_GEOMETRY = """\
# Parametric Geometry Generation

## Archimedean Spiral (Twisted Wire)

Generate spiral curves from parametric equations and sweep a cross-section along them:

```python
import cubit, math

cubit.cmd('reset')

# Spiral parameters
dtheta = 0.1    # angle increment (radians)
a = 0.1         # spiral growth rate: r = a * theta
offset_theta = 2 * math.pi  # start offset

# Generate vertices along spiral
for i in range(1000):
    theta = dtheta * i + offset_theta
    r = a * theta
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    cubit.cmd(f"create vertex {x} {y} 0")

# Create spline through all vertices
cubit.cmd("create curve spline location vertex all")

# Create cross-section and sweep
wire_w = 0.9 * (2 * a * math.pi)
wire_h = 1.0
cubit.cmd(f"create surface rectangle width {wire_w} height {wire_h} yplane")
cubit.cmd("move Surface 1 location vertex 1")
cubit.cmd("sweep surface 1 along curve 1")
```

## Icosahedron on Sphere (Regular Surface Mesh)

Generate a regular triangular mesh on a sphere using icosahedron vertices
and great-circle arcs. Uses APREPRO variables for golden ratio construction:

```python
import cubit, os

for n in [2, 4, 8, 16]:  # mesh refinement levels
    cubit.cmd('reset')

    # Golden ratio icosahedron vertices via APREPRO
    cubit.cmd('#{p = (1+sqrt(5))/2}')    # golden ratio
    cubit.cmd('#{q = sqrt(p+2)}')        # normalization
    cubit.cmd('#{r0 = 100}')             # sphere radius
    cubit.cmd('#{r1 = r0/q}')
    cubit.cmd('#{r2 = r1*p}')

    cubit.cmd('create sphere radius {r0}')

    # Create 12 icosahedron vertices
    cubit.cmd('create vertex location  0,{ r1},{ r2}')   # vertex 1
    cubit.cmd('create vertex location  0,{-r1},{ r2}')   # vertex 2
    # ... (10 more vertices for complete icosahedron)

    # Create great-circle arcs between adjacent vertices
    cubit.cmd('create vertex location 0, 0, 0')  # center (vertex 13)
    cubit.cmd('create curve arc center vertex 13 1 2')
    cubit.cmd('create curve arc center vertex 13 1 3')
    # ... (30 edges total for icosahedron)

    # Imprint arcs onto sphere surface
    cubit.cmd('imprint volume all with curve all')
    cubit.cmd('merge all')
    cubit.cmd('merge vertex all force')

    # Mesh with parametric refinement
    cubit.cmd(f'surface all interval {n}')
    cubit.cmd('Surface all scheme triadvance')
    cubit.cmd('mesh surf all')

    # Smooth for quality
    cubit.cmd('surface all smooth scheme smart laplacian')
    cubit.cmd('smooth surface all')

    # Export
    cubit.cmd('block 1 add surface all')
    folder = f"n={n}"
    os.makedirs(folder, exist_ok=True)
    cubit.cmd(f'export nastran "{folder}/sphere.nas" dimension 3 overwrite large')
```

## Key Techniques

| Technique | Command |
|-----------|---------|
| Parametric vertex | `create vertex {x} {y} {z}` in loop |
| Spline through vertices | `create curve spline location vertex all` |
| Great-circle arc | `create curve arc center vertex CENTER V1 V2` |
| Sweep along curve | `sweep surface ID along curve ID` |
| APREPRO variable | `#{var = expression}` then `{var}` in commands |
| Imprint curves on surface | `imprint volume all with curve all` |
"""

CUBIT_BATCH_PROCESSING = """\
# Batch Processing & Parametric Studies

## Multi-Case Mesh Generation

Generate meshes at multiple refinement levels in separate folders:

```python
import cubit, os

for mesh_size in [5.0, 2.0, 1.0, 0.5]:
    folder = f'mesh={mesh_size:.1f}'
    os.makedirs(folder, exist_ok=True)

    cubit.init([''])
    cubit.cmd('reset')
    cubit.cmd('create surface circle radius 3 zplane')
    cubit.cmd('surface all scheme circle interval')
    cubit.cmd(f'surface all size {{0.2*{mesh_size}}}')
    cubit.cmd('mesh surface all')
    cubit.cmd(f'export nastran "{folder}/mesh.bdf" dimension 2 overwrite large')
    cubit.cmd('exit')
```

## Batch CMD Generation for External Solvers

Generate batch scripts for running multiple solver cases:

```python
import os

cases = []
for current in range(1000, 21000, 2000):
    for element in ['HEX8', 'TET4']:
        for mesh in ['coarse', 'medium', 'fine']:
            cases.append((current, element, mesh))

with open('batch_run.cmd', 'w') as f:
    for current, element, mesh in cases:
        folder = f'{current}AT/{element}_{mesh}'
        f.write(f'cd /d {os.path.join(base_dir, folder)}\\n')
        f.write(f'solver.exe input.dat\\n')

os.system('batch_run.cmd')
```

## Directory Structure for Parametric Studies

```
project/
  mesh=5.0/mesh.bdf
  mesh=2.0/mesh.bdf
  mesh=1.0/mesh.bdf
  mesh=0.5/mesh.bdf
  batch_run.cmd
  results/
    summary.csv
```

## Best Practices

1. **Use `os.makedirs(folder, exist_ok=True)`** - avoids errors if folder exists
2. **Call `cubit.cmd('reset')` between cases** - clean state for each mesh
3. **Use `cubit.init([''])` per case** if running standalone (not in GUI)
4. **Format filenames with parameters** - enables post-processing automation
5. **Graphics off for batch**: `cubit.cmd('set echo off')` and `-nographics -batch` flags
"""

CUBIT_FORMAT_CONVERSION = """\
# Format Conversion Pipelines

## Cubit -> Nastran -> VTU (via meshio)

Convert Cubit mesh to VTU format for ParaView visualization using meshio:

```python
import cubit, meshio

# Create and mesh geometry in Cubit
cubit.cmd('reset')
cubit.cmd('brick x 10')
cubit.cmd('create Cylinder height 10 radius 10')
cubit.cmd('imprint all')
cubit.cmd('merge all')
cubit.cmd('mesh vol all')

# Register blocks (required for Nastran export)
cubit.cmd('block 1 add vol 1')
cubit.cmd('block 2 add vol 2')
cubit.cmd('export nastran "mesh.nas" dimension 3 overwrite everything')

# Convert to VTU using meshio
mesh = meshio.read("mesh.nas")
meshio.write("output.vtu", mesh, file_format="vtu", binary=True)
```

## Cubit -> Patran -> Coordinate Transform -> Reimport

Export to Patran format, apply coordinate transformation, reimport as geometry:

```python
import numpy as np
import re

def transform_patran_mesh(patran_file, transform_func):
    \"\"\"Read Patran file, transform node coordinates, write back.\"\"\"
    with open(patran_file, 'r') as f:
        lines = f.readlines()

    for n in range(len(lines)):
        fields = re.split(r"\\s+|\\n$", lines[n])
        if fields[0] == '01':  # Node record
            coord_line = re.split(r"\\s+|\\n$", lines[n+1])
            xyz = np.array([float(coord_line[1]),
                            float(coord_line[2]),
                            float(coord_line[3])])
            new_xyz = transform_func(xyz)
            lines[n+1] = f' {new_xyz[0]:15e} {new_xyz[1]:15e} {new_xyz[2]:15e}\\n'

    with open(patran_file, 'w') as f:
        f.writelines(lines)

# Example: Kelvin inversion (maps exterior to interior of sphere)
a = 30.0  # inversion radius
def kelvin_inversion(xyz):
    r = np.linalg.norm(xyz)
    return (a / r)**2 * xyz

cubit.cmd('block 1 add surface all')
cubit.cmd(f'export patran "temp.pat" block 1 dimension 3 overwrite')
transform_patran_mesh("temp.pat", kelvin_inversion)
cubit.cmd('reset')
cubit.cmd(f'import patran mesh geometry "temp.pat" feature_angle 135.00')

# Reconstruct surfaces from imported mesh
for sid in range(1, 7):
    cubit.cmd(f'create surface net from mapped surface {sid} heal')
cubit.cmd('delete body 1 to 6')
cubit.cmd('create volume surface 7 to 12 heal')
```

## Format Comparison for Conversion

| Source | Target | Method | Use Case |
|--------|--------|--------|----------|
| Cubit | VTU/VTK | meshio | ParaView visualization |
| Cubit | Netgen | cubit_mesh_export | NGSolve FEM (part of Radia) |
| Cubit | Gmsh | cubit_mesh_export | GMSH post-processing |
| Cubit | Patran | built-in export | Coordinate transforms |
| Gmsh GEO | Cubit | geo2jou converter | Geometry migration |

## meshio Installation

```bash
pip install meshio
# Supports: Nastran, VTK, VTU, Exodus, Gmsh, XDMF, and 30+ formats
```
"""

CUBIT_KELVIN_TRANSFORM = """\
# Kelvin Sphere Transformation in Cubit

The Kelvin (inversion) transform maps the exterior of a sphere to its interior,
enabling FEM solution of open-boundary problems on a bounded domain.

## Transform Formula

```
p' = (a/r)^2 * p
```

where `a` is the inversion sphere radius and `r = |p|` is the distance from origin.

## Complete Workflow

```python
import cubit
import numpy as np

a = 30.0   # inversion sphere radius
up = 2.1 * a  # vertical offset for Kelvin domain

# 1. Create and mesh the source geometry
cubit.cmd('reset')
cubit.cmd('brick x 60 y 60 z 5')
cubit.cmd('move Volume 1 z -22.5 include_merged')
cubit.cmd('volume all size 2.5')
cubit.cmd('mesh surf all')

# 2. Export to Patran, apply Kelvin transform, reimport
cubit.cmd('block 1 add surface all')
cubit.cmd('export patran "temp.pat" block 1 dimension 3 overwrite')
# (apply coordinate transform to temp.pat - see Format Conversion topic)

cubit.cmd('reset')
cubit.cmd('import patran mesh geometry "temp.pat" feature_angle 135.00')

# 3. Reconstruct smooth surfaces from imported mesh
for n in range(1, 7):
    cubit.cmd(f'create surface net from mapped surface {n} heal')
cubit.cmd('delete body 1 to 6')
cubit.cmd('create volume surface 7 to 12 heal')

# 4. Intersect with inversion sphere
cubit.cmd(f'create sphere radius {a}')
vid = cubit.get_last_id("volume")
cubit.cmd(f'move Volume {vid} z {up} include_merged')
cubit.cmd(f'intersect volume {vid-1} {vid}')

# 5. Create complementary region (sphere minus Kelvin volume)
cubit.cmd(f'create sphere radius {a}')
vid = cubit.get_last_id("volume")
cubit.cmd(f'move Volume {vid} z {up} include_merged')
cubit.cmd(f'subtract volume {vid-1} from volume {vid} imprint keep_tool')

# 6. Imprint and merge for conformal mesh
cubit.cmd('imprint vol all')
cubit.cmd('merge vol all')
```

## Tips

- **feature_angle 135**: Important for Patran reimport to detect surface boundaries
- **`create surface net from mapped surface`**: Reconstructs smooth geometry from mesh
- **heal**: Repairs small gaps in reconstructed geometry
- Use with NGSolve Kelvin CoefficientFunction for material property mapping
"""


CUBIT_NGSOLVE_PANEL = """\
# Cubit Panel Development: External Python + NGSolve

## Architecture: Two-Process Model

Cubit panels that use NGSolve must run computation in an **external Python
subprocess** because Cubit's bundled Python (3.10) does not have NGSolve.

```
Cubit GUI (bundled Python 3.10, has PyQt5)
  register_toolbar.py
    1. Qt dialog (modeless, non-blocking)
    2. cubit.cmd('save cub5 "temp.cub5" overwrite')
    3. QProcess.start(external_python, [calc_*.py, ...])
    4. QProcess.finished -> read JSON result file (--output)
    5. Update Qt dialog with results

External Python (system Python 3.12, has NGSolve)
  calc_volume.py / calc_surface.py / calc_inductance.py
    1. Import NGSolve FIRST (before cubit)
    2. Import cubit (suppress banner, -noinit)
    3. cubit.cmd('open "temp.cub5"')
    4. Compute (export_NGSolveCurvedMesh, BEM, etc.)
    5. Write result to --output JSON file
```

## Critical Rules

1. **NGSolve before cubit**: `from ngsolve import ...` MUST come before
   `import cubit`. Cubit adds its site-packages (Python 3.10 numpy) to
   sys.path, which conflicts with system Python 3.12 numpy.

2. **QProcess, not subprocess.run()**: subprocess.run() blocks the GUI thread,
   causing Cubit to freeze/crash. QProcess is non-blocking.

3. **Modeless dialogs**: Use `dialog.show()` not `dialog.exec()`. Modal
   exec() blocks Cubit's main window. Keep a reference on main_window to
   prevent garbage collection:
   ```python
   def _show_dialog():
       if not hasattr(main_window, '_dlg') or main_window._dlg is None:
           main_window._dlg = MyDialog(main_window)
       main_window._dlg.show()
       main_window._dlg.raise_()
   action.triggered.connect(_show_dialog)
   ```

4. **-noinit flag**: `cubit.init([..., "-noinit"])` prevents .cubit startup
   file from being replayed when opening cub5 in external Python.

5. **Suppress cubit banner**: `cubit.init()` outputs a large banner to
   C-level stderr. Use `contextlib.redirect_stderr(io.StringIO())`.

6. **JSON result via --output file**: Do NOT rely on stdout parsing.
   cubit.cmd() pollutes C-level stdout. Use `--output result.json`:
   ```python
   # calc_inductance.py
   parser.add_argument("--output", default="")
   if args.output:
       with open(args.output, "w") as f:
           json.dump(result, f)

   # register_toolbar.py
   args = [..., "--output", json_path]
   # In _on_extract_finished:
   with open(json_path, "r") as f:
       data = json.load(f)
   ```

7. **Cubit site-packages cleanup**: After `import cubit`, remove Cubit's
   bundled site-packages from sys.path to avoid numpy/scipy conflicts.

8. **startup.py encoding**: Use `exec(open(..., encoding="utf-8").read())`
   because Windows Japanese locale defaults to cp932. Never use non-ASCII
   characters (em dash, smart quotes, etc.) in register_toolbar.py.

## BEM Inductance Panel (DC, Source/Sink EFIE)

### Method

Saddle point EFIE with source/sink current injection:
```
[SL  D^T] [J] = [0]
[D   0  ] [p] = [g]
```
- SL = LaplaceSL (ngsolve.bem, Laplace kernel 1/(4*pi*r))
- D = divergence (HDivSurface -> SurfaceL2 order=0)
- g = +1/A_source on source faces, -1/A_sink on sink faces
- L = mu_0 * J^T @ SL @ J

### Block Registration for Source/Sink

CRITICAL: Register source/sink blocks AFTER conductor, BEFORE boundary.
In export_NGSolveCurvedMesh, last block wins for bcname. Source/sink
(subset) must have higher block ID than boundary (superset).

```python
# CORRECT: source/sink before boundary
cubit.cmd("set duplicate block elements on")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "conductor"')
cubit.cmd("block 2 add face in surface 1")   # gap face
cubit.cmd('block 2 name "source"')
cubit.cmd("block 3 add face in surface 3")   # gap face
cubit.cmd('block 3 name "sink"')
cubit.cmd("block 4 add face in surface all")  # superset last
cubit.cmd('block 4 name "boundary"')

# WRONG: boundary first -> source/sink labels overridden
cubit.cmd("block 1 add face all")           # boundary first
cubit.cmd('block 1 name "boundary"')
cubit.cmd("block 2 add face in surface 1")  # last block wins
cubit.cmd('block 2 name "source"')          # boundary overrides!
```

### Hex Mesh for BEM (Recommended)

Create gapped torus by revolve (NOT subtract). Subtract creates complex
topology that Cubit's sweeper cannot handle.

```
# Revolve approach (sweepable)
create surface circle radius 0.005 yplane
move surface 1 x 0.05 include_merged
sweep surface 1 zaxis angle 355
volume 1 scheme sweep source surface 1 target surface 3
surface 1 size 0.005
surface 1 scheme pave
mesh surface 1
mesh volume 1
```

Surface IDs from revolve are deterministic: 1=start, 2=body, 3=end.

### Panel UI Settings

| Setting | Values | Default | Description |
|---------|--------|---------|-------------|
| Curve order | 1-2 | 2 | Mesh geometry (2=curved elements) |
| FES order | 0-2 | 0 | HDivSurface basis (0=RWG) |

- Curve(2): improves geometric accuracy (integration points on curved surface)
- FES order=0: sufficient for smooth current distributions (torus)
- FES order=1: useful for sharp features (corners, edges)
- Constraint (SurfaceL2) is always order=0 (element-wise current conservation)

### Performance (source/sink EFIE, ToDense)

| DOFs | SL assembly | SL extract | Solve | Total |
|------|------------|------------|-------|-------|
| ~1000 | 0.6s | 1.6s | 0.1s | ~2s |
| ~2000 | 4.5s | 26s | 0.6s | ~31s |
| ~5000 | 19s | 143s | 2.5s | ~165s |

Bottleneck: SL extract (ToDense) is O(N^2). For larger problems, need
iterative solver + FMM matvec.

D matrix: use mat.ToDense().NumPy() (0.01s), NOT column-by-column extraction
(13s for 5000 DOFs).

### Cubit Surface Mesh Size Control

`volume all size X` only affects interior tets, NOT surface mesh.
Cubit's trimesher uses geometry-based sizing by default.

To control surface mesh density:
```
set trimesher geometry sizing off   # disable curvature-based sizing
surface all size 0.005              # explicit surface element size
```

Without this, all mesh size values give the same surface mesh count.

### GMSH Visualization

Open GMSH button launches GMSH GUI via external Python:
```python
ext_python -c "import sys, gmsh; gmsh.initialize(sys.argv, run=True); gmsh.finalize()" file.msh
```

CRITICAL: Do NOT use gmsh.bat or `python -m gmsh` from Cubit's Python.
- gmsh.bat is a Python wrapper that calls system python, but Cubit's Python
  intercepts the call
- `python -m gmsh` does NOT launch the GUI
- Only `gmsh.initialize(sys.argv, run=True)` launches the GMSH GUI window

Result .msh v2.2 contains two NodeData fields:
- |J| (ncomp=1): scalar current density magnitude, colormap display
- J (ncomp=3): vector current density, arrow display on curved surface

IMPORTANT: Vector data must be passed as 2D numpy array shape (nv, 3),
NOT flatten(). flatten() creates 1D array causing arr[ni,c] IndexError
in GmshPostExport.write_v22().

```python
# CORRECT
post.add_field("J", node_vec, ncomp=3)     # shape (nv, 3)

# WRONG - causes empty NodeData
post.add_field("J", node_vec.flatten(), ncomp=3)  # shape (nv*3,)
```

GMSH display: |J| shows current crowding on inner radius of torus bend
(physically correct: shorter path = higher current density at constant
voltage across source/sink ports).

### COMSOL Import Workflow

.msh v2.2 with NodeData is directly importable by COMSOL:

1. File > Import > GMSH (.msh v2.2)
2. Physical Groups ("source", "sink", "boundary") -> COMSOL Selections
3. NodeData (|J|, J) -> Results > External dataset
4. Heat Transfer: Q = |J|^2 / sigma as boundary heat source
5. Lorentz force: F = J x B (vector field from NodeData)

### Panel Code Reload (No Cubit Restart)

To reload panel code after editing register_toolbar.py, run in Cubit:
```
play "S:/Radia/01_GitHub/src/radia/panels/startup.py"
```
This re-exec's register_toolbar.py, updating all class definitions.
The menu shows "Radia menu already registered." (duplicate check), but
dialog classes (InductanceDialog, etc.) are reloaded with new code.

After reload, close and reopen the dialog to use updated code.
Existing dialog instances use the OLD class definition.
"""


CUBIT_RESULTS_VISUALIZATION = """\
# Cubit Results Visualization

## set_nodal_variable: Display Field Data in Cubit

Cubit can display per-node scalar fields using the results view mode
(rainbow icon in toolbar). Use `set_nodal_variable` to set values:

```python
import cubit

# Get all node IDs
node_ids = list(cubit.get_entities('node'))

# Compute values (e.g., from BEM solve)
values = [compute_field(nid) for nid in node_ids]

# Set nodal variable
cubit.set_nodal_variable(node_ids, "J_magnitude", values)

# Switch to results view mode to see the contour plot
# (click rainbow icon in Cubit toolbar, or use Display menu)
```

## Limitations of set_nodal_variable

**IMPORTANT**: `set_nodal_variable` is limited:
- Display only — values are NOT preserved on save/reload of .cub5
- No export — Exodus export does NOT include these variables automatically
  (despite claims in older Cubit docs). You must write a separate Exodus file.
- No high-order interpolation — Cubit renders piecewise-linear on mesh facets
- No vector field support — scalar only, one component at a time

For production post-processing, use GMSH .msh output instead (see below).

## Recommended: GMSH Post-Processing (.msh)

For high-order element visualization and persistent results, use GmshPostExport:
```python
from radia.gmsh_post_export import GmshPostExport
post = GmshPostExport(mesh, boundary=True)  # boundary=True for BND from volume mesh
post.add_field("|J|", node_J, ncomp=1)       # per-node scalar
post.write("results.msh")
# Open in GMSH GUI for high-order interpolation
```

For curved meshes (mesh.Curve(p)), GmshPostExport automatically outputs
high-order elements (Tri6/Tri10/Tri15/Tri21) with nodes on the curved surface.
Node coordinates are extracted via GetTrafo evaluated at GMSH reference points.

**GMSH display setting**: High-order elements require `Mesh.NumSubEdges = 4`
(or higher) in GMSH to render curved surfaces. Default is 1 (straight lines).
Set via Tools -> Command Line: `Mesh.NumSubEdges = 4;`

**Advantages over Cubit set_nodal_variable**:
- Persistent: .msh file saved alongside model
- High-order: GMSH renders curved elements correctly
- Vector fields: supports ncomp=3 for vector visualization
- Linked to Optuna DB: each trial's gmsh_file recorded in user_attrs

## Workflow: BEM Current Density Visualization

```
1. BEM solve (external Python) -> per-DOF current vector x
2. Convert HDivSurface DOFs -> per-node scalar (J magnitude)
3. Write .msh file via GmshPostExport (persistent)
4. Return result JSON (inductance, n_dofs, gmsh_file path)
5. Record trial in Optuna SQLite DB (params + result + gmsh_file)
6. Panel auto-launches GMSH with .msh file
7. optuna-dashboard shows parameter sweep with linked .msh files
```

## Cubit Import Formats (NO VTK)

Cubit does NOT support VTK import. Supported mesh import formats:
- Exodus II (.exo, .e, .g) - with nodal_var support
- Abaqus, Nastran, Patran, STEP, IGES, Fluent, LSDyna
"""


def get_cubit_documentation(topic: str = "all") -> str:
	"""Return Cubit scripting documentation by topic."""
	topics = {
		"overview": CUBIT_OVERVIEW,
		"blocks": CUBIT_BLOCKS,
		"blocks_only_policy": CUBIT_BLOCKS_ONLY_POLICY,
		"element_order": CUBIT_ELEMENT_ORDER,
		"mesh_schemes": CUBIT_MESH_SCHEMES,
		"step_exchange": CUBIT_STEP_EXCHANGE,
		"initialization": CUBIT_INITIALIZATION,
		"mixed_elements": CUBIT_MIXED_ELEMENTS,
		"common_mistakes": CUBIT_COMMON_MISTAKES,
		"hex_workflow": CUBIT_HEX_WORKFLOW,
		"tet_workflow": CUBIT_TET_WORKFLOW,
		"2d_mesh": CUBIT_2D_MESH_WORKFLOW,
		"troubleshooting": CUBIT_TROUBLESHOOTING,
		"design_philosophy": CUBIT_DESIGN_PHILOSOPHY,
		"aprepro_journal": CUBIT_APREPRO_JOURNAL,
		"mesh_to_geometry": CUBIT_MESH_TO_GEOMETRY,
		"parametric_geometry": CUBIT_PARAMETRIC_GEOMETRY,
		"batch_processing": CUBIT_BATCH_PROCESSING,
		"format_conversion": CUBIT_FORMAT_CONVERSION,
		"kelvin_transform": CUBIT_KELVIN_TRANSFORM,
		"ngsolve_panel": CUBIT_NGSOLVE_PANEL,
		"optuna": CUBIT_NGSOLVE_PANEL,  # Optuna integration documented in panel section
		"results_visualization": CUBIT_RESULTS_VISUALIZATION,
		"post_processing": CUBIT_RESULTS_VISUALIZATION,
	}

	topic = topic.lower().strip()
	if topic == "all":
		return "\n\n".join(topics.values())
	elif topic in topics:
		return topics[topic]
	else:
		return (
			f"Unknown topic: '{topic}'. "
			f"Available: all, {', '.join(topics.keys())}"
		)
