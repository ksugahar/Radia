"""
GMSH knowledge base for Radia MCP server.

Covers: visualization, post-processing, .msh file format, .geo scripting,
command-line usage, options system, and high-order element display.

IMPORTANT: In the Radia project, GMSH is used for VISUALIZATION AND
POST-PROCESSING ONLY, NOT for mesh generation. Mesh generation uses
Netgen (tet) or Cubit (hex). See GMSH_RADIA_POLICY for details.
"""

GMSH_RADIA_POLICY = """
# GMSH Usage Policy in Radia

**GMSH is used for visualization and post-processing only, NOT mesh generation.**

## Allowed Uses
- Opening and visualizing .msh files (GMSH GUI)
- Post-processing field data (NodeData, ElementData)
- Merging STEP geometry + .msh field data for overlay visualization
- .geo companion files for display settings
- Reading .msh format via `gmsh_mesh_import.py` (pure file reader)

## NOT Allowed
- GMSH Python API (`gmsh.model.occ.*`) for geometry creation
- GMSH Python API for mesh generation
- `import gmsh` in computation scripts
- `from radia.gmsh_builder import GmshBuilder` (removed)

## Mesh Generation: 2-Path Only
1. **STEP -> Netgen** (NGSolve OCC): tet meshes with `mesh.Curve(order)`
2. **STEP -> Cubit** (Coreform Cubit): structured hex meshes

## Visualization Workflow
```
GMSH GUI:
  Merge "coil.step"        <- CoilBuilder.write_step()
  Merge "magnet.step"      <- OCC shape -> STEP
  Merge "field.msh"        <- NGSolve -> GmshPostExport (.msh)
  -> Geometry + field overlay visualization
```

## Output Format
- Default: .msh v2.2 (NGSolve ReadGmsh() compatible)
- Optional: .msh v4.1 (large-scale, structured Physical Groups)
- `GmshPostExport.write()` defaults to v2.2
"""

GMSH_OVERVIEW = """
# GMSH Overview

GMSH is an open-source 3D finite element mesh generator with a built-in
CAD engine and post-processor. Version 4.15.2 (March 2026).

- **License**: GNU GPL v2+
- **Authors**: Christophe Geuzaine, Jean-Francois Remacle
- **Website**: https://gmsh.info/
- **Citation**: Geuzaine & Remacle, IJNME 79(11), pp. 1309-1331, 2009

## Installation (Radia project)
GMSH is installed as a standalone executable (not pip).
- Lab machines: `C:\\gmsh.exe`
- .msh files are associated with gmsh.exe via registry

## Modules
1. **Geometry** - CAD kernel (built-in or OpenCASCADE)
2. **Mesh** - Automatic meshing (1D/2D/3D)
3. **Solver** - ONELAB/GetDP integration
4. **Post-processing** - Field visualization with views
"""

GMSH_COMMAND_LINE = """
# GMSH Command-Line Reference

## Opening Files
```bash
gmsh file.msh                    # Open .msh in GUI
gmsh file.step                   # Open STEP geometry
gmsh file.step -merge field.msh  # Overlay STEP + field data
gmsh file.geo                    # Run .geo script
```

## Key Options
```bash
-n                  # Hide all meshes and views on startup
-merge              # Merge next file (instead of opening)
-o file             # Specify output file name
-format string      # Output format (msh2, msh4, vtk, stl, bdf...)
-order int          # Set mesh element order
-bin                # Binary output
-save               # Save and exit
-save_all           # Save all elements (even without physical groups)
-v int              # Verbosity level (0-99, default 5)
-nt int             # Number of threads
-numsubedges int    # Subdivisions for high-order display
-string "cmd"       # Parse command string at startup
-setnumber name val # Set option value
-option file        # Parse option file at startup
-nopopup            # No dialog popups in scripts
```

## Mesh Generation (batch)
```bash
gmsh file.geo -3              # Generate 3D mesh
gmsh file.geo -2 -o out.msh  # Generate 2D mesh, save
gmsh file.geo -3 -order 2    # 2nd order 3D mesh
gmsh file.geo -3 -format msh22 -save_all  # MSH v2.2 output
```

## Post-processing
```bash
gmsh field.msh                          # View field data
gmsh -merge file1.msh -merge file2.msh  # Multiple data files
gmsh display.geo                        # .geo with Merge + options
```

## Useful Combinations
```bash
# Open .msh with curved elements properly displayed
gmsh file.msh -numsubedges 4

# Merge STEP geometry with field results
gmsh coil.step -merge field.msh -numsubedges 4
```
"""

GMSH_KEYBOARD_SHORTCUTS = """
# GMSH Keyboard Shortcuts

## Module Switching
| Key | Action |
|-----|--------|
| `g` | Geometry module |
| `m` | Mesh module |
| `p` | Post-processing module |
| `s` | Solver module |

## Meshing
| Key | Action |
|-----|--------|
| `1` / `F1` | Mesh 1D (lines) |
| `2` / `F2` | Mesh 2D (surfaces) |
| `3` / `F3` | Mesh 3D (volumes) |

## Navigation
| Key | Action |
|-----|--------|
| Left/Right arrow | Previous/next time step |
| Up/Down arrow | Previous/next view |
| `0` | Reload geometry |
| `Ctrl+0` or `9` | Reload full project |

## File Operations
| Key | Action |
|-----|--------|
| `Ctrl+O` | Open file |
| `Ctrl+S` | Save mesh |
| `Ctrl+E` | Export |
| `Ctrl+N` | New model |
| `Shift+Ctrl+O` | Merge file(s) |

## Options Dialogs
| Key | Action |
|-----|--------|
| `Shift+O` | General options |
| `Shift+G` | Geometry options |
| `Shift+M` | Mesh options |
| `Shift+P` | Post-processing options |
| `Shift+W` | View options |

## Display Toggles (Alt+key)
| Key | Action |
|-----|--------|
| `Alt+L` | Geometry lines |
| `Alt+P` | Geometry points |
| `Alt+S` | Geometry surfaces |
| `Alt+M` | Toggle all mesh entities |
| `Alt+Shift+P` | Mesh nodes |
| `Alt+Shift+S` | Mesh surface edges |
| `Alt+Shift+D` | Mesh surface faces |
| `Alt+Shift+V` | Mesh volume edges |
| `Alt+H` | Hide/show all views |
| `Alt+I` | Hide/show all scales |
| `Alt+W` | Enable/disable lighting |
| `Alt+E` | Show element outlines for views |

## View Orientation
| Key | Action |
|-----|--------|
| `Alt+X/Y/Z` | View along +X/+Y/+Z axis |
| `Alt+Shift+X/Y/Z` | View along -X/-Y/-Z axis |
| `Alt+1` | 1:1 scale |
| `Alt+O` | Toggle orthographic/perspective |

## Other
| Key | Action |
|-----|--------|
| `Ctrl+L` | Message console |
| `Ctrl+F` | Full screen |
| `Ctrl+Q` | Quit |
| `Alt+C` | Cycle color schemes |
| `Alt+T` | Cycle interval modes for views |
"""

GMSH_OPTIONS = """
# GMSH Options Reference

Options can be set in:
- .geo files: `Mesh.NumSubEdges = 4;`
- Command line: `gmsh -setnumber Mesh.NumSubEdges 4`
- GUI: Tools -> Options
- Configuration: `~/.gmsh` or `%APPDATA%/gmsh.conf`

## General Options
| Option | Default | Description |
|--------|---------|-------------|
| `General.NumThreads` | 1 | Parallel threads |
| `General.Verbosity` | 5 | Verbosity (0-99) |
| `General.Terminal` | 0 | Terminal output |
| `General.Orthographic` | 1 | Orthographic projection |
| `General.ConfirmQuit` | 1 | Ask before quit (set 0 to disable) |
| `General.SmallAxes` | 1 | Show axes indicator |

## Mesh Options (Most Important for Radia)
| Option | Default | Description |
|--------|---------|-------------|
| `Mesh.NumSubEdges` | 2 | **Subdivisions for high-order display (set 4 for curved)** |
| `Mesh.ElementOrder` | 1 | Element polynomial order |
| `Mesh.MshFileVersion` | 4.1 | MSH format version |
| `Mesh.Binary` | 0 | Binary output |
| `Mesh.SaveAll` | 0 | Save all elements |
| `Mesh.Algorithm` | 6 | 2D algorithm (6=Frontal-Delaunay) |
| `Mesh.Algorithm3D` | 1 | 3D algorithm (1=Delaunay, 10=HXT) |
| `Mesh.MeshSizeFactor` | 1.0 | Global size scaling |
| `Mesh.MeshSizeMin` | 0 | Minimum element size |
| `Mesh.MeshSizeMax` | 1e22 | Maximum element size |
| `Mesh.ColorCarousel` | 1 | Coloring mode |
| `Mesh.Lines` | 0 | Show 1D elements |
| `Mesh.SurfaceEdges` | 1 | Show surface edges |
| `Mesh.SurfaceFaces` | 0 | Show surface faces |
| `Mesh.VolumeEdges` | 1 | Show volume edges |
| `Mesh.VolumeFaces` | 0 | Show volume faces |

## View Options (View[n].*)
| Option | Default | Description |
|--------|---------|-------------|
| `View.Visible` | 1 | Show/hide view |
| `View.IntervalsType` | 2 | Display (1=iso, 2=continuous, 3=discrete, 4=numeric) |
| `View.NbIso` | 10 | Number of intervals |
| `View.RangeType` | 1 | Scale (1=default, 2=custom, 3=per step) |
| `View.CustomMin` | 0 | Custom minimum value |
| `View.CustomMax` | 0 | Custom maximum value |
| `View.ShowScale` | 1 | Show color scale bar |
| `View.ShowElement` | 0 | Show element boundaries |
| `View.VectorType` | 4 | Vector display (1=segment, 2=arrow, 4=3D arrow, 5=displacement) |
| `View.GlyphLocation` | 1 | Glyph location (1=centroid, 2=node) |
| `View.Light` | 1 | Enable lighting |
| `View.SmoothNormals` | 0 | Smooth normals |
| `View.ArrowScale` | 20 | Arrow scale factor |

## Post-Processing Options
| Option | Default | Description |
|--------|---------|-------------|
| `PostProcessing.Link` | 0 | Link mode between views |
| `PostProcessing.Binary` | 0 | Binary post-processing files |
| `PostProcessing.ForceNodeData` | 0 | Force NodeData format |
| `PostProcessing.SaveMesh` | 1 | Save mesh when exporting |
"""

GMSH_MSH_FORMAT = """
# GMSH .msh File Format

## Version 2.2 (Standard for Radia/NGSolve)

```
$MeshFormat
2.2 0 8
$EndMeshFormat

$PhysicalNames
numNames
dimension physicalTag "name"
...
$EndPhysicalNames

$Nodes
numNodes
nodeNumber x y z
...
$EndNodes

$Elements
numElements
elmNumber elmType numTags tag1 tag2 ... nodeList
...
$EndElements

$NodeData
numStringTags
"viewName"
numRealTags
timeValue
numIntegerTags
timeStepIndex
numComponents        (1=scalar, 3=vector, 9=tensor)
numNodes
partitionIndex       (0=none)
nodeNumber value...
...
$EndNodeData
```

## Version 4.1 (Optional Large-Scale Output)

```
$MeshFormat
4.1 0 8
$EndMeshFormat

$PhysicalNames
numPhysicalNames
dimension physicalTag "name"
...
$EndPhysicalNames

$Entities
numPoints numCurves numSurfaces numVolumes
pointTag X Y Z numPhysicalTags physicalTag...
curveTag minX minY minZ maxX maxY maxZ numPhysicalTags physicalTag... numBoundingPoints pointTag...
surfaceTag minX minY minZ maxX maxY maxZ numPhysicalTags physicalTag... numBoundingCurves curveTag...
volumeTag minX minY minZ maxX maxY maxZ numPhysicalTags physicalTag... numBoundingSurfaces surfaceTag...
$EndEntities

$Nodes
numEntityBlocks numNodes minNodeTag maxNodeTag
entityDim entityTag parametric numNodesInBlock
  nodeTag
  ...
  x y z
  ...
...
$EndNodes

$Elements
numEntityBlocks numElements minElementTag maxElementTag
entityDim entityTag elementType numElementsInBlock
  elementTag nodeTag...
  ...
...
$EndElements
```

NodeData/ElementData sections are identical in v2.2 and v4.1.

## Element Type Codes

### Linear (Order 1)
| Code | Type | Nodes |
|------|------|-------|
| 1 | Line2 | 2 |
| 2 | Tri3 | 3 |
| 3 | Quad4 | 4 |
| 4 | Tet4 | 4 |
| 5 | Hex8 | 8 |
| 6 | Prism6 | 6 |
| 7 | Pyr5 | 5 |
| 15 | Point1 | 1 |

### Quadratic (Order 2)
| Code | Type | Nodes |
|------|------|-------|
| 8 | Line3 | 3 |
| 9 | Tri6 | 6 |
| 10 | Quad9 | 9 |
| 11 | Tet10 | 10 |
| 12 | Hex27 | 27 |
| 13 | Prism18 | 18 |
| 14 | Pyr14 | 14 |
| 16 | Quad8 | 8 (serendipity) |
| 17 | Hex20 | 20 (serendipity) |

### High-Order Triangles (Used by GmshPostExport)
| Code | Type | Nodes | Order | Nodes/edge | Interior |
|------|------|-------|-------|------------|----------|
| 2 | Tri3 | 3 | 1 | 0 | 0 |
| 9 | Tri6 | 6 | 2 | 1 | 0 |
| 21 | Tri10 | 10 | 3 | 2 | 1 |
| 23 | Tri15 | 15 | 4 | 3 | 3 |
| 25 | Tri21 | 21 | 5 | 4 | 6 |

### High-Order Tetrahedra
| Code | Type | Nodes | Order |
|------|------|-------|-------|
| 4 | Tet4 | 4 | 1 |
| 11 | Tet10 | 10 | 2 |
| 29 | Tet20 | 20 | 3 |
| 30 | Tet35 | 35 | 4 |
| 31 | Tet56 | 56 | 5 |

### High-Order Hexahedra
| Code | Type | Nodes | Order |
|------|------|-------|-------|
| 5 | Hex8 | 8 | 1 |
| 12 | Hex27 | 27 | 2 |
| 92 | Hex64 | 64 | 3 |
| 93 | Hex125 | 125 | 4 |

## Node Ordering Convention (High-Order)
1. Corner vertices (same as linear element)
2. Internal nodes for each edge (equispaced, low-to-high vertex index)
3. Internal nodes for each face (recursive)
4. Volume internal nodes

## Element Tags (v2.2)
- Tag 1: physical entity tag (Physical Group)
- Tag 2: elementary entity tag (geometry entity)
- Optional tag 3+: partition information

## Physical Groups
- Define material regions or boundary conditions
- By default, only elements in physical groups are saved
- Use `Mesh.SaveAll = 1` to save all elements
"""

GMSH_GEO_SCRIPTING = """
# GMSH .geo Scripting Language

## Purpose in Radia
.geo files are used as **companion files** for display settings when
opening .msh results. They are NOT used for mesh generation.

## Companion .geo File Pattern
```
// display_settings.geo - companion for results.msh
Merge "results.msh";

// Required for high-order curved elements
Mesh.NumSubEdges = 4;

// Clean display
Mesh.VolumeEdges = 0;
Mesh.SurfaceEdges = 0;

// View settings
View[0].IntervalsType = 2;  // Continuous
View[0].ShowScale = 1;
View[0].VectorType = 4;     // 3D arrows
```

## Multi-File Overlay
```
// overlay.geo - STEP geometry + field results
Merge "coil.step";
Merge "magnet.step";
Merge "field.msh";

Mesh.NumSubEdges = 4;
Mesh.VolumeEdges = 0;

// Field view settings
View[0].IntervalsType = 2;
View[0].ShowScale = 1;
```

## Basic Syntax Reference
```
// Variables
lc = 0.01;
pts[] = {1, 2, 3};

// Geometry primitives
Point(tag) = {x, y, z, meshSize};
Line(tag) = {startPt, endPt};
Circle(tag) = {startPt, centerPt, endPt};
Curve Loop(tag) = {curve_list};  // signs = orientation
Plane Surface(tag) = {loop_list};
Surface Loop(tag) = {surface_list};
Volume(tag) = {surfaceLoop_list};

// Physical groups
Physical Point("name", tag) = {point_list};
Physical Curve("name", tag) = {curve_list};
Physical Surface("name", tag) = {surface_list};
Physical Volume("name", tag) = {volume_list};

// Transformations
Translate {dx,dy,dz} { Point{1}; }
Rotate {{ax,ay,az}, {px,py,pz}, angle} { Volume{1}; }

// Extrusion
Extrude {dx,dy,dz} { Surface{1}; Layers{n}; Recombine; }

// Transfinite (structured)
Transfinite Curve {list} = n Using Progression p;
Transfinite Surface {list};
Transfinite Volume {list};
Recombine Surface {list};

// Control flow
If (expr) ... ElseIf (expr) ... Else ... EndIf
For i In {start:end:step} ... EndFor

// Include and Merge
Include "other.geo";
Merge "file.msh";
Merge "file.step";

// Mesh commands
Mesh 2;          // Generate 2D mesh
SetOrder 2;      // Set element order
Save "file.msh"; // Save mesh

// Options
General.Verbosity = 5;
Mesh.NumSubEdges = 4;
View[0].IntervalsType = 2;

// Output
Printf("format %g", value);
```

## Parsed Post-Processing Views
```
View "B field" {
  // Scalar point
  SP(x,y,z){val};

  // Vector point
  VP(x,y,z){vx,vy,vz};

  // Scalar triangle (3 vertices, 3 values)
  ST(x1,y1,z1, x2,y2,z2, x3,y3,z3){v1,v2,v3};

  // Vector triangle
  VT(x1,y1,z1, x2,y2,z2, x3,y3,z3){vx1,vy1,vz1, vx2,vy2,vz2, vx3,vy3,vz3};

  // Element type prefixes: S=scalar, V=vector, T=tensor
  // Geometry suffixes: P=point, L=line, T=tri, Q=quad,
  //                    S=tet, H=hex, I=prism, Y=pyramid

  TIME {t1, t2, ...};
};
```

Note: Parsed views are inefficient for large datasets.
Use .msh NodeData/ElementData instead.
"""

GMSH_HIGH_ORDER = """
# High-Order Element Display in GMSH

## Critical Setting: Mesh.NumSubEdges

**MUST set `Mesh.NumSubEdges = 4`** to correctly render curved high-order
elements. Default value (2) draws nearly straight edges.

### How to set:
1. **Command line**: `gmsh file.msh -numsubedges 4`
2. **In .geo file**: `Mesh.NumSubEdges = 4;`
3. **GMSH console**: Type `Mesh.NumSubEdges = 4;`
4. **GUI**: Tools -> Options -> Mesh -> NumSubEdges

## GmshPostExport (Radia -> GMSH)

`GmshPostExport` writes high-order curved mesh + field data to .msh.

```python
from radia.gmsh_post_export import GmshPostExport

# BEM/FEM surface visualization
post = GmshPostExport(mesh, boundary=True)  # BND from volume mesh
post.add_field("|J|", node_J, ncomp=1)      # per-vertex scalar
post.add_vector_field("J", gf_J)            # vector field
post.write("results.msh")                   # default v2.2
post.write("results.msh", version="4.1")    # optional v4.1
```

## Supported Triangle Orders
| mesh.Curve(p) | GMSH Type | Nodes | Code |
|---------------|-----------|-------|------|
| p=1 | Tri3 | 3 | 2 |
| p=2 | Tri6 | 6 | 9 |
| p=3 | Tri10 | 10 | 21 |
| p=4 | Tri15 | 15 | 23 |
| p=5 | Tri21 | 21 | 25 |

## Companion .geo File Template
```
// Automatically generated companion file
Merge "results.msh";
Mesh.NumSubEdges = 4;
Mesh.VolumeEdges = 0;
View[0].IntervalsType = 2;
```

## High-Order Node Extraction
Node positions are extracted via `mesh.GetTrafo(el)` evaluated at GMSH
reference coordinates. H1 GridFunction approach is unreliable for p>=4.
"""

GMSH_RADIA_WORKFLOW = """
# GMSH Workflow in Radia Project

## 1. Field Visualization (Primary Use)

```
NGSolve FEM solve -> GridFunction -> GmshPostExport -> .msh -> GMSH GUI
```

### Steps:
1. Solve FEM problem in NGSolve
2. Use `GmshPostExport` to write .msh with field data
3. Open in GMSH with `Mesh.NumSubEdges = 4`

## 2. Geometry + Field Overlay

```
CoilBuilder -> .step (coil geometry)
OCC shapes  -> .step (magnet/iron geometry)
NGSolve     -> .msh  (field data)
-> GMSH merges all for combined visualization
```

### .geo companion:
```
Merge "coil.step";
Merge "magnet.step";
Merge "field.msh";
Mesh.NumSubEdges = 4;
Mesh.VolumeEdges = 0;
```

## 3. Cubit Panel Integration

The Cubit panel opens GMSH for results display:
- Output: `examples/cubit_panels/<solver>/results/`
- Uses `pythonw.exe` to launch GMSH (no console window)
- .geo companion file auto-generated with display settings

## 4. BEM Inductance Results

Unified .msh output containing:
- Volume B field (NodeData)
- Surface J (NodeData)
- Coil wireframe (1D elements)

Node ID separation:
- Volume nodes: 1..nv_vol
- Surface nodes: nv_vol+1..end

Physical groups: "air", "coil_surface", "coil_wire"

## 5. GMSH .msh as Input to NGSolve

Cubit exports .msh v2.2 -> NGSolve `ReadGmsh()` reads it.
This is a pure file format path -- no GMSH dependency needed.

```python
from ngsolve import *
mesh = Mesh(unit_square.GenerateMesh())  # or:
# mesh = Mesh("cubit_export.msh")  # ReadGmsh() for .msh v2.2
```
"""

GMSH_ONELAB = """
# ONELAB Distribution (S:\\ONELAB)

GMSH is distributed as part of the ONELAB suite.

## Structure
```
S:\\ONELAB\\
  00_installer/
    gmsh.exe          (4.15.2, standalone executable)
    getdp.exe         (FEM solver, not used by Radia)
    onelab.py          (ONELAB Python API)
    tutorials/         (60 .geo tutorial files, t1-t21)
    examples/          (78 .geo examples)
    models/            (18 physics domains, 142 .geo files)
      ElectricMachines/  (PMSM, IM, SRM, WFSM models)
      Magnets/           (permanent magnet models)
      Inductor/          (inductor analysis)
      ...
    templates/         (GetDP problem templates)
```

## Tutorial Files (t1-t21)
Key tutorials for Radia users:
- t1.geo: Basic geometry and mesh
- t2.geo: Transformations, ruled surfaces
- t5.geo: Mesh sizes, attractors
- t6.geo: Transfinite meshing
- t8.geo: Post-processing, scalar/vector views
- t10.geo: Mesh size fields
- t13.geo: Remeshing
- t16.geo: OCC geometry kernel
- t20.geo: STEP/IGES import
- t21.geo: Physical groups

## Electric Machine Models
Complete models with .geo + .pro (GetDP problem definition):
- PMSM (8-pole, GRUCAD): pmsm.geo, pmsm.pro
- Induction Machine (4-pole, 3kW): im.geo, im.pro
- SRM (4-pole): srm.geo, srm.pro
- WFSM (4-pole): wfsm.geo, wfsm.pro
"""

GMSH_PITFALLS = """
# GMSH Common Pitfalls

## 1. Curved Elements Look Straight
**Problem**: High-order elements display as straight-edged polygons.
**Fix**: Set `Mesh.NumSubEdges = 4` (default is 2, too coarse).

## 2. Missing Elements in Output
**Problem**: .msh file has no elements or missing elements.
**Fix**: Either define Physical Groups, or use `Mesh.SaveAll = 1`.
By default, only elements in Physical Groups are saved.

## 3. Quit Confirmation Dialog
**Problem**: "Do you really want to quit?" dialog appears every time.
**Fix**: Set `General.ConfirmQuit = 0` in `%APPDATA%/gmsh.conf`.

## 4. Console Encoding on Windows (cp932)
**Problem**: Unicode characters cause errors in Japanese Windows.
**Fix**: Use ASCII only in GMSH scripts and output.

## 5. .msh v2.2 vs v4.1 Confusion
**Problem**: NGSolve ReadGmsh() expects v2.2 but GMSH defaults to v4.1.
**Fix**: Use `Mesh.MshFileVersion = 2.2` or `-format msh22`.

## 6. NodeData Field Not Displayed
**Problem**: .msh file has NodeData section but no view appears.
**Fix**: Ensure node IDs in NodeData match the $Nodes section.
Check numComponents (1, 3, or 9).

## 7. STEP File Shows No Mesh
**Problem**: Merged STEP file doesn't show mesh/field data.
**Fix**: STEP is geometry only. Merge .msh file separately for fields.
Use `Merge "file.step"; Merge "field.msh";` in .geo.

## 8. gmsh.bat vs gmsh.exe
**Problem**: `pip install gmsh` creates gmsh.bat (Python wrapper).
**Fix**: Radia uses standalone gmsh.exe (`C:\\gmsh.exe`).
The Python gmsh package is NOT installed.
"""


def get_gmsh_documentation(topic: str = "all") -> str:
    """Return GMSH usage documentation by topic.

    Args:
        topic: One of: all, policy, overview, cli, shortcuts, options,
               msh_format, geo, high_order, workflow, onelab, pitfalls

    Returns:
        Documentation string for the requested topic.
    """
    topics = {
        "policy": GMSH_RADIA_POLICY,
        "overview": GMSH_OVERVIEW,
        "cli": GMSH_COMMAND_LINE,
        "shortcuts": GMSH_KEYBOARD_SHORTCUTS,
        "options": GMSH_OPTIONS,
        "msh_format": GMSH_MSH_FORMAT,
        "geo": GMSH_GEO_SCRIPTING,
        "high_order": GMSH_HIGH_ORDER,
        "workflow": GMSH_RADIA_WORKFLOW,
        "onelab": GMSH_ONELAB,
        "pitfalls": GMSH_PITFALLS,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        available = ", ".join(topics.keys())
        return f"Unknown topic: '{topic}'. Available: all, {available}"
