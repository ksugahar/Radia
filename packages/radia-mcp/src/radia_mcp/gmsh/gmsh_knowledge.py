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
- Reading .msh format for field data (GmshPostExport output)

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
- Default: .msh v4.1 (GMSH visualization)
- `GmshPostExport.write()` and `vol2msh()` emit .msh v4.1
- Post-processing launch artifact: `case.geo` that merges the .msh/.step data
- Exact option sidecars: `case.geo.opt` for normal launch and `case.msh.opt`
  for raw mesh/data inspection. A plain `case.opt` is not auto-loaded.
- `cubit-mesh-export` `export gmsh "case.msh"` should attach `case.geo`,
  `case.geo.opt`, and `case.msh.opt`; open `case.geo` for normal review.
- For acoustic/FEM-BEM post artifacts, use the shared
  `gmsh_post_display_contract` / `write_gmsh_post_launch_artifact` pattern:
  one MSH v4.1 data file, a `.geo` launch file, exact `.geo.opt` and
  `.msh.opt` sidecars, named views, Z-up camera metadata, and optional
  cut-plane metadata (`General.Clip0A/B/C/D`).  This is the cross-language
  contract shared with the MATLAB/Gypsilab readable acoustic lane.
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
The canonical launcher is `gmsh` on PATH. On current LAB installs this is the
Python wrapper `C:\\Program Files\\Python312\\Scripts\\gmsh.bat`, backed by the
`gmsh` Python package (`gmsh.py`) at version 4.15.2. Do not hard-code
`C:\\gmsh.exe`; use `gmsh` / `shutil.which("gmsh")` / the registered file
association.
- `.geo` is the primary file association for Radia post-processing launch
- `.msh` association is optional raw mesh/data inspection; do not make it the
  user-facing post-processing contract
- Python scripts may import `gmsh` for post-processing inspection, but Radia
  computation scripts must not use `gmsh.model.*` for geometry or mesh creation

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
gmsh display.geo                        # auto-loads display.geo.opt if present
```

## Windows PowerShell `-string` quoting
When passing a Gmsh script snippet with a Windows path through PowerShell,
prefer single quotes inside the Gmsh snippet:

```powershell
gmsh display.geo -string "Print 'C:/temp/frame.png'; Exit;"
```

Do not pass backslash-escaped double quotes such as `\"C:/temp/frame.png\"`
through `gmsh.bat`; the wrapper can leave Gmsh seeing an unquoted `W:` token,
which fails as `Unknown variable 'W'`. `Print` also needs a graphical
OpenGL/FLTK context to create PNG/JPEG images; for unattended animation export,
use the Python API path in the `animation` topic.

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
- .msh.opt file (auto-loaded alongside .msh, see GMSH_OPT_FILE)

## General Options
| Option | Default | Description |
|--------|---------|-------------|
| `General.NumThreads` | 1 | Parallel threads |
| `General.Verbosity` | 5 | Verbosity (0-99) |
| `General.Terminal` | 0 | Terminal output |
| `General.Orthographic` | 1 | Orthographic projection (1=ortho, 0=perspective) |
| `General.ConfirmQuit` | 1 | Ask before quit (set 0 to disable) |
| `General.SmallAxes` | 1 | Show small axes indicator in corner |
| `General.Axes` | 0 | Show full axes (0=none, 1=simple axes, 2=box) |
| `General.AxesMikado` | 0 | Mikado-style axes |
| `General.AlphaBlending` | 1 | Enable alpha (transparency) blending |
| `General.Antialiasing` | 0 | Antialiasing |
| `General.BackgroundGradient` | 0 | Background gradient (0=none, 1=vertical, 2=horizontal, 3=spherical) |
| `General.ColorScheme` | 1 | Color scheme (0=dark, 1=light) |
| `General.GraphicsFont` | "Helvetica" | Font name |
| `General.GraphicsFontSize` | 15 | Font size for labels |
| `General.GraphicsFontSizeTitle` | 18 | Font size for titles |

## General.Light Options
| Option | Default | Description |
|--------|---------|-------------|
| `General.Light0` | 1 | Enable light 0 |
| `General.Light0X` | 0.65 | Light 0 X direction |
| `General.Light0Y` | 0.65 | Light 0 Y direction |
| `General.Light0Z` | 1 | Light 0 Z direction |
| `General.Light0W` | 0 | Light 0 W (0=directional, 1=positional) |
| `General.Light1..5` | 0 | Additional lights (disabled by default) |

## General.Color Options
| Option | Description |
|--------|-------------|
| `General.Color.Background` | Background color e.g. `{255,255,255}` (white) |
| `General.Color.BackgroundGradient` | Gradient end color |
| `General.Color.Foreground` | Foreground (axes, text borders) color |
| `General.Color.Text` | Text color |

## Camera / Trackball Options (saved in .msh.opt)
| Option | Description |
|--------|-------------|
| `General.Trackball` | 1=trackball rotation mode |
| `General.TrackballQuaternion0..3` | Rotation as unit quaternion (w,x,y,z) |
| `General.RotationX/Y/Z` | Rotation angles in degrees |
| `General.RotationCenterGravity` | 1=rotate around model center of gravity |
| `General.ScaleX/Y/Z` | Zoom scale factors |
| `General.TranslationX/Y/Z` | Pan translation |

## Clipping Plane Options
| Option | Description |
|--------|-------------|
| `General.Clip0A` | Clip plane 0 normal X (equation: Ax+By+Cz+D=0) |
| `General.Clip0B` | Clip plane 0 normal Y |
| `General.Clip0C` | Clip plane 0 normal Z |
| `General.Clip0D` | Clip plane 0 offset D (positive = shift along normal) |
| `General.ClipFactor` | Clip box size factor (default 5) |
| `General.ClipWholeElements` | 1=clip whole elements (don't cut through) |
| `General.ClipOnlyVolume` | 1=clip only volume elements (not surfaces) |
| `General.ClipOnlyDrawIntersectingVolume` | 0=draw all, 1=only intersecting |
| `Mesh.Clip` | 1=enable mesh clipping |

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
| `Mesh.ColorCarousel` | 1 | Coloring mode (0=solid, 1=by element type, 2=by physical group, 3=by partition) |
| `Mesh.Light` | 1 | Enable lighting on mesh |
| `Mesh.LightLines` | 2 | Light on lines (0=off, 1=on, 2=two-side) |
| `Mesh.LightTwoSide` | 1 | Two-sided lighting |
| `Mesh.Lines` | 0 | Show 1D elements |
| `Mesh.SurfaceEdges` | 1 | Show surface edges |
| `Mesh.SurfaceFaces` | 0 | Show surface faces |
| `Mesh.VolumeEdges` | 1 | Show volume edges |
| `Mesh.VolumeFaces` | 0 | Show volume faces |
| `Mesh.SmoothNormals` | 0 | Smooth normals for rendering |
| `Mesh.Nodes` | 0 | Show mesh nodes |
| `Mesh.Normals` | 0 | Show normal vectors |

## Mesh.Color Options (for ColorCarousel=2, physical group coloring)
Colors are assigned by physical group index (Zero=group 0, One=group 1, ...):
```
Mesh.Color.Nodes       = {0,0,255};
Mesh.Color.Lines       = {0,0,0};
Mesh.Color.Triangles   = {160,150,255};
Mesh.Color.Quadrangles = {130,120,225};
Mesh.Color.Tetrahedra  = {160,150,255};
Mesh.Color.Hexahedra   = {130,120,225};
Mesh.Color.Prisms      = {232,210,23};
Mesh.Color.Pyramids    = {217,113,38};
Mesh.Color.Zero        = {255,120,0};
Mesh.Color.One         = {204,38,38};
// ... up to Mesh.Color.Nineteen
```

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


GMSH_OPT_FILE = """
# GMSH .msh.opt File

## Overview

A `.msh.opt` file is a GMSH options file placed alongside a `.msh` file
with the same base name (e.g. `model.msh` -> `model.msh.opt`).
GMSH automatically loads it when the .msh file is opened, restoring the
exact camera angle, clipping, colors, and visibility from the last session.

## Auto-Load Mechanism

```
model.msh        <- mesh/field data
model.msh.opt    <- automatically loaded when model.msh is opened
```

GMSH searches for `<filename>.opt` in the same directory as `<filename>`.
This works for any file type: `.msh.opt`, `.geo.opt`, `.step.opt`, etc.

Radia post-processing should prefer `case.geo` as the launch target. Keep
`case.msh` as the raw mesh/field container, and attach durable display state to
`case.geo` either inline or through the exact sidecar `case.geo.opt`.

Important naming detail: `Gmsh` appends `.opt` to the exact filename that is
opened. Opening `case.msh` looks for `case.msh.opt`; opening `case.geo` looks
for `case.geo.opt`. A sidecar named only `case.opt` is not an auto-load
contract for either file. For double-click workflows, artifact writers should
emit `case.geo.opt` next to `case.geo`, or copy the critical options into the
`.geo` itself. A launcher may mirror old `case.opt` sidecars to `case.geo.opt`
for backward compatibility, but the durable contract is the exact auto-load
name.

## Generating a .msh.opt File

From the GMSH GUI:
1. Open the .msh file and adjust the view interactively
2. **File -> Save Options As Default** writes `~/.gmsh` (global)
3. **File -> Save Session State** saves a `.opt` file next to the current file
4. Or: **Tools -> Options -> General -> Save** exports current settings

From command line / .geo script:
```
// Save current options to file
Save "model.msh.opt";
```

## Structure of a .msh.opt File

A `.msh.opt` file uses the same syntax as `.geo` scripts.
It is a sequence of option assignments and GMSH commands:

```
// --- Projection ---
General.Orthographic = 1;       // 1=orthographic, 0=perspective
General.AlphaBlending = 1;
General.BackgroundGradient = 0; // 0=flat background

// --- Background / foreground colors ---
General.Color.Background = {255,255,255};  // white
General.Color.Foreground = {85,85,85};

// --- Axes ---
General.Axes = 0;        // hide full axes
General.SmallAxes = 0;   // hide corner axes indicator

// --- Light ---
General.Light0 = 1;
General.Light0X = 0.65;
General.Light0Y = 0.65;
General.Light0Z = 1;

// --- Camera (trackball) ---
General.Trackball = 1;
General.TrackballQuaternion0 = 0.2043;  // rotation quaternion
General.TrackballQuaternion1 = 0.0937;
General.TrackballQuaternion2 = 0.1295;
General.TrackballQuaternion3 = 0.9658;
General.ScaleX = 4.0;
General.ScaleY = 4.0;
General.ScaleZ = 4.0;
General.TranslationX = 0;
General.TranslationY = 0;

// --- Z-up x-z plane post view ---
// Useful for acoustic/BEM post artifacts with an x-z pressure plane
// and a 3-D drum/body surface. Z points upward in the corner axis triad.
General.Trackball = 0;
General.RotationX = -68;
General.RotationY = 0;
General.RotationZ = 0;
General.RotationCenterGravity = 1;

// --- Hide entities by physical group tag ---
Hide { Volume{5, 6}; }  // hide air_gap (5) and air (6)

// --- Clipping plane (x=-0.02565 plane, clip left side) ---
General.Clip0A = -1;    // normal direction (-X)
General.Clip0B = 0;
General.Clip0C = 0;
General.Clip0D = 0.02565;  // offset (distance from origin)
General.ClipFactor = 5;
General.ClipWholeElements = 1;    // don't cut through elements
General.ClipOnlyVolume = 1;       // clip only volume elements
General.ClipOnlyDrawIntersectingVolume = 0;
Mesh.Clip = 1;           // activate clipping on mesh

// --- Mesh display ---
Mesh.VolumeEdges = 1;
Mesh.VolumeFaces = 1;
Mesh.SurfaceEdges = 0;
Mesh.SurfaceFaces = 0;
Mesh.SmoothNormals = 1;
Mesh.Light = 1;
Mesh.LightTwoSide = 1;
Mesh.ColorCarousel = 2;  // color by physical group

// --- Mesh colors (by physical group index) ---
Mesh.Color.Tetrahedra = {160,150,255};
Mesh.Color.Zero = {255,120,0};
Mesh.Color.One  = {204,38,38};
// ... up to Nineteen
```

For a strict front-on x-z view along the y-axis, use
`General.RotationX = -90`, but 3-D bodies will look flatter. Keep the same
camera in `case.msh.opt` for raw mesh inspection, or users will see a different
orientation after opening the raw `.msh`.

## Shared Acoustic/FEM-BEM Display Contract

Radia acoustic and MATLAB/Gypsilab post artifacts should use the same durable
launch structure:

```
case.msh          # single Gmsh MSH v4.1 result/data container
case.geo          # launch target, Merge "case.msh"
case.geo.opt      # exact autoload sidecar for post display
case.msh.opt      # exact autoload sidecar for raw mesh/data inspection
case.display.json # manifest with camera, cut-plane, views, and schema
```

Use `gmsh_post_display_contract` to plan or validate this structure from MCP,
and `write_gmsh_post_launch_artifact` from Python helper code when Radia writes
the files directly.  The MATLAB/Gypsilab equivalent is
`writeGmshPostLaunchArtifact`.  Both record named scalar/vector/displacement
views, the Z-up x-z camera preset (`General.RotationX = -68`,
`General.RotationZ = 0`), and optional cut planes via
`General.Clip0A/B/C/D` plus `Mesh.Clip = 1`.

## Hide / Show Commands

```
// Hide specific entities
Hide { Volume{5, 6}; }       // hide volumes with tags 5, 6
Hide { Surface{3}; }          // hide surface
Hide "*";                     // hide all entities

// Show specific entities
Show { Volume{1, 2, 3, 4}; }
Show "*";                     // show all entities
```

These commands are persistent in a `.opt` file — GMSH replays them on load.

## Clipping Plane Formula

The clip plane equation is: `A*x + B*y + C*z + D >= 0` (visible side).

Examples:
```
// Clip at x = 0.02565 (show x > 0.02565)
General.Clip0A = -1;  General.Clip0D = 0.02565;  // -x + 0.02565 >= 0 -> x <= 0.02565
// Actually: visible when -x + D >= 0, i.e. x <= D

// Clip at z = 0 (show z > 0, i.e. upper half)
General.Clip0A = 0;  General.Clip0B = 0;  General.Clip0C = 1;  General.Clip0D = 0;
```

`General.ClipWholeElements = 1` is recommended for mesh display — it avoids
partially-cut tetrahedra which look bad. Only whole elements on the visible
side are drawn.

## Radia Project Usage

In the Radia project, `.msh.opt` files are used to set up a reproducible
visualization state for mesh review:

- **White background** (`General.Color.Background = {255,255,255}`)
- **Orthographic projection** (`General.Orthographic = 1`)
- **Air volumes hidden** (`Hide { Volume{5, 6}; }` for air_gap/air)
- **X-plane clipping** to show internal cross-section
- **ColorCarousel = 2** (color by physical group) for material identification
- **No axes** for clean screenshots

This state is saved per-.msh file so each mesh file remembers its own view.

## Programmatic Generation (Python)

```python
def write_msh_opt(msh_path: str, hidden_volumes: list[int],
                  clip_x: float | None = None) -> None:
    opt_path = msh_path + ".opt"
    lines = [
        "General.Orthographic = 1;",
        "General.Color.Background = {255,255,255};",
        "General.Color.Foreground = {85,85,85};",
        "General.Axes = 0;",
        "General.SmallAxes = 0;",
        "General.BackgroundGradient = 0;",
        "Mesh.Light = 1;",
        "Mesh.LightTwoSide = 1;",
        "Mesh.VolumeEdges = 1;",
        "Mesh.VolumeFaces = 1;",
        "Mesh.SurfaceEdges = 0;",
        "Mesh.SurfaceFaces = 0;",
        "Mesh.SmoothNormals = 1;",
        "Mesh.ColorCarousel = 2;",
    ]
    if hidden_volumes:
        tags = ", ".join(str(v) for v in hidden_volumes)
        lines.append(f"Hide {{ Volume{{{tags}}}; }}")
    if clip_x is not None:
        lines += [
            f"General.Clip0A = -1;",
            f"General.Clip0B = 0;",
            f"General.Clip0C = 0;",
            f"General.Clip0D = {clip_x};",
            "General.ClipFactor = 5;",
            "General.ClipWholeElements = 1;",
            "General.ClipOnlyVolume = 1;",
            "General.ClipOnlyDrawIntersectingVolume = 0;",
            "Mesh.Clip = 1;",
        ]
    with open(opt_path, "w") as f:
        f.write("\\n".join(lines) + "\\n")
```
"""

GMSH_MSH_FORMAT = """
# GMSH .msh File Format

## Version 4.1 (Lab-wide standard, 2026-04)

v4.1 is the only supported format across the Radia repository
(Cubit plugin output, NGSolve post-processing, all panels).
netgen I/O is always via `.vol` (never `.msh`).

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
.geo files are the **standard Radia post-processing launch artifact**.
Post-processing exporters should write a `case.geo` recipe next to the `.msh`
created by Radia's existing exporters (`GmshPostExport.write()`, `vol2msh()`,
or combined .msh writers). The `.geo` file should `Merge` the raw `.msh`,
`.step`, and field files, then carry the display policy needed for Explorer
double-click, LLM/headless review, and reproducible screenshots. The raw `.msh`
remains the mesh/field data container; the user-facing Open GMSH target is
`case.geo`.

They are NOT used for mesh generation.

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

## .geo + .opt Launch Contract

Do not rely on a plain `display.opt` being auto-loaded when a user double-clicks
`display.geo`. Gmsh auto-loads `display.geo.opt` for `display.geo`, not
`display.opt`. For reproducible lab artifacts, use one of these two contracts:

1. Put all critical display options directly in the `.geo` file after the
   `Merge` lines. This is the safest double-click path.
2. If you also keep a plain `display.opt` sidecar, the artifact writer should
   emit the exact auto-load twin `display.geo.opt` at the same time.

For clipped acoustic/FEM-BEM fields, the `.geo` should contain the clip options
itself if it is meant to be opened from Explorer:

```
Merge "field.msh";
Mesh.Clip = 0;       // keep mesh/surface geometry whole
General.Clip0A = 0;
General.Clip0B = -1;  // visible side y < 0 for the lab convention
General.Clip0C = 0;
General.Clip0D = 0;
General.ClipOnlyVolume = 1;  // keep separate surface/drum geometry whole
General.ClipWholeElements = 0;
View[0].Visible = 1;
View[0].Clip = 1;   // clip the post-processing acoustic pressure view
View[1].Visible = 0;
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

For directories with many related high-order examples, keep one shared
`_gmsh_display.geo` companion in that directory:

```
// Shared GMSH display companion
Mesh.NumSubEdges = 4;
// Merge "<result>.msh";
```

## GmshPostExport (Radia -> GMSH)

`GmshPostExport` writes high-order curved mesh + field data to .msh.

```python
from radia.gmsh_post_export import GmshPostExport

# BEM/FEM surface visualization
post = GmshPostExport(mesh, boundary=True)  # BND from volume mesh
post.add_field("|J|", node_J, ncomp=1)      # per-vertex scalar
post.add_vector_field("J", gf_J)            # vector field
post.write("results.msh")                   # v4.1 lab standard
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

## 2. Netgen/Cubit Mesh, GMSH Display

The mesh generator owns the mesh. GMSH owns the view.

```
Netgen/Cubit -> tri/tet or hex mesh -> NGSolve/Radia validation
             -> GMSH .msh/.pos data + .geo launch recipe -> standalone GMSH
```

For display examples, write a GMSH v4.1 `.msh` from the existing Radia/Cubit
exporters. That is a **mesh export**, not GMSH mesh generation. Keep the script
free of `import gmsh` and `gmsh.model.*`; write `.vol` for the solver path and
`.msh`/`.pos`/`.geo` only for display. Historical v2.2 snippets are legacy
references, not the preferred Radia output.
For post-processing, `case.geo` is the standard artifact to open in GMSH; it
should `Merge` the exported data and contain the critical display options or
the exact `case.geo.opt` twin.

### .vol and GMSH

Do not plan a Radia workflow around GMSH directly opening Netgen `.vol` files,
and do not add a GMSH-side `.vol` reader/plugin for normal Radia operation.
Radia already owns the `.msh` output path: use `GmshPostExport.write()` for an
NGSolve mesh with fields, `vol2msh()` for `.vol` plus field payloads, and
combined v4.1 writers for multi-mesh post-processing. Cubit-side mesh display
uses `export gmsh "case.msh" order N` / `cubit_mesh_export.export_Gmsh_ver4`.
On the lab GMSH 4.15.2 package, opening existing `.vol` files fails as a syntax
parse of a GMSH script, so `.vol` is not a portable GMSH input contract. The
supported path is:

```
Netgen/Cubit .vol -> NGSolve/Radia reads the mesh -> existing Radia .msh exporter -> case.geo
```

Reference workflow:
`docs/visualization/MESH_GUIDE.md`, plus `docs/visualization/_gmsh_display.geo`
for the minimal standalone display companion.

## 3. Geometry + Field Overlay

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

## 4. Cubit Panel Integration

The Cubit panel opens GMSH for results display:
- Output: panel or validation run directories; do not use the retired
  examples staging tree (for example `runs/radia_*`, `panels/samples/...`, or
  `validation_test/...` depending on owner)
- Uses `pythonw.exe` to launch GMSH (no console window)
- .geo companion file auto-generated with display settings

## 5. BEM Inductance Results

Unified .msh output containing:
- Volume B field (NodeData)
- Surface J (NodeData)
- Coil wireframe (1D elements)

Node ID separation:
- Volume nodes: 1..nv_vol
- Surface nodes: nv_vol+1..end

Physical groups: "air", "coil_surface", "coil_wire"

## 5. GMSH .msh as Input to NGSolve

Cubit exports .vol via `export netgen` -> NGSolve reads it directly.
The .vol file is the standard interface between Cubit and NGSolve.

```python
from ngsolve import *
mesh = Mesh("cubit_export.vol")  # .vol is the sole interface
```
"""

GMSH_ONELAB = """
# ONELAB Distribution (public-safe curated corpus)

GMSH is distributed as part of the ONELAB suite.

## Structure
```
public-safe curated corpus
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
**Problem**: .msh v2.2 vs v4.1 confusion when viewing in GMSH.
**Fix**: Radia, `GmshPostExport`, `vol2msh`, and `cubit-mesh-export` standardize
on `.msh v4.1`. Treat v2.2 snippets as legacy reference material only; do not
downgrade new Radia post-processing output to v2.2.

## 6. NodeData Field Not Displayed
**Problem**: .msh file has NodeData section but no view appears.
**Fix**: Ensure node IDs in NodeData match the $Nodes section.
Check numComponents (1, 3, or 9).

## 7. STEP File Shows No Mesh
**Problem**: Merged STEP file doesn't show mesh/field data.
**Fix**: STEP is geometry only. Merge .msh file separately for fields.
Use `Merge "file.step"; Merge "field.msh";` in .geo.

## 8. gmsh.bat vs gmsh.exe
**Problem**: A script assumes `C:\\gmsh.exe` and fails even though `gmsh` works
from the shell.
**Fix**: The Python-wrapper launcher (`gmsh.bat` backed by `gmsh.py`) is the
canonical current install. Resolve `gmsh` through PATH or `shutil.which("gmsh")`
instead of hard-coding an executable path.

## 9. .geo Opens but the .opt Settings Are Missing
**Problem**: A user double-clicks `display.geo` and sees default colors,
missing clip planes, wrong view visibility, or no y<0 cut even though
`display.opt` exists next to it.

**Cause**: Gmsh auto-loads an option file by appending `.opt` to the exact file
being opened. `display.geo` can auto-load `display.geo.opt`, but not
`display.opt`. Likewise `field.msh` can auto-load `field.msh.opt`, but not
`field.opt`.

**Fix**: For Explorer/double-click artifacts, either put the critical options
directly in `display.geo`, or write the exact Gmsh sidecar
`display.geo.opt` alongside it. For post-processing NodeData fields, remember
that `Mesh.Clip` clips mesh entities while `View[0].Clip` clips the displayed
field view. Use `Mesh.Clip = 0` plus `View[0].Clip = 1` when the geometry/drum
surface must remain whole but the acoustic pressure view should show a y<0
section. Keep `General.ClipOnlyVolume`, view visibility, colormap, and scale
settings in the `.geo` if the artifact must be portable without a sidecar.
Use a split viewer contract for CAE artifacts: `case.geo` is the post-processing
display recipe and should have `case.geo.opt`; `case.msh` is the raw mesh/data
inspection entry point and should have `case.msh.opt` with post views hidden and
mesh faces/edges visible. Radia post-processing exporters should therefore emit
`.geo` by default. A plain `case.opt` is only a compatibility mirror, not the
Explorer auto-load contract for either file.
Avoid relying on UserChoice; use HKLM/HKCR ProgID class associations that invoke
the `gmsh` command when Windows Explorer double-click behavior matters. Register
`.geo` as the primary Radia post-processing launch association. Keep `.msh`
association optional for raw mesh/data inspection only; panels and docs should
open `case.geo`.

## 10. GMSH window appears "invisible" -- it's on an off-screen monitor

**Symptom**: After running a viewer script (e.g. `python view_*_gmsh.py` or
`open_gmsh.py file.msh`), the user reports "GMSHが見えない" / "GMSH is
invisible" but `tasklist` shows python.exe is running.

**Root cause**: GMSH restores the LAST-USED window position from
`%APPDATA%/gmsh-options` (or its in-process cache). If the user EVER had
a second monitor attached at a coordinate like x=-2560, every subsequent
GMSH launch silently opens at THOSE coordinates -- which now correspond
to no physical display. The window is real (`IsWindowVisible=True`,
`IsIconic=False`), just at a coordinate region you can't see.

**Diagnostic** (PowerShell):
```powershell
# Step 1: confirm the window exists with title "Gmsh - <file>".
Get-Process python | Select-Object Id, MainWindowTitle, SessionId

# Step 2: P/Invoke GetWindowRect to inspect coordinates. Negative x = off-screen.
Add-Type -Namespace W -Name U -MemberDefinition @'
[System.Runtime.InteropServices.DllImport("user32.dll")]
public static extern bool GetWindowRect(System.IntPtr h, out RECT r);
public struct RECT { public int L; public int T; public int R; public int B; }
'@
$h = (Get-Process python | Where-Object MainWindowTitle -Match Gmsh).MainWindowHandle
$r = New-Object W.U+RECT; [W.U]::GetWindowRect($h, [ref]$r) | Out-Null
"$($r.L),$($r.T) - $($r.R),$($r.B)"
```

**Reactive rescue** (move the existing window onto the primary monitor):
```powershell
Add-Type -Namespace W -Name U -MemberDefinition @'
[System.Runtime.InteropServices.DllImport("user32.dll")]
public static extern bool MoveWindow(System.IntPtr h, int X, int Y, int W, int H, bool R);
'@
$h = (Get-Process python | Where-Object MainWindowTitle -Match Gmsh).MainWindowHandle
[W.U]::MoveWindow($h, 100, 100, 1280, 800, $true)
```

**Preventive (bake into every viewer script)**: force window geometry
BEFORE `gmsh.fltk.run()` so the stale saved position can't take effect.
```python
import gmsh
gmsh.initialize(["-noconfig"])  # skip %APPDATA%/gmsh-options stale geom
gmsh.option.setNumber("General.GraphicsPositionX", 100)
gmsh.option.setNumber("General.GraphicsPositionY", 100)
gmsh.option.setNumber("General.GraphicsWidth", 1280)   # NOT GraphicsSizeX
gmsh.option.setNumber("General.GraphicsHeight", 800)   # NOT GraphicsSizeY
gmsh.option.setNumber("General.MenuPositionX", 100)
gmsh.option.setNumber("General.MenuPositionY", 100)
gmsh.merge(msh_path)
gmsh.fltk.run()
gmsh.finalize()
```

The option names are Width/Height in GMSH 4.x. SizeX/SizeY were never
valid and raise "Could not set option" immediately; do NOT copy-paste
that variant from older gmsh tutorials.

**Wrong hypotheses to skip** when triaging "GMSH invisible":
1. "Claude Code background process lacks GUI access" -- WRONG. The python
   process IS on Console session 1, MainWindowTitle is set normally.
2. "Need `gmsh.initialize(sys.argv, run=True)` instead of `gmsh.fltk.run()`"
   -- that distinction matters ONLY in Cubit's embedded Python, NOT here.
3. "Use Start-Process to detach the launch" -- detaching does not help
   when the window restores to off-screen coordinates regardless of
   process owner.

**ALWAYS check window geometry FIRST** before re-launching, re-coding the
viewer, or hypothesising about launch paths. Reference incident
2026-05-30 (third recurrence): user observed Gx fingerprint coil viewer
producing python.exe with empty stdout; GetWindowRect returned
(-2560, 497) - (-1166, 1957), MoveWindow rescued it. Documented in
`memory/feedback_gmsh_gui_invisible_from_background.md`.
"""


GMSH_ANIMATION = """
# GMSH Animation: Displacement View for Moving Bodies

## Overview

GMSH can animate mesh displacement using $NodeData with VectorType=5
(displacement mode). Each time step specifies a displacement vector per node.
Combined with STEP geometry (static), this creates stator + moving body animations.

Runnable artifact-inspection notebook:
`docs/gmsh_animation/gmsh_animation.ipynb`.  Its synchronized JSON sidecars are
`docs/gmsh_animation/gmsh_animation_results.json` (domain values) and
`docs/gmsh_animation/gmsh_animation_result.json` (notebook-output sync).  The
docs-local artifact inspected there is `docs/gmsh_animation/`: MSH v4.1,
2430 nodes, 1003 elements, 21 vector NodeData frames, final displacement 0.15 m.

Runnable export notebook:
`docs/gmsh_animation/gmsh_animation_export.ipynb`.  It opens the docs-local
`animation.geo`, relies on `animation.geo.opt`, synchronizes all visible view
time steps, exports PNG frames, writes GIF/MP4 movies, and records which
companions (`.geo`, `.geo.opt`, `.msh.opt`) were used.  This is the small
teaching example for students and the MCP regression reference for animation
export behavior.

## .msh File Structure (v4.1 with time-stepped displacement)

```
$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
1
3 1 "mover"
$EndPhysicalNames
$Entities
...
$EndEntities
$Nodes
...  (TET10 nodes: vertices + mid-edge)
$EndNodes
$Elements
...  (type 11 = TET10)
$EndElements
$NodeData
1
"Displacement"
1
0.000000         <- time value (step 0)
3
0                <- step index
3                <- 3 components (vector)
2430             <- number of nodes
1 -1.0e-01 0.0 0.0   <- node 1: displacement (x, y, z)
2 -1.0e-01 0.0 0.0
...
$EndNodeData
$NodeData
1
"Displacement"
1
0.025000         <- time value (step 1)
3
1
3
2430
1 -9.5e-02 0.0 0.0
...
$EndNodeData
...  (repeat for each time step)
```

## Key View Options for Displacement Animation

From GMSH source analysis (PViewVertexArrays.cpp):
- VectorType=5 calls addScalarElement -> va_triangles (GL_FILL)
- IntervalsType MUST be 2 (Continuous) or 3 (Discrete) for solid fill
  IntervalsType=1 (Iso) draws iso-surface slices, NOT element faces
- TET10 (10 nodes) > PVIEW_NMAX (8) -> silently SKIPPED unless
  AdaptVisualizationGrid=1

```python
# Via gmsh Python API:

# Displacement mode
gmsh.option.setNumber("View[0].VectorType", 5)
gmsh.option.setNumber("View[0].DisplacementFactor", 1.0)

# CRITICAL: solid fill (NOT Iso which is the default)
gmsh.option.setNumber("View[0].IntervalsType", 2)     # 2=Continuous

# Exterior faces only
gmsh.option.setNumber("View[0].DrawSkinOnly", 1)
gmsh.option.setNumber("View[0].ShowElement", 0)        # 0=no wireframe

# CRITICAL for TET10/HEX20: adaptive subdivision
gmsh.option.setNumber("View[0].AdaptVisualizationGrid", 1)
gmsh.option.setNumber("View[0].MaxRecursionLevel", 1)  # 1-2

# Lighting
gmsh.option.setNumber("View[0].Light", 1)
gmsh.option.setNumber("View[0].LightTwoSide", 1)
gmsh.option.setNumber("View[0].SmoothNormals", 1)

# Ensure vectors/tets are drawn
gmsh.option.setNumber("View[0].DrawVectors", 1)
gmsh.option.setNumber("View[0].DrawTetrahedra", 1)

# Auto range (not custom)
gmsh.option.setNumber("View[0].RangeType", 1)
gmsh.option.setNumber("View[0].ShowScale", 0)
```

## Hiding Static Mesh (Original Position)

When using displacement View, GMSH also shows the mesh at its original
position. To hide it and show only the displaced version:

```
// Hide all mesh display
Mesh.SurfaceFaces = 0;
Mesh.VolumeEdges = 0;
Mesh.VolumeFaces = 0;
Mesh.SurfaceEdges = 0;
Mesh.Points = 0;
```

**Note**: `Mesh.*` options hide ALL meshes globally. Geometry (STEP)
surfaces are controlled separately via `Geometry.Surfaces`.

## STEP Geometry (Static Background)

```
// Show STEP as solid surface
Geometry.Surfaces = 1;
Geometry.SurfaceType = 2;        // filled surface (not wireframe)
```

## High-Order Elements in Displacement View

`Mesh.NumSubEdges` only affects MESH display, NOT View display.
For curved elements in displacement Views, use:

```
View[0].AdaptVisualizationGrid = 1;  // enable adaptive subdivision
View[0].MaxRecursionLevel = 2;       // subdivision depth (2-3 is good)
View[0].TargetError = 0.0001;        // adaptive error threshold
```

Without these, TET10/HEX20 elements render as straight-edged in the View.

## Solid (Filled) Display of Displaced Elements

VectorType=5 renders displaced elements as SOLID filled surfaces
colored by displacement magnitude. Key options:

```
View[0].VectorType = 5;              // displacement mode
View[0].ShowElement = 0;             // 0=no wireframe overlay, 1=with edges
View[0].DrawSkinOnly = 1;            // only external faces (faster)
View[0].IntervalsType = 2;           // 2=continuous colormap (filled)
View[0].Light = 1;                   // enable lighting for 3D shading
View[0].SmoothNormals = 1;           // smooth normals across elements
View[0].Boundary = 0;               // 0=volume, 1=surface only
View[0].Explode = 1.0;              // 1.0=no shrinking
```

Note: IntervalsType must be 2 (continuous) or 3 (discrete) for filled
rendering. Type 1 (iso-values) shows isolines, type 4 (numeric) shows
values as text.

## External View for Coloring Displaced Mesh

```
View[0].ExternalView = 1;  // color displaced mesh by View[1] scalar data
```
Enables showing stress/temperature on deformed shape. View[0] provides
displacement, View[1] provides scalar field for coloring.

## Transparency

GMSH has limited transparency support:

```
View[0].ColormapAlpha = 0.3;       // global alpha (0=transparent, 1=opaque)
View[0].ColormapAlphaPower = 0.0;  // nonlinear alpha mapping
View[0].FakeTransparency = 0;      // 0=real (sorted), 1=additive (faster)
```

- **Real transparency**: back-to-front triangle sorting per frame. Correct
  but slow, and only works within a single View.
- **Fake transparency**: additive blending. Fast but incorrect for
  overlapping geometry.
- **No per-element/per-material alpha**. Alpha is global per View.
- **STEP geometry surfaces have NO alpha support** at all.

For proper transparency, use ParaView (VTK) or Blender.

## .opt Companion File

GMSH auto-loads `<filename>.opt` when opening `<filename>.msh`.
Place display settings in the .opt file for reproducible visualization:

```
// file.msh.opt (auto-loaded with file.msh)
Mesh.SurfaceFaces = 1;
Mesh.VolumeEdges = 0;
Mesh.SurfaceEdges = 0;
Mesh.NumSubEdges = 4;
```

## Python API Example

```python
import gmsh
gmsh.initialize()

# Static geometry
gmsh.merge("stator.step")
gmsh.option.setNumber("Geometry.Surfaces", 1)
gmsh.option.setNumber("Geometry.SurfaceType", 2)

# Animated mesh with displacement NodeData
gmsh.merge("mover.msh")
gmsh.option.setNumber("Mesh.SurfaceFaces", 0)  # hide static mesh

# Displacement view
gmsh.option.setNumber("View[0].VectorType", 5)
gmsh.option.setNumber("View[0].DisplacementFactor", 1.0)
gmsh.option.setNumber("View[0].ShowElement", 1)
gmsh.option.setNumber("View[0].DrawSkinOnly", 1)
gmsh.option.setNumber("View[0].ShowScale", 0)
gmsh.option.setNumber("Mesh.NumSubEdges", 4)

gmsh.fltk.run()
gmsh.finalize()
```

## Programmatic PNG/GIF Export

Gmsh can export post-processing animation frames reliably through the Python
API when an FLTK/OpenGL context is initialized. A plain command-line
`Print 'frame.png'` can parse correctly but still fail with
`requires a graphical interface context`.

Use this pattern for `.geo` artifacts that merge a time-stepped `.msh` and
auto-load `case.geo.opt`:

```python
from pathlib import Path
import gmsh
from PIL import Image

geo = Path("case.geo")
out_gif = Path("case_zup.gif")
frames = []

gmsh.initialize(["-noconfig"])
gmsh.option.setNumber("General.GraphicsWidth", 800)
gmsh.option.setNumber("General.GraphicsHeight", 800)
gmsh.open(str(geo))             # auto-loads case.geo.opt
gmsh.fltk.initialize()          # required graphics context

for step in range(num_steps):
    gmsh.option.setNumber("View[0].TimeStep", step)
    gmsh.option.setNumber("View[1].TimeStep", step)
    gmsh.fltk.update()
    frame = Path("C:/temp") / f"frame_{step:03d}.png"
    gmsh.write(str(frame))
    frames.append(frame)

gmsh.fltk.finalize()
gmsh.finalize()

images = [Image.open(p).convert("P", palette=Image.Palette.ADAPTIVE)
          for p in frames]
images[0].save(out_gif, save_all=True, append_images=images[1:],
               duration=40, loop=0, disposal=2)
```

For multi-view animations, set `PostProcessing.Link = 1`, set each
`View[i].TimeStep` explicitly, and keep `PostProcessing.AnimationCycle = 0`.
`AnimationCycle = 1` cycles visible views instead of synchronizing them, which
causes flicker when a pressure view and a displacement view are meant to be
shown simultaneously.

On Windows FLTK builds, the exported PNG width can be smaller than
`General.GraphicsWidth` because the GUI sidebar consumes part of the window.
If an exact 800x800 exported frame is required and `GraphicsWidth=800` produces
600x800, set `General.GraphicsWidth=1000` and `General.GraphicsHeight=800`.

## Coordinate Transformation Pipeline (applied in order)

1. **Explode** (`View.Explode=1.0`) -- shrink toward barycenter
2. **Transform** (`View.TransformXX..ZZ`) -- 3x3 matrix
3. **Offset** (`View.OffsetX/Y/Z=0`) -- constant translation
4. **Raise** (`View.RaiseX/Y/Z=0`) -- elevation proportional to scalar
5. **NormalRaise** (`View.NormalRaise=0`) -- elevation along element normal
6. **Displacement** (VectorType=5) -- `xyz += DisplacementFactor * val`
7. **GeneralizedRaise** (`View.UseGeneralizedRaise=0`) -- formula-based

## GeneralizedRaise (Custom Coordinate Formulas)

```
View[0].UseGeneralizedRaise = 1;
View[0].GeneralizedRaiseFactor = 1.0;
View[0].GeneralizedRaiseX = "v0";   // variables: x,y,z,v0..v8,s,t
View[0].GeneralizedRaiseY = "v1";
View[0].GeneralizedRaiseZ = "v2";
```

Enables arbitrary coordinate transformations using field values.
14 input variables: xyz[0..2], val[3..11], step[12], time[13].

## Rotation + Translation Animation via Displacement

For combined rotation + translation, compute per-node displacement:
```python
# CW rotation around Z through (cx, cy) + X translation
theta = -2 * pi * n_rotations * t
cos_t, sin_t = cos(theta), sin(theta)
rx, ry = x - cx, y - cy
dx = (cos_t*rx - sin_t*ry + cx + x_disp) - x
dy = (sin_t*rx + cos_t*ry + cy) - y
# Write as NodeData displacement vector per time step
```

**NGSolve TET reference coordinate mapping** (important for correct curving):
```
ref(0,0,0) -> el.vertices[3]   (NOT el.vertices[0])
ref(1,0,0) -> el.vertices[0]
ref(0,1,0) -> el.vertices[1]
ref(0,0,1) -> el.vertices[2]
```
"""


def get_gmsh_documentation(topic: str = "all") -> str:
    """Return GMSH usage documentation by topic.

    Args:
        topic: One of: all, policy, overview, cli, shortcuts, options, opt_file,
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
        "opt_file": GMSH_OPT_FILE,
        "msh_format": GMSH_MSH_FORMAT,
        "geo": GMSH_GEO_SCRIPTING,
        "high_order": GMSH_HIGH_ORDER,
        "workflow": GMSH_RADIA_WORKFLOW,
        "onelab": GMSH_ONELAB,
        "pitfalls": GMSH_PITFALLS,
        "animation": GMSH_ANIMATION,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        available = ", ".join(topics.keys())
        return f"Unknown topic: '{topic}'. Available: all, {available}"
