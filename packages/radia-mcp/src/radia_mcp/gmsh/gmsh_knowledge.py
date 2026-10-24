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
- Headless inspection/validation/rendering through the radia-mcp tools
  (`gmsh_inspect_msh`, `gmsh_validate_msh`, `gmsh_validate_geo`,
  `gmsh_render`, `gmsh_export_animation`,
  `gmsh_write_post_launch_artifact`) -- these run the gmsh API only in
  a subprocess and only on existing files, never for mesh generation

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
- Default: .msh v4.2 (GMSH visualization)
- `GmshPostExport.write()` and `vol2msh()` emit .msh v4.2
- Post-processing launch artifact: `case.geo` that merges the .msh/.step data
- Exact option sidecars: `case.geo.opt` for normal launch and `case.msh.opt`
  for raw mesh/data inspection. A plain `case.opt` is not auto-loaded.
- `cubit-mesh-export` `export gmsh "case.msh"` should attach `case.geo`,
  `case.geo.opt`, and `case.msh.opt`; open `case.geo` for normal review.
- For acoustic/FEM-BEM post artifacts, use the shared
  `gmsh_post_display_contract` / `write_gmsh_post_launch_artifact` pattern:
  one MSH v4.2 data file, a `.geo` launch file, exact `.geo.opt` and
  `.msh.opt` sidecars, named views, Z-up camera metadata, and optional
  cut-plane metadata (`General.Clip0A/B/C/D`).  This is the cross-language
  contract shared with the MATLAB/Gypsilab readable acoustic lane.

## Dumped Figures: Axis Equal for Spatial Plots (POLICY)

Every figure dumped from the gmsh post lane whose BOTH axes are
spatial lengths -- contour / flux-line plots, streamline plots, cut
sections, mesh pictures -- MUST be axis equal (1:1:1).  A stretched
axis misrepresents exactly what these figures exist to show: field-
line density, the orthogonality of field lines and equipotentials,
and the geometry itself.

- gmsh renders are axis equal by default: `General.ScaleX`,
  `General.ScaleY`, `General.ScaleZ` all stay at 1.0.  Do NOT move
  them for contour/streamline/section figures.  `gmsh_render` /
  `gmsh_export_animation` return a warning note whenever a Scale
  option is overridden away from 1.0.
- Deliberate exaggeration (e.g. Z-scaling a thin plate or a warp
  display) is allowed ONLY when the caption states the scale factor
  -- it is a labeled exception, not a default.
- matplotlib figures: any spatial x-y plot (streamline coordinates,
  section outlines, exported contour data) must set
  `ax.set_aspect("equal")`.  Value-vs-parameter plots (line/curve
  profiles, histograms, point histories) are NOT spatial figures and
  keep the auto aspect -- axis equal would be wrong there.
"""

GMSH_OVERVIEW = """
# Gmsh 5 post-processing runtime

Radia uses qualified official Gmsh 5 with ASCII MSH 4.2. Solver meshes are
checked Cubit/Netgen `.vol`; Gmsh is post-processing only.

Installation and wheel/provenance commands have one maintained source:
[the Gmsh runtime guide](https://github.com/ksugahar/Radia/blob/main/packages/radia-mcp/docs/operations/gmsh-development-runtime.md).
Dependency requirements are declared only in the package's pyproject extras.
Do not copy dated pins or a second install recipe into scripts or knowledge.

Native tools require `RADIA_GMSH_PROVENANCE`, bind Python/native payload
hashes to the adopted official wheel, and probe actual MSH 4.2 field I/O.
Every run records the exact distribution/native version and wheel SHA256.
A new development or released build must pass writer/field roundtrips and
relevant post-processing semantics checks before adoption. Tools never run
pip or substitute a different runtime. Pure-Python readers stay exact-4.2
without requiring native Gmsh.

Gmsh is GPL software by Christophe Geuzaine and Jean-Francois Remacle.
`gmsh_detect()` reports available capabilities. Inspection and CSV/statistics
can run without native Gmsh; rendering and native plugins require it. See
`gmsh_usage("workflow")` for the complete result-to-picture workflow.
"""

GMSH_COMMAND_LINE = """
# Opening existing post-processing results

## Windows PowerShell `-string` quoting

Use single quotes inside a Gmsh snippet passed through PowerShell:

```powershell
gmsh display.geo -string "Print 'frame.png'; Exit;"
```

Backslash-escaped double quotes passed through a batch wrapper can leave
an unquoted drive token and produce `Unknown variable 'W'`. Image export
also requires an initialized FLTK/OpenGL context. Prefer `gmsh_render`
or `gmsh_export_animation` for unattended image output.


```sh
gmsh case.geo
gmsh result.msh -numsubedges 4
```

Open the generated `.geo` launch artifact to retain its exact `.geo.opt`
settings. Opening raw `.msh` uses `.msh.opt`. A plain `case.opt` is not an
autoload companion. Use `gmsh_render` and `gmsh_export_animation` for automated
output and recorded settings. Do not invoke Gmsh to create a solver mesh.
Solver meshes originate in Cubit or Netgen and cross the solver boundary as
checked `.vol` plus the versioned label contract.
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
| `Mesh.MshFileVersion` | 4.2 | MSH format version |
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
| `View.ArrowSizeMax` | 20 | Max arrow size in pixels (`ArrowScale` was the pre-4.x name) |
| `View.ArrowSizeMin` | 20 | Min arrow size in pixels (set = Max for fixed-size arrows) |

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
case.msh          # single Gmsh MSH v4.2 result/data container
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
# MSH 4.2 result contract

MSH 4.2 is the sole supported mesh/data container. The supported pure-Python
inspection boundary is ASCII, with `$MeshFormat` declaring `4.2 0 8`.
Unsupported versions, binary files and malformed sections fail explicitly.
Use `gmsh_inspect_msh` and `gmsh_validate_msh` before native processing.

Native model-based views use `NodeData`, `ElementData`, or `ElementNodeData`.
Tags are size_t. Data storage is blocked: after string/real/integer tags, the
data header is `numBlocks numEntries minTag maxTag`; a node/element block
stores its entry count, all tags, then component values. Element-node blocks
start with `elementType numValuesPerElement numElementsInBlock`. Embedded
entities belong in `$Entities`/`$PartitionedEntities`; `$GhostElements` is also
blocked. Changing a file's version header does not convert its data layout.

The pure-Python reader currently accepts ordinary unpartitioned meshes and
embedded entities. PartitionedEntities, GhostElements, Polytopes, Edges and
Faces sections fail explicitly; describing their format does not imply reader
support. There is no partial parse or fallback.

Use the existing result exporters or Gmsh 5 view APIs; do not hand-write a
second serializer. Validate every producer's output by opening it in Gmsh 5,
reading back prescribed values, writing the view, and comparing it again.
Record view name, scalar/vector/tensor component count, time index and time in
seconds, tags, SI units, and exact runtime/wheel identity.

Before every native field hand-off, node and element field tags greater than
2147483647 are rejected explicitly. There is no remapping or fallback.
Pure-Python inspection and mesh-only large tags retain their own support;
that does not imply native field support at those tags.

Solver inputs are checked Cubit/Netgen `.vol` meshes. MSH is exclusively
derived visualization/post-processing, including conductor display objects.
Display tessellation in `raster.py` is transient visualization data, never
written as a solver mesh or fed to a solver.
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
5. **MCP**: `gmsh_render` defaults to `numsubedges=4` (verified on a
   Curve(3) Tri10 sphere: numsubedges=1 renders a faceted polyhedron,
   4 renders the smooth sphere)

## Verify High-Order Meshes Before Trusting the Picture

Display looks fine even when the element node ordering is wrong -- the
GUI never evaluates Jacobians. Gate every high-order export with
`gmsh_validate_msh(path, check_jacobians=True)`: zero negative
determinants and the integrated per-type volume matching CAD are the
repo acceptance criteria (see "GMSH API Node Ordering Verification
Policy"). NumSubEdges only affects MESH rendering; post views need
`View[i].AdaptVisualizationGrid = 1` (>8-node elements are silently
skipped otherwise), which `gmsh_render` / `gmsh_export_animation` set
by default.

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
post.write("results.msh")                   # v4.2 lab standard
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
# Gmsh 5 / MSH 4.2 post-processing workflow

Start with an existing result file exported by Radia or cubit-mesh-export.
Keep the original checked solver `.vol`, its label contract and result
provenance. Gmsh does not create or import the solver mesh in this workflow.
The calls below are MCP tools, using absolute paths for files.

## 1. Inspect the result before native access

```python
gmsh_detect()
gmsh_inspect_msh(msh_path="/results/field.msh")
gmsh_validate_msh(msh_path="/results/field.msh")
gmsh_field_stats(msh_path="/results/field.msh")
```

Check the returned errors, view names, component counts, step times, bounding
box and finite ranges. Use the actual view name in subsequent calls. A vector
view is not three unrelated scalar views. Mesh-only validation does not
certify a physical solution; Jacobian checks can be requested separately.

## 2. Create the normal launch artifact

```python
gmsh_write_post_launch_artifact(
    msh_path="/results/field.msh", output_base="/results/field",
    title="Magnetic field", camera_preset="z_up_xz_from_positive_y")
```

This writes `field.geo`, `field.geo.opt`, `field.msh.opt` and
`field.display.json`. Open `field.geo` in Gmsh for interactive viewing.
The existing MSH data are not modified. STEP overlays are display geometry;
their coordinates and units must match the result. They are not solver input.

## 3. Render a view or a cut

```python
gmsh_render(path="/results/field.geo", png_out="/results/B.png",
            view="B", time_step=0, camera_preset="iso")
gmsh_render(path="/results/field.msh", png_out="/results/B_cut.png",
            view="B", time_step=0,
            cut_plane={"enabled": True, "normal": [0, 1, 0], "offset": 0})
```

Use a shared explicit color range when comparing frames/cases. Preserve equal
spatial axes and state any intentional displacement exaggeration. High-order
edge subdivision and adaptive view display are enabled by the renderer;
display refinement does not change the solver discretization. Native Gmsh
graphics is the normal interactive path. Performance must be measured on the
adopted build; version choice alone is not a speedup measurement.

## 4. Extract values using the same view and time convention

```python
gmsh_probe(msh_path="/results/field.msh", view_name="B",
           points=[[0, 0, 0]], step=0)
gmsh_line_profile(msh_path="/results/field.msh", view_name="B",
                  start=[0, 0, 0], end=[0, 0, 0.1], n=101)
gmsh_export_csv(msh_path="/results/field.msh", view_name="B",
                csv_out="/results/B.csv")
```

An outside probe reports `found=false`; it is not a zero field. Integrals,
flux/force operations and line samples need their declared domain, normal,
units and interpolation semantics. They do not replace numerical validation.
Use `gmsh_reference("python_api_postproc")` for native view/plugin details
and `gmsh_usage("paraview")` for filter correspondences and limitations.

## 5. Export time-dependent results

```python
gmsh_export_animation(path="/results/field.geo", view="B",
                      gif_out="/results/B.gif", delay_ms=80,
                      link_views=True)
```

Verify real step times in the MSH data before choosing playback delay.
Synchronize relevant views; displacement views need the displacement display
mode and documented scale. Keep camera and color limits fixed for comparison.
See `gmsh_usage("animation")` and `gmsh_usage("magnetic_movies")`.

## 6. Save evidence and respect output types

Retain input/output hashes, exact Gmsh wheel/version, tool arguments, view/time
selection, companion files, returned warnings and actual image dimensions.
Model-based view exports must remain MSH 4.2 and be re-inspected after write.
List-based plugin views use their native `.pos` output; a `.msh` target is
rejected rather than converted implicitly. There is no solver import route.
"""

GMSH_ONELAB = """
# ONELAB and post-processing

ONELAB parameterized solver/mesh-generation models are not a Radia solver
input route. Use the result display companions and MCP post-processing tools
on existing MSH 4.2 results. Physical parameters and solver execution remain
with the owning Radia application, using checked Cubit/Netgen `.vol` inputs.
"""

GMSH_PITFALLS = """
# GMSH Common Pitfalls

## 1. Curved Elements Look Straight
**Problem**: High-order elements display as straight-edged polygons.
**Fix**: Set `Mesh.NumSubEdges = 4` (default is 2, too coarse). The MCP
`gmsh_render` tool applies this by default. If the mesh STILL looks
faceted at NumSubEdges=4, the mid-edge nodes probably lie on straight
lines (a "linear high-order" export) -- confirm with `gmsh_inspect_msh`
(element order) and re-export with real curving (`mesh.Curve(p)` /
Cubit `export ... order N`).

## 1b. High-Order Mesh Displays Fine but Is Actually Broken
**Problem**: The GUI shows a plausible mesh, yet solvers or volume
integrals misbehave; nobody notices for months.
**Cause**: Rendering never evaluates Jacobians, so systematically
inverted node ordering (negative det everywhere) is invisible.
**Fix**: Gate exports with `gmsh_validate_msh(path, check_jacobians=True)`
-- 0 negative determinants + per-type volume vs CAD. Real case 2026-08:
all four docs/gmsh_animation TET10 meshes were fully inverted while
animating perfectly in the GUI.

## 1c. gmsh.open() Silently Changes Mesh Display Options
**Problem**: A script sets display options, opens a view-bearing .msh,
and the picture shows filled surfaces that were never requested (or a
"hide everything" recipe still draws the mesh).
**Cause**: `gmsh.open()` on a post-processing file adjusts mesh display
options itself -- observed on gmsh 4.15.2: opening a .msh with NodeData
flips `Mesh.SurfaceFaces` to 1.
**Fix**: Set display options AFTER `gmsh.open()`, and override
explicitly (`Mesh.SurfaceFaces = 0` etc.) when you need full control.
The MCP `gmsh_render` applies its options after open; its `blank_check`
result exposes the content fraction so an unexpectedly full/empty frame
is visible programmatically.

## 2. Missing Elements in Output
**Problem**: .msh file has no elements or missing elements.
**Fix**: Either define Physical Groups, or use `Mesh.SaveAll = 1`.
By default, only elements in Physical Groups are saved.

## 3. Stale Option Names from Older GMSH Versions
**Problem**: Recipes copied from old tutorials set options that no
longer exist (e.g. `General.ConfirmQuit`, `View.ArrowScale`) -- gmsh
warns "unknown option" or a .geo merge errors out.
**Fix**: Verify names against the installed gmsh with the MCP
`gmsh_probe_options` tool (or `gmsh_validate_geo(check_options=True)`
for a whole .geo). Known renames: `View.ArrowScale` ->
`View.ArrowSizeMin`/`ArrowSizeMax`; `General.ConfirmQuit` was removed.

## 4. Console Encoding on Windows (cp932)
**Problem**: Unicode characters cause errors in Japanese Windows.
**Fix**: Use ASCII only in GMSH scripts and output.


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
(-2560, 497) - (-1166, 1957), MoveWindow rescued it.

Note: the MCP `gmsh_render` / `gmsh_export_animation` tools bake this
prevention in (subprocess with `-noconfig` + explicit
GraphicsPositionX/Y + Width/Height), so headless screenshots never
inherit a stale off-screen geometry.
"""


GMSH_ANIMATION = """
# GMSH Animation: Displacement View for Moving Bodies

## Overview

GMSH can animate mesh displacement using $NodeData with VectorType=5
(displacement mode). Each time step specifies a displacement vector per node.
Combined with STEP geometry (static), this creates stator + moving body animations.

Use newly validated MSH 4.2 animation artifacts. Saved historical notebook
output does not qualify a current runtime or regenerate its input data.

Runnable export notebook:
`docs/gmsh_animation/gmsh_animation_export.ipynb`.  It opens the docs-local
`animation.geo`, relies on `animation.geo.opt`, synchronizes all visible view
time steps, exports PNG frames, writes GIF/MP4 movies, and records which
companions (`.geo`, `.geo.opt`, `.msh.opt`) were used.  This is the small
teaching example for students and the MCP regression reference for animation
export behavior.

## .msh File Structure (v4.2 with time-stepped displacement)

```
$MeshFormat
4.2 0 8
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
1 2430 1 2430    <- blocks, count, min/max node tag
2430             <- nodes in this block
1
2
...              <- all 2430 node tags first
-1.0e-01 0.0 0.0 <- node 1 displacement
-1.0e-01 0.0 0.0 <- node 2 displacement
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
1 2430 1 2430
2430
1
...              <- all node tags first
-9.5e-02 0.0 0.0
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
```

- **Real transparency**: back-to-front triangle sorting per frame. Correct
  but slow, and only works within a single View.
- **Fake transparency**: additive blending. Fast but incorrect for
  overlapping geometry.
- **No per-element/per-material alpha**. Alpha is global per View.
- **STEP geometry surfaces have NO alpha support** at all.

For see-through rendering, isolate the region with per-view visibility or a
cut plane, or export a surface (STL/glTF) to Blender.

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

The MCP tool `gmsh_export_animation` packages this whole recipe (linked
views, AnimationCycle=0, per-view TimeStep stepping, PNG frames, Pillow
GIF assembly, subprocess isolation); `gmsh_render` does the single-frame
PNG case. Reach for the raw pattern below only when a notebook or script
needs custom per-frame logic.

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


GMSH_PARAVIEW_PARITY = """
# ParaView -> gmsh correspondence (radia-mcp post verbs)

Everything the lab used ParaView for maps onto the radia-mcp gmsh
tools below.  Semantics were MEASURED on gmsh 4.15.2 and are locked by
tests (test_gmsh_paraview_parity.py, test_gmsh_post_process.py).

## Filter correspondence matrix

| ParaView filter | radia-mcp tool | Notes (measured semantics) |
|---|---|---|
| Calculator | gmsh_math_eval | NODAL application: result interpolates f(node values) |
| Gradient / Curl / Divergence | gmsh_derived_field | exact derivative of the P1 interpolant, element-wise constant |
| Eigenvalues (tensor) | gmsh_derived_field(operation="eigenvalues") | writes Min/Mid/Max as 3 views in one file |
| Contour (isosurface) | gmsh_isosurface | recur_level > 0 = adaptive extraction on order-2 data (0.21 -> 0.008 measured) |
| Slice | gmsh_cut_plane_extract | plane A x + B y + C z + D = 0 |
| Clip (visual) | gmsh_render(cut_plane=...) | display-side clipping, data untouched |
| Threshold | gmsh_threshold | selects on the ELEMENT MEAN of the scalar |
| Extract Surface | gmsh_extract_skin | boundary skin with the field interpolated |
| Reflect | gmsh_mirror_expand | parity-aware: B/H are PSEUDOVECTORS (v' = det(M) M v) |
| Transform | gmsh_transform_view | coordinates only -- data rewrite via value_expressions |
| Warp By Vector | gmsh_warp | Plugin(Warp) MOVES THE MODEL NODES (in place) |
| Cell Data to Point Data | gmsh_smooth_to_nodes | node = mean of adjacent elements (10/20 -> 15) |
| Stream Tracer | gmsh_streamlines | adaptive RK4 + CLOSED-LOOP detection + termination reasons (Plugin StreamLines is broken on this build) |
| Evenly Spaced Streamlines 2D | gmsh_streamlines_2d | Jobard-Lefer on ANY plane slice of 3D data (ParaView: native-2D datasets only) |
| Particle Tracer / ParticlePath | gmsh_particle_trace | BEYOND parity: ParaView advects massless tracers along a velocity field; this integrates the relativistic LORENTZ DYNAMICS dp/dt = q(E + v x B) (Boris pusher), which ParaView has no native filter for |
| (no ParaView filter) | gmsh_poincare | BEYOND parity: plane-piercing section of the orbits, interpolated between samples |
| Integrate Variables (surface) | gmsh_flux_integral | ParaView integrates a dataset over ITS OWN cells; this integrates n.B over a patch YOU specify (rect/disc) |
| (no ParaView filter) | gmsh_line_integral | BEYOND parity: circulation with the analytic tangent in parameter space + an Ampere-law NI check |
| (no ParaView filter) | gmsh_maxwell_force | BEYOND parity: Maxwell stress force/torque on a box -- ParaView has no EM stress tensor |
| (no ParaView filter) | gmsh_gap_harmonics | BEYOND parity: space-harmonic FFT around the machine air gap |
| (no ParaView filter) | gmsh_compare_fields | BEYOND parity: CROSS-MESH solution comparison (ParaView's Resample With Dataset is the nearest, without the norms) |
| (FEMM-style 2D flux plot) | gmsh_flux_lines | contours of A_z / psi = EXACT equal-flux field lines, no integration |
| Glyph (vector arrows) | gmsh_render options | View.VectorType=4 (3D arrow), View.GlyphLocation, View.ArrowSizeMax/ArrowSizeMin |
| Plot Over Line | gmsh_line_profile | straight line + PNG graph |
| Plot over custom curve | gmsh_curve_profile | parametric x(u),y(u),z(u); air-gap B(theta) in one call |
| Plot Data Over Time | gmsh_point_history | per-step values + recorded step times |
| Resample To Image | gmsh_resample_grid | CutBox regular grid, W varies fastest |
| Spreadsheet / Save Data | gmsh_export_csv | pure-Python on the MSH parser; nodes or element centroids |
| Histogram | gmsh_field_histogram | value / |field| distribution + PNG |
| Find Data (max/min location) | gmsh_view_min_max | argmin/argmax coordinates + values |
| Integrate Variables | gmsh_integrate | P1 accuracy on high-order views; exact FE integrals live in NGSolve |
| Temporal shift/harmonics | gmsh_harmonic_to_time / gmsh_modulus_phase | AC phasor -> rotating field / amplitude+phase |
| Comparative views | gmsh_render_montage | labeled PNG grid |
| Camera orbit / fly-around | gmsh_export_animation(orbit_axis=...) | camera sweeps, data fixed at time_step |
| CAD + data overlay (Merge) | gmsh_render(merge_files=[...]) | coil STEP + filament .msh + field .msh in ONE scene; CAD merge auto-enables shaded faces |
| Warp-like displacement display | View.DisplacementFactor | vector view drawn as displaced (display only) |
| Color map controls | gmsh_render options | View.RangeType, View.CustomMin/CustomMax, View.SaturateValues, View.IntervalsType, View.NbIso |

## Field lines: pick the right tool (the ParaView pain point, solved)

Stream tracing is where ParaView disappoints (manual seeds, bunched
density, loops that overdraw or stop mid-circle).  Three tools, by
decreasing exactness -- always prefer the highest one that applies:

1. **2D / axisymmetric with a potential available ->
   gmsh_flux_lines.**  Contours of A_z (or psi = r*A_theta) ARE the
   field lines, with EQUAL FLUX between adjacent lines: exact
   geometry, exact physical density, closed curves by construction.
   No seeding, no integration error.  This is how FEMM draws motor
   flux plots; if the solve can export A_z, use this.
2. **Plane slice of a 3D field -> gmsh_streamlines_2d.**
   Jobard-Lefer evenly spaced placement: seeds spawn automatically
   d_sep from accepted lines, lines stop at d_sep/2 from neighbors.
   Uniform visual density (line spacing does NOT encode |B|; keep
   the |B| information in the line color).  EXACT field lines when
   the plane is a symmetry plane (B.n = 0); otherwise the standard
   projected-field portrait -- say so in captions.
3. **True 3D lines -> gmsh_streamlines.**  Curvature-adaptive RK4
   (step halves above max_turn_deg per step), closed loops detected
   and closed exactly, per-line termination reasons.  A line that
   reports "stagnation" INSIDE the domain marks a field zero (or a
   numerically incoherent region) -- a real diagnostic, not a
   rendering artifact.
4. **Particle ORBITS are NOT field lines -> gmsh_particle_trace.**
   A charged particle GYRATES around a field line (gyroradius
   r = gamma m v_perp / |q| B) instead of following it -- never
   depict a beam with streamlines.  The relativistic Boris pusher
   integrates dp/dt = q(E + v x B) with the time step defaulting to
   1/64 of the LOCAL gyration period at each seed (B = 0 at a seed
   demands an explicit dt_s -- fail-loud, no guessed time scale).
   Color by time to read direction of travel, by energy when an
   e_view_name accelerates the particle.  Physics self-checks built
   in: in a pure magnetic field the Boris rotation conserves speed
   EXACTLY, so the reported per-track speed_change_rel (~1e-15) is
   an integrator health gate, and the seed gyroradius lets you check
   the orbit scale against the mesh by hand.  For tracks computed on
   the radia side (rad.Fld exact open-boundary fields via
   radia.particle_tracking.track_lorentz_ivp), write them with
   radia.gmsh_post_export.export_particle_tracks_msh and render the
   .msh together with the coil STEP (merge_files) -- same views
   (time / speed / kinetic energy / v arrows), solver-grade fields.

   **Flying beam (`animation_frames=N`, both sides).**  Step k of the
   extra `beam` view carries the real quantity only where the particle
   has already been and a far out-of-range SENTINEL everywhere else.
   MEASURED mechanism: with a FIXED colour range and
   `SaturateValues = 0` gmsh DOES NOT DRAW an element whose value falls
   outside the range, so the trail grows while the colour scale stays
   put -- strictly better than sweeping the range itself, which
   rescales the colours every frame.  Two measured traps: the sentinel
   must be FINITE (gmsh writes a NaN as the literal `nan` and its own
   parser then rejects the file: "Unknown variable 'nan'"), and it must
   sit far out on the scale of the VALUES, not of their spread -- a
   monoenergetic beam coloured by energy has a spread of ~1e-13, so
   "one spread above the maximum" would carry the beam's own colour.
   Always render the beam view with the returned
   `animation["render_hint"]["color"]`; an autoscaled colour bar folds
   the sentinel back into the picture.  Within ONE file the camera does
   NOT drift between frames (the bounding box comes from the mesh
   nodes, not from the drawn subset), unlike the cross-file case that
   `render_panels` has to fix with a hidden frame.

   **A beamline is the axis-equal exception.**  A beam is millimetres
   across and metres long, so a 1:1 render of a beamline is literally a
   line (measured: a 0.85 m quadrupole track with a 12 mm aperture
   renders as one stripe, magnet CAD included).  Write the `.msh` for
   interactive 3D and PRINT an envelope plot with the transverse axis
   magnified AND the magnification stated -- see
   docs/gmsh_post/em_particle_orbits.ipynb.

## EM-native reductions: the number, not the picture (all measured)

Everything above turns a field into an IMAGE.  These five turn it into
the NUMBER an electromagnetics engineer actually reports, each locked
by a closed form:

| Verb | Golden (MEASURED on the 6-tet Kuhn cube unless noted) |
|---|---|
| `gmsh_flux_integral` | uniform B = 2 T through a 0.4 x 0.4 rect tilted 30 deg: `\\|B\\| A cos30 = 0.27712812921102`, rel 2.4e-14; disc r = 0.2 normal z: `2 pi (0.04) = 0.2513274122871834`, rel 1.1e-15 |
| `gmsh_line_integral` | H = (-y, x, 0), circle r = 0.3: `2 pi r^2 = 0.565486677646163`, rel 3.9e-16 on the parametric route |
| `gmsh_maxwell_force` | uniform B = 2 T, box half 0.4: EVERY face `B^2/(2 mu0) A = 1.018592e6 N` (rel 1.0e-14) and the total EXACTLY zero -- a real cancellation of megaNewton terms. Nonzero case: `B = (a x, 0, 0)` has `div B = a`, so `F_x = a^2 c_x V / mu0`, rel 9.1e-16 |
| `gmsh_gap_harmonics` | scalar `S = x + 2y` on r = 0.25, N = 64: `a_1 = 0.25`, `b_1 = 0.5`, leakage over every other bin 7.2e-17 |
| `gmsh_compare_fields` | two DIFFERENT meshings of one cube: a LINEAR field agrees to 5.6e-16 (P1 is exact on both), a QUADRATIC one shows real interpolation error that falls 0.245 -> 0.0622 under refinement, ratio 3.95 vs the O(h^2) 4 |

Traps these verbs are built around, each MEASURED:

1. **A partially-sampled locus is refused, not averaged.**  All of them
   return `ok: False` with `n_outside` when a sample point leaves the
   mesh.  A "flux" quoted over the part of a patch that happened to be
   inside is exactly the unauditable number this lane exists to
   prevent.
2. **Read `per_face` before believing a small force.**  In a uniform
   field each face of the box carries ~1e6 N and they cancel to
   machine zero.  A near-zero TOTAL is evidence of correct
   cancellation only if the per-face terms are large.
3. **The air-gap circle must EXCLUDE its endpoint.**  Sampling
   `theta` over a closed `[0, 2pi]` double-counts `theta = 0` and
   contaminates every coefficient; `gmsh_curve_profile` samples the
   closed interval, which is why reusing it for harmonics is wrong.
4. **THD can legitimately be `None`.**  It is referenced to order
   n = 1, and a p-pole-pair machine sampled over a full revolution has
   its fundamental at n = p.  A THD computed against a roundoff-level
   n = 1 (MEASURED: 2.86 for a pure DC field) is noise dressed as a
   result, so it is refused instead.
5. **`gmsh_diff_msh` cannot compare across meshes.**  It compares
   counts and per-view statistics: MEASURED, two meshings of one cube
   carrying a bit-identical field come back
   `fields_match=False, count mismatch=True`.  Correct for its
   same-mesh regression contract, useless for FEM-vs-HDiv-VIM --
   that is `gmsh_compare_fields`.

## Beautiful isosurfaces: the four levers (all measured)

A raw Plugin(Isosurface) render looks faceted and flat.  Four levers
turn it into the publication picture; the radia-mcp tools apply the
first two for you:

1. **Adaptive extraction on high-order data**:
   ``gmsh_isosurface(recur_level=3..4)`` / ``gmsh_flux_lines(...)``
   (accepted range 0..6).
   Plugin(Isosurface) honors RecurLevel ONLY when the source view has
   adaptive visualization data -- the tools enable
   ``View.AdaptVisualizationGrid`` + ``View.MaxRecursionLevel`` on
   the source view before running the plugin (without that the
   option is SILENTLY ignored).  Measured on a TET10 quadratic
   field: radial error 0.21 -> 0.008, 1 -> 143 triangles at level 4.
   Order-2 GmshPostExport output is exactly what this feeds on.
2. **Smooth shading**: ``View.SmoothNormals = 1`` (gmsh default is
   OFF -- the main cause of the faceted look).  gmsh_render /
   gmsh_export_animation set it on every view by default
   (``smooth_normals=False`` opts out); ``View.AngleSmoothNormals``
   (default 30 deg) keeps genuine creases sharp.
3. **Nesting with transparency works** --
   ``color={"alpha": 0.3..0.5}`` on a stack of levels
   (``View.LightTwoSide`` already defaults to 1).  **CORRECTION
   (2026-08-07, measured):** an earlier note here blamed "hairline
   T-junction cracks from per-element adaptive subdivision".  That was
   WRONG.  Counting background pixels enclosed by the surface
   silhouette on a curved p2 field gives **0 crack pixels at every
   recursion level 0..3 and at alpha 1.0 and 0.35** (235 -> 18653 iso
   elements), for single and nested shells alike.  Adaptive
   subdivision does not crack the surface.
4. **What DOES look like cracks: an OPEN shell.**  Where the level set
   crosses the outer domain boundary the surface is cut open, and a
   semi-transparent open shell is see-through by construction.
   Measured on the same field: a shell closing inside the domain gives
   0%, one crossing the box faces gives 19-29% see-through, and the
   figure is **independent of recursion level** (261 vs 5184 elements
   -> same 22%) -- which is what rules the subdivision explanation out.
   ``gmsh_isosurface`` now reports ``open_surface`` /
   ``touches_outer_boundary`` / ``boundary_vertices`` and a note, so the
   artefact is a stated fact
   rather than a mystery.  Fixes: choose a level whose shell closes
   inside the domain, enlarge the domain, or clip the view.  This is
   explicitly an outer-bounding-box contact check; Gmsh emits
   element-local isosurface polygons without shared topology, so internal
   openings are not classified by this field.
5. **Clip + opaque** stays the cleanest nesting recipe when the inner
   levels matter: ``cut_plane={"enabled": True, "normal": [0,-1,0],
   "offset": 0}``.  Opaque clipped surfaces keep full colormap
   saturation.

Direct display alternative (no extraction): rendering the volume
scalar with ``View[0].IntervalsType = 1`` + ``View[0].NbIso`` draws
isosurfaces on the fly (RangeType=2 + CustomMin/CustomMax pin the
levels) -- handy for quick looks.

## Time series across FILES (one .msh per step)

gmsh's own time steps live inside one view; a transient solver writes
one mesh per step.  ``gmsh_time_series`` treats an ordered file list as
the time axis:

- per-tag ``min``/``max``/``mean``/``std``/``rms``/``ptp`` **and
  ``argmax_time``/``argmin_time``** written as views into one .msh --
  "where is the peak, and WHEN" becomes a picture;
- per-step global aggregates (min/max/mean/rms), i.e. the plot-over-
  time series, plus optional interpolated point histories (a real gmsh
  probe per file, not a nearest-node lookup) and a matplotlib PNG.

The files must share one node/element numbering; a series whose mesh
changed is not one time series, and that is checked rather than
averaged over silently.

## STEP overlay + PEEC filament post (the Merge workflow, headless)

The lab's classic GUI recipe (Merge coil.step + filaments.msh +
field.msh) is one headless call:

    gmsh_render("field.msh",
                merge_files=["coil.step", "filaments.msh"])

- A merged CAD file (.step/.brep) auto-enables shaded solid faces
  (Geometry.Surfaces=1, SurfaceType=2, vertex dots off); override
  with geometry_display=False or Geometry.* options.
- UNITS (measured): radia/netgen ``WriteStep`` writes METER
  coordinates, so the STEP overlays radia field data 1:1 with no
  conversion.  External CAD STEP files are commonly in mm -- check
  the overlay bbox before trusting a composite figure.

PEEC filaments are first-class data via
``radia.gmsh_post_export.export_filaments_msh`` (per-filament
physical groups, element tags offset so sibling merges never drop
lines):

- ``currents=I`` writes ``|I| [A]`` per line element (color).
- ``direction_view=True`` adds ``I direction [A]``: unit tangent x
  SIGNED Re(I) arrows -- the winding's current flow direction, with
  reverse-carrying strands pointing backwards.
- ``complex_steps=True`` (AC) adds ``I_complex [A]`` with step 0 =
  Re I, step 1 = Im I: feed it to gmsh_harmonic_to_time +
  gmsh_export_animation for the rotating-phasor current animation
  (outer-filament phase lead = visible skin effect), or to
  gmsh_modulus_phase for amplitude/phase maps.

A SOLVED PEEC circuit goes straight to a picture -- no hand-assembled
polylines -- through ``export_peec_topology_msh`` (the
``build_topology()`` dict + branch currents) or its one-hop wrapper
``PEECCircuitSolver.export_gmsh(msh, freq, port_currents)``:

    solver.export_gmsh("coil_I.msh", 50e3, [1.0],
                       direction_view=True, complex_steps=True)
    gmsh_render("coil_I.msh", "coil_I.png", merge_files=["coil.step"])

It writes ``|I| [A]``, ``|J| [A/m^2]`` (= |I| / (width*height): the
current-DENSITY map, which is what actually exposes skin/proximity
when strands have unequal cross-sections) and ``P [W]`` (per-segment
ohmic loss; ``current_convention="amplitude"`` default gives
0.5|I|^2 R, ``"rms"`` gives |I|^2 R -- the choice is echoed in the
returned summary so the reported power is never ambiguous).

REVERSE CURRENTS -- ``|I|`` HIDES THEM (measured, 4 mm square bar at
100 MHz): in the inductance-limited regime interior strands run
BACKWARDS -- 1 reverse strand at nwinc=nhinc=3 (Re I = -0.091 A of
1 A total), 4 at 4x4, 8 at 5x5, the high-frequency limit reproducing
L^-1 1 exactly.  A |I| colormap paints them as ordinary positive
current.  Read the reversal from the signed ``I direction`` arrows or
from step 0 of ``I_complex`` (= Re I, with a diverging colormap);
``n_reverse_segments`` in the returned summary counts them.

CoilBuilder closes the loop end-to-end: ``write_step()`` for the CAD
solid, ``to_filaments(nw, nh, frequency=...)`` for the Tier-A current
distribution, ``to_radia()`` for the Biot-Savart field -- the three
layers of one overlay figure come from one builder.

Executable showcase: ``docs/gmsh_post/em_fieldlines.ipynb`` runs the
whole battery on two analytic-field cases (circular coil with exact
psi contours + cross-check; opposed-PM gap with mid-plane
streamlines, gap profile, nested adaptive isosurfaces), with the key
numbers locked in its synchronized result JSON.

## Figure controls: named values, no gmsh option strings

gmsh HAS every display knob ParaView has (verified by probing this
build: camera, orthographic, 6 clip planes, labelled axes, colorbar
range/format/intervals/log, glyph decimation, 2D text).  What used to be
missing was reaching them without knowing their option names, so
gmsh_render / gmsh_export_animation take them as named values:

- ``camera_preset="+x"|"-x"|"+y"|"-y"|"+z"|"-z"|"iso"`` -- the named
  axis points AT the camera, so ``"+y"`` shows the x-z plane face-on.
  MEASURED with a marker rig, not copied from another package: the
  minus views are a HORIZONTAL mirror with up preserved.  ``(180,0,0)``
  is NOT ``-z`` -- it flips vertically while still looking from +z.
  Guessing raw ``rotation`` angles is how a plane gets rendered
  edge-on as a single line.
- ``color={"range": [lo, hi] | "shared", "log":, "intervals":,
  "style": "continuous|iso|discrete|numeric", "format":, "colormap":,
  "alpha":, "show_scale":, "saturate":, "views": [i]}``.
  **The range is the publication-critical one**: gmsh autoscales EVERY
  view to its own extrema, so two panels of the same quantity are NOT
  comparable until they share a scale.  ``"shared"`` unifies the views
  of one render; for a cross-FILE comparison read the union range from
  gmsh_field_stats and pass it explicitly to each render.
- ``glyphs={"type": "arrow3d", "sampling": n, "size_max":, "center":,
  "location": "cog|vertex"}`` -- ``sampling`` draws every n-th element,
  the difference between a readable arrow field and a solid mat.
- ``clip=[{"normal": [nx,ny,nz], "offset": d, "apply_to":
  ["views","mesh","geometry"], "whole_elements": bool}]`` (max 6);
  keeps the ``n . x + d >= 0`` half-space.
- ``axes=True`` or ``{"mode": "box|frame|open|full|open_grid",
  "labels": ["x [m]", ...], "format":, "tics":}`` -- labels carry the
  units for a publication figure.
- ``annotations=["text"]`` or ``[{"text":, "x":, "y":, "align":,
  "size":}]`` in window pixels (negative counts from the far edge --
  keep clear of the colorbar strip at the bottom).

Unknown names raise with the valid list; raw ``options={}`` still wins
over the structured form, so an option gmsh gains later is reachable
without waiting for a wrapper.

## The five ParaView gaps, and what closes each

- **Volume rendering**: gmsh has no ray-caster.  ``gmsh_volume_render``
  composites N semi-transparent cut planes with a value-dependent
  opacity (``ColormapAlphaPower``: alpha grows as value**power, so low
  values fade instead of fogging).  Named for what it does.  Limits:
  per-slice not per-ray compositing, slices seen edge-on read as
  stripes (keep ``axis`` near the view direction), one CutPlane pass
  per slice.
- **Surface LIC**: not available.  ``gmsh_flow_texture`` packs
  Jobard-Lefer evenly spaced streamlines densely enough to read as a
  texture (``density`` = spacings across the plane diagonal; 60 is a
  texture, 15-20 stays countable).  Better in one way -- every curve is
  a real trajectory, so it stays probe-able -- worse in another: it
  does not fill every pixel.
- **Multi-view with a shared camera AND scale**: ``gmsh_render_panels``.
  ``gmsh_render_montage`` pastes independent renders, so each panel
  auto-fits its own scene and autoscales its own colour bar -- the
  panels LOOK comparable while encoding different scales.  Panels share
  the zoom via a hidden 8-point frame spanning the union bounding box
  (gmsh refits on every draw and IGNORES ``General.Min*/Max*`` and
  ``ZoomFactor`` -- measured -- so a common bounding box is the only
  mechanism), and share the colour range from ``gmsh_field_range``.
  Sharing a range across DIFFERENT quantities is REFUSED (T and A/m^2
  on one bar means nothing): pass ``view=``, an explicit range, or
  ``share_color=False``.
- **Cross-file colour range**: ``gmsh_field_range`` unions min/max over
  files/views with the pure-Python reader (no gmsh launch); feed it to
  ``color={"range": [...]}``.
- **Compound selection (Find Data)**: ``gmsh_select`` evaluates a
  boolean expression per element over ``x, y, z``, ``v0, v1, ...`` and
  the view names (``B`` -> ``b``), so a query can mix fields with each
  other and with position.  ``carry`` rides the chosen view's VALUES
  into the extraction -- extracting the bare 1/0 mask gives a flat blob
  whose colour bar reads "1".  Unknown names raise with the list.

## Honest gaps (do not fake these)

- VOLUME RENDERING: gmsh itself has no ray-caster, and that has not
  changed -- but the lane now carries a REAL one: gmsh_volume_raycast
  resamples the field by probing gmsh, then marches orthographic rays
  front-to-back (emission-absorption, per-ray occlusion; for a uniform
  medium the returned transmittance_min equals (1-alpha)^n exactly --
  tested).  It is the same resample-to-image approach ParaView's GPU
  volume mode uses.  ``step_files`` depth-composites STEP solids as
  opaque shaded surfaces INSIDE the march (a near plate hides the
  volume, a far plate is hidden by it -- tested both ways), so the
  coil sits IN the glowing field, not pasted on it.  Trade-off:
  standalone labelled PNG, no gmsh interactivity -- keep
  gmsh_volume_render (slice stack) when the figure must stay inside
  gmsh's scene.
- Surface LIC: gmsh still has none, but gmsh_lic is real line integral
  convolution (RK2 advection, box kernel, |v|-modulated colour) over a
  gmsh-probed plane -- every pixel carries flow direction, and the
  smear axis provably rotates with the field (tested).  Trade-off:
  direction texture, not trajectories -- gmsh_flow_texture keeps the
  probe-able-curve property.  ``step_files`` overlays the CAD section
  outline (triangle-plane intersection) in black -- the conductor
  cross-section on the texture; same standalone-PNG limitation.
- Plugin(CutSphere): returns an EMPTY view on this build -- not
  exposed.  Use gmsh_cut_plane_extract / gmsh_threshold instead.
- Plugin(Summation): absent from this build ("Unknown plugin"); a
  bare run() call even fabricates a misleading view.  Sum two views
  with gmsh_math_eval(other_view=..., "v0+w0").
- Plugin(StreamLines): only re-emits seed points -- gmsh_streamlines
  implements its own RK4 tracer instead.

## Measured pitfalls behind the tools

- Transform / Warp / Smooth / ModulusPhase run IN PLACE on the input
  view (Warp even moves the model nodes).  The tools materialize
  copies where the original must survive (mirror/transform) and
  write the view directly where the displacement IS the result
  (warp).
- gmsh.view.combine("elements", "name") RENAMES the merged view to
  "<name>_Combine".
- gmsh.view.probe on VECTOR list views returns the NEAREST value at
  ANY distance (scalar views return empty outside).  Always gate
  "found" on the returned distance -- gmsh_probe and the streamline
  tracer do this internally.
- A wrong count in a $Nodes/$Elements header header does not raise in
  gmsh: it CRASHES the process with heap corruption (0xC0000374).
  Run gmsh_validate_msh before feeding hand-written MSH to anything.
"""


GMSH_MAGNETIC_MOVIES = """
# Magnetic field lines and transient movies (TEAM28 recipe)

## Physical input contract
Use checked Cubit -> .vol -> NGSolve results; Gmsh is post-processing only.
Trace the signed instantaneous vector B(x,t), never |B| or the norm of a
complex phasor. For exp(+i*omega*t), reconstruct Re(B_hat*exp(i*omega*t)).
Field lines are instantaneous integral curves, not particle trajectories.
Reversing B can leave the curve geometry unchanged: show direction arrows
and a signed B probe alongside magnitude colours. Do not claim that moving
geometry alone establishes a coupled electromagnetic/mechanical solution.

## Choose the existing MCP tool
- gmsh_streamlines traces a vector view using adaptive arc-length RK4,
  bidirectional integration and loop handling. Select the intended view
  and time step explicitly; do not rely on Plugin(StreamLines).
- gmsh_flux_lines contours a scalar potential/flux function. In axisymmetry
  use psi=r*A_theta; equal psi intervals encode flux, not equal pixel spacing.
- gmsh_particle_trace is for charged-particle trajectories, not magnetic lines.
- Inspect tool schemas before execution. Validate input MSH with
  gmsh_validate_msh; use gmsh_render and gmsh_export_animation for existing
  post-views. The TEAM28 renderer below additionally reconstructs moving
  geometry and field lines at every saved solver step.

## Seeding and three-dimensional coverage
For an axisymmetric field, weight seed radii on a transverse plane by
2*pi*r*abs(Bz), integrate the crossing flux, then take equal-flux quantiles.
Repeat over the full azimuth, not just the XZ plane. TEAM28 uses 10 radii
and 12 azimuths (120 seeds), midpoint RK2 along normalized signed B,
0.8 mm steps and at most 400 steps in each direction. These are example
display parameters, not universal convergence settings. Stop at domain
boundaries and near nulls; check step refinement and field interpolation.
Equal crossing flux does not give uniform screen density. Projection,
occlusion and field strength affect density; a closed line can cross the
seed plane twice, so 120 seeds need not mean 120 distinct flux tubes.
For general 3D fields, distribute seeds over a surface and trace the full
vector field. Revolve a meridian only when the physical solution is axisymmetric.

## Time resolution and moving meshes
Save and render every required physical step. At dt=0.2 ms there are 100
samples per 50 Hz cycle and 50 per 100 Hz cycle. Slowing playback cannot
recover omitted samples. Magnitude and force may vary at 100 Hz while signed
B varies at 50 Hz. Synchronize current, force, signed B, geometry and lines
with the same physical timestamp; rebuild lines for each frame.
For ALE, distinguish reference mesh lookup coordinates from physical
coordinates. Evaluate with the matching deformation, then map/resample
onto the physical tracing grid. Check mapping monotonicity/Jacobians.
Use fixed colour ranges, equal spatial scaling and visible coil geometry.
Render the plate side/bottom at physical thickness; keep surface-field
offsets separate from thickness. Record camera and any exaggerated scale.

## Gmsh 5 rendering traps and acceptance
Initialize FLTK after loading the first complete scene. Clear previous
views explicitly: selecting one view's TimeStep does not synchronize
separate line views. Use gmsh.fltk.wait(), not fltk.check(). Do not use
`Delete View;` as a scene reset. Use plain native .pos files unless gzip
loading has been verified. Preserve entity colours in the launch .geo or
reapply them through the API; MSH alone does not preserve the display setup.
View.IntervalsType=2 is continuous; vector arrow sizes are pixels.
Use numeric frame order (e.g. frame-%05d), no implicit interpolation or
periodic reset. Check unique element tags, complete/finite field data,
physical geometry, all source steps exactly once, and decoded movie frame
counts. Keep JSON evidence with solver identity, versions, units, dt,
camera, seeding, source hashes and rendering limitations. Native replay
keys may be sparse: report their stride separately from the full movie.

## Executable reference in the TEAM28 development checkout
The independently committed TEAM28 development work (commit bd256fbb1)
contains validation_test/maglev/team28_fixed/README.md and the scripts below.
This knowledge change does not integrate those numerical samples into main.
Check that the files exist in the selected checkout before offering execution;
if absent, use the existing MCP tools above or obtain that reviewed source.
Run the scripts with --help first, then follow that README in order:
transient.py -> moving_movie.py -> verify_movie.py -> encode_motion.py.
The renderer supports disjoint --start-frame/--end-frame chunks and
--merge-chunks; verify complete coverage before encoding. encode_motion.py
creates the full slow movie and a three-cycle 50/100 Hz waveform zoom.
These validation scripts are checkout resources, not installed-wheel tools.
The packaged knowledge is available without NGSolve/Gmsh; execution needs
the documented solver inputs and optional runtimes. Render on a 100/102
desktop, preserve caller-owned sessions, and close only operation-owned
windows. This recipe does not certify an experimental TEAM28 trajectory.
"""


def get_gmsh_documentation(topic: str = "all") -> str:
    """Return GMSH usage documentation by topic.

    Args:
        topic: One of: all, policy, overview, cli, shortcuts, options, opt_file,
               msh_format, geo, high_order, workflow, onelab, pitfalls,
               animation, magnetic_movies, paraview

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
        "magnetic_movies": GMSH_MAGNETIC_MOVIES,
        "paraview": GMSH_PARAVIEW_PARITY,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        available = ", ".join(topics.keys())
        return f"Unknown topic: '{topic}'. Available: all, {available}"
