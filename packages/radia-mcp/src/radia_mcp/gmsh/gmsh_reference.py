"""Native Gmsh post-processing reference; solver meshing is owned by Netgen/Cubit.

API recipes require the qualified Gmsh 5 runtime and an existing MSH 4.2 result.
Read the workflow and format contract with gmsh_usage before native access.
Historical observations in retained API notes are provenance, not new-build
qualification: run the current semantics and round-trip checks on adoption.
"""

GMSH_PLUGINS = """
# GMSH Built-in Plugins (Selection)

## Field Computation
| Plugin | Description |
|--------|-------------|
| Curl | Compute curl of vector field |
| Divergence | Compute divergence |
| Gradient | Compute gradient |
| Eigenvalues | Tensor eigenvalues |
| Eigenvectors | Tensor eigenvectors |
| Lambda2 | Vortex criterion |

## Extraction / Slicing
| Plugin | Description |
|--------|-------------|
| CutPlane | Extract on plane |
| CutSphere | Extract on sphere |
| CutBox | Extract in box region |
| CutGrid | Extract on regular grid |
| CutParametric | Extract on parametric curve |
| Isosurface | Extract isosurface |
| Skin | Extract mesh skin |

## Analysis
| Plugin | Description |
|--------|-------------|
| Integrate | Integrate field over domain |
| MinMax | Find min/max values |
| MeshVolume | Compute mesh volume |
| AnalyseMeshQuality | Jacobian, IGE, ICN quality |
| Probe | Probe field at point |

## Modification
| Plugin | Description |
|--------|-------------|
| MathEval | Evaluate expression on view |
| ModifyComponents | Modify view components |
| Smooth | Smooth view data |
| Warp | Warp by vector field |
| Transform | Geometric transformation |
| Summation | Sum views |

## Conversion
| Plugin | Description |
|--------|-------------|
| HarmonicToTime | Harmonic to time domain |
| ModulusPhase | Compute modulus and phase |
| Scal2Vec | Scalar to vector |
| Scal2Tens | Scalar to tensor |

## Visualization
| Plugin | Description |
|--------|-------------|
| Annotate | Add text annotations |
| StreamLines | Compute streamlines |
| Particles | Particle tracing |
"""

GMSH_VIEW_OPTIONS = """
# GMSH View Options (View[n].*)

## Display Type
| Option | Default | Description |
|--------|---------|-------------|
| View.Visible | 1 | Show/hide view |
| View.Type | 1 | 1: 3D, 2: 2D space, 3: 2D time, 4: 2D |
| View.IntervalsType | 2 | 1: iso, 2: continuous, 3: discrete, 4: numeric |
| View.NbIso | 10 | Number of intervals |
| View.DrawSkinOnly | 0 | Draw only 3D skin |
| View.ShowElement | 0 | Show element boundaries |
| View.Boundary | 0 | Draw N-b dimensional boundary |

## Scale & Range
| Option | Default | Description |
|--------|---------|-------------|
| View.RangeType | 1 | 1: default, 2: custom, 3: per step |
| View.CustomMin | 0 | Custom minimum |
| View.CustomMax | 0 | Custom maximum |
| View.ScaleType | 1 | 1: linear, 2: log, 3: double log |
| View.ShowScale | 1 | Show color bar |
| View.SaturateValues | 0 | Clamp to custom range |

## Vector Display
| Option | Default | Description |
|--------|---------|-------------|
| View.VectorType | 4 | 1: segment, 2: arrow, 3: pyramid, 4: 3D arrow, 5: displacement, 6: comet |
| View.GlyphLocation | 1 | 1: centroid, 2: node |
| View.ArrowSizeMax | 60 | Max arrow size (pixels) |
| View.ArrowSizeMin | 0 | Min arrow size |
| View.DisplacementFactor | 1 | Displacement amplification |

## Tensor Display
| Option | Default | Description |
|--------|---------|-------------|
| View.TensorType | 1 | 1: Von-Mises, 2: max eigen, 3: min eigen, 4: eigenvectors, 5: ellipse, 6: ellipsoid |

## Coloring
| Option | Default | Description |
|--------|---------|-------------|
| View.ColormapNumber | 2 | 0: black, 1: vis5d, 2: jet, 3: lucie, 4: rainbow, 5: emc2000, 6: incandescent, 7: hot, 8: pink, 9: grayscale, 10: french, 11: hsv |
| View.ColormapInvert | 0 | Invert colormap |

## Lighting
| Option | Default | Description |
|--------|---------|-------------|
| View.Light | 1 | Enable lighting |
| View.SmoothNormals | 0 | Smooth normals |

## Time Steps
| Option | Default | Description |
|--------|---------|-------------|
| View.TimeStep | 0 | Current time step |
| View.ShowTime | 3 | 0: none, 1: time, 2: harmonic, 3: auto |

## Adaptive Visualization
| Option | Default | Description |
|--------|---------|-------------|
| View.AdaptVisualizationGrid | 0 | Adaptive refinement for high-order |
| View.MaxRecursionLevel | 0 | Max recursion depth |
| View.TargetError | 0.0001 | Target error |
"""

GMSH_GEOMETRY_OPTIONS = """
# GMSH Geometry Options (Selection)

## OpenCASCADE (OCC) Options

| Option | Default | Description |
|--------|---------|-------------|
| Geometry.OCCTargetUnit | "" | Unit for STEP import ("M" for meters) |
| Geometry.OCCAutoEmbed | 1 | Auto-embed internal entities |
| Geometry.OCCAutoFix | 1 | Auto-fix orientations |
| Geometry.OCCBooleanPreserveNumbering | 1 | Preserve numbering |
| Geometry.OCCBooleanSimplify | 1 | Simplify results |
| Geometry.OCCBooleanNonDestructive | 0 | Keep original shapes |
| Geometry.OCCFixDegenerated | 0 | Fix degenerated edges |
| Geometry.OCCFixSmallEdges | 0 | Fix small edges |
| Geometry.OCCImportLabels | 1 | Import STEP labels |
| Geometry.OCCMakeSolids | 0 | Fix shells to solids |
| Geometry.OCCParallel | 0 | Parallel booleans |
| Geometry.OCCScaling | 1 | Import scaling |
| Geometry.OCCSewFaces | 0 | Sew faces |
| Geometry.Tolerance | 1e-8 | Geometrical tolerance |
| Geometry.ScalingFactor | 1 | Global geometry scaling |

## Display

| Option | Default | Description |
|--------|---------|-------------|
| Geometry.Points | 1 | Show points |
| Geometry.Curves | 1 | Show curves |
| Geometry.Surfaces | 0 | Show surfaces |
| Geometry.SurfaceType | 0 | 0: cross, 1: wireframe, 2: solid |
| Geometry.NumSubEdges | 100 | Curve display subdivisions |
"""

GMSH_PYTHON_API_POSTPROC = '''
# GMSH Python API — Post-Processing (`gmsh.view.*`, `gmsh.plugin.*`, I/O)

Gmsh 5 provides interactive GUI display and Python-driven post-processing.
Use `.geo` launch companions for repeatable GUI settings and the MCP render
tools for automated output. Cubit/Netgen own solver meshing; Gmsh consumes
existing MSH 4.2 results. STEP overlays are display geometry only.

API names/signatures are checked against the adopted official wheel; native
plugin semantics are qualified by the separately recorded reference tests.

## Top-level file I/O

```python
import gmsh
gmsh.initialize()

gmsh.open(fileName)    # Open .msh / .geo / .brep / .step etc. into current model
gmsh.merge(fileName)   # Merge file into current model (useful for .pos post-view)
gmsh.write(fileName)   # Write current model/views. Format inferred from extension
                       # (.msh / .pos / .step / .vtk / .med / .x3d / ...)
gmsh.clear()           # Clear all models and views
gmsh.finalize()
```

## `gmsh.view.*` — view management (18 functions)

| Function | Signature | Purpose |
|---|---|---|
| add | `add(name, tag=-1) -> int` | Create a new post-processing view, returns tag |
| remove | `remove(tag)` | Remove view by tag |
| getIndex | `getIndex(tag) -> int` | Get dynamic list index (for view option access) |
| getTags | `getTags() -> list[int]` | Tags of all loaded views |
| addModelData | `addModelData(tag, step, modelName, dataType, tags, data, time=0., numComponents=-1, partition=0)` | Attach scalar/vector/tensor data to model nodes or elements |
| addHomogeneousModelData | `addHomogeneousModelData(...)` | Same, flattened data vector (for uniform per-tag sizes) |
| getModelData | `getModelData(tag, step) -> (dataType, tags, data, time, numComponents)` | Retrieve model data |
| getHomogeneousModelData | `getHomogeneousModelData(tag, step) -> (...)` | Retrieve flattened data |
| addListData | `addListData(tag, dataType, numEle, data)` | Attach "list data" (geometry + values, not referencing model) |
| getListData | `getListData(tag) -> (dataType, numElements, data)` | Retrieve list data |
| addListDataString | `addListDataString(tag, coord, data, style=[])` | Add text annotations at coords |
| getListDataStrings | `getListDataStrings(tag, dim) -> (coord, data, style)` | Retrieve text annotations |
| setInterpolationMatrices | `setInterpolationMatrices(tag, type, d, coef, exp, dGeo=0, coefGeo=[], expGeo=[])` | High-order interpolation matrices |
| addAlias | `addAlias(refTag, copyOptions=False, tag=-1) -> int` | Create an alias view (shared data, separate options) |
| combine | `combine(what, how, remove=True, copyOptions=True)` | Merge views: what="elements"/"steps"/"views", how="all"/"visible"/"name" |
| probe | `probe(tag, x, y, z, step=-1, numComp=-1, gradient=False, distanceMax=0., xElemCoord=[], yElemCoord=[], zElemCoord=[], dim=-1) -> (values, distance)` | Sample a view at a point |
| write | `write(tag, fileName, append=False)` | Export a single view to file (.pos / .msh / .txt) |
| setVisibilityPerWindow | `setVisibilityPerWindow(tag, value, windowIndex=0)` | GUI-only visibility control |

Snake_case aliases exist: `add_model_data`, `get_tags`, etc.

## `gmsh.view.option.*` — per-view display options (7 functions)

```python
gmsh.view.option.setNumber(tag, optionName, value)      # numeric option
gmsh.view.option.getNumber(tag, optionName) -> float
gmsh.view.option.setString(tag, optionName, value)
gmsh.view.option.getString(tag, optionName) -> str
gmsh.view.option.setColor(tag, optionName, r, g, b, a=255)
gmsh.view.option.getColor(tag, optionName) -> (r, g, b, a)
gmsh.view.option.copy(refTag, tag)                       # Copy all options from refTag
```

Common view option names: `Visible`, `IntervalsType`, `NbIso`, `Min`, `Max`,
`CustomMin`, `CustomMax`, `ColormapNumber`, `RangeType`, `ScaleType`,
`TimeStep`, `VectorType`, `GlyphLocation`, `ShowElement`, `ShowScale`,
`PointType`, `LineType`, `LineWidth`, `Boundary`, `DrawLines`, `DrawSkinOnly`.
Full list: upstream `gmsh.html` section "View options".

## `gmsh.plugin.*` — post-processing plugins (3 functions)

```python
gmsh.plugin.setNumber(name, option, value)   # Set plugin numeric option
gmsh.plugin.setString(name, option, value)   # Set plugin string option
gmsh.plugin.run(name) -> int                 # Run plugin, returns created view tag (or -1)
```

Useful plugins for post-processing:

| Plugin | Purpose | Key options |
|---|---|---|
| `CutPlane` | Extract planar cut of a 3D view | A, B, C, D (plane coefficients), View |
| `CutSphere` | Spherical cut | Xc, Yc, Zc, R, View |
| `CutBox` | Box cut | X0..Z0, X1..Z1, X2..Z2, X3..Z3, View |
| `Isosurface` | Generate iso-contours on 3D view | Value, View |
| `Integrate` | Integrate field over elements | View, OverTime |
| `Gradient` | Compute gradient of scalar view | View |
| `Curl` | Compute curl of vector view | View |
| `Divergence` | Compute divergence | View |
| `ModulusPhase` | Split complex view into magnitude/phase | View |
| `MathEval` | Arbitrary expression on existing views | Expression0..7, View |
| `Transform` | Affine transform view coords | A11..A33, Tx..Tz, View |
| `Skin` | Extract skin (outer boundary) | View |
| `Smooth` | Smooth a view | View |
| `Triangulate` | Triangulate 2D point cloud | View |
| `Warp` | Warp mesh by vector field | Factor, TimeStep, View |
| `StreamLines` | Stream line integration | X0..DZ, NumPointsU, MaxIter, View |

Plugins always operate on existing views (by view tag); `run` creates a
new view and returns its tag (or -1 on error).

## Canonical post-processing workflow (lab-style)

```python
import gmsh
gmsh.initialize()

# 1) Load a .msh produced by Netgen or Cubit
gmsh.open("build123d_netgen_output.msh")

# 2) Attach computed field data (e.g. from NGSolve, Radia, or solver output)
view_tag = gmsh.view.add("B_field")
gmsh.view.addModelData(
    tag=view_tag,
    step=0,
    modelName=gmsh.model.getCurrent(),
    dataType="NodeData",
    tags=node_tags,       # list of node IDs
    data=B_vectors,       # list of [Bx, By, Bz] per node
    numComponents=3,
)

# 3) View options (display control)
gmsh.view.option.setNumber(view_tag, "VectorType", 5)    # 3D arrow
gmsh.view.option.setNumber(view_tag, "ShowScale", 1)
gmsh.view.option.setString(view_tag, "Name", "B [T]")

# 4) Derived views via plugins — no manual post-math needed
gmsh.plugin.setString("CutPlane", "View", str(view_tag))
gmsh.plugin.setNumber("CutPlane", "A", 1); gmsh.plugin.setNumber("CutPlane", "B", 0)
gmsh.plugin.setNumber("CutPlane", "C", 0); gmsh.plugin.setNumber("CutPlane", "D", 0)
cut_tag = gmsh.plugin.run("CutPlane")

# 5) Point probe
values, dist = gmsh.view.probe(view_tag, 0.01, 0.0, 0.05)

# 6) Export — .msh bundles mesh + view data (lab visualization format)
gmsh.view.write(view_tag, "B_field.pos")
gmsh.write("full_postproc.msh")

gmsh.finalize()
```

## Non-GUI batch mode

Useful for CI and headless runs:

```python
gmsh.initialize(argv=["gmsh", "-nopopup"])
# ... load + post-process ...
gmsh.finalize()
# No fltk.run() means no GUI window opens; everything goes through files.
```

## Pitfalls

- **modelName must match**: `gmsh.view.addModelData(..., modelName=...)` needs the
  currently loaded model's name (use `gmsh.model.getCurrent()`).
- **Tag ordering**: `tags` in addModelData must match the loaded mesh's node/element
  IDs exactly; `.msh` writer preserves IDs, Netgen→Gmsh conversion should too.
- **dataType choice**: "NodeData" for nodal fields, "ElementData" for cell-centered
  (one value per element), "ElementNodeData" for DG/high-order nodal per element.
- **numComponents**: pass explicitly (1 scalar / 3 vector / 9 tensor) to avoid
  inference errors on non-standard sizes.
- **Plugin View indexing**: `CutPlane` / `Isosurface` etc. take `View` as a
  **string index** (tag encoded as decimal string), not an int.
- **Color options**: use `setColor(r, g, b, a)` — NOT a 4-tuple argument.
- **.pos vs .msh**: `.pos` = gmsh parsed view format (ASCII, self-contained);
  `.msh` = model-based MSH 4.2 mesh + views; use ASCII for repository readers.
'''

GMSH_POSTPROC_RECIPES = """
# Post-processing recipes through the existing native Gmsh tools

Use gmsh_usage("workflow") first: inspect an existing MSH 4.2 result and
select its actual view name. The examples below assume T is scalar and B
is a three-component vector, coordinates are meters and field units are
declared by the producer. These are MCP calls, not a second postprocessor.

## Scalar colormap and vector glyphs

gmsh_render selects an existing view and accepts color/glyph/option controls.
Use a common fixed range for comparisons. Native View.VectorType=4 draws
3D arrows; 5 is displacement. Mesh.NumSubEdges and adaptive visualization
refine display only. Save the exact launch companions and returned settings.

## Cut and iso-surface of a scalar field

```python
gmsh_cut_plane_extract(msh_path="/results/field.msh", view_name="T",
                       normal=[0, 0, 1], offset=0,
                       out_file="/results/T_cut.pos")
gmsh_isosurface(msh_path="/results/field.msh", view_name="T", value=300,
                out_file="/results/T_iso.pos")
```

The extracted views are native list-data .pos files. Merge them with the
original result for display. Do not give them a .msh suffix or feed them to
the mesh/NodeData-only inspection and statistics reader. Adaptive extraction
is explicit and needs its declared resolution/error settings.

## Magnitude of a vector field

```python
gmsh_math_eval(msh_path="/results/field.msh", view_name="B",
               expressions=["Sqrt(v0^2+v1^2+v2^2)"], result_name="B_magnitude",
               out_file="/results/B_magnitude.pos")
```

MathEval applies the expression at the input sample values; interpolation
and nonlinear expressions do not commute. Probe the derived .pos with the
native probe tool, or display it. Never interpret a missing probe as zero.

## Integral and line sample

```python
gmsh_integrate(msh_path="/results/field.msh", view_name="T", dimension=3)
gmsh_line_profile(msh_path="/results/field.msh", view_name="B",
                  start=[0, 0, -0.05], end=[0, 0, 0.05], n=101)
```

The integral is a post-processing integral on the selected view/domain;
its physical unit depends on the input and integration dimension. It is
not an exact FE energy/force validation. Exact solver integrals stay with
NGSolve on the checked solver mesh. A vector line profile retains components;
choose the physical component or magnitude deliberately.

## Region display and moving results

Inspect physical names rather than hard-coding region tags. Model entity
visibility does not automatically hide field-view samples in that region.
Prefer producer-created per-region views when independent visibility is
needed. For time-dependent or moving results use gmsh_usage("animation")
and gmsh_usage("magnetic_movies"); preserve real step times, fixed camera
and comparison range. The GUI opens the generated .geo launch artifact.

All native operations require the qualified Gmsh 5 runtime. Retain the exact
wheel/version, inputs, output hashes and tool arguments. A list-view export
does not qualify as a model-based MSH 4.2 result or as any solver mesh.
"""

def get_gmsh_reference(topic: str = "all") -> str:
    """Return post-processing API, plugins, display options and recipes."""
    topics = {
        "python_api_postproc": GMSH_PYTHON_API_POSTPROC,
        "postproc_recipes": GMSH_POSTPROC_RECIPES,
        "plugins": GMSH_PLUGINS,
        "view_options": GMSH_VIEW_OPTIONS,
        "geometry_options": GMSH_GEOMETRY_OPTIONS,
    }
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    if topic in topics:
        return topics[topic]
    return f"Unknown topic: '{topic}'. Available: all, {', '.join(topics)}"
