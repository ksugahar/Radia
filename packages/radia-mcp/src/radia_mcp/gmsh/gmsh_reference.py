"""
GMSH comprehensive reference: all options, fields, plugins, formats.

Extracted from https://gmsh.info/doc/texinfo/gmsh.html (v4.15.2).
"""

GMSH_MESH_ALGORITHMS = """
# GMSH Mesh Algorithms

## 2D Algorithms (Mesh.Algorithm)

| Code | Name | Best For |
|------|------|----------|
| 1 | MeshAdapt | Complex size fields, adaptive |
| 2 | Automatic | Auto-select |
| 3 | Initial mesh only | No refinement |
| 5 | Delaunay | Large size gradients |
| 6 | Frontal-Delaunay | **Default**. Highest quality |
| 7 | BAMG | Anisotropic 2D |
| 8 | Frontal-Delaunay for Quads | Right-triangles for quad recombination |
| 9 | Packing of Parallelograms | |
| 11 | Quasi-structured Quad | |

## 3D Algorithms (Mesh.Algorithm3D)

| Code | Name | Best For |
|------|------|----------|
| 1 | Delaunay | **Default** |
| 3 | Initial mesh only | |
| 4 | Frontal | |
| 7 | MMG3D | Remeshing |
| 9 | R-tree | |
| 10 | HXT | Parallel Delaunay |
"""

GMSH_OUTPUT_FORMATS = """
# GMSH Output Formats

## Mesh Formats (Mesh.Format)

| Code | Format | Extension | Description |
|------|--------|-----------|-------------|
| 1 | msh | .msh | GMSH native |
| 2 | unv | .unv | Universal (IDEAS) |
| 10 | auto | * | **Default**. Auto from extension |
| 16 | vtk | .vtk | VTK |
| 19 | vrml | .wrl | VRML |
| 27 | stl | .stl | STL |
| 30 | mesh | .mesh | MEDIT |
| 31 | bdf | .bdf | Nastran BDF |
| 32 | cgns | .cgns | CGNS |
| 33 | med | .med | MED |
| 39 | inp | .inp | Abaqus INP |
| 42 | su2 | .su2 | SU2 |
| 49 | neu | .neu | Gambit Neutral |
| 50 | matlab | .m | MATLAB |

## GMSH API: High-Order Node Ordering Verification

**POLICY**: Verify HO node ordering via `gmsh.model.mesh.getElementProperties()`.
This returns reference node positions that definitively establish edge ordering.

```python
import gmsh
gmsh.initialize()

# Get reference node positions for TET10 (type 11)
name, dim, order, nn, ref_pts, _ = gmsh.model.mesh.getElementProperties(11)
# ref_pts: [x0,y0,z0, x1,y1,z1, ...] for all 10 nodes
# Nodes 0-3: vertices, Nodes 4-9: edge midpoints
# Midpoint of edge (a,b) at 0.5*(ref[a] + ref[b])
```

**Verified GMSH edge ordering** (differs from Nastran!):

| Element | GMSH edge order (from API) |
|---------|---------------------------|
| TET10   | (0,1),(1,2),(0,2),(0,3),(2,3),(1,3) |
| HEX20   | (0,1),(0,3),(0,4),(1,2),(1,5),(2,3),(2,6),(3,7),(4,5),(4,7),(5,6),(6,7) |
| PRISM15 | (0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,4),(3,5),(4,5) |
| PYRAMID13 | (0,1),(0,3),(0,4),(1,2),(1,4),(2,3),(2,4),(3,4) |
| TRI6    | (0,1),(1,2),(0,2) |
| QUAD8   | (0,1),(1,2),(2,3),(0,3) |

**Volume verification via getJacobians**:

```python
gmsh.open("exported.msh")
etypes, etags, ntags = gmsh.model.mesh.getElements(dim=3)
for et in etypes:
    local_coords, weights = gmsh.model.mesh.getIntegrationPoints(int(et), "Gauss4")
    jac, det, pts = gmsh.model.mesh.getJacobians(int(et), local_coords)
    # Signed volume = sum(det[i*n_gp+j] * weights[j])
    # Negative det at some Gauss points is expected for coarse curved meshes
    # (isoparametric self-intersection, NOT a node ordering bug)
```

## PostProcessing Formats (PostProcessing.Format)

| Code | Format |
|------|--------|
| 0 | ASCII view |
| 1 | Binary view |
| 2 | Parsed view |
| 5 | Gmsh mesh |
| 10 | Automatic (**default**) |
"""

GMSH_SIZE_FIELDS = """
# GMSH Mesh Size Fields

Size fields define spatially-varying element sizes. Chain them with
Min/Max fields and set as `Background Field`.

## Usage Pattern
```
// Disable other size sources
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

// Define fields
Field[1] = Distance;
Field[1].CurvesList = {1, 2};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.01;
Field[2].SizeMax = 0.1;
Field[2].DistMin = 0.05;
Field[2].DistMax = 0.5;

Field[3] = Min;
Field[3].FieldsList = {2};

Background Field = 3;
```

## Field Types

### Distance
Distance to points, curves, or surfaces.
```
Field[n] = Distance;
Field[n].CurvesList = {curve_tags};
Field[n].PointsList = {point_tags};
Field[n].SurfacesList = {surface_tags};
Field[n].Sampling = 20;  // points per curve
```

### Threshold
Size ramp based on distance field (most common).
```
Field[n] = Threshold;
Field[n].InField = distance_field_id;
Field[n].SizeMin = 0.01;   // size at DistMin
Field[n].SizeMax = 0.1;    // size at DistMax
Field[n].DistMin = 0.05;   // start ramp
Field[n].DistMax = 0.5;    // end ramp
Field[n].Sigmoid = 0;      // 0: linear, 1: sigmoid
Field[n].StopAtDistMax = 0; // 1: keep SizeMax beyond DistMax
```

### MathEval
Arbitrary mathematical expression of x, y, z.
```
Field[n] = MathEval;
Field[n].F = "0.01 + 0.1*Sqrt(x*x + y*y)";
```

### Box
Step-change inside a box region.
```
Field[n] = Box;
Field[n].VIn = 0.01;   // size inside
Field[n].VOut = 0.1;    // size outside
Field[n].XMin = -1; Field[n].XMax = 1;
Field[n].YMin = -1; Field[n].YMax = 1;
Field[n].ZMin = -1; Field[n].ZMax = 1;
Field[n].Thickness = 0.1;  // transition zone
```

### Ball
Spherical region.
```
Field[n] = Ball;
Field[n].VIn = 0.01;
Field[n].VOut = 0.1;
Field[n].XCenter = 0; Field[n].YCenter = 0; Field[n].ZCenter = 0;
Field[n].Radius = 0.5;
Field[n].Thickness = 0.1;
```

### Cylinder
Cylindrical region.
```
Field[n] = Cylinder;
Field[n].VIn = 0.01;
Field[n].VOut = 0.1;
Field[n].XCenter = 0; Field[n].YCenter = 0; Field[n].ZCenter = 0;
Field[n].XAxis = 0; Field[n].YAxis = 0; Field[n].ZAxis = 1;
Field[n].Radius = 0.5;
```

### Constant
Constant size for specific entities.
```
Field[n] = Constant;
Field[n].VIn = 0.01;
Field[n].VOut = 0.1;
Field[n].SurfacesList = {1, 2};
Field[n].VolumesList = {1};
```

### BoundaryLayer
Boundary layer mesh insertion.
```
Field[n] = BoundaryLayer;
Field[n].CurvesList = {1, 2, 3};
Field[n].Size = 0.001;       // first layer height
Field[n].Ratio = 1.2;        // growth ratio
Field[n].NbLayers = 10;      // number of layers
Field[n].Thickness = 0.01;   // total thickness
Field[n].Quads = 1;           // 1: quad layers
Field[n].AnisoMax = 1e10;
```

### AutomaticMeshSizeField
Automatic sizing from curvature and features.
```
Field[n] = AutomaticMeshSizeField;
Field[n].features = 1;
Field[n].gradation = 1.1;
Field[n].nPointsPerCircle = 20;
Field[n].nPointsPerGap = 0;
```

### Combination Fields
```
Field[n] = Min;     // Minimum of fields (finest mesh wins)
Field[n].FieldsList = {1, 2, 3};

Field[n] = Max;     // Maximum of fields
Field[n].FieldsList = {1, 2};

Field[n] = Restrict; // Apply field only to specific entities
Field[n].InField = source_id;
Field[n].SurfacesList = {1};
Field[n].VolumesList = {1};
```

### PostView
Use a post-processing view as size field.
```
Merge "bgmesh.pos";
Field[n] = PostView;
Field[n].ViewIndex = 0;
```

## Size Determination Priority
1. Model bounding box size
2. Mesh.MeshSizeFromPoints (geometry point lc values)
3. Mesh.MeshSizeFromCurvature (automatic)
4. Background mesh size field
5. Per-entity constraints

Result clamped to [Mesh.MeshSizeMin, Mesh.MeshSizeMax] * Mesh.MeshSizeFactor.
"""

GMSH_ALL_MESH_OPTIONS = """
# GMSH Mesh Options (Complete)

## Algorithm & Control

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.Algorithm | 6 | 2D algorithm |
| Mesh.Algorithm3D | 1 | 3D algorithm |
| Mesh.AlgorithmSwitchOnFailure | 1 | Switch on failure |
| Mesh.ElementOrder | 1 | Element order |
| Mesh.SecondOrderIncomplete | 0 | Incomplete 2nd order (serendipity) |
| Mesh.SecondOrderLinear | 0 | Linear interpolation only |
| Mesh.SubdivisionAlgorithm | 0 | 0: none, 1: quads, 2: hexes, 3: barycentric |

## Size Control

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.MeshSizeFactor | 1 | Global scaling |
| Mesh.MeshSizeMin | 0 | Minimum size |
| Mesh.MeshSizeMax | 1e22 | Maximum size |
| Mesh.MeshSizeFromPoints | 1 | Use geometry point sizes |
| Mesh.MeshSizeFromCurvature | 0 | Elements per 2*pi |
| Mesh.MeshSizeExtendFromBoundary | 1 | Extend from boundary |
| Mesh.AnisoMax | 1e33 | Max anisotropy |
| Mesh.SmoothRatio | 1.8 | Size ratio (BAMG) |

## High-Order Optimization

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.HighOrderOptimize | 0 | 0: none, 1: opt, 2: elastic+opt, 3: elastic, 4: fast |
| Mesh.HighOrderIterMax | 100 | Max iterations |
| Mesh.HighOrderPassMax | 25 | Max passes |
| Mesh.HighOrderThresholdMin | 0.1 | Min quality threshold |
| Mesh.HighOrderThresholdMax | 2 | Max quality threshold |
| Mesh.HighOrderPoissonRatio | 0.33 | Poisson ratio for elastic |
| Mesh.HighOrderDistCAD | 0 | Optimize CAD distance |

## Optimization

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.Optimize | 1 | Optimize tet quality |
| Mesh.OptimizeThreshold | 0.3 | Quality threshold |
| Mesh.OptimizeNetgen | 0 | Use Netgen optimizer |
| Mesh.Smoothing | 1 | Smoothing steps |

## Recombination (Quad/Hex)

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.RecombinationAlgorithm | 1 | 0: simple, 1: blossom, 2: full-quad, 3: blossom full |
| Mesh.RecombineAll | 0 | Recombine all surfaces |
| Mesh.RecombineOptimizeTopology | 5 | Topology opt passes |
| Mesh.RecombineMinimumQuality | 0.01 | Min quad quality |
| Mesh.Recombine3DAll | 0 | Recombine all volumes |

## Transfinite

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.FlexibleTransfinite | 0 | Allow modification |
| Mesh.TransfiniteTri | 0 | 3-sided surface support |
| Mesh.QuasiTransfinite | 0 | Non-matching sides |

## File Format

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.Format | 10 | Output format (10=auto) |
| Mesh.MshFileVersion | 4.1 | MSH version |
| Mesh.Binary | 0 | Binary output |
| Mesh.SaveAll | 0 | Save all elements |
| Mesh.SaveParametric | 0 | Save parametric coords |
| Mesh.ScalingFactor | 1 | Global scaling on save |
| Mesh.BdfFieldFormat | 1 | Nastran format |

## Partitioning (METIS)

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.NbPartitions | 0 | Number of partitions |
| Mesh.MetisAlgorithm | 1 | 1: Recursive, 2: K-way |
| Mesh.PartitionCreateTopology | 1 | Create partition BRep |
| Mesh.PartitionCreatePhysicals | 1 | Physical groups per partition |
| Mesh.PartitionCreateGhostCells | 0 | Ghost cells |
| Mesh.PartitionSplitMeshFiles | 0 | One file per partition |

## Display

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.ColorCarousel | 1 | Color mode |
| Mesh.SurfaceEdges | 1 | Show surface edges |
| Mesh.SurfaceFaces | 0 | Show surface faces |
| Mesh.VolumeEdges | 1 | Show volume edges |
| Mesh.VolumeFaces | 0 | Show volume faces |
| Mesh.Nodes | 0 | Show nodes |
| Mesh.NodeSize | 4 | Node pixel size |
| Mesh.NumSubEdges | 2 | **High-order display subdivisions** |
| Mesh.Explode | 1 | Element shrink factor |
| Mesh.Light | 1 | Lighting |
| Mesh.DrawSkinOnly | 0 | 3D skin only |

## Misc

| Option | Default | Description |
|--------|---------|-------------|
| Mesh.Renumber | 1 | Renumber continuously |
| Mesh.FirstNodeTag | 1 | First node tag |
| Mesh.FirstElementTag | 1 | First element tag |
| Mesh.MinCircleNodes | 7 | Min nodes on circles |
| Mesh.MeshOnlyVisible | 0 | Mesh visible only |
"""

GMSH_TRANSFINITE = """
# Transfinite (Structured) Meshing

## Curves
```
Transfinite Curve {tag_list} = N;
Transfinite Curve {tag_list} = N Using Progression P;
Transfinite Curve {tag_list} = N Using Bump B;
```
- N: number of nodes (N-1 elements)
- Progression P: geometric ratio (P=2 -> each segment 2x previous)
- Bump B: refinement at both ends (symmetric)

## Surfaces
```
Transfinite Surface {tag_list};
Transfinite Surface {tag_list} = {p1, p2, p3, p4};
Transfinite Surface {tag_list} Left|Right|Alternate;
```
- 3 or 4 corner points (auto-detected if omitted)
- 3-sided: set `Mesh.TransfiniteTri = 1`
- Orientation: triangle diagonal direction

## Volumes
```
Transfinite Volume {tag_list};
Transfinite Volume {tag_list} = {p1, p2, p3, p4, p5, p6};  // prism
Transfinite Volume {tag_list} = {p1, ..., p8};               // hex
```
- 5 or 6 face volumes only

## Recombine (Quad/Hex)
```
Recombine Surface {tag_list};
Recombine Surface {tag_list} = angle;  // max angle for recombination
```

## Requirements
- All curves bounding a transfinite surface must have transfinite constraints
- Opposite edges must have the same number of nodes
"""

GMSH_BOOLEAN_OPS = """
# Boolean Operations (OpenCASCADE Only)

Requires `SetFactory("OpenCASCADE");`

## Operations
```
BooleanIntersection { Volume{1}; Delete; } { Volume{2}; Delete; }
BooleanUnion { Volume{1}; Delete; } { Volume{2}; Delete; }
BooleanDifference { Volume{1}; Delete; } { Volume{2}; Delete; }
BooleanFragments { Volume{1}; Delete; } { Volume{2}; Delete; }
```
- First block = object, second block = tool
- `Delete` removes original after operation
- Works with `[Physical] Point|Curve|Surface|Volume`

## BooleanFragments (Key for FEM)
Computes all intersections, creating conformal interfaces.
Essential for multi-body meshing (no duplicate surfaces at shared boundaries).
```
BooleanFragments { Volume{1}; Delete; } { Volume{2,3,4}; Delete; }
// All volumes now share conformal interfaces
```

## Assign to tag
```
BooleanDifference(100) = { Volume{1}; Delete; } { Volume{2}; Delete; };
// Result stored as Volume 100
```

## OCC Primitives
```
Box(tag) = {x, y, z, dx, dy, dz};
Sphere(tag) = {xc, yc, zc, radius};
Sphere(tag) = {xc, yc, zc, radius, angle1, angle2, angle3};
Cylinder(tag) = {x, y, z, dx, dy, dz, radius};
Cylinder(tag) = {x, y, z, dx, dy, dz, radius, angle};
Torus(tag) = {xc, yc, zc, r1, r2};
Cone(tag) = {x, y, z, dx, dy, dz, r1, r2};
Wedge(tag) = {x, y, z, dx, dy, dz};
Wedge(tag) = {x, y, z, dx, dy, dz, ltx};
Rectangle(tag) = {x, y, z, dx, dy};
Rectangle(tag) = {x, y, z, dx, dy, roundedRadius};
Disk(tag) = {xc, yc, zc, rx};
Disk(tag) = {xc, yc, zc, rx, ry};
```

## Fillet / Chamfer
```
Fillet { Volume{v}; } { Curve{c_list}; } { radius_list }
Chamfer { Volume{v}; } { Curve{c_list}; } { Surface{s_list}; } { distance_list }
```
"""

GMSH_EXTRUSION = """
# Extrusion

## Translation
```
out[] = Extrude {dx, dy, dz} { Surface{1}; };
// out[0] = top surface, out[1] = volume, out[2:] = lateral surfaces
```

## Rotation
```
out[] = Extrude { {ax,ay,az}, {px,py,pz}, angle } { Surface{1}; };
// Rotate around axis (ax,ay,az) through point (px,py,pz) by angle (radians)
```

## Twist (Translation + Rotation)
```
out[] = Extrude { {dx,dy,dz}, {ax,ay,az}, {px,py,pz}, angle } { Surface{1}; };
```

## Along Wire (OCC only)
```
Wire(1) = {curve_tags};
out[] = Extrude { Surface{1}; } Using Wire {1};
```

## Structured Layers
```
Extrude {0,0,1} { Surface{1}; Layers{10}; }              // 10 uniform layers
Extrude {0,0,1} { Surface{1}; Layers{ {3,7}, {0.3,1} }; } // 3 layers to 30%, 7 to 100%
Extrude {0,0,1} { Surface{1}; Layers{10}; Recombine; }   // Hex/prism layers
```

## ThruSections (OCC only)
```
ThruSections(tag) = {curve_loop_list};
Ruled ThruSections(tag) = {curve_loop_list};
```
"""

GMSH_PERIODIC = """
# Periodic Mesh

## Periodic Curves
```
Periodic Curve {slave} = {master};
```

## Periodic Surfaces
```
Periodic Surface {slave} = {master};
Periodic Surface slave_tag { slave_edges } = master_tag { master_edges };
```

## With Transformation
```
Periodic Curve {slave} = {master} Translate {dx, dy, dz};
Periodic Curve {slave} = {master} Rotate { {ax,ay,az}, {px,py,pz}, angle };
Periodic Surface {slave} = {master} Affine { a11,a12,a13,a14, a21,..., a44 };
```
"""

GMSH_MESH_COMMANDS = """
# GMSH Mesh Scripting Commands

## Generation
```
Mesh 1;                   // 1D mesh
Mesh 2;                   // 2D mesh
Mesh 3;                   // 3D mesh
```

## Modification
```
SetOrder N;               // Set element order (1, 2, ...)
RefineMesh;               // Uniform refinement
RecombineMesh;            // Recombine to quads/hexes
PartitionMesh N;          // Partition into N parts
```

## Optimization
```
OptimizeMesh "Gmsh";              // Default tet optimizer
OptimizeMesh "Netgen";            // Netgen optimizer
OptimizeMesh "HighOrder";         // High-order optimization
OptimizeMesh "HighOrderElastic";  // Elastic smoother
OptimizeMesh "HighOrderFastCurving"; // Fast curving
OptimizeMesh "Laplace2D";         // Laplace smoothing
```

## Embedding
```
Point {p_tags} In Surface {s_tag};    // Embed points in surface
Curve {c_tags} In Surface {s_tag};    // Embed curves in surface
Point {p_tags} In Volume {v_tag};     // Embed points in volume
Curve {c_tags} In Volume {v_tag};     // Embed curves in volume
Surface {s_tags} In Volume {v_tag};   // Embed surfaces in volume
```

## Save
```
Save "output.msh";
```
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


GMSH_PYTHON_API_FIELD = '''
# GMSH Python API — `gmsh.model.mesh.field.*` (reference only)

**Lab policy (2026-04-19)**: **This API is NOT used in the lab workflow.**
GMSH は研究室方針として **ポスト処理専用**であり、mesh 生成は **Netgen (tet)
/ Cubit (hex)** の責任範囲。mesh size field は mesh 生成の道具なので、本節
は **他文献/外部コードを読む際の参照情報**として保持する。

- 実操作のメッシュ制御は Netgen / Cubit 側で行う（Cubit は `.jou`、Netgen は
  NGSolve Python API）
- ポスト処理で gmsh を駆動する場合の API は `gmsh_reference("python_api_postproc")` を見る
- 既存 lint rule（`gmsh.model.mesh.generate()` 等を CRITICAL で禁止）はこの
  方針に沿って意図的に維持されている

Source: `api/gmsh.py` in upstream gmsh (verified against 4.15.x).

## Function list (12 functions)

| Function | Signature | Purpose |
|---|---|---|
| add | `add(fieldType, tag=-1) -> int` | Create a field; returns field tag |
| remove | `remove(tag)` | Delete field by tag |
| list | `list() -> vector<int>` | List all existing field tags |
| getType | `getType(tag) -> str` | Get field type name |
| setNumber | `setNumber(tag, option, value)` | Set double-valued option |
| getNumber | `getNumber(tag, option) -> float` | Read double-valued option |
| setNumbers | `setNumbers(tag, option, values)` | Set list-of-doubles option |
| getNumbers | `getNumbers(tag, option) -> list` | Read list-of-doubles option |
| setString | `setString(tag, option, value)` | Set string-valued option |
| getString | `getString(tag, option) -> str` | Read string-valued option |
| setAsBackgroundMesh | `setAsBackgroundMesh(tag)` | Use this field as background mesh size |
| setAsBoundaryLayer | `setAsBoundaryLayer(tag)` | Use this field as boundary layer |

Snake_case aliases exist for every function:
`set_number = setNumber`, `get_type = getType`, etc.

**Note**: there is **no** `setStrings` / `getStrings` — only singular string
options are supported. Multi-value options are numeric (PointsList,
CurvesList, SurfacesList, FieldsList etc.) and use `setNumbers`.

## Canonical usage pattern (from upstream tutorial t10.py)

```python
import gmsh
gmsh.initialize()

# ... build geometry via gmsh.model.occ.* or gmsh.model.geo.*
gmsh.model.geo.synchronize()

# 1. Distance field on point/curve
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "PointsList", [5])
gmsh.model.mesh.field.setNumbers(1, "CurvesList", [2])
gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

# 2. Threshold field using field 1 as input
gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", lc/30)
gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15)
gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5)

# 3. Math eval (spatial function)
gmsh.model.mesh.field.add("MathEval", 3)
gmsh.model.mesh.field.setString(3, "F",
    "cos(4*3.14*x) * sin(4*3.14*y) / 10 + 0.101")

# 4. MathEval referencing other field values by name F{tag}
gmsh.model.mesh.field.add("Distance", 4)
gmsh.model.mesh.field.setNumbers(4, "PointsList", [1])
gmsh.model.mesh.field.add("MathEval", 5)
gmsh.model.mesh.field.setString(5, "F", f"F4^3 + {lc/100}")

# 5. Box field (step change inside region)
gmsh.model.mesh.field.add("Box", 6)
gmsh.model.mesh.field.setNumber(6, "VIn", lc/15)
gmsh.model.mesh.field.setNumber(6, "VOut", lc)
gmsh.model.mesh.field.setNumber(6, "XMin", 0.3); gmsh.model.mesh.field.setNumber(6, "XMax", 0.6)
gmsh.model.mesh.field.setNumber(6, "YMin", 0.3); gmsh.model.mesh.field.setNumber(6, "YMax", 0.6)
gmsh.model.mesh.field.setNumber(6, "Thickness", 0.3)

# 6. Combine: minimum of chosen fields
gmsh.model.mesh.field.add("Min", 7)
gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2, 3, 5, 6])

# 7. Activate as background size field
gmsh.model.mesh.field.setAsBackgroundMesh(7)

# When using a size field, disable competing size sources
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# Delaunay algorithm (=5) handles large size gradients better than
# Frontal-Delaunay (=6, default).
gmsh.option.setNumber("Mesh.Algorithm", 5)

gmsh.model.mesh.generate(2)
gmsh.write("out.msh")
gmsh.finalize()
```

## Optional: size callback (per-element adjustment)

```python
def mesh_size_cb(dim, tag, x, y, z, lc):
    return min(lc, 0.02 * x + 0.01)
gmsh.model.mesh.setSizeCallback(mesh_size_cb)
```

## Field option reference by field type

Each `field.add(<Type>, tag)` accepts different options. Use `setNumber`,
`setNumbers`, or `setString` depending on the option's value type below.
See `gmsh_reference("fields")` for the full option list (documented in .geo
syntax — options names are identical in the Python API).

### Distance
- `PointsList` (numbers), `CurvesList` (numbers), `SurfacesList` (numbers),
  `Sampling` (number) — points per curve
- Type mapping: setNumbers / setNumbers / setNumbers / setNumber

### Threshold
- `InField` (number, tag of source field), `SizeMin`, `SizeMax`, `DistMin`,
  `DistMax`, `Sigmoid` (0 or 1), `StopAtDistMax` (0 or 1)
- Type mapping: all setNumber

### MathEval
- `F` (string) — expression of x, y, z and `F<tag>` referencing other fields
- Type mapping: setString

### Box / Ball / Cylinder / Frustum
- `VIn`, `VOut`, `XMin/XMax/YMin/YMax/ZMin/ZMax` or `XCenter/YCenter/ZCenter/Radius`,
  `Thickness` — transition width
- Type mapping: setNumber

### BoundaryLayer
- `CurvesList` (numbers), `PointsList`, `Size`, `Ratio`, `NbLayers`, `Thickness`,
  `Quads`, `AnisoMax`
- Options mix setNumber / setNumbers; additionally activate with
  `setAsBoundaryLayer(tag)` not `setAsBackgroundMesh(tag)`

### Min / Max / Restrict
- `FieldsList` (numbers for Min/Max), `InField` (number) + entity lists for Restrict

### PostView
- `ViewIndex` or `ViewTag` (number) — reference to a view loaded with `Merge`

### AutomaticMeshSizeField
- `features`, `gradation`, `nPointsPerCircle`, `nPointsPerGap`, `hBulk`, `hMin`,
  `hMax` — all setNumber

## Common pitfalls

- **Forgot synchronize**: call `gmsh.model.geo.synchronize()` or
  `gmsh.model.occ.synchronize()` before defining fields that reference entities.
- **Tag reuse**: `field.add(..., tag)` with an existing tag will overwrite. Use
  `-1` for auto-assigned tags.
- **Background vs boundary layer**: `setAsBackgroundMesh(tag)` for scalar size
  fields; `setAsBoundaryLayer(tag)` only for `BoundaryLayer` type.
- **Multiple background**: only one background field active at a time — combine
  multiple fields with `Min` or `Max` before calling setAsBackgroundMesh.
- **Disable other size sources** when using a field-driven mesh (see tutorial
  code above), otherwise point-based sizes will dominate.
'''


GMSH_PYTHON_API_POSTPROC = '''
# GMSH Python API — Post-Processing (`gmsh.view.*`, `gmsh.plugin.*`, I/O)

**Lab policy (2026-04-19)**: gmsh は **Python API 経由で駆動する軽量ポスト
処理ツール**として使う。
- **メッシュ生成には使わない** — mesh 生成は Netgen (tet) / Cubit (hex) 側
- **Netgen はポストに向かない** — Tcl 実装、mesh ビューワを改造した造り、
  Python 駆動性が弱い → gmsh で Python 駆動のポスト処理を担う
- **.geo / GUI は使わない** — すべて `gmsh` Python モジュール経由
- 既存 lint rule はこの方針に沿って維持（mesh 生成系 API の使用を検出）

Source: `api/gmsh.py` upstream (4.15.x).

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
  `.msh` = mesh + embedded views (binary preferred for large data).
'''


GMSH_POSTPROC_RECIPES = '''
# GMSH post-processing recipes (lab-specific common patterns)

Every recipe below assumes `gmsh.initialize()` was called and a mesh
(`.msh`) was opened with field data attached via
`gmsh.view.addModelData(...)`. See `python_api_postproc` for the base
API. These are *cookbook* compositions aimed at the research-lab
application domains.

## R1. Scalar field colormap (B-magnitude, temperature, etc.)

```python
import gmsh, math

# Assume mesh already loaded and a scalar view exists.
# View contains |B| per node (1 component).
view = some_view_tag

# Use the jet colormap and fix range so multi-run comparisons align.
gmsh.view.option.setNumber(view, "ColormapNumber", 2)    # 2 = jet
gmsh.view.option.setNumber(view, "RangeType", 2)         # 2 = custom
gmsh.view.option.setNumber(view, "CustomMin", 0.0)
gmsh.view.option.setNumber(view, "CustomMax", 1.5)       # Tesla
gmsh.view.option.setNumber(view, "NbIso", 20)            # color bands
gmsh.view.option.setNumber(view, "ShowScale", 1)
gmsh.view.option.setString(view, "Name", "|B| [T]")

gmsh.write("B_mag_colormap.msh")
```

Common colormaps: 0=default, 1=spectrum, 2=jet, 5=blue-red, 8=grayscale.

## R2. Vector field (B vectors) with arrow glyphs

```python
# Vector view (3 components per node)
gmsh.view.option.setNumber(view, "VectorType", 5)   # 5 = 3D arrow
gmsh.view.option.setNumber(view, "GlyphLocation", 1)  # 1 = barycenter
gmsh.view.option.setNumber(view, "ArrowSizeMax", 40)
gmsh.view.option.setNumber(view, "ArrowSizeMin", 10)
# Sample subset for clarity (avoid one-arrow-per-element clutter)
gmsh.view.option.setNumber(view, "Sampling", 5)
gmsh.view.option.setNumber(view, "NbIso", 0)          # hide iso-bands
```

## R3. Cut plane through a 3D field

```python
# Extract a plane cut at z=0 of a 3D view to visualize interior
gmsh.plugin.setString("CutPlane", "View", str(view))
gmsh.plugin.setNumber("CutPlane", "A", 0.0)    # plane: A*x + B*y + C*z + D = 0
gmsh.plugin.setNumber("CutPlane", "B", 0.0)
gmsh.plugin.setNumber("CutPlane", "C", 1.0)
gmsh.plugin.setNumber("CutPlane", "D", 0.0)
cut_view = gmsh.plugin.run("CutPlane")

# Style the cut view separately
gmsh.view.option.setNumber(cut_view, "IntervalsType", 2)  # filled
gmsh.view.option.setString(cut_view, "Name", "B on z=0")
```

Multiple cuts (e.g., three orthogonal planes for a dipole inspection):
run `CutPlane` three times with different A/B/C/D, collect the three
`cut_view` tags.

## R4. Iso-surface (B = 0.5 T contours)

```python
gmsh.plugin.setString("Isosurface", "View", str(view))
gmsh.plugin.setNumber("Isosurface", "Value", 0.5)   # Tesla
gmsh.plugin.setNumber("Isosurface", "OtherView", -1)
iso_view = gmsh.plugin.run("Isosurface")
gmsh.view.option.setString(iso_view, "Name", "|B| = 0.5 T")
```

Multiple levels: loop over values, collect tags, use `view.combine` or
`view.addAlias` if you want shared display options.

## R5. Line probe along an arbitrary path (e.g., beam axis)

```python
# Probe |B| along the z-axis from z=-50 to z=+50 (beam trajectory)
import numpy as np
zs = np.linspace(-50, 50, 201)
b_along_z = []
for z in zs:
    vals, dist = gmsh.view.probe(view, 0.0, 0.0, float(z))
    # vals is empty list if outside mesh; guard for that
    b_along_z.append(vals[0] if vals else float("nan"))

# Outside: plot with matplotlib / export CSV, etc.
import csv
with open("Bz_along_beam.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["z_mm", "B_T"])
    for z, b in zip(zs, b_along_z): w.writerow([z, b])
```

For a vector view, probe returns 3 components — pick the relevant one
(e.g., `vals[1]` for By). For tensor/high-order views, check
`len(vals)`.

## R6. Field integration over a region (Integrate plugin)

```python
# Integrate a scalar view (e.g., Joule heating density) over its domain
# Useful for: total losses in a coil, total flux through a pole face
gmsh.plugin.setString("Integrate", "View", str(view))
total_view = gmsh.plugin.run("Integrate")
# The resulting view holds a single scalar. Read it back:
values, _ = gmsh.view.probe(total_view, 0, 0, 0)
total = values[0] if values else None
print(f"integrated scalar = {total}")
```

## R7. Stream lines on a vector field

```python
# Start seeds from a 2D patch (u,v grid) and integrate
gmsh.plugin.setString("StreamLines", "View", str(view))
gmsh.plugin.setNumber("StreamLines", "X0", -40.0)
gmsh.plugin.setNumber("StreamLines", "Y0", -40.0)
gmsh.plugin.setNumber("StreamLines", "Z0", 0.0)
gmsh.plugin.setNumber("StreamLines", "DX1", 80.0)   # first direction
gmsh.plugin.setNumber("StreamLines", "DY1", 0.0)
gmsh.plugin.setNumber("StreamLines", "DX2", 0.0)     # second direction
gmsh.plugin.setNumber("StreamLines", "DY2", 80.0)
gmsh.plugin.setNumber("StreamLines", "NumPointsU", 10)
gmsh.plugin.setNumber("StreamLines", "NumPointsV", 10)
gmsh.plugin.setNumber("StreamLines", "MaxIter", 500)
gmsh.plugin.setNumber("StreamLines", "DT", 0.2)
sl_view = gmsh.plugin.run("StreamLines")
```

Use: visualize eddy current paths, flux lines in yokes.

## R8. Compute derived fields on the fly (MathEval plugin)

```python
# Given a vector B view, create a new view of |B|
gmsh.plugin.setString("MathEval", "Expression0",
    "Sqrt(v0*v0 + v1*v1 + v2*v2)")
gmsh.plugin.setString("MathEval", "View", str(vec_view))
gmsh.plugin.setString("MathEval", "OtherView", "-1")
mag_view = gmsh.plugin.run("MathEval")
gmsh.view.option.setString(mag_view, "Name", "|B| derived")
```

Variables available in expressions:
- `v0..v8` — components of the input view
- `x, y, z` — spatial coordinates
- `w0..w8` — components of an OtherView (if set)

Composes with R3 / R4: run MathEval to get scalar magnitude, then
CutPlane / Isosurface on that magnitude view.

## R9. Multi-region selective display (masking by physical group)

When using `run_pipeline_multi` output, each region has a named
physical group.  To hide the surrounding "air" and show only
workpiece + coil:

```python
air_entity_tags = gmsh.model.getEntitiesForPhysicalGroup(3, 100003)
for ent in air_entity_tags:
    # gmsh uses model entity visibility flags
    gmsh.model.setVisibility([(3, int(ent))], 0, recursive=True)

# Still need to hide view data in the air region — but since
# addModelData is tied to element tags, per-element view visibility
# is controlled via gmsh.view.option.setNumber(view, "ShowElement", 0)
# combined with element-level filter.  Simpler: write per-region views.
```

Tip: the cleanest approach is to **attach one view per region**
(separate `gmsh.view.add` + `addModelData` calls with only the
region's element tags), then toggle view visibility with
`gmsh.view.option.setNumber(view_tag, "Visible", 0 or 1)`.

## R10. Batch snapshot (PNG/PDF) for figures

```python
gmsh.initialize(argv=["gmsh", "-nopopup"])
gmsh.open("result.msh")
# ... configure view options as in R1/R2/R3 ...
gmsh.option.setNumber("General.BackgroundGradient", 0)
gmsh.option.setNumber("General.GraphicsWidth", 1200)
gmsh.option.setNumber("General.GraphicsHeight", 900)
gmsh.fltk.initialize()             # requires a display on Linux; OK headless on Windows
gmsh.write("figure.png")           # PNG
gmsh.write("figure.pdf")           # vector
gmsh.fltk.finalize()
gmsh.finalize()
```

For CI without display, use `gmsh.graphics.draw()` + PostScript export
or run under Xvfb on Linux. On Windows this works off-screen natively
in most GMSH builds.

## Recipe → research application mapping

| Recipe | Applies to |
|---|---|
| R1 colormap | IH temperature map, |B| map in magnets |
| R2 vector glyphs | B-field direction in dipole/quad |
| R3 cut plane | Dipole mid-plane field view, IH axial cut |
| R4 iso-surface | Good-field region of a magnet |
| R5 line probe | Beam trajectory integrand, axial B profile |
| R6 integrate | Joule heating total, flux through surface |
| R7 stream lines | Eddy current visualization |
| R8 MathEval | |B| from vector B, derived energy densities |
| R9 multi-region | Hide air to declutter coil/yoke plots |
| R10 batch PNG/PDF | Paper/thesis figures, CI artifacts |
'''


def get_gmsh_reference(topic: str = "all") -> str:
    """Return GMSH reference documentation by topic.

    Args:
        topic: algorithms, formats, fields, python_api_field, python_api_postproc,
               mesh_options, transfinite, boolean, extrusion, periodic,
               mesh_commands, plugins, view_options, geometry_options

    Returns:
        Documentation string.
    """
    topics = {
        "algorithms": GMSH_MESH_ALGORITHMS,
        "formats": GMSH_OUTPUT_FORMATS,
        "fields": GMSH_SIZE_FIELDS,
        "python_api_field": GMSH_PYTHON_API_FIELD,
        "python_api_postproc": GMSH_PYTHON_API_POSTPROC,
        "postproc_recipes": GMSH_POSTPROC_RECIPES,
        "mesh_options": GMSH_ALL_MESH_OPTIONS,
        "transfinite": GMSH_TRANSFINITE,
        "boolean": GMSH_BOOLEAN_OPS,
        "extrusion": GMSH_EXTRUSION,
        "periodic": GMSH_PERIODIC,
        "mesh_commands": GMSH_MESH_COMMANDS,
        "plugins": GMSH_PLUGINS,
        "view_options": GMSH_VIEW_OPTIONS,
        "geometry_options": GMSH_GEOMETRY_OPTIONS,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        available = ", ".join(topics.keys())
        return f"Unknown topic: '{topic}'. Available: all, {available}"
