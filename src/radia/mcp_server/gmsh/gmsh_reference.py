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


def get_gmsh_reference(topic: str = "all") -> str:
    """Return GMSH reference documentation by topic.

    Args:
        topic: algorithms, formats, fields, mesh_options, transfinite,
               boolean, extrusion, periodic, mesh_commands, plugins,
               view_options, geometry_options

    Returns:
        Documentation string.
    """
    topics = {
        "algorithms": GMSH_MESH_ALGORITHMS,
        "formats": GMSH_OUTPUT_FORMATS,
        "fields": GMSH_SIZE_FIELDS,
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
