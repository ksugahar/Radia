"""
GMSH examples and tutorials knowledge base.

Content from public-safe curated corpus tutorials, examples, and electric machine models.
"""

GMSH_TUTORIALS = """
# GMSH Tutorial Summary (public-safe curated corpus)

## t1.geo -- Geometry Basics, Physical Groups
Fundamental building blocks: Point, Line, Curve Loop, Plane Surface.
Physical Groups for tagging mesh entities.
```
Point(1) = {0, 0, 0, lc};        // x, y, z, mesh_size
Line(1) = {1, 2};                  // start, end point
Curve Loop(1) = {1, 2, 3, 4};     // ordered curves (sign = direction)
Plane Surface(1) = {1};            // bounded by curve loop
Physical Curve("boundary") = {1, 2, 3, 4};
Physical Surface("domain") = {1};
```

## t2.geo -- Transformations, Extrusion
```
Translate {dx,dy,dz} { Point{1}; }
Rotate { {ax,ay,az}, {px,py,pz}, angle } { Surface{1}; }
new_entities[] = Duplicata { Surface{1}; };  // copy + return tags
out[] = Extrude {0,0,1} { Surface{1}; };
// out[0] = top surface, out[1] = volume, out[2:] = lateral surfaces
```

## t3.geo -- Extruded Meshes, ONELAB Parameters
Structured extrusion with layer control:
```
Extrude {0,0,h} {
  Surface{1}; Layers{ {8,2}, {0.5,1} }; Recombine;
}
// 8 layers to 50% height, 2 layers to 100%, quad/hex elements
```
ONELAB interactive parameters:
```
DefineConstant[ h = {1, Min 0.1, Max 10, Step 0.1, Name "Geometry/Height"} ];
```

## t5.geo -- Macros, Loops, Holes
User-defined macros and automatic tag assignment:
```
Macro CheeseHole
  p1 = newp; Point(p1) = {x, y, z, lc};
  ...
Return
For t In {1:5}
  Call CheeseHole;
EndFor
```
Volume with holes: `Volume(1) = {exterior_loop, hole1_loop, hole2_loop};`

## t6.geo -- Transfinite (Structured) Meshes
```
Transfinite Curve {1,3} = 10;                    // 10 nodes
Transfinite Curve {2,4} = 20 Using Progression 1.2; // graded
Transfinite Surface {1} = {1,2,3,4};             // corner points
Recombine Surface {1};                             // -> quads
Mesh.Smoothing = 100;                              // elliptic smoothing
```
All 4 edges must have matching constraints. Opposite edges same node count.

## t8.geo -- Post-Processing Views
Loading and configuring views:
```
Include "view1.pos";
View[0].IntervalsType = 2;   // continuous
View[0].NbIso = 25;
View[0].ShowScale = 1;
View[0].SmoothNormals = 1;
```
Animation:
```
For step In {0:View[0].NbTimeStep-1}
  View[0].TimeStep = step;
  Draw;
  Print Sprintf("frame_%03g.png", step);
  Sleep 0.1;
EndFor
```

## t10.geo -- Mesh Size Fields
Chain Distance -> Threshold -> Min for localized refinement:
```
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

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

## t13.geo -- STL Remeshing
Remesh an STL file without CAD:
```
Merge "file.stl";
ClassifySurfaces {angle, includeBoundary, forceParametrizablePatches};
CreateGeometry;
Surface Loop(1) = Surface{:};
Volume(1) = {1};
```

## t16.geo -- OpenCASCADE CSG
```
SetFactory("OpenCASCADE");
Box(1) = {0,0,0, 1,1,1};
Sphere(2) = {0.5,0.5,0.5, 0.6};
BooleanDifference(3) = { Volume{1}; Delete; } { Volume{2}; Delete; };
```
`BooleanFragments` for conformal multi-body meshing (no duplicate interfaces).

## t20.geo -- STEP Import
```
SetFactory("OpenCASCADE");
v() = ShapeFromFile("part.step");
// Geometry.OCCTargetUnit = "M" for meters
```

## t21.geo -- Mesh Partitioning
```
Mesh 2;
PartitionMesh 4;  // Partition into 4 parts (METIS)
// Or: Plugin(SimplePartition) for chessboard-style
```
"""

GMSH_ELECTRIC_MACHINES = """
# Electric Machine Models (public-safe curated corpus)

## Overview
Complete 2D FEM models for rotating electrical machines using GMSH + GetDP.
Developed by ULiege-ULB (Sabariego, Gyselinck, Geuzaine).

## Available Models

| Model | Type | Poles | Description |
|-------|------|-------|-------------|
| pmsm | PMSM | 8 | Permanent magnet synchronous (GRUCAD) |
| pmsm_cbmag | PMSM | 8 | Complete bridge magnet variant |
| lomonova | PMSM | 8 | In-wheel motor (Lomonova) |
| im | IM | 4 | Induction motor, 18.5kW (WEG) |
| team30a | IM | 4 | TEAM Problem 30a |
| srm | SRM | 4 | Switched reluctance |
| wfsm | WFSM | 4 | Wound-field synchronous |

## File Structure (per model)
```
<model>_data.geo    # Parametric dimensions (DefineConstant)
<model>.geo         # Main geometry (includes data + rotor + stator)
<model>_rotor.geo   # Rotor geometry
<model>_stator.geo  # Stator geometry
<model>.pro         # GetDP problem definition
<model>_circuit.pro # Circuit coupling
```

## Physical Region Convention
```
ROTOR_FE     = 1000  # Rotor iron
ROTOR_AIR    = 1001  # Rotor air
ROTOR_MAGNET = 1010  # Permanent magnets
ROTOR_BAR_nn = 10nn  # Induction motor bars (individually numbered)

STATOR_FE    = 2000  # Stator iron
STATOR_IND   = 2300  # Stator windings (AP=A+, AM=A-, BP=B+, etc.)
STATOR_AIR   = 2001  # Stator air

AIR_EXT      = 3000  # Exterior air
AIR_INF      = 3100  # Infinite shell

MOVING_BAND  = 9999  # Rotor-stator coupling
```

## Symmetry Exploitation
```
NbrPolesTot = 8;
NbrPolesInModel = 2;  // Quarter model
SymmetryFactor = NbrPolesTot / NbrPolesInModel;  // = 4
```
Supported: 1 pole, 2 poles, 4 poles (quarter), 8 poles (full).

## PMSM Key Parameters (pmsm_data.geo)
- Remanent induction: 1.2 T (NdFeB)
- Relative permeability: mu_r = 1000 (iron)
- Air gap: ~1mm
- All dimensions in mm, converted to meters

## Induction Machine Key Parameters (im_data.geo)
- Power: 18.5 kW, 60 Hz
- Air gap: 0.47 mm
- Rotor: 40 bars, sigma = 28 MS/m (aluminum)
- Stator: 48 slots, 3-phase, 128 turns/phase
- Moving band at 1/3 and 2/3 of air gap
- Each rotor bar: individual Physical Region (ROTOR_BAR01..40)
"""

GMSH_EXAMPLE_PATTERNS = """
# Key GMSH Patterns and Techniques

## 1. Automatic Tag Assignment
```
p1 = newp; Point(p1) = {x, y, z, lc};
c1 = newc; Line(c1) = {p1, p2};
s1 = news; Plane Surface(s1) = {cl1};
v1 = newv; Volume(v1) = {sl1};
```

## 2. Entity Queries
```
BoundingBox Volume{1};         // Get bounding box
PointsOf { Volume{:}; }       // All points in all volumes
Surface In BoundingBox {x1,y1,z1, x2,y2,z2}  // Spatial query
CombinedBoundary { Volume{1}; }  // Outer boundary surfaces
Boundary { Surface{1}; }      // Boundary curves of surface
```

## 3. Array Arithmetic
```
pts[] = {1, 2, 3, 4, 5};
n = #pts[];           // length = 5
pts2[] = pts[{0:2}];  // slice = {1, 2, 3}
Printf("%g", pts[0]); // first element
```

## 4. Bulk Mesh Size
```
MeshSize { PointsOf{ Volume{:}; } } = 0.1;
// or per-point:
MeshSize {1, 2, 3} = 0.01;
```

## 5. Entity Coloring
```
Recursive Color Red { Volume{1}; }
Recursive Color Blue { Volume{2}; }
```

## 6. Memory-Aware Sizing
```
If (GetValue("General.TotalMemory", 2048) <= 2048)
  lc = 0.01;  // coarser for low memory
Else
  lc = 0.005;
EndIf
```

## 7. ONELAB Interactive Parameters
```
DefineConstant[
  freq = {50, Min 1, Max 1000, Step 1,
          Name "Parameters/Frequency [Hz]",
          Highlight "AliceBlue"},
  flag = {0, Choices {0, 1},
          Name "Parameters/Enable Feature"}
];
```

## 8. Include Hierarchy
```
Include "data.geo";      // Parameters
Include "rotor.geo";     // Rotor geometry
Include "stator.geo";    // Stator geometry
// Main file assembles and meshes
```

## 9. Companion .geo for Display
```
// display.geo - open this in GMSH
Merge "coil.step";       // Geometry
Merge "results.msh";     // Field data
Mesh.NumSubEdges = 4;    // Curved elements
Mesh.VolumeEdges = 0;    // Clean display
View[0].IntervalsType = 2;
```
"""

GMSH_MAGNETS_INDUCTOR = """
# Magnet and Inductor Models (public-safe curated corpus)

## Permanent Magnet Force Computation (models/Magnets)
3D Maxwell stress tensor force between multiple magnets.

### Geometry Pattern
- N magnets: cylinders (revolve profile) or boxes (extrude)
- Air domain: auto-sized from BoundingBox + margin
- Volume with holes: `Volume(1) = {exterior_sl, mag1_sl, mag2_sl, ...}`
- Each magnet: Physical Volume + Physical Surface (skin)

### Key Technique: Auto-Sizing Air Domain
```
BoundingBox; // triggers computation
General.MinX, General.MaxX, ... // access extents
air_xmin = General.MinX - margin;
// Create air box from computed bounds
```

## EI-Core Inductor (models/Inductor)
2D/3D inductor with E-core and I-core.

### Features
- Supports 2D and 3D, OCC and built-in kernel
- Half/quarter/full symmetry via `Flag_Symmetry`
- Open/closed core (air gap) via `Flag_OpenCore`
- Shell transformation to infinity (`AIRINF` region)
- Parametric mesh density via single `md` parameter

### Symmetry Pattern
```
Flag_Symmetry = 2; // quarter model
SymmetryFactor = (Flag_Symmetry) ? 2*Flag_Symmetry : 1;
```

### Coil Construction (3D)
Combined extrusion for 3D coils with bends:
1. Extrude straight section
2. Rotational extrusion for 90-degree bend (Pi/2)
3. Extrude next straight section
4. Repeat for full coil path
```
out[] = Extrude {0,dl,0} { Surface{coil_face}; Layers{n}; };
out[] = Extrude { {0,0,1}, {xbend,ybend,0}, Pi/2 } { Surface{out[0]}; Layers{n}; };
```

### Physical Regions
```
ECORE  = 1000  // E-core iron
ICORE  = 1100  // I-core iron
COIL   = 2000  // Coil winding
AIR    = 3000  // Air domain
AIRINF = 3100  // Infinite shell
AIRGAP = 3200  // Air gap
```
"""


def get_gmsh_examples(topic: str = "all") -> str:
    """Return GMSH examples and tutorial documentation.

    Args:
        topic: tutorials, machines, patterns, magnets

    Returns:
        Documentation string.
    """
    topics = {
        "tutorials": GMSH_TUTORIALS,
        "machines": GMSH_ELECTRIC_MACHINES,
        "patterns": GMSH_EXAMPLE_PATTERNS,
        "magnets": GMSH_MAGNETS_INDUCTOR,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        available = ", ".join(topics.keys())
        return f"Unknown topic: '{topic}'. Available: all, {available}"
