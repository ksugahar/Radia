"""
Knowledge base for build123d MCP server.

Focused on CAE (Computer-Aided Engineering) workflows:
geometry creation for FEM/BEM meshing, not 3D printing.
"""

# ============================================================
# Overview
# ============================================================

OVERVIEW = """\
# build123d for CAE — Overview

build123d is a Python-native parametric CAD library built on the
Open Cascade (OCCT) kernel. It provides two programming styles:

- **Algebra mode**: stateless, operator-driven (`Box() - Cylinder()`)
- **Builder mode**: stateful, context-managed (`with BuildPart(): ...`)

## Why build123d for CAE

- Python-native: loops, functions, parametric studies
- OCCT kernel: industry-proven BREP geometry
- STEP/BREP export: direct pipeline to Netgen, Cubit, GMSH meshers
- Label/color attributes: map to material regions and boundary conditions
- OSS (Apache-2.0): no license restrictions

## CAE Pipeline

```
build123d (geometry) -> .step/.brep -> Mesher (Netgen/Cubit/GMSH) -> Solver (NGSolve/Radia)
```

## CAE-Safe Subset (IMPORTANT)

For CAE meshing, prefer **primitives + boolean** operations:
- Box, Cylinder, Cone, Sphere, Torus, Wedge + operator `-`, `+`, `&`
- extrude, revolve, sweep, loft from clean sketches

**Avoid** unless necessary:
- fillet/chamfer with very small radii (creates micro-edges)
- Complex Sketch constraints (over-constraint -> topology instability)
- Mixing boolean operations across different CAD kernels (ACIS vs OCCT)

The guiding principle: **Cubit-like simplicity in build123d syntax**.
Keep geometry clean and meshable.
"""

# ============================================================
# 3D Primitives
# ============================================================

PRIMITIVES_3D = """\
# 3D Primitives

All primitives support `rotation`, `align`, and `mode` parameters.
Default alignment is CENTER on all axes.

## Box
```python
Box(length, width, height)
Box(10, 20, 5)  # 10mm x 20mm x 5mm box centered at origin
```

## Cylinder
```python
Cylinder(radius, height, arc_size=360)
Cylinder(25, 50)  # r=25mm, h=50mm full cylinder
Cylinder(25, 50, arc_size=180)  # half cylinder
```

## Cone
```python
Cone(bottom_radius, top_radius, height, arc_size=360)
Cone(25, 10, 50)  # truncated cone
Cone(25, 0, 50)   # pointed cone
```

## Sphere
```python
Sphere(radius, arc_size1=-90, arc_size2=90, arc_size3=360)
Sphere(30)  # full sphere r=30mm
Sphere(30, arc_size1=0)  # hemisphere (upper half)
```

## Torus
```python
Torus(major_radius, minor_radius, minor_start_angle=0,
      minor_end_angle=360, major_angle=360)
Torus(50, 10)  # full torus, R=50mm, r=10mm
```

## Wedge
```python
Wedge(xsize, ysize, zsize, xmin, zmin, xmax, zmax)
# Near face: xsize x zsize, far face: (xmin..xmax) x (zmin..zmax)
```

## Boolean Operations (Algebra Mode)
```python
# Union (fuse)
result = Box(10, 10, 10) + Cylinder(3, 20)

# Subtraction (cut)
result = Box(10, 10, 10) - Cylinder(3, 20)

# Intersection
result = Box(10, 10, 10) & Cylinder(8, 20)
```

## Positioning
```python
from build123d import Pos, Rot

# Translate
part = Pos(5, 0, 0) * Cylinder(3, 10)

# Rotate (degrees)
part = Rot(0, 0, 45) * Box(10, 10, 5)

# Combine
result = Box(20, 20, 10) - Pos(5, 0, 0) * Cylinder(3, 20)
```
"""

# ============================================================
# 2D Sketch Primitives
# ============================================================

PRIMITIVES_2D = """\
# 2D Sketch Primitives

Sketch objects create 2D faces on a workplane. Used as input to
extrude, revolve, sweep, and loft.

## Circle
```python
Circle(radius)
Circle(25)  # circle r=25mm
```

## Rectangle
```python
Rectangle(width, height)
Rectangle(10, 20)  # 10mm x 20mm rectangle
```

## Ellipse
```python
Ellipse(x_radius, y_radius)
```

## RegularPolygon
```python
RegularPolygon(radius, side_count, major_radius=True)
RegularPolygon(10, 6)  # hexagon inscribed in r=10 circle
```

## Polygon (arbitrary vertices)
```python
Polygon((0, 0), (10, 0), (10, 5), (5, 10), (0, 10))
```

## Trapezoid
```python
Trapezoid(width, height, left_side_angle, right_side_angle=None)
Trapezoid(20, 10, 70)  # symmetric trapezoid
```

## RectangleRounded
```python
RectangleRounded(width, height, radius)
# Caution for CAE: small radius creates short edges
```

## Boolean on Sketches
```python
cross = Rectangle(2, 10) + Rectangle(10, 2)  # cross shape
ring = Circle(10) - Circle(8)  # annular ring
```
"""

# ============================================================
# Operations
# ============================================================

OPERATIONS = """\
# Operations (3D from 2D)

## extrude — Linear extrusion
```python
# Algebra mode
part = extrude(Rectangle(10, 20), amount=5)

# Builder mode
with BuildPart() as bp:
    with BuildSketch():
        Rectangle(10, 20)
    extrude(amount=5)
result = bp.part
```

## revolve — Axisymmetric rotation
```python
# Revolve rectangle around Y axis -> cylindrical shell
sketch = Rectangle(5, 20, align=(Align.MIN, Align.MIN))
part = revolve(Pos(10, 0) * sketch, axis=Axis.Y)

# Partial revolve
part = revolve(sketch, axis=Axis.Y, revolution_arc=180)
```

## sweep — Along a path
```python
# Sweep circle along a helix -> coil
path = Helix(pitch=10, height=50, radius=25)
part = sweep(Circle(2), path=path)
```

## loft — Between cross-sections
```python
with BuildPart() as bp:
    with BuildSketch(Plane.XY):
        Circle(10)
    with BuildSketch(Plane.XY.offset(20)):
        Rectangle(15, 15)
    loft()
```

## offset — Shell/thicken
```python
# Create hollow box (shell)
box = Box(10, 10, 10)
top_face = box.faces().sort_by(Axis.Z)[-1]
shelled = offset(box, amount=-1, openings=[top_face])
```

## thicken — Non-planar face to solid
```python
part = thicken(some_face, amount=2)
```

## mirror
```python
part = mirror(Box(10, 10, 5), about=Plane.YZ)
```

## split
```python
top_half = split(Sphere(10), bisect_by=Plane.XY, keep=Keep.TOP)
```

## section — Cross-section extraction
```python
cross = section(my_part, section_by=Plane.XY)
# Returns a Sketch (2D face)
```

## fillet / chamfer (use with care for CAE)
```python
box = Box(10, 10, 10)
# Fillet all edges
filleted = fillet(box.edges(), radius=1)

# Chamfer specific edges
top_edges = box.edges().filter_by(Axis.Z).sort_by(Axis.Z)[-4:]
chamfered = chamfer(top_edges, length=0.5)
```
**CAE warning**: fillet/chamfer with radius < 0.1 * characteristic_length
will create micro-edges that degrade mesh quality. Prefer defeaturing
(removing small features) over adding them.

## scale
```python
scaled = scale(my_part, by=2.0)           # uniform
scaled = scale(my_part, by=(1, 1, 0.5))   # non-uniform (changes geom type)
```
"""

# ============================================================
# Curves
# ============================================================

CURVES = """\
# Curve Objects

Curves (Edge/Wire) are used as sweep paths, sketch boundaries, etc.

## Line
```python
Line((0, 0), (10, 5))
```

## Polyline
```python
Polyline((0, 0), (10, 0), (10, 10), (0, 10), close=True)
```

## Spline
```python
Spline((0, 0), (5, 3), (10, 0), tangents=((0, 1), (0, -1)))
```

## Arc types
```python
CenterArc(center=(0, 0), radius=10, start_angle=0, arc_size=90)
RadiusArc(start_point=(0, 0), end_point=(10, 0), radius=15)
SagittaArc(start_point=(0, 0), end_point=(10, 0), sagitta=3)
ThreePointArc((0, 0), (5, 3), (10, 0))
TangentArc((0, 0), (10, 5), tangent=(1, 0))
```

## Helix (for coils)
```python
Helix(pitch=10, height=100, radius=25)
Helix(pitch=10, height=100, radius=25, cone_angle=5)  # conical helix
```

## Bezier
```python
Bezier((0, 0), (5, 10), (10, 0))  # quadratic
Bezier((0, 0), (3, 10), (7, 10), (10, 0))  # cubic
```

## BuildLine (Builder mode)
```python
with BuildLine() as bl:
    Line((0, 0), (10, 0))
    ThreePointArc((10, 0), (15, 5), (10, 10))
    Line((10, 10), (0, 10))
path = bl.line
```
"""

# ============================================================
# Export / Import
# ============================================================

EXPORT_IMPORT = """\
# Export / Import

## STEP (primary CAE interchange format)
```python
from build123d import export_step, import_step

# Export
export_step(my_part, "model.step")
export_step(my_part, "model.step", unit=Unit.MM)

# Import
compound = import_step("external_model.step")
```

## BREP (OCCT native, lossless for Netgen)
```python
from build123d import export_brep, import_brep

export_brep(my_part, "model.brep")
shape = import_brep("model.brep")
```
BREP is preferred when staying within the OCCT ecosystem
(build123d -> Netgen/NGSolve) because no STEP translation
tolerance issues arise.

## STL (triangulated surface, for visualization)
```python
from build123d import export_stl, import_stl

export_stl(my_part, "model.stl", tolerance=0.001)
face = import_stl("model.stl")  # returns Face (reference only)
```

## glTF (web visualization)
```python
from build123d import export_gltf

export_gltf(my_part, "model.glb", binary=True)
```

## SVG import (2D cross-sections)
```python
from build123d import import_svg

wires_and_faces = import_svg("cross_section.svg")
```

## CAE Pipeline Examples

### build123d -> Netgen (OSS, direct BREP)
```python
import netgen.occ as occ

export_brep(my_part, "model.brep")
geo = occ.OCCGeometry("model.brep")
mesh = geo.GenerateMesh(maxh=0.01)
```

### build123d -> Cubit (hex mesh via STEP)
```python
export_step(my_part, "model.step")
# In Cubit .jou:
# import step "model.step" heal
# volume all scheme tetmesh
# mesh volume all
```
"""

# ============================================================
# Topology Queries
# ============================================================

TOPOLOGY = """\
# Topology Queries and Selectors

## Basic queries
```python
part = Box(10, 20, 30)

part.faces()      # ShapeList of all 6 faces
part.edges()      # ShapeList of all 12 edges
part.vertices()   # ShapeList of all 8 vertices
part.wires()      # ShapeList of wires
part.solids()     # ShapeList of solids
part.shells()     # ShapeList of shells
```

## Properties
```python
part.volume        # volume of solid
part.area          # total surface area
face.area          # area of single face
edge.length        # length of edge
part.center()      # center of mass (CenterOf.MASS by default)
part.bounding_box()  # BoundBox object
```

## Filtering (ShapeList methods)
```python
# By axis direction
top_face = part.faces().sort_by(Axis.Z)[-1]     # highest Z face
bottom_face = part.faces().sort_by(Axis.Z)[0]   # lowest Z face
z_edges = part.edges().filter_by(Axis.Z)         # edges parallel to Z

# By geometry type
planar = part.faces().filter_by(GeomType.PLANE)
curved = part.faces().filter_by(GeomType.CYLINDER)

# By position
right_faces = part.faces().filter_by(
    lambda f: f.center().X > 5
)

# Group by
face_groups = part.faces().group_by(Axis.Z)  # group by Z position
```

## Labels and Colors (for CAE region assignment)
```python
# Assign labels
coil = Cylinder(25, 50)
coil.label = "coil"
coil.color = Color("copper")

core = Box(100, 100, 50)
core.label = "iron_core"

# Assembly with labeled parts
assembly = Compound(children=[coil, core])
assembly.label = "electromagnet"

# Export with labels (preserved in STEP)
export_step(assembly, "model.step")
```

## CAE Quality Checks
```python
# Minimum edge length (detect micro-edges)
min_edge = min(e.length for e in part.edges())
print(f"Minimum edge length: {min_edge}")

# Check validity
print(f"Valid: {part.is_valid}")

# Face count (complexity indicator)
print(f"Faces: {len(part.faces())}")

# Bounding box dimensions
bb = part.bounding_box()
print(f"Size: {bb.size}")
```
"""

# ============================================================
# CAE Guidelines
# ============================================================

CAE_GUIDELINES = """\
# CAE-Specific Guidelines for build123d

## 1. Prefer Primitives + Boolean

For CAE meshing (FEM/BEM), geometry must be **clean and meshable**.
The safest approach mirrors what Cubit users do naturally:

```python
# GOOD: primitives + boolean (clean topology)
workpiece = Cylinder(radius=25, height=50)
coil = Pos(0, 0, 30) * Cylinder(radius=30, height=10) - \\
       Pos(0, 0, 30) * Cylinder(radius=20, height=10)
model = workpiece + coil

# RISKY: complex sketch with constraints
# (may produce micro-edges or degenerate faces)
```

## 2. Minimum Edge Length Rule

Before exporting for meshing, check:
```python
min_edge = min(e.length for e in part.edges())
char_length = max(part.bounding_box().size)
ratio = min_edge / char_length

if ratio < 0.001:
    print(f"WARNING: micro-edge detected ({min_edge:.6f})")
    print("This may cause mesh quality issues")
```

Target: min_edge / characteristic_length > 0.01

## 3. BREP vs STEP for Netgen Pipeline

- **BREP** (export_brep): lossless within OCCT ecosystem.
  Use when meshing with Netgen directly.
- **STEP**: universal interchange. Use when meshing with Cubit.
  May introduce tolerance differences at import.

```python
# For Netgen: prefer BREP
export_brep(part, "model.brep")

# For Cubit: must use STEP
export_step(part, "model.step")
```

## 4. Label Convention for CAE Regions

Use labels to define material regions and boundary conditions:
```python
# Material regions (-> NGSolve/Cubit blocks)
iron.label = "iron"
air.label = "air"
coil.label = "coil"

# Boundary faces (-> NGSolve/Cubit sidesets)
# Access specific faces and label them
top = part.faces().sort_by(Axis.Z)[-1]
top.label = "symmetry_plane"
```

## 5. Axisymmetric 2D Models

For 2D axisymmetric problems, create the cross-section as a Sketch:
```python
# Cross-section in XZ plane (r-z plane)
with BuildSketch(Plane.XZ) as sk:
    with BuildLine():
        Polyline((10, 0), (25, 0), (25, 50), (10, 50), close=True)
    make_face()
cross_section = sk.sketch

# Export as BREP for 2D meshing
export_brep(cross_section, "cross_section.brep")
```

## 6. Assembly Structure for Multi-Region Models

```python
# Create individual bodies
coil = Cylinder(30, 10) - Cylinder(20, 10)
coil.label = "coil"
coil.color = Color("orange")

workpiece = Cylinder(25, 50)
workpiece.label = "workpiece"
workpiece.color = Color("gray")

air = Cylinder(100, 100) - coil - workpiece
air.label = "air"

# Combine as assembly
model = Compound(children=[coil, workpiece, air])
export_step(model, "ih_model.step")
```

## 7. Avoiding Dirty Geometry

Common pitfalls that create unmeshable geometry:
- **fillet(radius < 0.1 * min_dimension)**: creates micro-faces
- **boolean with tangent surfaces**: OCCT may succeed but produce
  degenerate edges at tangent contact
- **non-manifold results**: boolean of coincident faces
- **import from other kernels**: STEP from ACIS (Cubit) may have
  slightly different trimming curves

Mitigation:
```python
# After boolean operations, validate
result = box - cylinder
assert result.is_valid, "Boolean produced invalid geometry"
min_e = min(e.length for e in result.edges())
assert min_e > 1e-6, f"Micro-edge: {min_e}"
```
"""

# ============================================================
# Examples
# ============================================================

EXAMPLES = """\
# CAE Examples with build123d

## Example 1: Induction Heating (IH) Coil + Workpiece

```python
from build123d import *

# Workpiece (steel cylinder)
workpiece = Cylinder(radius=0.025, height=0.025)
workpiece.label = "workpiece"

# Coil (hollow cylinder, single turn)
coil = (
    Pos(0, 0, 0.005)
    * (Cylinder(radius=0.035, height=0.005)
       - Cylinder(radius=0.030, height=0.005))
)
coil.label = "coil"

# Air domain
air = Cylinder(radius=0.1, height=0.1) - workpiece - coil
air.label = "air"

# Export
model = Compound(children=[workpiece, coil, air])
export_step(model, "ih_model.step")
export_brep(model, "ih_model.brep")

# Quality check
for child in [workpiece, coil, air]:
    min_e = min(e.length for e in child.edges())
    print(f"{child.label}: {len(child.faces())} faces, "
          f"min edge = {min_e:.6f}")
```

## Example 2: Axisymmetric Cross-Section (2D)

```python
from build123d import *

# Define cross-section in r-z plane (Plane.XZ)
# r = x-axis, z = z-axis

with BuildSketch(Plane.XZ) as coil_cs:
    with BuildLine():
        Polyline(
            (0.030, 0.000),
            (0.035, 0.000),
            (0.035, 0.005),
            (0.030, 0.005),
            close=True,
        )
    make_face()
coil_face = coil_cs.sketch
coil_face.label = "coil"
export_brep(coil_face, "coil_2d.brep")
```

## Example 3: Parametric E-Core Transformer

```python
from build123d import *

def e_core(W=10, H=15, D=5, t=3):
    \"\"\"E-shaped transformer core.

    Args:
        W: total width
        H: total height
        D: depth (extrusion)
        t: leg/yoke thickness
    \"\"\"
    # E-shape cross-section
    with BuildSketch() as sk:
        # Bottom yoke
        Rectangle(W, t, align=(Align.CENTER, Align.MIN))
        # Left leg
        Pos(-W/2 + t/2, t) * Rectangle(t, H - t,
            align=(Align.CENTER, Align.MIN))
        # Center leg
        Pos(0, t) * Rectangle(t, H - t,
            align=(Align.CENTER, Align.MIN))
        # Right leg
        Pos(W/2 - t/2, t) * Rectangle(t, H - t,
            align=(Align.CENTER, Align.MIN))

    core = extrude(sk.sketch, amount=D)
    core.label = "iron_core"
    return core

# Parametric study
for t in [2, 3, 4]:
    core = e_core(t=t)
    export_step(core, f"e_core_t{t}.step")
    print(f"t={t}: volume={core.volume:.1f}")
```

## Example 4: Accelerator Dipole Yoke (Quarter Model)

```python
from build123d import *

# Quarter yoke with symmetry
R_bore = 30   # bore radius [mm]
R_outer = 150  # outer radius [mm]
L_half = 200   # half length [mm]

# Yoke cross-section (quarter, XZ plane)
with BuildSketch(Plane.XZ) as yoke_sk:
    # Quarter annulus
    with BuildLine():
        Line((R_bore, 0), (R_outer, 0))
        CenterArc((0, 0), R_outer, 0, 90)
        Line((0, R_outer), (0, R_bore))
        CenterArc((0, 0), R_bore, 90, -90)
    make_face()

yoke = extrude(yoke_sk.sketch, amount=L_half)
yoke.label = "iron_yoke"
export_step(yoke, "dipole_quarter_yoke.step")
```
"""

# ============================================================
# Topic lookup
# ============================================================

COIL_MODELING = """
# Coil Modeling for PEEC Filament Extraction

## Variable Cross-Section Coils (Tapered, Shaped)

Use `loft()` between rectangular cross-sections placed along a helix
path. Each section can have different (w, h), creating a smooth taper.

```python
from build123d import *
import math

R, pitch, n_turns = 50, 10, 5   # mm
w_bot, h_bot = 4, 4              # bottom cross-section
w_top, h_top = 2, 2              # top cross-section
n_sec = n_turns * 12 + 1         # 12 sections per turn

sections = []
for i in range(n_sec):
    t = i / (n_sec - 1)
    angle = 2 * math.pi * n_turns * t
    z = pitch * n_turns * t
    cx, cy = R * math.cos(angle), R * math.sin(angle)
    w = w_bot + (w_top - w_bot) * t
    h = h_bot + (h_top - h_bot) * t
    # tangent direction
    dx = -R * math.sin(angle) * 2 * math.pi * n_turns
    dy = R * math.cos(angle) * 2 * math.pi * n_turns
    dz = pitch * n_turns
    norm = math.sqrt(dx**2 + dy**2 + dz**2)
    tangent = (dx/norm, dy/norm, dz/norm)
    plane = Plane(origin=(cx, cy, z),
                  x_dir=(tangent[1], -tangent[0], 0),
                  z_dir=tangent)
    with BuildSketch(plane) as sk:
        Rectangle(w, h)
    sections.append(sk.sketch)

coil = loft(sections)
```

## MCP Tool Shortcut

Instead of writing the loft code manually, use:

    generate_helix_coil(radius=50, pitch=10, n_turns=5,
                        w_start=4, h_start=4, w_end=2, h_end=2)

This returns the solid + the helix path points for section_along_path.

## PEEC Filament Extraction Pipeline

1. **Generate solid**: `generate_helix_coil(...)` or custom `execute_build123d(...)`
2. **Export STEP**: set `export_dir` to save the solid
3. **Section along path**: `section_along_path(step_file, path_json)` returns
   per-segment (area, w_est, h_est) in CAD units
4. **Build filaments**: feed (w, h) into `PEECBuilder.add_connected_segment()`
   with per-segment local dimensions (convert CAD units to meters)
5. **Solve**: `PEECCircuitSolver(topo, use_hacapk=True, outer_method="saddle")`

Python-side helper: `from radia.coil_from_cad import filaments_from_step`
wraps steps 3-5 in one call.

## Best Practices

- **sections_per_turn >= 8** for smooth loft (12 recommended)
- **Label the solid** (`coil.label = "coil"`) for STEP export filename
- **Square cross-section**: w_est = h_est = sqrt(area)
- **Rectangular**: pass aspect_ratio to section_along_path (TODO)
- **Multi-turn proximity**: section_along_path uses centroid-proximity
  filtering so one cut plane hitting multiple turns picks the right face
- **nwinc/nhinc >= 3** for skin/proximity at f > 10 kHz in PEEC
- **Units**: build123d works in mm, PEEC in meters. scale = 1000.
"""

_TOPICS = {
    "overview": OVERVIEW,
    "primitives_3d": PRIMITIVES_3D,
    "primitives_2d": PRIMITIVES_2D,
    "operations": OPERATIONS,
    "curves": CURVES,
    "export_import": EXPORT_IMPORT,
    "topology": TOPOLOGY,
    "cae_guidelines": CAE_GUIDELINES,
    "examples": EXAMPLES,
    "coil_modeling": COIL_MODELING,
}


def get_build123d_documentation(topic: str = "overview") -> str:
    """Return knowledge base content for requested topic."""
    key = topic.strip().lower()
    if key == "all":
        return "\n\n".join(_TOPICS.values())
    if key in _TOPICS:
        return _TOPICS[key]
    available = ", ".join(sorted(_TOPICS.keys()))
    return f"Unknown topic '{topic}'. Available: {available}, all"
