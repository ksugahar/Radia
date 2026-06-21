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

## Footgun: `is_valid` is an attribute, not a method

In build123d 0.10+, `Shape.is_valid` is a **bool attribute**, not a
method. Older training-data patterns that call `part.is_valid()` will
raise `TypeError: 'bool' object is not callable`. Always write:
```python
assert part.is_valid           # correct
# NOT: assert part.is_valid()  # wrong in 0.10+
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

COIL_FROM_POLYLINE = """
# Coil Construction from a Centerline Polyline (build123d loft)

## When to Use

The centerline is given as **discrete 3D points** `r[0..N-1]` (measured
wire path, optimization output, FEA tool export, scan data) — NOT an
analytical helix and NOT extracted from existing CAD.

Counterparts:
- `coil_modeling` (this file) — analytical helix with `loft()` over
  parametric sections.
- `radia_mcp.cubit` `coil_from_polyline` topic — same algorithm, Cubit
  flavour (creates one body per station, lofts pairwise).
- `radia.coil_from_cad` — centerline EXTRACTION from existing STEP CAD.

This topic is the **forward** problem: discrete points → solid CAD.

## Algorithm (frame at each station + global loft)

```
Input:  r[0], r[1], ..., r[N-1]   (ordered 3D centerline points, mm)
        R                         (cross-section radius / dims)

For each interior point n = 1 .. N-2:
    # 1. Tangent at r[n] via centered finite difference
    w = (r[n+1] - r[n-1]) / |r[n+1] - r[n-1]|

    # 2. "Up" vector — smallest-component trick to avoid up // w
    k  = argmin |w[i]|     (the tangent's smallest absolute component)
    up = e_k

    # 3. Orthonormal frame at r[n]
    u = normalize(w x up)        # in-plane x-axis
    v = w x u                    # in-plane y-axis (= z_dir x x_dir)

    # 4. build123d Plane: x_dir=u, z_dir=w => y_dir is auto = w x u = v
    plane = Plane(origin=r[n], x_dir=u, z_dir=w)
    section_n = Circle(R)  applied on plane     # or Rectangle, etc.

# 5. Single global loft over all sections in order
coil = loft(sections)                           # build123d Solid
```

`Plane(origin, x_dir, z_dir)` automatically derives the in-plane y-axis
as `z_dir x x_dir`. With our conventions this gives **exactly** the
columns of `M = [u | v | w]` from the Cubit flavour — so the two
flavours produce geometrically identical cross-section placements.

## Why the Smallest-Component Up Vector

`up = e_k` with `k = argmin|w_i|` GUARANTEES `up` is not parallel to
`w`, so `w x up` is well-conditioned. A naive `up = (0, 0, 1)` fails
when the centerline is vertical.

It is **NOT** a parallel-transport / rotation-minimizing frame (RMF):

| Frame                   | Twist between stations | Cost                |
|-------------------------|------------------------|---------------------|
| Smallest-component up   | uncontrolled (random)  | O(1) per station    |
| Parallel transport      | minimal (Bishop / RMF) | O(N) sequential     |
| Frenet (binormal)       | ill-defined at straights | requires curvature |

For **circular** cross-section the twist is invisible (rotational
symmetry), so this cheap frame is fine. For **rectangle / racetrack /
ellipse**, the twist between adjacent stations may produce a visibly
twisted lofted solid — switch to RMF (double-reflection / Bishop).

## Reference Implementation (algebra API, recommended)

```python
import numpy as np
from build123d import Plane, Circle, loft, export_step

# r is a list of np.array([x, y, z]) centerline points (mm), length N
R_cross = 1.0   # cross-section radius (mm)

sections = []
for n in range(1, len(r) - 1):
    # 1. tangent
    w = r[n + 1] - r[n - 1]
    w = w / np.linalg.norm(w)

    # 2. smallest-component up
    up = np.zeros(3)
    up[int(np.argmin(np.abs(w)))] = 1.0

    # 3. orthonormal frame
    u = np.cross(w, up); u /= np.linalg.norm(u)
    # v = np.cross(w, u)  # implied by Plane (z_dir x x_dir)

    # 4. build123d plane and section
    plane = Plane(origin=tuple(r[n]), x_dir=tuple(u), z_dir=tuple(w))
    sections.append(plane * Circle(R_cross))   # algebra: place sketch on plane

# 5. global loft
coil = loft(sections)
coil.label = "coil"
export_step(coil, "coil_from_polyline.step")
```

`plane * Circle(R)` is the algebra-API idiom: the `Plane.__mul__`
operator places the sketch on that plane (equivalent to
`Plane.from_local_coords(sketch)`).

## Builder API Variant

```python
from build123d import BuildPart, BuildSketch, Plane, Circle, loft

with BuildPart() as part:
    sections = []
    for n in range(1, len(r) - 1):
        w  = r[n + 1] - r[n - 1]; w  /= np.linalg.norm(w)
        up = np.zeros(3); up[int(np.argmin(np.abs(w)))] = 1.0
        u  = np.cross(w, up); u  /= np.linalg.norm(u)
        plane = Plane(origin=tuple(r[n]), x_dir=tuple(u), z_dir=tuple(w))
        with BuildSketch(plane) as sk:
            Circle(R_cross)
        sections.append(sk.sketch)
    loft(sections)         # appended to BuildPart context

coil = part.part
```

## Pairwise (Chain) Loft Variant — for Wiggly Polylines

`loft([s_0, s_1, ..., s_{N-2}])` is a SINGLE OCC global loft. For
densely sampled, wiggly polylines this can fail with self-intersection
or surface tolerance errors. The Cubit flavour avoids this by lofting
**pairwise** (one volume per adjacent pair, then `imprint+merge+unite`).
The build123d analog is:

```python
from build123d import Compound

pieces = []
for n in range(len(sections) - 1):
    pieces.append(loft([sections[n], sections[n + 1]]))
coil = Compound(children=pieces)
# Optional fuse: from build123d import Part; coil = Part() + pieces
```

`Compound` keeps each segment as a separate solid (good if you want
to tag per-segment labels for PEEC). `Part() + pieces` fuses them into
one body with shared internal faces removed (analogous to Cubit's
`unite volume all` after `imprint+merge+compress`).

## Caveats and Failure Modes

1. **Endpoints are dropped.** The centered finite difference needs
   both neighbours, so the first and last input points are NOT used as
   stations — the resulting solid is shorter than the input polyline.
   To preserve length, prepend / append ghost points or use one-sided
   differences at the ends.

2. **Closed loops require a wrap-around section.** For a closed coil
   `r[N-1] approx r[0]`, append `sections[0]` to the end of the list
   before the final `loft(...)` so the global loft wraps around — and
   compute the centered tangent at the seam with periodic neighbours.

3. **Twist visible for non-circular cross-section.** The
   smallest-component frame can flip discretely when the dominant
   tangent component changes between adjacent stations (e.g. tangent
   rotates from +x-dominant to +y-dominant, k jumps from 1 to 2).
   Adjacent sections then differ by a 90 deg roll → loft self-intersects
   or twists. Fix: parallel-transport (Bishop / double-reflection) frame.

4. **Polyline density matters.** Too sparse → faceted polygonal solid;
   too dense → many input sections (slow loft, may fail OCC tolerance).
   For wire paths from physical measurement, smooth and resample to a
   uniform arc-length spacing before invoking this algorithm.

5. **Units.** build123d default unit is mm. Make sure `r` is in mm
   (or scale before/after). For PEEC the downstream pipeline scales
   by 1/1000 to meters when calling `PEECBuilder.add_connected_segment`.

6. **Global vs pairwise loft choice.** Default to `loft(sections)`
   (smoother solid, single body); fall back to pairwise + Compound if
   OCC raises self-intersection or surface tolerance errors on a
   wiggly polyline.

## Comparison to Existing Coil Topics in This Knowledge Base

| Method                    | Input                  | Best For                        |
|---------------------------|------------------------|----------------------------------|
| `coil_modeling`           | analytical helix       | Regular cylindrical solenoid     |
| `coil_from_polyline` (this) | discrete 3D points    | Measured / optimized free-form  |
| `radia.coil_from_cad`     | existing STEP CAD      | Reverse-engineering filament from CAD |
| `radia.coil_builder`      | analytical primitives  | Beam-optics design (straight + arc) |

Pipeline downstream is identical: STEP export →
`coil_from_cad.extract_centerline_from_step` → filaments →
`PEECBuilder` → `PEECCircuitSolver`.

## Cubit Flavour — Same Algorithm, Different Engine

Identical math; different engine. See the `coil_from_polyline` topic
in `radia_mcp.cubit` for the Cubit version. Key engine differences:

| Aspect              | build123d                    | Cubit                       |
|---------------------|------------------------------|-----------------------------|
| Frame placement     | `Plane(origin, x_dir, z_dir)` | `rotate Surface ... about origin` then `move` |
| Loft strategy       | `loft([s0..sM])` (global)     | `create volume loft surface n n+1` (pairwise) |
| Cleanup             | `Part() + pieces` to fuse     | `imprint all; merge all; compress; unite volume all` |
| Section primitive   | `Circle(R)`, `Rectangle(w,h)` | `create surface circle ... zplane` |
| STEP export         | `export_step(coil, ...)`      | `export step "..." overwrite` |
| Use it when ...     | Scripted, no Cubit license needed | Already in a Cubit meshing pipeline |

## Downstream Consumers (after STEP export)

Once the lofted coil is exported to STEP, it feeds two production
paths in this codebase. The path FROM the build123d solid is identical
to the Cubit case — engine choice does not affect downstream usage.

### A. Radia PEEC (filament-based circuit extraction)

```python
from build123d import export_step
export_step(coil, "coil.step")

from radia.coil_from_cad import extract_centerline_from_step
from radia.peec_topology import PEECBuilder, PEECCircuitSolver

path_m, w_m, h_m = extract_centerline_from_step("coil.step")  # mm -> m
b = PEECBuilder()
prev = b.add_node_at(*path_m[0])
for p, w, h in zip(path_m[1:], w_m[1:], h_m[1:]):
    nxt = b.add_node_at(*p)
    b.add_connected_segment(prev, nxt, w, h, sigma=5.8e7,
                            nwinc=3, nhinc=3)
    prev = nxt
b.add_port(b.nodes[0], b.nodes[-1])
solver = PEECCircuitSolver(b.build_topology())
Z = solver.compute_port_impedance(freq=1e6)
```

`extract_centerline_from_step` 5-Path dispatch handles the topology
produced by this algorithm (single lofted tube, possibly closed).

### B. BEM-A (port-driven self-inductance via ngsolve.bem)

BEM-A consumes a **meshed** coil with `source` / `sink` BND labels on
the two terminal faces. build123d alone cannot mesh — pipe the STEP
through Cubit (preferred for hex) or Netgen (tet):

```python
# Option 1: Cubit hex mesh
#   import "coil.step"; mesh volume all
#   sideset 1 add surface <source_id>; sideset 1 name "source"
#   sideset 2 add surface <sink_id>;   sideset 2 name "sink"
#   export netgen "coil.vol" overwrite

# Option 2: Netgen tet mesh from build123d directly
from build123d import export_step
from ngsolve import Mesh
from netgen.occ import OCCGeometry
# Tag the two end faces in build123d BEFORE export:
coil.faces().sort_by(Axis.Z)[0].label = "source"   # bottom cap
coil.faces().sort_by(Axis.Z)[-1].label = "sink"    # top cap
export_step(coil, "coil.step")
geo = OCCGeometry("coil.step")
mesh = Mesh(geo.GenerateMesh(maxh=2.0))

from radia.bem.coil_inductance_ngsolve import compute_inductance_source_sink
res = compute_inductance_source_sink(
    mesh, source_label="source", sink_label="sink",
    Z_s_re=1.0/(5.8e7 * 9.3e-5),  # optional AC SIBC at 50 kHz Cu
)
print("L =", res["L"], "H,  R =", res["R"], "Ω")
```

**Critical: BEM-A requires terminal faces.**
- *Open-ended path* (`r[0] != r[N-1]`): the two end-cap faces of the
  lofted solid are natural source/sink. Tag them before export.
- *Closed loop* (`r[N-1] ≈ r[0]`): the solid has no end faces. Either
  (a) leave a small gap in the polyline (skip the wrap station) so two
  caps appear, or (b) cut a thin slice perpendicular to the path
  post-loft and tag the two new faces. PEEC has no equivalent
  constraint — it can drive a closed loop via a port between any two
  interior nodes.

| Pipeline | Input from STEP                   | Output                | Use when                   |
|----------|-----------------------------------|-----------------------|----------------------------|
| PEEC     | centerline + (w, h)               | R + jωL, SPICE / Z(f) | DC-1 MHz, multi-port, MOR  |
| BEM-A    | surface mesh + source/sink labels | L (+ optional AC R)   | high-fidelity self-L on free-form paths |
"""

EXAMPLES_INTRO = """\
# build123d Introductory Examples (36 walk-through examples)

Source: build123d `docs/introductory_examples.rst`. Each example is shown in
both **Builder mode** (context-manager based) and **Algebra mode** (operator-
based). For CAE work, either mode is acceptable — algebra mode is often more
concise.

## Ex 1 — Simple Rectangular Plate (Box)
```python
ex1 = Box(100, 100, 10)
```

## Ex 2 — Plate with Hole (Subtract)
```python
ex2 = Box(100, 100, 10) - Cylinder(radius=20, height=10)
```

## Ex 3 — Extruded Prismatic Solid (2D sketch + extrude)
```python
ex3_sk = Circle(40) - Rectangle(20, 20)
ex3 = extrude(ex3_sk, amount=10)
```

## Ex 4 — Profile from Lines/Arcs + make_face + extrude
```python
l1 = Line((0,0), (100,0)); l2 = Line((100,0), (100,100))
l3 = Line((100,100), (0,100)); l4 = Line((0,100), (0,0))
profile = make_face(Curve() + [l1, l2, l3, l4])
ex4 = extrude(profile, amount=10)
```

## Ex 5 — Moving the Working Point (Pos)
```python
box = Box(20, 20, 20)
ex5 = Pos(0,0,0)*box + Pos(50,0,0)*box
```

## Ex 6 — Point Lists (Locations)
```python
cyl = Cylinder(radius=10, height=20)
ex6 = Curve() + [Pos(p)*cyl for p in [(0,0,0),(50,0,0),(100,0,0)]]
```

## Ex 7 — Polygons at multiple locations
```python
poly = RegularPolygon(radius=20, side_count=6)
sketch = Sketch() + [Pos(x,0)*poly for x in (0, 50, 100)]
ex7 = extrude(sketch, amount=10)
```

## Ex 8 — Polylines + mirror
```python
poly = Polyline((0,0),(20,0),(20,10),(40,10),(40,20),(0,20))
profile = make_face(poly)
ex8 = extrude(profile + mirror(profile, about=Plane.YZ), amount=10)
```

## Ex 9 — Selectors + fillet
```python
box = Box(100, 100, 20)
top_edges = box.edges().group_by(Axis.Z)[-1]
ex9 = fillet(top_edges, radius=5)
```

## Ex 10 — Select LAST + Hole feature
```python
box = Box(100, 100, 20)
filled = fillet(box.edges(), radius=5)
ex10 = filled - Hole(radius=10, depth=20)
```

## Ex 11 — Face as plane + GridLocations
```python
box = Box(100, 100, 20)
topf = box.faces().sort_by(Axis.Z)[-1]
grid = GridLocations(x_spacing=30, y_spacing=30, x_count=2, y_count=2)
pent = RegularPolygon(radius=10, side_count=5)
ex11 = box - extrude(topf * [loc*pent for loc in grid], amount=-5)
```

## Ex 12 — Spline curve
```python
pts = [(0,0),(10,20),(20,10),(30,30)]
spline = Spline(*pts)       # -> used inside BuildLine
```

## Ex 13 — CounterBoreHole / CounterSinkHole + PolarLocations
```python
box = Box(100, 100, 20)
locs = PolarLocations(radius=40, count=4)
ex13 = box - [loc * CounterBoreHole(pilot_diameter=8, bore_diameter=12, bore_depth=10) for loc in locs]
```

## Ex 14 — '@'/'%' operators + sweep
```python
l1 = Line((0,0),(100,0)); l2 = Line((100,0),(100,100))
ex14 = sweep(Circle(radius=5), path=l1 + l2)
```

## Ex 15 — Mirror 2D geometry
```python
base = Curve() + [l1, l2, l3]
curve = make_face(base + mirror(base, about=Plane.YZ))
ex15 = extrude(curve, amount=10)
```

## Ex 16 — Mirror 3D object
```python
box = Box(100,100,20)
ex16 = box + mirror(box, about=Plane.XY.offset(10))
```

## Ex 17 — Mirror from face
```python
box = Box(100,100,20)
mir = Plane(box.faces().sort_by(Axis.Y)[-1])
ex17 = box + mirror(box, about=mir)
```

## Ex 18 — Workplane on face + subtract-extrude
```python
box = Box(100,100,20)
topf = box.faces().sort_by(Axis.Z)[-1]
ex18 = box - extrude(topf * Rectangle(50,50), amount=-10)
```

## Ex 19 — Workplane on vertex
```python
box = Box(100,100,20)
topf = box.faces().sort_by(Axis.Z)[-1]
vtx = topf.vertices().group_by(Axis.X)[-1]
ex19 = box - [Pos(v.X, v.Y, 0) * Cylinder(radius=10, height=20) for v in vtx]
```

## Ex 20 — Offset sketch workplane
```python
box = Box(100,100,20)
plane = Plane(box.faces().sort_by(Axis.X)[0]).offset(10)
ex20 = box + extrude(plane * Circle(radius=20), amount=10)
```

## Ex 21 — Plane in center of cylinder
```python
cyl = Cylinder(radius=20, height=100)
plane = Plane(origin=cyl.origin + (0,0,50), z_dir=cyl.z_dir)
ex21 = cyl + extrude(plane * Circle(radius=10), amount=50)
```

## Ex 22 — Rotated workplane + GridLocations
```python
box = Box(100,100,20)
plane = Plane.XY.rotated(Axis.X, 45)
grid = GridLocations(x_spacing=30, y_spacing=30, x_count=2, y_count=2)
ex22 = box + extrude(plane * [loc*Circle(radius=5) for loc in grid], amount=10)
```

## Ex 23 — Revolve + split
```python
poly = Polyline((0,0),(20,0),(20,20),(0,20))
profile = split(make_face(poly - Rectangle(10,10)), Plane.ZY, keep=Keep.FORWARD)
ex23 = revolve(profile, axis=Axis.Z)
```

## Ex 24 — Loft between circle and rectangle
```python
circle = Circle(radius=30)
rect = Pos(0,0,50) * Rectangle(40, 40)
ex24 = loft([circle, rect])
```

## Ex 25 — Offset sketch (2D)
```python
rect = Rectangle(100, 100)
ex25 = extrude(offset(rect, amount=10), amount=10)
```

## Ex 26 — 3D offset / shell
```python
box = Box(100,100,100)
ex26 = offset(box, amount=5)
```

## Ex 27 — Split object
```python
box = Box(100,100,20)
plane = Plane(box.faces().sort_by(Axis.X)[0]).offset(50)
ex27 = split(box, plane, keep=Keep.FORWARD)
```

## Ex 28 — Locate features using faces
```python
prism = extrude(RegularPolygon(radius=30, side_count=3), amount=20)
sphere = Sphere(radius=50)
ex28 = sphere - [Plane(f) * Cylinder(radius=10, height=30) for f in prism.faces()]
```

## Ex 29 — Classic OCC Bottle
```python
profile = make_face(Spline((0,0),(10,5),(10,30),(5,40),(5,70)))
bottle = revolve(profile, axis=Axis.Z)
ex29 = offset(bottle, amount=-2)
```

## Ex 30 — Bezier curve with weights
```python
pts = [(0,0),(20,10),(40,0),(50,20)]; wts = [1,1,1,1]
curve = Polyline(*pts) + Bezier(*pts, weights=wts)
ex30 = extrude(make_face(curve), amount=10)
```

## Ex 31 — Nested locations (PolarLocations)
```python
hex_solid = extrude(RegularPolygon(radius=10, side_count=6), amount=10)
ex31 = [loc * hex_solid for loc in PolarLocations(radius=50, count=6)]
```

## Ex 32 — Python for-loop over GridLocations
```python
circles = Circle(radius=10)
grid = GridLocations(x_spacing=30, y_spacing=30, x_count=2, y_count=2)
ex32 = [extrude(loc*circles, amount=(i+1)*5) for i, loc in enumerate(grid)]
```

## Ex 33 — Function + for-loop
```python
def make_square(size): return Rectangle(size, size)
ex33 = [Pos(0,0,i*10) * extrude(make_square(20+i*10), amount=10) for i in range(5)]
```

## Ex 34 — Embossed / debossed text
```python
box = Box(100,100,20)
topf = box.faces().sort_by(Axis.Z)[-1]
ex34 = box + extrude(topf * Text("Hello", font_size=20), amount=5) \\
           - extrude(topf * Text("World", font_size=20), amount=5)
```

## Ex 35 — Slots (SlotCenterToCenter, SlotArc)
```python
slot = SlotCenterToCenter(center_point_distance=50, width=10, height=20)
arc = RadiusArc((0,0),(10,10), radius=15)
ex35 = extrude(slot + SlotArc(arc=arc, width=10), amount=10)
```

## Ex 36 — extrude until face (Until.LAST)
```python
box = Box(100,100,50)
ex36 = extrude(Circle(radius=20), amount=Until.LAST)
```

## Key primitives/operations index

| Category | Names |
|---|---|
| 3D primitives | Box, Cylinder, Sphere, Cone, Torus |
| 2D primitives | Rectangle, Circle, RegularPolygon, SlotCenterToCenter, SlotArc, Text |
| Curves | Line, Polyline, Spline, Bezier, RadiusArc, CenterArc |
| Operations | extrude, revolve, sweep, loft, offset, split, fillet, chamfer, mirror |
| Placement | Pos, Plane, Locations, GridLocations, PolarLocations, HexLocations |
| Selectors | faces().sort_by(Axis.Z), edges().group_by(Axis.X), Select.LAST |
| Modes | Mode.ADD, Mode.SUBTRACT, Mode.INTERSECT, Mode.PRIVATE |
| Termination | Until.NEXT, Until.LAST, Until.FIRST |
"""

EXAMPLES_GALLERY = """\
# build123d Gallery Examples (21 complete examples)

Source: build123d `docs/examples_1.rst`. Real-world parametric designs
showing idiomatic use of build123d. For each: title, what it demonstrates,
and key primitives/operations.

## 1. Benchy — STL import + BREP replacement
- Imports STL model as Solid via Mesher class, modifies by replacing chimney
  with BREP version. Shows STL/BREP interop.
- Key: Mesher class, STL import, model editing

## 2. Bicycle Tire — wrap_faces for 2D->3D projection
- Realistic tire with patterned tread by projecting 2D planar geometry onto
  curved 3D surface.
- Key: wrap_faces, pattern projection

## 3. Bracelet — Sweep + Gordon surfaces + text projection
- Doubly-curved bracelet with embossed label. OCCT stress test with freeform
  surfaces, swept sections, Gordon surfaces.
- Key: sweep, Gordon surfaces, curve/text projection, location_at

## 4. Canadian Flag in Wind — Face projection to non-planar
- Projects planar faces onto non-planar face.
- Key: Face projection, complex line snapping

## 5. Cast Bearing Unit — Draft operations
- Castable flanged bearing housing with draft angles for mold release.
- Key: draft, flanged geometry, mold design

## 6. Circuit Board With Holes — Locations / product / range
- Placing holes around a part.
- Key: Locations, product, range

## 7. Clock Face — PolarLocations for radial features
- 3D clock face with minute indicators, arcs, lines, fillets, hour labels, slots.
- Key: PolarLocations, arcs, fillets, extrude

## 8. Fast Grid Holes — Face-with-holes pattern (625 holes)
- Efficient approach: 2D Face construction with hole wires, then extrude in
  a single operation (vs 3D boolean subtraction).
- Key: HexLocations, Face constructor with hole wires, perf optimization

## 9. Handle — Multisection sweep
- Drawer handle by sweeping through multiple profiles.
- Key: multisection sweep, complex profiles

## 10. Heat Exchanger — HexLocations + fillet
- Parametric heat exchanger core, tube positions from HexLocations,
  circular end caps, filleted tube ends.
- Key: HexLocations, parametric design, fillet, circular constraints

## 11. Key Cap — extrude with taper + extrude_until_next
- Cherry MX key cap.
- Key: extrude(..., taper=...), extrude until next

## 12. Build123d Logo — Text and lines
- Former build123d logo; builder mode generates SVG output.
- Key: Text, lines, SVG export

## 13. Maker Coin — DoubleTangentArc + text emboss
- Smooth transitions via DoubleTangentArc, embossed text by projection
  (even-depth text).
- Key: DoubleTangentArc, text projection, embossing

## 14. Multi-Sketch Loft — loft + type selection + shell
- Key: loft, sort_by type, shell

## 15. Peg Board J Hook — Complex path + flattening
- J-shaped pegboard hook: sweep complex path, flatten sides for 3D printing.
- Key: sweep, complex paths, flattening

## 16. Platonic Solids — Custom Part object (PlatonicSolid)
- Five highly symmetrical 3D shapes with congruent regular polygon faces.
- Key: custom Part class, symmetrical geometry, polyhedra

## 17. Playing Cards — Custom Sketch objects + Bezier
- Club, Spade, Heart, Diamond, PlayingCard; two-part box with suit cutouts.
- Key: Bezier curves, custom Sketch objects, SVG import

## 18. Stud Wall — Custom Part + RigidJoints assembly
- Custom Part objects, assemblies with snap-point RigidJoints.
- Key: custom Part, assemblies, RigidJoints, object copy

## 19. Tea Cup — Revolve + sweep + offset + fillet
- Complex, non-flat geometry: revolve body, sweep handle, shell, fillet.
- Key: revolve, sweep (2x), offset/shell, fillet

## 20. Toy Truck — Sketch reuse + symmetry + tapered extrusions
- Detailed body, cab, grill, bumper with joints for assembly.
- Key: sketch reuse, symmetry, taper, selective fillet, joints

## 21. Vase — Revolve + shell + selective fillet
- Revolving, shelling, edge fillet by position and type.
- Key: revolve, shell, edge selection by position/type, fillet

## Distilled patterns for CAE (applicable ones)

| Pattern from gallery | CAE relevance |
|---|---|
| Fast Grid Holes (Ex8): Face-with-holes + single extrude | Efficient boolean for perforated plates / cooling fins |
| Heat Exchanger (Ex10): HexLocations + fillet | IH coil arrays, hex-packed tube fields |
| Multi-Sketch Loft (Ex14): Loft + shell | Variable cross-section coils, taper windings |
| Clock Face (Ex7): PolarLocations | Stator tooth arrangements, Halbach segment placement |
| Tea Cup (Ex19): Revolve + shell | Axisymmetric coil bobbins, pressure vessels |
| Key Cap (Ex11): extrude with taper | Drafted pole pieces, wedged magnets |
| Vase (Ex21): Edge selection by position/type | Selective fillet on magnet poles / armature teeth |
"""

LAB_POLICY = """\
# Sugawara Lab (菅原研) CAD / meshing stance (2026-04-20, reaffirmed)

**Official lab position: build123d + Cubit are the primary pair.**
These two cover the vast majority of lab work — build123d for CAD
authoring (AI-writable, Python-native, OCCT), Cubit for
industrial-grade hex meshing and visualization.

## Stance toward neighbouring CAD tools

| Tool | Lab stance |
|---|---|
| **build123d** | **主力 (push)** — 新規 lab work はこれで執筆 |
| **Cubit** | **主力 (push)** — hex mesh + visualization の正面 |
| **FreeCAD** | **応援 (friendly)** — 菅原研は FreeCAD コミュニティを支持する。lab では主力ではないが、assembly/設計共有は FreeCAD が強いシーンがあるので interop を大切にする。`freecad_to_cubit_hex` 経由で pipeline に乗る |
| CadQuery | interop / compat — OCCT sibling、既存資産の移行 |
| OpenSCAD | interop / compat — 既存スクリプト対応 |

FreeCAD は「前面には出さないが味方」のスタンス。コミュニティを応援
し、貶める気は無い。ただし**研究室の運用 (lab-internal workflow) は
build123d**であり、新規モデリング・反復・パラメータスイープ・CAE
パイプラインの全てで build123d を使う。FreeCAD を lab 内部で使うこ
とはない。`freecad_exec_safely` は外部 / community 向け interop。

## Lab iteration pattern (build123d 派)

build123d は stateless (毎 `execute_build123d` が fresh Python) なので、
Cubit/FreeCAD 流の `*_exec_safely` (live state を checkpoint→ batch
dry-run → commit) は不要。代わりに lab は:

- **`build123d_try`** — 1 変種を subprocess 隔離試行
- **`build123d_try_race(scripts, max_concurrent)`** — N 変種を並列
  試行、`valid_largest_volume` / `first_valid` / `first_ok` 規則で
  winner 自動選定。fillet 半径や寸法のパラメータスイープに最適
- **`build123d_to_cubit_hex(script)`** — winner を STEP → Cubit hex
  まで一発で

source は git で管理 → revert は `git restore`、安全網は VCS 任せ。
build123d パイプラインに `*_exec_safely` 不要、というのが lab 設計。

Use this topic when deciding which tool to use or when translating
between them.

## Role split

| Tool | Role | When to use |
|---|---|---|
| **build123d** (Python / OCCT) | **tet path, Python-driven** | Default CAD. VSCode-native authoring. build123d → Netgen (tet) → Radia / Gmsh pipeline. AI-writable. |
| **Cubit** (`.jou`, `cubit.cmd(...)`) | **hex path + CAD fallback** | (1) Radia/ELF requires hex — **only Cubit can supply**. (2) Shapes build123d cannot express cleanly (netgen.occ has limits). (3) Legacy `.jou` assets. |
| (not used) | | FreeCAD (reproducibility), raw OCCT Python (verbose) |

**This is a complementary split, not an either/or.** Both tools live in
the lab toolchain permanently. Cubit is **not** being deprecated.

## When you are asked to "translate a .jou to build123d"

- **tet-only geometry** (boolean primitives, sweep, revolve, fillet):
  translate directly — see `cubit_rosetta` topic for the mapping table.
- **Cubit-only primitives** (`imprint all`, `merge all`, hex-sweep
  schemes, `block`, `sideset`, `nodeset`): do **not** invent build123d
  equivalents — these concepts don't exist here. Flag them and keep
  the `.jou` for the hex path instead.
- **Mesh directives** (`mesh volume N`, `scheme tetmesh`, size
  settings): in the build123d pipeline, meshing is **Netgen's job**.
  Drop these commands; pass `maxh` to `OCCGeometry.GenerateMesh()`.
- **Block/sideset labels**: in the build123d pipeline, use
  `part.label = "name"` + `Compound(children=[...])` + the
  `run_pipeline_multi` helper in `examples/build123d_netgen_gmsh_flow/`
  which bridges labels to Gmsh physical groups.

## Related topics
- `cubit_rosetta` — verb-by-verb mapping Cubit ↔ build123d
- `examples_intro` — 36 build123d patterns from upstream docs
- `cae_guidelines` — clean-geometry rules for meshing

The companion mcp-server-cubit has the matching
`scripting_lab_policy` / `scripting_build123d_crossref` topics.
"""

CUBIT_ROSETTA = """\
# Cubit ↔ build123d Rosetta Stone (CAD subset)

Side-by-side mapping for the most common Cubit CAD verbs and their
build123d equivalents.  Use this when porting `.jou` scripts or when
writing build123d by a Cubit-trained mental model.

**Scope**: geometry only. Mesh directives, blocks, sidesets, imprint,
merge, and hex-sweep schemes are **intentionally absent** — they have
no direct build123d equivalent (see `lab_policy`).

## 3D primitives

| Cubit journal | build123d |
|---|---|
| `create brick x 10 y 20 z 5` | `Box(10, 20, 5)` |
| `create cylinder radius 3 height 10` | `Cylinder(radius=3, height=10)` |
| `create sphere radius 5` | `Sphere(radius=5)` |
| `create prism height 10 sides 6 radius 5` | `extrude(RegularPolygon(radius=5, side_count=6), amount=10)` |
| `create torus major radius 10 minor radius 2` | `Torus(major_radius=10, minor_radius=2)` |
| `create cone radius 3 radius2 1 height 10` | `Cone(bottom_radius=3, top_radius=1, height=10)` |

## 2D sketch primitives (prefixed with `BuildSketch` or algebra mode)

| Cubit journal | build123d |
|---|---|
| `create surface rectangle width 10 height 5` | `Rectangle(10, 5)` |
| `create surface circle radius 5` | `Circle(5)` |
| `create surface polygon sides 6 radius 5` | `RegularPolygon(5, side_count=6)` |

## Boolean operations

| Cubit journal | build123d |
|---|---|
| `subtract volume 2 from volume 1` | `v1 - v2` |
| `unite volume 1 2` | `v1 + v2` |
| `intersect volume 1 2` | `v1 & v2` |
| `chop volume 1 with plane zplane offset 0` | `split(v1, Plane.XY, keep=Keep.BOTH)` |

## Transforms

| Cubit journal | build123d |
|---|---|
| `move volume 1 x 10 y 0 z 0` | `Pos(10, 0, 0) * v1` |
| `rotate volume 1 angle 45 about z` | `Rot(0, 0, 45) * v1` or `v1.rotate(Axis.Z, 45)` |
| `scale volume 1 2` | `Pos((0,0,0), (0,0,0)) * v1.scale(2)` |
| `reflect volume 1 about yplane` | `mirror(v1, about=Plane.YZ)` |

## Sweep / revolve / loft

| Cubit journal | build123d |
|---|---|
| `sweep surface 1 along curve 1` | `sweep(sketch, path=line)` |
| `sweep surface 1 perpendicular distance 10` | `extrude(sketch, amount=10)` |
| `revolve surface 1 about z angle 360` | `revolve(sketch_on_plane_containing_axis, axis=Axis.Z, revolution_arc=360)` |
| `create volume loft surface 1 2 3` | `loft([sk1, sk2, sk3])` |

### Revolve idiom (common Cubit → build123d remap)

Cubit pattern — sketch on `zplane`, move out to radius, revolve about Z:
```
create surface rectangle width 5 height 10 zplane
move surface 1 x 20 y 0 z 0
revolve surface 1 about z angle 360
```
build123d equivalent — sketch placed on `Plane.XZ` (which contains the
Z axis) at offset via Pos, then `revolve`:
```python
profile = Plane.XZ * Pos(20, 0) * Rectangle(5, 10)
ring = revolve(profile, axis=Axis.Z)
```
Key swap: Cubit's "sketch in zplane then move" ↔ build123d's
"Plane.XZ * Pos(...) * sketch". The revolve axis must lie **in** the
sketch plane, so `Plane.XZ` or `Plane.YZ` are the common choices for
`axis=Axis.Z`.

## Fillet / chamfer

| Cubit journal | build123d |
|---|---|
| `modify curve 1 blend radius 0.5` | `fillet(edge_selection, radius=0.5)` |
| `modify curve in volume 1 blend radius 0.5` | `fillet(part.edges(), radius=0.5)` (all edges) |
| `chamfer curve 1 length 0.5` | `chamfer(edge_selection, length=0.5)` |

Edge selection in build123d uses selectors, e.g.
`part.edges().sort_by(Axis.Z)[-1]` for top edge,
`part.edges().group_by(Axis.Z)[-1]` for top-face ring,
`part.edges()` for **all edges** (the "modify curve in volume N" case).

## IDs vs variables

| Cubit model | build123d model |
|---|---|
| Implicit numbered IDs: `volume 1`, `surface 3`, `curve 5` | Explicit Python variables: `part`, `sk`, `edges` |
| `rename volume 1 "coil"` | `part.label = "coil"` |
| `group "coils" add volume 1 2 3` | `coils = Compound(children=[v1, v2, v3])`; `coils.label = "coils"` |
| `reset` | (no equivalent — build123d scripts start fresh each run, nothing to reset) |

### Compound vs fused Part

- `Compound(children=[v1, v2, v3])` — **tree-structured assembly**. Each
  child keeps its own `.label` and geometry; good for multi-region
  meshing via the `run_pipeline_multi` bridge (labels → physical groups).
- `v1 + v2 + v3` — **fused single solid** (boolean union). Labels on the
  children are lost; the result has a single volume. Use when the
  downstream consumer treats the geometry as one part (e.g., single
  domain FEM).
- Heuristic: if the downstream pipeline will want **per-region material
  names**, build a Compound. If you want one contiguous solid for a
  single-domain mesh, fuse with `+`.

## What has NO build123d equivalent (keep `.jou`!)

- `imprint all` / `merge all` — build123d doesn't expose multi-body
  shared-topology operations; it uses OCCT's Compound model directly.
  If a design needs imprinted+merged faces for conformal meshing with
  hex, **stay in Cubit**.
- `mesh volume N`, `scheme tetmesh`, `size ...`, `interval N` —
  meshing is Netgen's job downstream. Pass `maxh=...` to
  `OCCGeometry(...).GenerateMesh(maxh=maxh)`.
- `block 1 volume all`, `sideset 1 surface X`, `nodeset 1 curve Y` —
  use build123d labels + Compound children + the
  `run_pipeline_multi` bridge that maps labels to Gmsh physical groups.
- `export mesh "file.e"` (Exodus) — not part of the build123d pipeline.
  Export via NGSolve Mesh: `mesh.Export("file.msh", "Gmsh Format")`.
- Cubit APREPRO variables (`{var = 10}`) — use plain Python variables.

## Decision help

If more than ~20% of a `.jou` file is in the "no equivalent" table
above, **don't translate it** — keep the `.jou`, run it through Cubit,
export hex to Radia/ELF. build123d is the wrong tool for that design.

If the `.jou` is purely CAD primitives + booleans + sweep/revolve,
translation to build123d is straightforward using the tables above,
and the result will plug into the Netgen tet pipeline.
"""


EXAMPLES_LAB_PATTERNS = """\
# Lab-domain CAD patterns (build123d)

Parametric archetypes for the research lab's three application
domains: induction heating (IH), accelerator magnets, and magnetic
levitation, plus PEEC filament CAD. Each pattern is tet-pipeline
friendly (build123d → Netgen → Gmsh). Hex-required designs (e.g.,
large-scale Radia/ELF hex workpieces) stay in Cubit — see
`lab_policy`.

## IH coil archetypes

### 1. Single-turn ring coil (simplest)

```python
from build123d import Cylinder, Pos

r_in, r_out, h = 22.0, 28.0, 10.0   # mm
coil = Cylinder(radius=r_out, height=h) - Cylinder(radius=r_in, height=h)
coil.label = "coil"
```
Use: low-frequency IH concept demos, circuit-extraction sanity checks.

### 2. Helical multi-turn coil (realistic IH / inductor)

Use the `generate_helix_coil` MCP tool directly; or inline:
```python
import math
from build123d import BuildSketch, Rectangle, Plane, loft

R, pitch, n_turns = 30.0, 6.0, 5.0    # mm
w, h = 3.0, 3.0                        # cross-section
n_sec = int(n_turns * 12) + 1
sketches = []
for i in range(n_sec):
    t = i / (n_sec - 1)
    angle = 2 * math.pi * n_turns * t
    z = pitch * n_turns * t
    cx, cy = R * math.cos(angle), R * math.sin(angle)
    dx, dy, dz = (
        -R * math.sin(angle) * 2 * math.pi * n_turns,
         R * math.cos(angle) * 2 * math.pi * n_turns,
         pitch * n_turns,
    )
    norm = math.sqrt(dx*dx + dy*dy + dz*dz)
    tangent = (dx/norm, dy/norm, dz/norm)
    x_dir = (tangent[1], -tangent[0], 0)
    x_norm = math.hypot(x_dir[0], x_dir[1]) or 1.0
    x_dir = (x_dir[0]/x_norm, x_dir[1]/x_norm, 0)
    plane = Plane(origin=(cx, cy, z), x_dir=x_dir, z_dir=tangent)
    with BuildSketch(plane) as sk:
        Rectangle(w, h)
    sketches.append(sk.sketch)

coil = loft(sketches); coil.label = "coil"
```
Use: IH inductor field solve, PEEC filament extraction (feed the
centerline path to `section_along_path` tool).

### 3. Pancake coil (flat spiral)

```python
import math
from build123d import BuildSketch, BuildLine, Spline, make_face, Plane, extrude

R0, dR_per_turn, n_turns, track_w, thickness = 5.0, 2.5, 6, 1.5, 1.0
n_pts = 400
pts = []
for i in range(n_pts):
    t = i / (n_pts - 1) * n_turns
    r = R0 + dR_per_turn * t
    ang = 2 * math.pi * t
    pts.append((r * math.cos(ang), r * math.sin(ang)))

with BuildLine() as center:
    Spline(*pts)
# Sweep a rectangular cross-section of width track_w along this path,
# then extrude the resulting 2D trace by the PCB thickness.
# (For CAE, a simpler extruded-strip approximation suffices.)
```
Use: PCB-style IH pancakes, wireless power TX/RX coils.

## Accelerator magnet archetypes

### 4. Dipole yoke quarter (BEM-friendly with symmetry)

```python
from build123d import BuildSketch, BuildLine, Line, CenterArc, make_face, Plane, extrude

R_bore, R_outer, L_half = 30.0, 150.0, 200.0   # mm
with BuildSketch(Plane.XZ) as yoke_sk:
    with BuildLine():
        Line((R_bore, 0), (R_outer, 0))
        CenterArc((0, 0), R_outer, 0, 90)
        Line((0, R_outer), (0, R_bore))
        CenterArc((0, 0), R_bore, 90, -90)
    make_face()
yoke = extrude(yoke_sk.sketch, amount=L_half)
yoke.label = "iron_yoke"
```
Use: Radia dipole field, beam trajectory calculation (scope upper bound
= beam tracking; MBD is out of current scope per `lab_policy`).

### 5. Quadrupole pole tip (quarter model)

```python
# Hyperbolic pole face (simplified by a parabolic approximation
# adequate for initial BEM modeling)
from build123d import BuildSketch, BuildLine, Polyline, make_face, Plane, extrude

aperture, pole_half_width, L_half = 25.0, 30.0, 300.0
with BuildSketch(Plane.XZ) as pole:
    with BuildLine():
        Polyline(
            (aperture, aperture),
            (aperture + pole_half_width, aperture),
            (aperture + pole_half_width, aperture + 50),
            (aperture, aperture + 50),
            close=True,
        )
    make_face()
pole_tip = extrude(pole.sketch, amount=L_half)
pole_tip.label = "pole_tip_q1"
```
Use: quadrupole field modeling, pole shape optimization.

## Halbach segmented PM array

### 6. Halbach ring (labels-per-segment for magnetization direction)

Use the tested generator (no hand-rolled wedge math):

```python
from radia_mcp.build123d.archetypes import halbach_ring, parse_magnetization

# 12-segment DIPOLE Halbach ring (pole_pairs=1 -> uniform transverse bore field)
hb = halbach_ring(r_in=40, r_out=55, h=20, n_segments=12, pole_pairs=1, name="hb")
# each child is a real annular wedge whose label encodes its easy-axis angle:
for seg in hb.children:
    m_deg = parse_magnetization(seg.label)     # e.g. "hb_00_M30.000" -> 30.0
    # solver side: M = Br * (cos(m_deg), sin(m_deg), 0)
```

Key idea: each segment's **label encodes the magnetization (easy-axis) angle** via the
``_M<deg>`` suffix, set by the Mallinson law ``alpha_k = (pole_pairs+1)*theta_k`` (so the dipole
array's easy axis advances at TWICE the mechanical angle, a quadrupole at 3x).  The solver (Radia)
recovers it with ``archetypes.parse_magnetization`` and sets the magnetization direction;
``run_pipeline_multi`` (`examples/build123d_netgen_gmsh_flow/`) carries the labels to Gmsh physical
groups, preserving the mapping.  The underlying wedge is
``radia_mcp.build123d.modeling.annular_segment`` (full annulus intersected with a pie slice -- robust
for any span); ``polar_array`` tiles it.  See the `parametric_library` topic for the full op +
archetype set.

## PEEC filament CAD

### 7. Variable-cross-section helix for PEEC

Use the `generate_helix_coil` MCP tool with a taper:
```python
generate_helix_coil(
    radius=50.0, pitch=10.0, n_turns=5.0,
    w_start=4.0, h_start=4.0,   # bottom cross-section
    w_end=2.0,   h_end=2.0,     # top cross-section (thinner)
    sections_per_turn=12,
    label="tapered_coil",
    export_dir="runs",
)
```
Then call `section_along_path(step_file, helix_path_json)` with the
centerline from the same tool, and feed per-segment (w, h) into
`PEECBuilder.add_connected_segment(...)`.

## Magnetic levitation (scope: CAD + SimMechanics connection)

### 8. EMS attractive-force electromagnet + rail (SimMechanics ready)

```python
from build123d import Box, Cylinder, Pos, Compound

# Electromagnet (C-core)
core_w, core_h, core_d, gap = 60.0, 40.0, 50.0, 5.0
core = (Box(core_w, core_h, core_d)
        - Pos(0, 0, 0) * Box(core_w - 20, core_h - 10, core_d + 2))
core.label = "magnet_core"

# Rail (levitated body)
rail = Pos(0, -(core_h/2 + gap + 10/2), 0) * Box(core_w * 3, 10, core_d)
rail.label = "rail"

levitation = Compound(children=[core, rail])
levitation.label = "maglev_system"
```
Use: Radia force calculation per vertical position → feed to
SimMechanics as a force-vs-displacement table. Scope per lab_policy:
stops at SimMechanics connection (control loop is out of scope).

## Pattern → pipeline fit

| Pattern | Pipeline | Notes |
|---|---|---|
| Single-turn ring | tet, multi-region | `run_pipeline_multi`, label per region |
| Helical coil | tet | single region enough for field; PEEC via section helper |
| Pancake coil | tet | thin; watch min-edge quality |
| Dipole / Quad yoke | tet with symmetry planes | Radia handles quarter model natively |
| Halbach ring | tet, **multi-region** | segment labels drive magnetization |
| Tapered coil (PEEC) | post-solve (no mesh) | PEEC operates on filaments, not FEM mesh |
| Maglev EMS | tet, multi-region | rail + core, fuse or compound |
"""


JOINTS_AND_MATES = """# Joints, Locations, and Assembly Mates

build123d has first-class joint objects for kinematic / assembly
relationships, mirroring FreeCAD / Fusion 360 mates.

## Core joint types

| class            | motion                        | typical use                |
|------------------|-------------------------------|----------------------------|
| `RigidJoint`     | none — glued pose             | bracket ↔ plate            |
| `RevoluteJoint`  | single rotation axis          | hinge, axle                |
| `LinearJoint`    | single translation axis       | slide, drawer              |
| `CylindricalJoint` | rotate + slide along axis   | piston                     |
| `BallJoint`      | 3-dof rotation                | knee, trackball            |

## Minimal example

```python
from build123d import BuildPart, Box, Location, RigidJoint, export_step

# Two parts
with BuildPart() as plate:
    Box(40, 40, 5)
plate_part = plate.part
plate_part.label = "plate"

with BuildPart() as block:
    Box(10, 10, 10)
block_part = block.part
block_part.label = "block"

# Add named joints to each
RigidJoint("top",  plate_part, Location((0, 0, 2.5)))
RigidJoint("base", block_part, Location((0, 0, -5)))

# Connect
plate_part.joints["top"].connect_to(block_part.joints["base"])

export_step(plate_part, "plate_with_block.step")
```

## Lab conventions

1. Name every joint. Anonymous joints survive the STEP boundary as
   anonymous locations, losing semantic meaning.
2. Use `RigidJoint` for everything unless you explicitly need DOFs —
   FEM consumers (Cubit, ngsolve) treat parts as rigidly glued
   anyway.
3. For radia-coil assemblies (coil + iron pole + air gap), rigid-
   join each part to a shared origin frame so the exported STEP has
   coincident coordinate systems (mesh imprint/merge picks up the
   interfaces automatically).
"""


ASSEMBLIES_AND_COMPOUNDS = """# Assemblies: Compound, naming, STEP round-trip

build123d represents assemblies as `Compound` objects — a tree of
`Shape`s with shared children.  The STEP exporter walks the tree,
preserving labels as block names (critical for Cubit downstream).

## Build a compound

```python
from build123d import (
    BuildPart, Box, Cylinder, Compound, Location, Mode, export_step,
)

with BuildPart() as magnet:
    Box(20, 20, 10)
    magnet.part.label = "magnet"

with BuildPart() as pole:
    Box(30, 30, 5, align=(0.5, 0.5, 0))
    pole.part.label = "pole"
pole.part.locate(Location((0, 0, 10)))

with BuildPart() as air:
    Box(100, 100, 100, align=(0.5, 0.5, 0))
    with BuildPart(mode=Mode.SUBTRACT) as _cut:
        Box(20, 20, 10)
    air.part.label = "air"

asm = Compound(children=[magnet.part, pole.part, air.part])
export_step(asm, "asm.step")
```

## Labels matter

STEP round-trip preserves `shape.label`.  After `import step "asm.step"`
in Cubit, `block N add volume all; block N name "<label>"` will show
"magnet", "pole", "air" — so physical-group naming is free.

## Reverse: read assembly from STEP

```python
from build123d import import_step
compound = import_step("asm.step")
for child in compound.children:
    print(child.label, child.volume)
```
"""


CAE_WORKFLOW_TIPS = """# CAE workflow tips — getting from build123d to clean hex mesh

## Pipeline (lab standard, 2026-04-19)

```
build123d (Python) → STEP → Cubit (hex mesh) → .msh v4.1 or .vol
                                              → ngsolve / radia_ngsolve
                                              → Abaqus / Nastran / other FEA
```

One-call equivalents (radia-mcp ≥ 0.17):
  * `build123d_to_cubit_hex(script, target_size=1.0)` — build123d →
    STEP → `cubit_mesh_auto` → live GUI replay.
  * `cadquery_to_cubit_hex(...)` — same from cadquery source.

## Hex-friendly geometry checklist

1. **Prismatic**: extruded 2D sketches sweep to pure hex trivially.
2. **Simple topology**: ≤ 8 surfaces per volume usually sweeps. The
   `cubit_mesh_auto` geometry-split rung (0.16+) auto-webcuts when
   `surfaces/volume > 7`.
3. **Clean edges**: no micro-edges (< 1e-3 of characteristic length);
   `lint_build123d_script` flags fillets/chamfers below threshold.
4. **Coincident faces**: when assembling, rigidly-join parts so
   shared faces are exactly coincident → Cubit's `imprint all;
   merge all` picks them up without tolerance drift.
5. **Named bodies**: set `part.label = "..."` for every region you
   want as a Cubit block.

## Unit convention

build123d default units are **mm**. STEP preserves units. Cubit's
mesh-size command takes the same unit scale, so `volume all size 1.0`
= 1 mm hex edge on a 30-mm-side block = ~30 hexes per edge.

For SI workflows (meters), scale the Cubit side: `volume all size
0.001` = 1 mm for a model drawn in meters. Getting this wrong is
the #1 reason "mesh is too fine / too coarse" calls come in.

## Why build123d over cadquery for this pipeline

1. Builder API is less error-prone in complex assemblies (explicit
   context > chained `Workplane`).
2. First-class `Compound` + labels survive STEP → Cubit better.
3. Same OCCT backend, so interop is lossless via `execute_cadquery`
   / `cadquery_to_cubit_hex` when users bring cadquery scripts.
"""


PLANE_AXIS_LOCATION_COOKBOOK = """# Plane / Axis / Location cookbook

Three small classes cause ≥50% of "why doesn't this work" moments in
build123d. Memorize the patterns below.

## Plane: a sketching / operation surface with orientation

### Built-in planes

```python
from build123d import Plane

Plane.XY    # z-up, normal = +Z, x_dir = +X, origin (0,0,0)
Plane.XZ    # y-up, normal = +Y, x_dir = +X
Plane.YZ    # x-up, normal = +X, x_dir = +Y
Plane.front # alias for XZ (visualisation term)
Plane.top   # alias for XY
Plane.right # alias for YZ
# ... also Plane.back / Plane.bottom / Plane.left  (the negatives)
```

### Offset a plane along its normal

```python
sketch_plane = Plane.XY.offset(10)          # 10 mm up in +Z
```

### Build a plane at arbitrary origin + orientation

```python
# Plane with origin (10, 0, 0), normal pointing in -X, sketch x-axis = +Y
p = Plane(origin=(10, 0, 0), x_dir=(0, 1, 0), z_dir=(-1, 0, 0))
```

`x_dir` fixes the in-plane "right" direction; `z_dir` fixes the
out-of-plane normal. Both must be unit-ish; build123d normalizes.

### Plane × Location (rotate / translate a plane)

```python
rotated = Plane.XY * Rotation(0, 0, 45)      # 45° about Z
moved   = Plane.XY * Location((0, 0, 20))    # 20 mm up
combined = Plane.XY * Location((5, 0, 0), (90, 0, 0))  # translate + rotate
```

### Use a plane to place a BuildSketch

```python
from build123d import BuildPart, BuildSketch, Rectangle, extrude, Plane

with BuildPart() as bp:
    with BuildSketch(Plane.XY.offset(10)) as sk:    # sketch lives 10 mm up
        Rectangle(10, 5)
    extrude(amount=5)
```

### Use a face as a plane

```python
top_face = bp.part.faces().sort_by(Axis.Z)[-1]
with BuildSketch(top_face) as sk:                   # accept Face where Plane expected
    Rectangle(3, 3)
extrude(amount=2)
```

## Axis: an infinite line used for rotation / selection

### Built-in axes

```python
from build123d import Axis
Axis.X   # line through origin in +X
Axis.Y   # line through origin in +Y
Axis.Z   # line through origin in +Z
```

### Custom axis

```python
# Axis through (1, 0, 0) in direction (0, 1, 1) (normalized internally)
a = Axis(origin=(1, 0, 0), direction=(0, 1, 1))
```

### Rotate a shape about an axis

```python
# Note: Shape.rotate(axis, angle) — axis first, angle second
rotated = part.rotate(Axis.Z, 45)
rotated = part.rotate(Axis((1, 0, 0), (0, 0, 1)), 30)
```

Common mistake: `part.rotate(45, axis=(0,0,1))` → `TypeError: got
multiple values for argument 'axis'`. The shape method takes `axis`
as the FIRST positional. Pass an `Axis`, not a tuple.

### Select faces along an axis

```python
from build123d import SortBy
top_face = part.faces().sort_by(Axis.Z)[-1]   # highest face along Z
bottom   = part.faces().sort_by(Axis.Z)[0]    # lowest
```

## Location: rigid-body pose (position + orientation)

### Just a translation

```python
from build123d import Location
loc = Location((10, 0, 0))
```

### Translation + Euler rotation (degrees)

```python
loc = Location((10, 0, 0), (0, 0, 45))       # translate + rotate 45° about Z
loc = Location((0, 0, 0),  (90, 0, 0))       # rotate 90° about X
```

### Rotation-only

```python
from build123d import Rotation
loc = Rotation(0, 0, 45)                     # pure 45° about Z
```

### Apply a Location to a shape

```python
# Absolute placement (replaces any existing location)
placed = part.locate(loc)
# Relative placement (multiplies current location)
placed = part.located(loc)
# Shorthand via operator
placed = part * loc
```

### Use a Location to place a BuildPart context

```python
# Open a BuildPart whose "identity" is already at a specific pose
from build123d import BuildPart, Box, Mode, Location

with BuildPart() as outer:
    Box(10, 10, 10)
    with BuildPart(Location((0, 0, 10)), mode=Mode.SUBTRACT):
        Box(3, 3, 3)          # subtract a 3-cube centered at (0,0,10)
```

## Common recipes

### Sketch on the top face of an existing part

```python
top = part.faces().sort_by(Axis.Z)[-1]
with BuildSketch(top) as sk:
    Circle(3)
extrude(amount=5)
```

### Sketch on a plane rotated about Z

```python
with BuildSketch(Plane.XY * Rotation(0, 0, 30)) as sk:
    Rectangle(10, 5)
```

### Place a part at polar coordinates

```python
import math
for i in range(6):
    deg = 60 * i
    r = 20
    pose = Location((r * math.cos(math.radians(deg)),
                     r * math.sin(math.radians(deg)), 0), (0, 0, deg))
    with BuildPart(pose):
        Box(3, 3, 3)
```

### Offset a face's own plane outward by some distance

```python
cap = part.faces().sort_by(Axis.X)[-1]
cap_plane = Plane(cap).offset(5)    # 5 mm along the face normal
```
"""


BUILDER_VS_ALGEBRA_ROSETTA = """# Builder API vs Algebra API — side-by-side rosetta

build123d ships two styles of script. Pick ONE per script; mixing
leads to confusion. The lab standard is **Builder API** (clearer scope,
explicit context, less clever-one-liners).

| Goal | Builder API (`with ...` context) | Algebra API (direct operators) |
|---|---|---|
| 3-D box | `with BuildPart() as bp:`<br>`    Box(10, 10, 10)`<br>`part = bp.part` | `part = Box(10, 10, 10)` |
| Extrude sketch | `with BuildPart() as bp:`<br>`    with BuildSketch() as sk:`<br>`        Rectangle(10, 5)`<br>`    extrude(amount=3)` | `sk = Rectangle(10, 5)`<br>`part = extrude(sk, amount=3)` |
| Sweep along path | `with BuildPart() as bp:`<br>`    with BuildLine() as ln:`<br>`        Line((0,0,0), (10,0,0))`<br>`    with BuildSketch() as sk:`<br>`        Rectangle(2, 2)`<br>`    sweep(sk.sketch, path=ln.line)` | `ln = Line((0,0,0), (10,0,0))`<br>`sk = Rectangle(2, 2)`<br>`part = sweep(sk, path=ln)` |
| Subtract hole | `with BuildPart() as bp:`<br>`    Box(20, 20, 10)`<br>`    with BuildPart(mode=Mode.SUBTRACT):`<br>`        Cylinder(radius=3, height=20)` | `base = Box(20, 20, 10)`<br>`hole = Cylinder(radius=3, height=20)`<br>`part = base - hole` |
| Union (add) | `with BuildPart() as bp:`<br>`    Box(10, 10, 10)`<br>`    Sphere(radius=4)  # default Mode.ADD` | `part = Box(10, 10, 10) + Sphere(radius=4)` |
| Intersect | `with BuildPart() as bp:`<br>`    Box(10, 10, 10)`<br>`    Sphere(radius=6, mode=Mode.INTERSECT)` | `part = Box(10, 10, 10) & Sphere(radius=6)` |
| Fillet all edges | `fillet(bp.part.edges(), radius=0.5)` | `part = fillet(part.edges(), radius=0.5)` |
| Place on rotated plane | `with BuildSketch(Plane.XY * Rotation(0, 0, 45)):`<br>`    Rectangle(5, 5)` | `sk = Rectangle(5, 5) * Rotation(0, 0, 45) * Location((0,0,0))` |
| Assembly | `with BuildPart() as a: Box(10, 10, 10)`<br>`with BuildPart() as b: Sphere(5)`<br>`asm = Compound(children=[a.part, b.part])` | `a = Box(10, 10, 10)`<br>`b = Sphere(5)`<br>`asm = Compound(children=[a, b])` |

## When to pick which

| Situation | Prefer |
|---|---|
| Multi-step composition, nested modes | Builder (explicit scope) |
| One-liner parametric formula | Algebra (terser) |
| Teaching / onboarding | Builder (verbose → easier to debug) |
| Porting from CadQuery | Either, but Builder maps more clearly |
| Functional / immutable style | Algebra |

## Operators (Algebra API)

| Operator | Meaning |
|---|---|
| `a + b` | union |
| `a - b` | subtract |
| `a & b` | intersect |
| `a * loc` | apply Location (translate / rotate) |
| `a * rot` | apply Rotation |

## Mode (Builder API, and mode kwarg on Algebra ops)

```python
from build123d import Mode
Mode.ADD         # union (default for most primitives inside BuildPart)
Mode.SUBTRACT    # boolean subtract
Mode.INTERSECT   # boolean intersect
Mode.REPLACE     # overwrite (useful in BuildSketch)
Mode.PRIVATE     # scratchpad; not merged into context
```

## Sketch-side operators

| Goal | Builder | Algebra |
|---|---|---|
| Rectangle | `Rectangle(10, 5)` inside `BuildSketch` | `sk = Rectangle(10, 5)` |
| Subtract hole | `Circle(2, mode=Mode.SUBTRACT)` inside `BuildSketch` | `sk = Rectangle(10, 5) - Circle(2)` |
| Union | `Rectangle(10, 5)` + another shape in same `BuildSketch` | `sk = Rectangle(10, 5) + Circle(2)` |

## Gotcha: the "pending" concept

Inside a `BuildPart` with a nested `BuildLine`, the line becomes a
**pending path** that subsequent `sweep(sketch)` can consume. The
Algebra API has no pending concept — you pass `path=` explicitly.

```python
# Builder: sweep picks up the pending line automatically
with BuildPart() as bp:
    with BuildLine() as ln:
        Line((0,0,0), (10,0,0))
    with BuildSketch() as sk:
        Rectangle(2, 2)
    sweep(sk.sketch)                     # <-- no path=, uses ln.line

# Algebra: always pass path=
sk = Rectangle(2, 2)
ln = Line((0,0,0), (10,0,0))
part = sweep(sk, path=ln)                # <-- explicit
```

`lint_build123d_script` flags `sweep(sketch)` positional calls
without `path=` because the pending-curve implicit mode is a frequent
source of `BRep_API: Make Shape failed` when forgotten.
"""


from .build123d_api import get_api_reference as _get_b3d_api_reference


PARAMETRIC_LIBRARY = r"""
# Parametric modeling library (radia_mcp.build123d.modeling + .archetypes)

A small, **tested** library of solid-modeller operations and EM-device archetypes -- the parametric
building blocks you would author in a parametric solid modeller, as clean Netgen-meshable build123d
solids with region labels.  Pure build123d (no radia_mcp dependency at runtime), so they import as a
library OR copy-paste into an `execute_build123d` subprocess.

## modeling -- generic operations (the building blocks)

| helper | what it makes |
|---|---|
| `annular_segment(r_in, r_out, h, start_angle, end_angle)` | a radial WEDGE of an annulus (robust for any span < 360: annulus INTERSECTED with a pie slice -- no self-intersecting sketch) |
| `tube(r_in, r_out, h)` | hollow cylinder |
| `racetrack_coil(length, width, band, h, corner_radius)` | rounded-rectangle (racetrack) coil outline, extruded |
| `polar_array(part, count, total_angle=360)` | rotate-with-copies -> labelled Compound (full ring or partial fan) |
| `linear_array(part, count, spacing, direction)` | translate-with-copies |
| `mirrored(part, about=Plane.XZ)` | symmetry completion (original + mirror) |
| `assembly(*parts)` | group labelled regions into ONE multi-region Compound (does NOT fuse) |

## archetypes -- EM devices (composed from the ops)

| helper | device |
|---|---|
| `cylindrical_magnet(radius, h, m_angle_deg)` / `block_magnet(length, width, h, m_angle_deg)` | a PM primitive, easy-axis angle in the label |
| `halbach_ring(r_in, r_out, h, n_segments, pole_pairs=1)` | segmented Halbach PM ring; per-segment easy axis = Mallinson `(pole_pairs+1)*theta` |
| `pole_tip(base_width, tip_width, height, depth)` | trapezoidal shaped pole piece |
| `c_core(width, height, depth, leg, gap)` | C-shaped electromagnet yoke with a pole gap |
| `multipole_yoke(n_poles, r_bore, pole_len, pole_width, yoke_thickness, depth)` | n-pole iron yoke (dipole/quad/sextupole) + return ring |
| `h_dipole(width, height, depth, leg, pole_width, gap)` | H-frame dipole yoke (window frame + 2 poles + gap) |
| `solenoid(r_in, r_out, h)` / `helmholtz_pair(r_in, r_out, h, separation)` | solenoid bundle / coaxial coil pair |
| `cos_theta_dipole(radius, conductor_w, conductor_h, length, n_per_half)` | cos-theta winding (arcsin-spaced axial bars -> pure dipole) |
| `magnetization_map(compound, Br)` | close the loop: PM region labels -> `{label: (Mx, My)}` for the solver |

The magnetization convention is verified end-to-end against physics in
`tests/test_build123d_halbach_field.py`: a `halbach_ring` -> `magnetization_map` -> 2D A_z PM solve
gives a UNIFORM transverse bore field of `Br*ln(r_out/r_in)` (dipole), and a quadrupole Halbach
(`pole_pairs=2`) a field that grows from a null centre.  (A swept saddle coil is future work -- the
closed-path sweep is fragile in the OCCT kernel; use `cos_theta_dipole` for transverse-field windings.)

## Region labels & magnetization (the solver hand-off)

build123d keeps region labels on a Compound's **`.children`** (NOT on the flattened `.solids()`, which
drops them).  PM regions encode their easy-axis angle in the label via `archetypes.magnetization_tag`
(`"{name}_{i:02d}_M{deg:.3f}"`); the solver recovers it with `archetypes.parse_magnetization` and sets
`M = Br*(cos, sin, 0)`.  `run_pipeline_multi` carries these labels to Gmsh physical groups.

```python
from radia_mcp.build123d.modeling import annular_segment, polar_array, assembly, tube
from radia_mcp.build123d.archetypes import halbach_ring, c_core, parse_magnetization

hb = halbach_ring(40, 55, 20, 12, pole_pairs=1, name="pm")     # dipole Halbach
yoke = c_core(120, 90, 60, 18, 26, name="yoke")
machine = assembly(*hb.children, yoke, tube(10, 18, 40, "coil"), label="device")
mags = {seg.label: parse_magnetization(seg.label) for seg in hb.children}
```

Validated: tests/test_build123d_modeling.py, tests/test_build123d_archetypes.py (analytic volumes,
OCCT-valid, region labels, Netgen-meshable, Mallinson easy axes).
"""


_TOPICS = {
    "overview": OVERVIEW,
    "parametric_library": PARAMETRIC_LIBRARY,
    "lab_policy": LAB_POLICY,
    "cubit_rosetta": CUBIT_ROSETTA,
    "examples_lab_patterns": EXAMPLES_LAB_PATTERNS,
    "primitives_3d": PRIMITIVES_3D,
    "primitives_2d": PRIMITIVES_2D,
    "operations": OPERATIONS,
    "curves": CURVES,
    "export_import": EXPORT_IMPORT,
    "topology": TOPOLOGY,
    "cae_guidelines": CAE_GUIDELINES,
    "examples": EXAMPLES,
    "examples_intro": EXAMPLES_INTRO,
    "examples_gallery": EXAMPLES_GALLERY,
    "coil_modeling": COIL_MODELING,
    "coil_from_polyline": COIL_FROM_POLYLINE,
    "polyline_coil": COIL_FROM_POLYLINE,  # alias
    "centerline_coil": COIL_FROM_POLYLINE,  # alias
    "coil_from_points": COIL_FROM_POLYLINE,  # alias
    "frame_sweep": COIL_FROM_POLYLINE,  # alias
    "loft_polyline": COIL_FROM_POLYLINE,  # alias
    "joints_and_mates": JOINTS_AND_MATES,
    "assemblies_and_compounds": ASSEMBLIES_AND_COMPOUNDS,
    "cae_workflow_tips": CAE_WORKFLOW_TIPS,
    "plane_axis_location_cookbook": PLANE_AXIS_LOCATION_COOKBOOK,
    "builder_vs_algebra_rosetta": BUILDER_VS_ALGEBRA_ROSETTA,
    "api_reference": _get_b3d_api_reference(),
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
