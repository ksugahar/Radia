# PEEC Conductor Modeling Guide

**Date**: 2026-05-16 (revised from 2026-04-17; original proposal 2026-02-13)
**Status**: Implemented — filament model via `coil_from_cad.py`

---

> ## v4.48.2 -> v4.55.0 STEP-loading subsystem overhaul (2026-05-15 / 2026-05-16)
>
> This guide's architectural overview is correct, but the implementation
> underwent an 8-release fail-fast hardening campaign that significantly
> changed the internal behaviour of `filaments_from_step`.  Summary of
> the new layered defense (all enforced in current code):
>
> 1. **Multi-solid entry guard** -- raise on STEP containing > 1 solid
>    (v4.49.0 Tier D).
> 2. **Classification-based single dispatch** in
>    `extract_centerline_from_step`: 5 positive-match predicates
>    (Loft cross-sections / Circle-edge stations /
>    Revolution-sweep / OPEN longest-edge / CLOSED topology-spine) --
>    NO `try/except` cascade, NO `path_points_m` JSON override
>    (v4.49.0 Tier A).
> 3. **Adaptive resampling upfront** in `_centerline_from_open_spine`:
>    midpoint section to estimate `wire_radius`, then cap `n_segments`
>    at `floor(spine_length / (1.10 * wire_radius))` so adjacent
>    stations are always >= 1.10 wire-radius apart (v4.53.0).
> 4. **Wang-Joe Rotation-Minimizing Frame** replaces Rodrigues
>    parallel-transport in `_parallel_transport_frame` -- provably
>    minimum accumulated twist on polylines with kinks (v4.54.0).
> 5. **Corner densification** in `_centerline_from_open_spine`:
>    after adaptive resampling, insert intermediate spine points ON
>    the OCC curve near sharp bends until per-step bend angle <= 20 deg
>    (v4.54.0).
> 6. **Cap-centroid endpoint anchoring** in
>    `_centerline_from_open_spine`: replace rim-spine endpoints with
>    cap-face centroids from `coil_topology.extract_coil_topology` --
>    fixes the rim-to-cap kink at vertex N-1 that produced asymmetric
>    |I| distributions (v4.55.0).
> 7. **CCW winding normalisation** in `_sample_face_perimeter_in_pt_frame`:
>    signed-area check + reverse if CW; fixes per-segment Cubit lofts
>    where shared cross-section faces have alternating orientation
>    (v4.53.0, keiko's patch verbatim).
> 8. **Three orthogonal positive proofs** on every centerline / topo:
>    `_check_filaments_cover_solid_bbox` (under-coverage),
>    `_check_centerline_inside_solid` (gross wrong-location, bbox+5%
>    slack), `_check_centerline_near_solid_surface` (per-point
>    distance to solid surface via `BRepExtrema_DistShapeShape`,
>    tolerance 1.10 * wire_radius, 20-point subsample).  ALL three
>    must pass; orthogonal failure modes, NOT a fallback chain
>    (v4.50.0 + v4.51.0).
> 9. **Pre-Ruehli singular-corner check** in `filaments_from_polyline`:
>    raise if `bend > 60 deg AND adj_min_seg_len < wire_radius` --
>    covers HACApK path too because it runs BEFORE solver assembly
>    (v4.49.0 Tier B).
> 10. **Post-assembly finite-L safety net** in
>     `peec_bundle._assert_solver_L_finite`: raise on any non-finite
>     entry in `solver.L`, with diagnostic naming the offending
>     filament/segment pairs (v4.48.2; mostly belt-and-suspenders
>     since the v4.49.0 corner detect catches the same condition
>     pre-assembly).
>
> The `filaments_from_step(step_path, sigma=..., n_peri=...)` API is
> unchanged from v4.48.1; `path_points_m` kwarg was REMOVED in v4.48.0.
>
> **For the runnable how-to layer** (Cubit + build123d recipes that
> produce STEPs the auto-detect dispatch can solve, anti-patterns,
> 10-line build123d probe script), query
> `radia-mcp peec_inductance(topic="step_authoring")` -- the MCP
> knowledge module is updated with every release and is the
> canonical source per CLAUDE.md "MCP Knowledge Placement Policy".

---

## 1. Overview

This guide consolidates the design and workflow for conductor modeling in the Radia PEEC framework. It covers:

- **Conductor shape interface**: The current API (`PEECBuilder`) and proposed extensions for additional geometric primitives and path-based shapes.
- **Mesh import workflow**: How to bring arbitrary CAD geometry into the PEEC solver via Coreform Cubit and GMSH 1D edge meshes.
- **Design philosophy**: Why mesh import is the preferred approach and how cross-section metadata and port definitions are handled.
- **API reference**: Full specifications for existing and proposed builder methods, helper functions, and configuration options.

### High-Level Workflow

```
CAD Model (STEP)
    |
coil_from_cad.filaments_from_step()    -- auto-extract centerline + cross-sections
    |
PEECBuilder topology -> PEECCircuitSolver
    |
Impedance, Inductance, AC loss
```

The original proposal (2026-02-13) planned a Cubit 1D edge mesh -> GMSH import path.
This was replaced by the **filament model** (`coil_from_cad.py`), which extracts
PEEC topology directly from a STEP solid without meshing.

### Current Interface (Working)

```python
from peec_matrices import PEECBuilder

builder = PEECBuilder()
builder.create_wire([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10)  # Straight wire
builder.create_loop([0,0,0], 0.05, [0,0,1], 2e-3, 2e-3, 36)  # Circle
L, R, P, M_LS = builder.build()
```

---

## 2. Conductor Shape Interface

### 2.1 Basic Shapes (High Priority)

| Shape | API | Use Case |
|-------|-----|----------|
| **Racetrack** | `create_racetrack()` | Elongated coils, MagLev |
| **Arc** | `create_arc()` | Partial loops, sector coils |
| **Helix** | `create_helix()` | Solenoids, helical inductors |
| **Spiral** | `create_spiral()` | Planar inductors, PCB coils |

#### API Examples

```python
# Racetrack (elongated loop)
builder.create_racetrack(
    center=[0, 0, 0],
    length=0.2,         # Straight section length
    radius=0.05,        # End cap radius
    normal=[0, 0, 1],   # Plane normal
    width=2e-3,
    height=2e-3,
    n_segments=72       # Total segments (including straight + arc)
)

# Arc (partial circle)
builder.create_arc(
    center=[0, 0, 0],
    radius=0.05,
    normal=[0, 0, 1],
    start_angle=0,      # degrees
    end_angle=180,      # degrees
    width=2e-3,
    height=2e-3,
    n_segments=18
)

# Helix (3D spiral)
builder.create_helix(
    center=[0, 0, 0],
    radius=0.05,
    pitch=0.01,         # Axial distance per turn
    n_turns=10,
    axis=[0, 0, 1],     # Helix axis
    width=2e-3,
    height=2e-3,
    n_segments_per_turn=36
)

# Planar spiral
builder.create_spiral(
    center=[0, 0, 0],
    r_inner=0.01,       # Inner radius
    r_outer=0.05,       # Outer radius
    n_turns=5,
    normal=[0, 0, 1],   # Spiral plane
    width=2e-3,
    height=2e-3,
    n_segments=180
)
```

### 2.2 Path-Based Shapes (Medium Priority)

For arbitrary conductor paths:

```python
# From list of points
points = [[0,0,0], [0.1,0,0], [0.1,0.1,0], [0,0.1,0]]
builder.create_polyline(points, width=2e-3, height=2e-3, n_segments=40)

# From parametric curve
def curve(t):
    # t in [0, 1]
    x = 0.1 * np.cos(2*np.pi*t)
    y = 0.1 * np.sin(2*np.pi*t)
    z = 0.05 * t
    return [x, y, z]

builder.create_parametric_curve(curve, t_samples=100,
                                width=2e-3, height=2e-3)
```

### 2.3 Composite Shapes (Low Priority)

For complex assemblies:

```python
# Group multiple conductors
coil1 = builder.create_loop([0,0,0], 0.05, [0,0,1], 2e-3, 2e-3, 36)
coil2 = builder.create_loop([0,0,0.1], 0.05, [0,0,1], 2e-3, 2e-3, 36)

# Series connection
builder.connect_series([coil1, coil2], port_start=0, port_end=-1)

# Parallel connection
builder.connect_parallel([coil1, coil2])
```

### 2.4 Implementation Priority

| Feature | Priority | Complexity | Applications |
|---------|----------|------------|--------------|
| **Racetrack** | HIGH | Low | MagLev, induction heating |
| **Arc** | HIGH | Low | Sector coils, partial loops |
| **Helix** | MEDIUM | Medium | Solenoids, RF coils |
| **Spiral** | MEDIUM | Medium | Planar inductors, WPT |
| **Polyline** | MEDIUM | Low | Arbitrary paths from CAD |
| **Parametric** | LOW | Medium | Research, custom geometries |
| **Composite** | LOW | High | Multi-coil systems |

---

## 3. Conductor Modeling (Filament Model)

> **Note**: The original 1D edge mesh import workflow (2026-02-13 proposal) has been
> **superseded** by the filament model in `coil_from_cad.py`. See Section 3.4-3.6 below.

### 3.1 Legacy: Cubit 1D Edge Mesh (Superseded)

The following workflow was the original proposal. It is no longer recommended:

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'coil.step'")

# 1D edge mesh
curve_id = cubit.get_last_id("curve")
cubit.cmd(f"curve {curve_id} interval 36")
cubit.cmd(f"mesh curve {curve_id}")
cubit.cmd(f"block 1 add curve {curve_id}")
cubit.cmd("block 1 name 'conductor'")
cubit.cmd('export gmsh "coil_mesh.msh" overwrite')
```

**Use `filaments_from_step()` instead** — it extracts topology directly from the STEP solid.

### 3.2 Step 2 -- Import to Radia PEEC

```python
import gmsh
from peec_matrices import PEECBuilder

# Load mesh
gmsh.initialize()
gmsh.open("coil_mesh.msh")

# Get nodes and edges
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coords = node_coords.reshape(-1, 3) * 1e-3  # mm to m

elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

# Extract edge elements
edges = []
for i, elem_type in enumerate(elem_types):
    if elem_type == 1:  # 2-node line
        node_tags_flat = elem_node_tags[i]
        n_elems = len(node_tags_flat) // 2
        for j in range(n_elems):
            n0 = node_tags_flat[j*2]
            n1 = node_tags_flat[j*2 + 1]
            edges.append((n0, n1))

# Create PEEC segments
builder = PEECBuilder()

# Cross-section parameters (NOT in mesh)
width = 4e-3   # 4mm
height = 4e-3  # 4mm
sigma = 5.8e7  # S/m (copper)

for n0, n1 in edges:
    idx0 = np.where(node_tags == n0)[0][0]
    idx1 = np.where(node_tags == n1)[0][0]
    p0 = coords[idx0]
    p1 = coords[idx1]

    builder.create_wire(p0, p1, width, height, 1, sigma)

# Build matrices
L, R, P, M_LS = builder.build()
```

### 3.3 Known Problems and Solutions

#### Problem 1: Manual Cross-Section Specification

`width` and `height` must be specified in the Python script, not in the mesh.

**Option A: Block Attributes (Recommended)**

Use Cubit block attributes to store cross-section info:

```python
# In Cubit:
cubit.cmd("block 1 attribute count 3")
cubit.cmd("block 1 attribute index 1 4.0")   # width [mm]
cubit.cmd("block 1 attribute index 2 4.0")   # height [mm]
cubit.cmd("block 1 attribute index 3 5.8e7") # sigma [S/m]

# In Python (auto-read from mesh):
block_attrs = get_block_attributes(mesh, block_id=1)
width = block_attrs['width'] * 1e-3   # mm to m
height = block_attrs['height'] * 1e-3
sigma = block_attrs['sigma']
```

**Option B: GMSH Physical Group Names**

Encode cross-section in the physical group name:

```python
# In Cubit:
cubit.cmd("block 1 name 'conductor_w4.0_h4.0_s5.8e7'")

# In Python (parse name):
import re
name = "conductor_w4.0_h4.0_s5.8e7"
match = re.match(r'conductor_w([\d.]+)_h([\d.]+)_s([\d.e+]+)', name)
width = float(match.group(1)) * 1e-3
height = float(match.group(2)) * 1e-3
sigma = float(match.group(3))
```

**Option C: Separate Configuration File**

```yaml
# peec_config.yaml
conductors:
  - block: 1
    name: "coil"
    width: 4.0e-3   # m
    height: 4.0e-3  # m
    sigma: 5.8e7    # S/m
```

**Recommendation**: Use **Option A (Block Attributes)** -- cleanest and most robust.

#### Problem 2: Manual Port Definition

Ports are currently defined by coordinate-based search, which is error-prone and tedious.

**Current approach**:

```python
# Find node closest to target position
port_positive_target = np.array([r_mean, 0, 0])
min_dist = float('inf')
for i, tag in enumerate(node_tags):
    dist = np.linalg.norm(coords[i] - port_positive_target)
    if dist < min_dist:
        port_positive_node = tag
```

**Proposed solution**: Use Cubit **nodesets** to mark port nodes.

```python
# In Cubit:
cubit.cmd("nodeset 1 add node <ID>")  # Port positive
cubit.cmd("nodeset 1 name 'port_positive'")
cubit.cmd("nodeset 2 add node <ID>")  # Port negative
cubit.cmd("nodeset 2 name 'port_negative'")

# In Python (auto-read from mesh):
port_positive_nodes = get_nodeset(mesh, "port_positive")
port_negative_nodes = get_nodeset(mesh, "port_negative")
```

**Note**: The radia Cubit plugin does NOT currently support nodeset export to GMSH format. Until that support is added, the coordinate-based search remains the workaround.

### 3.4 Filament Model (Implemented)

The 1D edge mesh approach was superseded by the **filament model** in
`coil_from_cad.py`. This module extracts PEEC topology directly from a
STEP solid by perpendicular cross-sectioning along the coil centerline.

#### Pipeline

```
STEP solid -> extract_centerline_from_step() -> path + per-segment (w, h)
           -> build_peec_from_path()          -> PEECBuilder topology
           -> PEECCircuitSolver               -> Z, L, R
```

#### Key Functions (`radia.coil_from_cad`)

| Function | Description |
|----------|-------------|
| `helix_path(radius, pitch, n_turns, n_points)` | Generate discrete helix centerline |
| `build_peec_from_path(path, widths, heights, sigma)` | Build PEECBuilder topology from path + cross-sections |
| `extract_centerline_from_step(step_path, n_segments)` | Auto-extract centerline (marching engine + CAD-feature fast paths) |
| `filaments_from_step(step_path, ...)` | End-to-end: STEP solid -> PEEC topology |

### 3.5 Example End-to-End Workflow

#### Automatic (no explicit path needed)

```python
from radia.coil_from_cad import filaments_from_step

# Auto-extract centerline + cross-sections from STEP solid
topo = filaments_from_step(
    "induction_coil.step",
    sigma=5.8e7,
    nwinc=3, nhinc=3,          # sub-filament subdivision
    cad_units_per_meter=1000,   # STEP file in mm
    n_slices=200,
)

# Solve
from radia.peec_topology import PEECCircuitSolver
solver = PEECCircuitSolver(topo, use_hacapk=True, outer_method="saddle")
Z = solver.solve_impedance(frequency=50e3)
print(f"Z @ 50 kHz: {Z:.4e} Ohm")
```

#### With explicit path (e.g. helix)

```python
from radia.coil_from_cad import helix_path, build_peec_from_path
import numpy as np

# Define helix path: R=50mm, pitch=10mm, 5 turns
path = helix_path(radius=0.05, pitch=0.01, n_turns=5, n_points=501)

n_seg = path.shape[0] - 1
widths = np.full(n_seg, 4e-3)    # 4mm width
heights = np.full(n_seg, 4e-3)   # 4mm height

topo = build_peec_from_path(path, widths, heights, sigma=5.8e7)
```

### 3.6 Implementation Status

| Feature | Status |
|---------|--------|
| Filament model (`coil_from_cad.py`) | Done |
| Auto centerline extraction from STEP | Done |
| Per-segment cross-section from STEP solid | Done |
| Sub-filament subdivision (nwinc, nhinc) | Done |
| PEECCircuitSolver + HACApK | Done |
| ~~1D edge mesh import from GMSH~~ | Superseded by filament model |
| ~~`create_peec_from_mesh()` helper~~ | Not needed (filaments_from_step replaces) |

---

## 4. Design Philosophy

### Key Decisions

- **Mesh import** is the primary approach (flexible, handles arbitrary geometry).
- **Rectangular cross-section only** (simple and practical for the current scope).
- **Coreform Cubit for CAD** (professional tool -- avoids the need for custom importers).
- Parametric shapes (racetrack, spiral, etc.) are proposed extensions but are considered lower priority than a robust mesh import pipeline.

### Design Questions

1. **Cross-section specification**:
   - Current: Rectangular (`width`, `height`)
   - Future: Circular? Trapezoidal? Profile function?

2. **Segmentation control**:
   - Current: `n_segments` (uniform)
   - Future: Adaptive? Grading?

3. **Port definition**:
   - Current: Manual (find nodes by coordinates)
   - Future: Automatic via nodesets? Named ports?

4. **Units**:
   - Radia always uses meters (SI)
   - Keep consistent with Radia convention

---

## 5. API Reference

### PEECBuilder -- Existing Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `create_wire` | `create_wire(p0, p1, width, height, n_segments, sigma=5.8e7)` | Straight wire segment between two points |
| `create_loop` | `create_loop(center, radius, normal, width, height, n_segments)` | Full circular loop |
| `build` | `build()` | Compute and return `(L, R, P, M_LS)` matrices |

### PEECBuilder -- Proposed Methods

| Method | Signature | Priority | Description |
|--------|-----------|----------|-------------|
| `create_racetrack` | `create_racetrack(center, length, radius, normal, width, height, n_segments)` | HIGH | Elongated loop with straight sections and semicircular end caps |
| `create_arc` | `create_arc(center, radius, normal, start_angle, end_angle, width, height, n_segments)` | HIGH | Partial circular arc |
| `create_helix` | `create_helix(center, radius, pitch, n_turns, axis, width, height, n_segments_per_turn)` | MEDIUM | Helical coil along an axis |
| `create_spiral` | `create_spiral(center, r_inner, r_outer, n_turns, normal, width, height, n_segments)` | MEDIUM | Planar spiral from inner to outer radius |
| `create_polyline` | `create_polyline(points, width, height, n_segments)` | MEDIUM | Conductor following a sequence of points |
| `create_parametric_curve` | `create_parametric_curve(curve_fn, t_samples, width, height)` | LOW | Conductor following an arbitrary parametric curve `f(t)` |
| `connect_series` | `connect_series(conductors, port_start, port_end)` | LOW | Series-connect multiple conductor groups |
| `connect_parallel` | `connect_parallel(conductors)` | LOW | Parallel-connect multiple conductor groups |

### Filament Model (`radia.coil_from_cad`)

| Function | Description |
|----------|-------------|
| `filaments_from_step(step_path, ...)` | End-to-end: STEP solid -> PEEC topology |
| `build_peec_from_path(path, widths, heights, sigma)` | Path + cross-sections -> topology |
| `extract_centerline_from_step(step_path, n_segments)` | Auto-extract centerline from STEP |
| `helix_path(radius, pitch, n_turns, n_points)` | Generate helix centerline |

### Build Output

`builder.build()` returns a tuple of four matrices:

| Matrix | Symbol | Description |
|--------|--------|-------------|
| `L` | Partial inductance | Segment self- and mutual inductances |
| `R` | Resistance | Segment resistances |
| `P` | Potential coefficients | Segment self- and mutual potential coefficients |
| `M_LS` | Inductance-segment coupling | Coupling between inductance and segment meshes |

---

## Next Steps

1. **Variable cross-section optimization** — LLM-driven coil design via build123d + filament PEEC
2. **Profile abstraction** — non-rectangular cross-sections (circular, trapezoidal)
3. **Multi-turn coils** — series/parallel connection of filament groups

## See Also

- [../kelvin/KELVIN_TRANSFORMATION.md](../kelvin/KELVIN_TRANSFORMATION.md) — Open boundary for FEM+PEEC coil source
- [../api/API_REFERENCE.md](../api/API_REFERENCE.md) — Full Python API reference
