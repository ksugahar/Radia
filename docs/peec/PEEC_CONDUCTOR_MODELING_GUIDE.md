# PEEC Conductor Modeling Guide

**Date**: 2026-02-13
**Status**: Proposal / Implementation Plan

---

## 1. Overview

This guide consolidates the design and workflow for conductor modeling in the Radia PEEC framework. It covers:

- **Conductor shape interface**: The current API (`PEECBuilder`) and proposed extensions for additional geometric primitives and path-based shapes.
- **Mesh import workflow**: How to bring arbitrary CAD geometry into the PEEC solver via Coreform Cubit and GMSH 1D edge meshes.
- **Design philosophy**: Why mesh import is the preferred approach and how cross-section metadata and port definitions are handled.
- **API reference**: Full specifications for existing and proposed builder methods, helper functions, and configuration options.

### High-Level Workflow

```
CAD Model (STEP/IGES)
    |
Coreform Cubit (geometry + meshing)
    |
1D Edge Mesh Export (GMSH format)
    |
Radia PEEC (segment creation + solve)
```

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

## 3. Mesh Import Workflow

### 3.1 Step 1 -- Cubit Mesh Generation

```python
import cubit
import cubit_mesh_export

cubit.init(['cubit', '-nojournal', '-batch'])

# Import CAD or create geometry
cubit.cmd("import step 'coil.step'")

# OR create directly in Cubit
cubit.cmd(f"create curve arc radius 50 center 0 0 0 normal 0 0 1 "
          f"start angle 0 stop angle 360")

# Mesh with 1D edge elements
curve_id = cubit.get_last_id("curve")
cubit.cmd(f"curve {curve_id} interval 36")
cubit.cmd(f"curve {curve_id} scheme equal")
cubit.cmd(f"mesh curve {curve_id}")

# Define block (physical group)
cubit.cmd(f"block 1 add curve {curve_id}")
cubit.cmd("block 1 name 'conductor'")

# Export to GMSH v2.2
cubit_mesh_export.export_Gmesh(cubit, "coil_mesh.msh")
```

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

**Note**: `cubit_mesh_export` does NOT currently support nodeset export to GMSH format. Until that support is added, the coordinate-based search remains the workaround.

### 3.4 Helper Function

A wrapper to simplify mesh import:

```python
def create_peec_from_mesh(mesh_file, cross_section_config):
    """
    Create PEEC model from GMSH 1D edge mesh.

    Parameters:
    -----------
    mesh_file : str
        Path to GMSH .msh file
    cross_section_config : dict
        {block_id: {'width': float, 'height': float, 'sigma': float}}

    Returns:
    --------
    builder : PEECBuilder
        Builder with segments loaded
    """
    import gmsh
    from peec_matrices import PEECBuilder

    gmsh.initialize()
    gmsh.open(mesh_file)

    # Get nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = node_coords.reshape(-1, 3) * 1e-3  # mm to m

    # Get elements by block
    builder = PEECBuilder()

    for block_id, params in cross_section_config.items():
        # Get edges in this block
        edges = get_edges_in_block(gmsh, block_id)

        # Create segments
        for n0, n1 in edges:
            p0 = get_node_coord(coords, node_tags, n0)
            p1 = get_node_coord(coords, node_tags, n1)

            builder.create_wire(p0, p1,
                              params['width'],
                              params['height'],
                              1,
                              params['sigma'])

    gmsh.finalize()
    return builder

# Usage:
builder = create_peec_from_mesh(
    "coil_mesh.msh",
    cross_section_config={
        1: {'width': 4e-3, 'height': 4e-3, 'sigma': 5.8e7}
    }
)

L, R, P, M_LS = builder.build()
```

### 3.5 Example End-to-End Workflow (Target)

#### Cubit Script

```python
import cubit
import cubit_mesh_export

cubit.init(['cubit', '-nojournal', '-batch'])

# Import CAD
cubit.cmd("import step 'induction_coil.step'")

# Mesh
cubit.cmd("curve all interval 50")
cubit.cmd("mesh curve all")

# Define conductor with cross-section
cubit.cmd("block 1 add curve all")
cubit.cmd("block 1 name 'primary_coil'")
cubit.cmd("block 1 attribute count 3")
cubit.cmd("block 1 attribute index 1 6.0")   # width [mm]
cubit.cmd("block 1 attribute index 2 6.0")   # height [mm]
cubit.cmd("block 1 attribute index 3 5.8e7") # sigma [S/m]

# Export
cubit_mesh_export.export_Gmesh(cubit, "induction_coil.msh")
```

#### Radia Python Script

```python
from peec_mesh_import import create_peec_from_mesh

# Auto-load with attributes from mesh
builder = create_peec_from_mesh("induction_coil.msh", auto_config=True)

# OR manual override:
builder = create_peec_from_mesh(
    "induction_coil.msh",
    cross_section_config={1: {'width': 6e-3, 'height': 6e-3, 'sigma': 5.8e7}}
)

# Build matrices
L, R, P, M_LS = builder.build()

# Port impedance at 50 kHz
I_port = define_port_excitation(builder, port_positive=(0.1, 0, 0),
                                          port_negative=(-0.1, 0, 0))
Z = builder.compute_impedance(50e3, I_port)
print(f"Z @ 50 kHz: {Z:.4e} Ohm")
```

### 3.6 Implementation Plan

#### Phase 1: Improve Mesh Import (Current)

- [x] Basic 1D edge mesh import from GMSH
- [x] Manual cross-section specification
- [x] Coordinate-based port search
- [ ] **Add `create_peec_from_mesh()` helper function**

#### Phase 2: Block Attributes (Optional)

- [ ] Implement block attribute reading from GMSH
- [ ] Auto-extract cross-section from mesh metadata
- [ ] Fallback to manual specification if not found

#### Phase 3: Port Handling (Future)

- [ ] Request nodeset support in `cubit_mesh_export`
- [ ] Implement port auto-detection from nodesets
- [ ] Fallback to coordinate-based search

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

### Mesh Import Helper

| Function | Signature | Description |
|----------|-----------|-------------|
| `create_peec_from_mesh` | `create_peec_from_mesh(mesh_file, cross_section_config=None, auto_config=False)` | Load a GMSH 1D edge mesh and create a `PEECBuilder` with segments for each block |

**Parameters for `cross_section_config`**:

```python
{
    block_id: {
        'width': float,   # Cross-section width [m]
        'height': float,  # Cross-section height [m]
        'sigma': float    # Conductivity [S/m]
    }
}
```

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

1. **Implement `create_peec_from_mesh()` helper function** -- simplify mesh import
2. **Implement racetrack shape** -- most requested for MagLev applications
3. **Implement arc shape** -- useful for partial loops
4. **Test with Cubit-generated meshes** -- verify end-to-end workflow
5. **Add block attribute support** (optional) -- auto cross-section from mesh metadata
6. **Add parametric curve support** -- for CAD integration
7. **Document workflow** -- write user guide

---

## Open Questions

1. Which shapes are highest priority for current projects?
2. Is rectangular cross-section sufficient, or are circular/trapezoidal profiles needed?
3. Is manual port definition acceptable, or is nodeset-based auto-detection critical?
4. Should block attribute support be prioritized over new parametric shapes?
5. Are there other improvements needed in the mesh import pipeline?
