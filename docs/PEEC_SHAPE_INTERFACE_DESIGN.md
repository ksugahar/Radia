# PEEC Conductor Shape Interface Design

**Date**: 2026-02-13
**Status**: Proposal

---

## Current Interface (Implemented)

```python
from peec_matrices import PEECBuilder

builder = PEECBuilder()
builder.create_wire([0,0,0], [0.1,0,0], 1e-3, 1e-3, 10)  # Straight wire
builder.create_loop([0,0,0], 0.05, [0,0,1], 2e-3, 2e-3, 36)  # Circle
L, R, P, M_LS = builder.build()
```

---

## Proposed Extensions

### 1. Basic Shapes (High Priority)

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

---

### 2. Path-Based Shapes (Medium Priority)

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

---

### 3. Mesh Import (Already Working)

```python
# From GMSH 1D edge mesh
import gmsh
gmsh.initialize()
gmsh.open("coil_mesh.msh")

# Extract edges and create segments
node_tags, coords, _ = gmsh.model.mesh.getNodes()
elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

for edge in edges:
    n0, n1 = edge
    p0 = coords[n0]
    p1 = coords[n1]
    builder.create_wire(p0, p1, width, height, n_segments=1)

L, R, P, M_LS = builder.build()
```

---

### 4. Composite Shapes (Low Priority)

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

---

## Implementation Priority

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

## Design Questions

1. **Cross-section specification**:
   - Current: Rectangular (width, height)
   - Future: Circular? Trapezoidal? Profile function?

2. **Segmentation control**:
   - Current: `n_segments` (uniform)
   - Future: Adaptive? Grading?

3. **Port definition**:
   - Current: Manual (find nodes by coordinates)
   - Future: Automatic? Named ports?

4. **Units**:
   - Current: SI meters (via `rad.FldUnits('m')`)
   - Keep consistent with Radia convention

---

## Next Steps

1. **User feedback**: Which shapes are most needed?
2. **Implement racetrack**: Most requested for MagLev applications
3. **Implement arc**: Useful for partial loops
4. **Add parametric curve**: For CAD integration

---

**Questions for User**:

1. どの形状が最優先ですか？ (Which shapes are highest priority?)
2. 断面形状は矩形のみで十分ですか？ (Is rectangular cross-section sufficient?)
3. ポート定義は手動で良いですか？ (Is manual port definition OK?)
4. メッシュインポート vs パラメトリック形状、どちらが重要？ (Mesh import vs parametric shapes - which is more important?)
