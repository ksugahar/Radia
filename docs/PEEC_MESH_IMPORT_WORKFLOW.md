# PEEC Mesh Import Workflow Design

**Date**: 2026-02-13
**Status**: Implementation Plan

---

## Design Philosophy

**Key Decisions**:
- ✅ **Mesh import** (flexible, handles arbitrary geometry)
- ✅ **Rectangular cross-section only** (simple, practical)
- ✅ **Coreform Cubit for CAD** (professional tool, no custom importers)
- ❌ No parametric shapes (racetrack, spiral, etc.) - unnecessary complexity

---

## Workflow

```
CAD Model (STEP/IGES)
    ↓
Coreform Cubit (geometry + meshing)
    ↓
1D Edge Mesh Export (GMSH format)
    ↓
Radia PEEC (segment creation + solve)
```

---

## Current Implementation (Working)

### Step 1: Cubit Mesh Generation

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
cubit_mesh_export.export_gmsh_v2(cubit, "coil_mesh.msh")
```

### Step 2: Import to Radia PEEC

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

---

## Problems with Current Workflow

### 1. Manual Cross-Section Specification

**Problem**: `width` and `height` must be specified in Python script, not in mesh.

**Solution Options**:

#### Option A: Block Attributes (Recommended)

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

#### Option B: GMSH Physical Group Names

Encode cross-section in physical group name:

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

#### Option C: Separate Configuration File

```yaml
# peec_config.yaml
conductors:
  - block: 1
    name: "coil"
    width: 4.0e-3   # m
    height: 4.0e-3  # m
    sigma: 5.8e7    # S/m
```

**Recommendation**: Use **Option A (Block Attributes)** - cleanest, most robust.

---

### 2. Manual Port Definition

**Problem**: Ports defined by coordinate-based search (error-prone, tedious).

**Current Approach**:
```python
# Find node closest to target position
port_positive_target = np.array([r_mean, 0, 0])
min_dist = float('inf')
for i, tag in enumerate(node_tags):
    dist = np.linalg.norm(coords[i] - port_positive_target)
    if dist < min_dist:
        port_positive_node = tag
```

**Proposed Solution**: Use Cubit **nodesets** to mark port nodes.

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

**Problem**: `cubit_mesh_export` does NOT support nodeset export to GMSH format.

**Workaround**: Use coordinate-based search (current method) until nodeset support is added.

---

### 3. Helper Function Needed

Simplify mesh import with wrapper function:

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

---

## Implementation Plan

### Phase 1: Improve Mesh Import (Current)

- [x] Basic 1D edge mesh import from GMSH
- [x] Manual cross-section specification
- [x] Coordinate-based port search
- [ ] **Add `create_peec_from_mesh()` helper function**

### Phase 2: Block Attributes (Optional)

- [ ] Implement block attribute reading from GMSH
- [ ] Auto-extract cross-section from mesh metadata
- [ ] Fallback to manual specification if not found

### Phase 3: Port Handling (Future)

- [ ] Request nodeset support in `cubit_mesh_export`
- [ ] Implement port auto-detection from nodesets
- [ ] Fallback to coordinate-based search

---

## Example Workflow (Target)

### Cubit Script

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
cubit_mesh_export.export_gmsh_v2(cubit, "induction_coil.msh")
```

### Radia Python Script

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

---

## Next Steps

1. **Implement `create_peec_from_mesh()` helper function** - simplify mesh import
2. **Test with Cubit-generated meshes** - verify workflow
3. **Add block attribute support** (optional) - auto cross-section
4. **Document workflow** - write user guide

---

**Questions**:

1. `create_peec_from_mesh()` ヘルパー関数を実装しますか？ (Implement helper function?)
2. Block属性サポートは必要ですか？ (Need block attribute support?)
3. 他に改善すべき点はありますか？ (Any other improvements needed?)
