# PEEC + Magnetic Core: Practical Workflow

## Quick Start: Wire + Ferrite Core

```python
import radia as rad
import numpy as np
from peec_matrices import PEECBuilder
from peec_coupled import CoupledPEECSolver

MU_0 = 4e-7 * np.pi

# Step 1: Clear Radia state
rad.UtiDelAll()

# Step 2: Define conductor (PEEC)
builder = PEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)          # 100mm wire
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)  # 1x1mm Cu
builder.add_port(n1, n2)
topo = builder.build_topology()

# Step 3: Define magnetic core (Radia hexahedra)
core_verts = [
    [0.02, 0.005, -0.005], [0.08, 0.005, -0.005],
    [0.08, 0.015, -0.005], [0.02, 0.015, -0.005],
    [0.02, 0.005,  0.005], [0.08, 0.005,  0.005],
    [0.08, 0.015,  0.005], [0.02, 0.015,  0.005],
]
core = rad.ObjHexahedron(core_verts, [0, 0, 0])
rad.MatApl(core, rad.MatLin(999))  # mu_r = 1000

# Step 4: Compute coupling
solver = CoupledPEECSolver(topo, [core])
solver.compute_coupling_matrix(solver_method=0, mu_r_real=1000)

# Step 5: Extract impedance
Z = solver.compute_port_impedance(1e6)  # at 1 MHz
L_total = np.imag(Z) / (2 * np.pi * 1e6)
print(f"L = {L_total*1e9:.2f} nH")
```

## Core Subdivision

A single hexahedral element is often too coarse. Subdivide for accuracy:

```python
# 3x1x1 subdivision
inp = """\
.Units mm
.default sigma=5.8e7
N1 x=0 y=0 z=0
N2 x=100 y=0 z=0
E1 N1 N2 w=1 h=1
.external N1 N2
.magnetic
  type=box
  center=50,10,0
  size=60,10,10
  divisions=3,1,1
  mu_r=1000
.endmagnetic
.freq fmin=1e3 fmax=10e6 ndec=10
.end
"""
from fasthenry_parser import FastHenryParser
parser = FastHenryParser()
parser.parse_string(inp)
result = parser.solve()
```

## Common Pitfalls and Solutions

### 1. Coordinates must be in meters

Radia always uses meters. A 60mm core is `0.06`, not `60`.

```python
# WRONG
core = rad.ObjHexahedron([[20,5,-5], ...], [0,0,0])  # millimeters!

# CORRECT
core = rad.ObjHexahedron([[0.02, 0.005, -0.005], ...], [0,0,0])  # meters
```

### 2. Call `rad.UtiDelAll()` before creating objects

Radia keeps global state. Previous objects persist and interfere.

```python
rad.UtiDelAll()  # ALWAYS call first
core = rad.ObjHexahedron(verts, [0,0,0])
```

### 3. Hexahedron vertex ordering matters

MMMBuilder face ordering:
```
Face 0: bottom (-Z)  v0,v3,v2,v1
Face 1: top    (+Z)  v4,v5,v6,v7
Face 2: front  (-Y)  v0,v1,v5,v4
Face 3: back   (+Y)  v2,v3,v7,v6
Face 4: left   (-X)  v0,v4,v7,v3
Face 5: right  (+X)  v1,v2,v6,v5
```

Standard hex vertex order: bottom face CCW (v0-v3), top face CCW (v4-v7).

### 4. NGSBEM: Use surface-only mesh for BEM

Volume meshes cause ill-conditioning in the BEM single-layer operator.

```python
from netgen.occ import Box, Pnt, Glue, OCCGeometry

wire = Box(Pnt(0, -0.5e-3, -0.5e-3), Pnt(0.1, 0.5e-3, 0.5e-3))
# WRONG: volume mesh
# geo = OCCGeometry(wire)

# CORRECT: surface-only mesh
geo = OCCGeometry(Glue(wire.faces))  # No volume, surface only
mesh = Mesh(geo.GenerateMesh(maxh=0.5e-3))
```

### 5. NGSBEM: Set maxh proportional to smallest dimension

For a 1mm x 1mm cross-section wire, use `maxh <= 0.5mm`:

```python
# WRONG: elements are 5mm x 1mm (aspect ratio 5:1)
geo.GenerateMesh(maxh=0.005)

# CORRECT: elements are ~0.5mm (nearly equilateral)
geo.GenerateMesh(maxh=0.0005)
```

### 6. MSC system matrix sign convention

The MSC equation is `(1/chi + N) sigma = H_ext`, NOT `(-1/chi - N) sigma = H_ext`.

```python
# System matrix
K_msc = np.diag(inv_chi) + N   # CORRECT
# K_msc = -np.diag(inv_chi) - N  # WRONG (old sign convention)
```

### 7. Yano-Sugahara evaluation point

For MSC hexahedra, evaluate at midpoint between face center and element center:

```python
eval_point = (face_center + element_center) / 2.0  # CORRECT
# eval_point = face_center  # WRONG (singular self-term)
```

### 8. Point charge correction for multi-element

Without point charge correction, multi-element (3x3x3) accuracy degrades from
0.03% to 650% error.

### 9. Loop port definition

For a closed loop, don't use `add_port(n1, n1)`. Instead, split the loop
at one point and use two separate nodes:

```python
n1 = builder.add_node_at(x, y, z)      # port node A
n1b = builder.add_node_at(x, y, z)     # port node B (same location!)
# ... connect segments: n1b -> n2 -> n3 -> n4 -> n1
builder.add_port(n1, n1b)  # port across the gap
```

## Solver Selection Flowchart

```
Is core conducting (sigma > 0)?
  ├─ No → Is core nonlinear?
  │    ├─ Yes → Radia MSC (core_model='radia')
  │    └─ No  → Radia MSC (simplest, fastest)
  └─ Yes → Is mu_r > 1?
       ├─ No  → Scalar FEM-BEM (core_model='fembem')
       └─ Yes → Vector FEM-BEM (core_model='vector_fembem')
```

## Validated Accuracy

| Test Case | MMMBuilder/Radia | DOF |
|-----------|-----------------|-----|
| 1 hex (6 DOF) | 1.000002 | 6 |
| 3 hex (18 DOF) | 1.000002 | 18 |
| 9 hex (54 DOF) | 1.000003 | 54 |
| 27 hex (162 DOF) | 0.999685 | 162 |

PEEC vs NGSBEM (air inductance): 4.5% systematic difference
(filament approximation vs surface current distribution).
