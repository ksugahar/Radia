# Mesh Guide: Cubit and Netgen Workflows for Radia

## Overview

This guide consolidates the mesh generation workflows for Radia. All mesh generation uses Coreform Cubit. The only input format for NGSolve/Radia is Netgen `.vol` (via `export netgen`). GMSH is used for visualization only.

```
+-----------------------------------------------------------------+
|                    CAD -> Mesh -> Radia Workflow                 |
+-----------------------------------------------------------------+
|                                                                  |
|  CAD (STEP/IGES) -> Cubit -> export netgen -> .vol -> Radia|
|                                                                  |
|  Mesh types by application:                                      |
|    - Magnetic materials:        Volume mesh (Tet4, Hex8)         |
|    - Conductors (PEEC):            Surface mesh only (Tri3, Quad4)|
|                                                                  |
|  Mesh file formats:                                              |
|    - GMSH:   .msh -> NGSolve -> Radia                            |
|    - Netgen:  .vol -> NGSolve -> Radia                            |
|    - Cubit:  export netgen -> NGSolve -> Radia       |
|                                                                  |
+-----------------------------------------------------------------+
```

---

## 1. Mesh Types

### 1.1 Volume Elements

Volume elements -- tetrahedra, hexahedra, wedges -- fill the interior of a 3-D domain. They are required for magnetic material modelling (permanent magnets, soft magnetic materials).

**Requirement**: Volume Mesh

| Element Type | GMSH Element | Radia API | Use Case |
|----------|---------|----------|------|
| Tetrahedron | Tet4 | `ObjTetrahedron()` | Complex shapes |
| Hexahedron | Hex8 | `ObjHexahedron()` | Structured grid |
| Wedge/Prism | Wedge6 | `ObjWedge()` | Transition elements |

**GMSH generation**:
```python
gmsh.model.mesh.generate(3)  # 3D volume mesh
```

### 1.2 Surface Elements

Surface elements -- triangles or quadrilaterals -- cover the **boundary** of a 3-D domain. They are the outer skin of a mesh and serve two purposes:

1. **PEEC conductors**: Surface mesh is all that is needed for surface-current modelling.
2. **Netgen GUI display**: The Netgen GUI renders surface elements, not volume elements directly.

| Element Type | GMSH Element | Use Case |
|----------|---------|------|
| Triangle | Tri3 | Surface current distribution / Boundary display |
| Quadrilateral | Quad4 | Surface current distribution / Boundary display |

**GMSH generation for PEEC**:
```python
gmsh.model.mesh.generate(2)  # 2D surface mesh only
```

**Important**: PEEC uses a surface current model, so **volume mesh is not required**

**Reason**:
- Skin effect: Handled by SIBC (Surface Impedance Boundary Condition)
- Conductor interior: Current density decays exponentially (represented by surface impedance)
- Computational efficiency: Surface only provides sufficient accuracy

### 1.3 Auto-Generation in Standard Workflows

In every standard mesh-generation workflow surface elements are created automatically:

| Workflow | Surface elements | Reason |
|----------|-----------------|--------|
| **Netgen direct** (`geo.GenerateMesh()`) | Auto | Boundary mesh generated automatically |
| **NGSolve `Mesh()`** | Auto | STEP/OCC import recognises boundaries |
| **Cubit -> `export netgen`** | Auto | Cubit sidesets are converted to boundary elements |
| **GMSH -> NGSolve** | Auto | `.msh` files include boundary elements |

**In short, normal mesh generation requires no extra steps.**

```
CAD (STEP) -> Netgen / Cubit / GMSH -> Mesh generation -> .vol / .msh file
                                                             |
                                                     Surface elements
                                                   generated automatically
                                                             |
                                                       Netgen GUI OK
```

### 1.4 NGSolve Sample Meshes

Every mesh shipped under `share/ngsolve/` already contains surface elements:

| File | Volume elements | Surface elements | Netgen GUI |
|------|----------------|-----------------|-----------|
| cube.vol | 756 | 338 (Triangle) | OK |
| coil.vol | 1709 | Present | OK |
| coilshield.vol | 1798 | 376 (Tri+Quad) | OK |
| beam.vol | 31 | Present | OK |
| shaft.vol | 1622 | Present | OK |
| chip.vol | 0 | Present (Surface-only) | OK |
| doubleglazing.vol | 0 | Present (Surface-only) | OK |
| square.vol | 0 | Present (Surface-only) | OK |

All NGSolve sample `.vol` files display correctly in the Netgen GUI.

---

## 2. Cubit Export Workflows

All mesh generation uses Coreform Cubit. The **only** input format for NGSolve/Radia is Netgen `.vol`.
GMSH `.msh` is used only for visualization output (not as input to NGSolve).

### 2.1 Workflow 1: Magnetic Material (Volume Mesh)

#### Cubit -> .vol -> Radia

```python
import cubit
from ngsolve import Mesh
from radia.netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

# Cubit initialization
cubit.init(['cubit', '-nojournal', '-batch'])

# CAD file loading
cubit.cmd('import step "core.step" heal')

# Mesh generation
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.005")
cubit.cmd("mesh volume all")

# Block definition (becomes material label in .vol)
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "core"')

# Export to .vol (order 3 curved)
cubit.cmd('export netgen "core.vol" order 3 overwrite')

# Load in NGSolve -> Radia
mesh = Mesh("core.vol")
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='core')

# Apply material and solve
mat = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(mag_obj, mat)
rad.Solve(mag_obj, 0.0001, 1000, 1)
```

### 2.2 Workflow 2: Conductor (Surface Mesh / PEEC)

For PEEC conductor analysis, Cubit generates the surface mesh and exports via `.vol`.
The BEM solver uses BND elements only; volume elements are ignored.

```python
# Cubit surface mesh for PEEC
cubit.cmd('import step "coil.step" heal')
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.001")
cubit.cmd("mesh volume 1")

cubit.cmd("block 1 add volume 1")
cubit.cmd('block 1 name "coil"')

# Terminal faces: sidesets -> boundary labels in .vol
cubit.cmd("sideset 1 add surface 3")
cubit.cmd('sideset 1 name "source"')
cubit.cmd("sideset 2 add surface 5")
cubit.cmd('sideset 2 name "sink"')

cubit.cmd('export netgen "coil.vol" order 2 overwrite')
```

For simple coils without surface mesh, use analytical current sources:

```python
import numpy as np

# Circular coil: R=50mm, 1mm cross-section, J=1e6 A/m^2
coil = rad.ObjArcCur([0, 0, 0], [0.0495, 0.0505],
                     [-np.pi, np.pi], 0.001, 100, 1e6)
```

### 2.3 Workflow 3: Combined Model (Magnetic Material + Coil)

#### Example: Electromagnet (Iron Core + Coil)

```python
import cubit
import numpy as np
from ngsolve import Mesh
from radia.netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

rad.UtiDelAll()

# ===============================
# 1. Iron Core (Volume Mesh via Cubit)
# ===============================
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd('import step "core.step" heal')
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.005")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")
cubit.cmd('block 1 name "core"')
cubit.cmd('export netgen "core.vol" order 3 overwrite')

mesh_core = Mesh("core.vol")
core_obj = netgen_mesh_to_radia(mesh_core,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m')
mat_iron = rad.MatLin(1000)
rad.MatApl(core_obj, mat_iron)

# ===============================
# 2. Coil (Analytical Source)
# ===============================
coil_obj = rad.ObjArcCur([0, 0, 0], [0.0495, 0.0505],
                         [-np.pi, np.pi], 0.001, 100, 1e6)

# ===============================
# 3. Combine and Solve
# ===============================
container = rad.ObjCnt([core_obj, coil_obj])
rad.Solve(container, 0.0001, 1000, 1)

# Field calculation
B = rad.Fld(container, 'b', [0, 0, 0.1])
print(f"Field at (0, 0, 0.1): {B} T")
```

---

## 3. SetGeomInfo API (High-Order Curving)

> **Source**: [ksugahar/ngsolve](https://github.com/ksugahar/ngsolve) fork with SetGeomInfo API (netgen PR [#232](https://github.com/NGSolve/netgen/pull/232)).

### 3.1 Problem Statement

When meshes are imported from external mesh generators (Gmsh, Cubit, etc.) without geometry, `mesh.Curve(order)` fails because the UV parametric coordinates (geominfo) are not set. The `SetGeomInfo` API enables setting geominfo programmatically.

### 3.2 API

```python
Element2d.SetGeomInfo(vertex_index, u, v, trignum=0)
```

**Parameters:**
- `vertex_index`: 0-based index of the vertex within the element
- `u`, `v`: Surface parametric coordinates from the OCC geometry
- `trignum`: Triangle number for STL meshing (default: 0)

### 3.3 Recommended Workflow: Coreform Cubit + Name-based Mapping

```
1. OCC: Create geometry and name faces (name_occ_faces)
2. OCC: Export to STEP (face names preserved)
3. Cubit: Import STEP, generate mesh
4. Export: Use name-based face mapping to Netgen mesh
5. SetGeomInfo: Compute UV parameters analytically
6. mesh.Curve(order): High-order curving works correctly!
```

### 3.4 Code Example

```python
import cubit
from ngsolve import Mesh

cubit.init(['cubit', '-nojournal', '-batch'])

# 1. Import geometry into Cubit
cubit.cmd('import step "geometry.step" heal')

# 2. Mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size 0.15")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# 3. Export with high-order curving (NetgenCurver + ACIS projection)
cubit.cmd('export netgen "mesh.vol" order 2 overwrite')

# 4. Load in NGSolve — high-order nodes embedded in .vol
mesh = Mesh("mesh.vol")
```

### 3.5 Automatic Curving

`export netgen` handles all high-order curving automatically via
NetgenCurver + ACIS CallbackGeometry. No SetGeomInfo, no STEP file, no
`mesh.Curve()` call needed:

```python
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
mesh = Mesh("mesh.vol")  # already curved to order 3
```

### 3.6 Accuracy Results

| Geometry | Curve(2) Error | Curve(3) Error |
|----------|----------------|----------------|
| Complex (Boolean ops) | **0.0021%** | **0.0004%** |
| Cylinder | 0.0027% | 0.0006% |
| Sphere | 0.0027% | 0.0004% |
| Torus | 0.0010% | 0.0003% |

All results achieve **Netgen-native accuracy** (<0.003%).

### 3.7 Requirements

- NGSolve: Build from `ksugahar/ngsolve` branch `feature/setgeominfo`
- Coreform Cubit 2025.12+
- `pip install "radia[cubit]" && cubit-plugin-install`

### 3.8 Examples and Links

- Netgen PR: [NGSolve/netgen#232](https://github.com/NGSolve/netgen/pull/232)
- Forum: [Feature Request - SetGeomInfo API](https://forum.ngsolve.org/t/feature-request-python-api-for-high-order-curving-of-externally-imported-meshes/3810)
- PyPI: [cubit-mesh-export](https://pypi.org/project/cubit-mesh-export/)

---

## 4. Surface Elements: Display and Workarounds

### 4.1 Netgen GUI Display Behaviour

The Netgen GUI renders **surface elements**, not volume elements directly:

- **Surface elements present** -- the mesh boundary is displayed normally.
- **Volume elements only (no surface elements)** -- nothing is displayed, or a clipping plane is required to see anything.

### 4.2 Viewer Selection Guide

| Mesh type | Netgen GUI | ParaView | PyVista | webgui |
|-----------|-----------|----------|---------|--------|
| **Surface elements present** | **Recommended** | Overkill | Overkill | Geometry only |
| **Volume elements only** | Cannot display | **Recommended** | **Recommended** | Needs GridFunction |
| **Geometry check** | **Best** | Requires meshing | Requires meshing | OCC direct |
| **Field visualisation** | Not supported | **Best quality** | Fast | Interactive |

### 4.3 When Netgen GUI is the right tool

- The mesh contains surface elements.
- You need to inspect geometry or mesh quality.
- You want a lightweight, fast viewer.
- You are following an integrated shape -> mesh -> review workflow.

### 4.4 When to use ParaView / PyVista instead

- The mesh contains only volume elements.
- You need to visualise internal structure via slicing or clipping.
- You need to visualise field data (B, H, etc.).
- You need publication-quality figures.


```python
# 1. Export to VTS
import radia as rad

# 2. Open in ParaView

# 3. Filters > Slice
#    - Origin: [0, 0, 0.05]
#    - Normal: [0, 0, 1]
#    - Apply

# 4. Filters > Clip
#    - Clip Type: Plane
#    - Normal: [0, 0, 1]
#    - Apply
```


```python

# Load VTS

# Create slice at z=0.05 m
slice_z = grid.slice(normal='z', origin=[0, 0, 0.05])
slice_z.plot(scalars='B_magnitude', cmap='coolwarm')

# Or clip half of the domain
clipped = grid.clip(normal='z', origin=[0, 0, 0])
clipped.plot(scalars='B_magnitude', cmap='viridis')
```

### 4.7 Workaround: Regenerate mesh with surface elements via NGSolve

```python
from ngsolve import *
from netgen.occ import Box, Pnt, OCCGeometry

# Create geometry
box = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
geo = OCCGeometry(box)

# Generate mesh (includes surface elements)
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Export to .vol (surface elements included)
mesh.ngmesh.Save('mesh_with_surface.vol')

# Open in Netgen GUI
from netgen.gui import StartGUI
StartGUI()
mesh.ngmesh.Draw()
```

### 4.8 Workaround: Add surface elements to an existing volume-only mesh

```python
from ngsolve import *

# Load volume-only mesh
mesh = Mesh('volume_only.vol')

# NGSolve recognises boundaries automatically
# (the original mesh must still contain boundary information)

# Re-export -- surface elements will be included
mesh.ngmesh.Save('with_surface.vol')
```

### 4.9 Cubit Meshes and Surface Elements

When you define a **sideset** in Cubit and export via `export netgen`, the sideset surfaces become surface elements:

```python
import cubit
import cubit
from ngsolve import Mesh
from netgen.gui import StartGUI

# Cubit mesh generation
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'model.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")
cubit.cmd("block 1 add volume all")

# Define sidesets (these become boundary labels)
cubit.cmd("sideset 1 surface all")
cubit.cmd("sideset 1 name 'boundary'")

# Export to Netgen .vol
cubit.cmd('export netgen "model.vol" order 2 overwrite')
mesh = Mesh("model.vol")

# Verify
print(f"Volume elements:  {mesh.ngmesh.ne}")
print(f"Surface elements: {mesh.ngmesh.nse}")  # should be > 0

# Display in Netgen GUI
StartGUI()
mesh.ngmesh.Draw()
```

**Key point**: defining a sideset in Cubit is what creates the corresponding surface elements on export.

---

## 5. Troubleshooting

### 5.1 How to check whether surface elements are present

```bash
python scripts/check_vol_surface_elements.py mesh.vol
```

Or programmatically:

```python
from netgen.meshing import Mesh

mesh = Mesh()
mesh.Load('mesh.vol')

print(f"Volume elements:  {mesh.ne}")
print(f"Surface elements: {mesh.nse}")

if mesh.nse == 0:
    print("Warning: No surface elements")
    print("Netgen GUI may not display this mesh")
```

**Example output** of the check script:

```
Analyzing: mesh.vol
============================================================

Mesh Statistics:
  Vertices:        228
  Volume elements: 756
  Surface elements: 338

  Volume element types:
    Tet: 756

  Surface element types:
    Triangle: 338

============================================================
Display Compatibility:
============================================================

  Netgen GUI: COMPATIBLE
   - Surface elements present: 338
   - Mesh will be displayed as surface
   - Recommended viewer: Netgen GUI
```

### 5.2 Netgen GUI shows nothing -- possible causes

| Cause | Diagnosis | Fix |
|-------|-----------|-----|
| **No surface elements** | `mesh.nse == 0` | Regenerate mesh or use ParaView/PyVista (see Section 4) |
| **Mesh too small / too large for viewport** | Elements exist but view is empty | Mouse-wheel zoom, or *View > Center* |
| **Corrupt `.vol` file** | Load raises an error or counts are unexpected | Regenerate the mesh from geometry |

To regenerate a simple test mesh and confirm the GUI works:

```python
from netgen.occ import Box, Pnt, OCCGeometry

geo = OCCGeometry(Box(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
mesh = geo.GenerateMesh(maxh=0.2)
mesh.Save('test.vol')
```

### 5.3 Quick-reference checklist

- [ ] Run `check_vol_surface_elements.py` to confirm surface element count.
- [ ] Identify the mesh source (Netgen / Cubit / GMSH).
- [ ] If Cubit, verify that sidesets are defined before export.
- [ ] Test with an NGSolve sample mesh to rule out environment issues.
- [ ] If none of the above helps, switch to ParaView or PyVista.

### 5.4 Practical summary

| Situation | Surface elements | Action |
|-----------|-----------------|--------|
| **Netgen-generated mesh** | Auto-generated | Nothing to do |
| **Cubit -> Netgen** | Auto-converted | Define sidesets |
| **GMSH -> NGSolve** | Auto-converted | Nothing to do |
| **NGSolve samples** | All included | Nothing to do |
| **Volume-only (rare)** | None | Use ParaView / PyVista |

---

## 6. Tool Comparison (GMSH vs Netgen vs Cubit)

| Aspect | GMSH | Netgen | Coreform Cubit |
|------|------|--------|----------------|
| **CAD Import** | STEP/IGES direct | STEP/OCC | STEP/IGES direct |
| **License** | Open source | Open source | Commercial |
| **NGSolve Integration** | Direct .msh import | Native | `export netgen` |
| **2D/Axisymmetric** | Supported | 3D only recommended | Supported |
| **Surface Mesh** | `generate(2)` | Auto-generated | Auto via sideset |
| **Volume Mesh** | Tet/Hex supported | Tet (Hex via external tools) | Tet/Hex supported |
| **Hexahedral Mesh** | Structured grid only | Not supported (external tools) | High quality (recommended) |
| **High-order curving** | Not supported (external processing) | Native | Via SetGeomInfo API |
| **Visualization** | GMSH GUI | Netgen GUI | Cubit GUI |

**Recommended**:
- **All mesh generation**: Coreform Cubit (tet, hex, wedge, pyramid, boundary layers)
- **High-order curving**: `export netgen` (order 1-5 via ACIS CallbackGeometry)
- **Simple test geometries**: Netgen OCC (code generation, automatic meshing)
- **GMSH**: Visualization only (not mesh generation). Use `export gmsh` for export.

### GMSH Role in Radia

GMSH is used **only** for visualization and post-processing:
- View mesh exported via `export gmsh "mesh.msh"`
- View field results exported via `GmshPostExport.write("results.msh")`
- GMSH is NOT used for mesh generation (GmshBuilder was removed)

---

## 7. FAQ

### Q1: Is volume mesh not needed for PEEC?

**A: Not required.** PEEC uses a surface current approximation.

**Reason**:
1. **Skin effect**: At high frequencies, current concentrates at the surface
2. **SIBC**: Surface impedance represents the current distribution inside the conductor
3. **Computational efficiency**: Surface mesh alone provides sufficient accuracy

**Applicable range**: When frequency x size is larger than the skin depth

### Q2: Can GMSH generate hexahedral meshes?

**A: GMSH is not used for mesh generation in Radia.** Use Coreform Cubit for all mesh generation, including hex meshing. Cubit supports `scheme map`, `scheme sweep`, `scheme tetmesh`, boundary layers, and webcut decomposition for complex geometries.

### Q3: Surface elements are always required?

**A: In practice, no extra steps are needed.** In every standard mesh-generation workflow (Netgen, Cubit, GMSH), surface elements are created automatically. The only edge case is a volume-only mesh, which is rare. See Section 5 for troubleshooting if this occurs.

### Q4: What if `mesh.Curve(order)` fails on an imported mesh?

**A:** Use `export netgen "mesh.vol" order N` which handles high-order curving automatically via NetgenCurver + ACIS CallbackGeometry. No `mesh.Curve()` call or SetGeomInfo API needed.

---

## 8. References

### Sample Scripts

| File | Description |
|---------|------|
| `docs/cubit_mesh_export/hex_sphere_highorder/hex_sphere_curved_ngsolve.py` | Curved high-order `.vol` mesh loaded through NGSolve |
| `docs/peec_integration/demos/demo_gmsh_surface_mesh.py` | Surface mesh (PEEC conductor) |
| `docs/visualization/_gmsh_display.geo` | Minimal GMSH display companion |

### Recommended Radia Workflows (Summary)

```
Magnetic materials (permanent magnets / iron cores):
  CAD -> Cubit -> export netgen "mesh.vol" -> NGSolve Mesh() -> netgen_mesh_to_radia()

Conductors (coils / shields):
  Analytical: rad.ObjArcCur(), rad.ObjRaceTrk(), rad.ObjFlmCur()
  PEEC: CAD -> Cubit -> export netgen "coil.vol" -> PEEC solver

Combined model (electromagnets, etc.):
  Iron mesh + Coil source -> rad.ObjCnt() -> rad.Solve()

High-order curving (Cubit):
  CAD -> Cubit -> export netgen "mesh.vol" order N -> Mesh("mesh.vol")
  (NetgenCurver + ACIS projection, no mesh.Curve() needed)
```

### Key Points

1. **Cubit**: All mesh generation (tet, hex, wedge, pyramid, BL)
2. **Mesh types**: Magnetic material = volume, Conductor = surface (PEEC) or analytical
3. **Input format**: `.vol` is the only input to NGSolve/Radia
4. **High-order curving**: Automatic via `export netgen` + NetgenCurver + ACIS
5. **GMSH**: Visualization only (view mesh via `export gmsh`, view fields via `GmshPostExport`)

### External Links

- Netgen PR #232: [NGSolve/netgen#232](https://github.com/NGSolve/netgen/pull/232)
- SetGeomInfo Forum: [Feature Request - SetGeomInfo API](https://forum.ngsolve.org/t/feature-request-python-api-for-high-order-curving-of-externally-imported-meshes/3810)
- radia Cubit plugin: see `src/cubit_plugin/` (replaces the old `cubit_mesh_export` PyPI package)
- [PEEC_INTEGRATION.md](PEEC_INTEGRATION.md) (to be created in the future)

---

**Created**: 2026-02-22
**Scope**: Radia Mesh Generation Workflows (GMSH, Netgen, Cubit)
