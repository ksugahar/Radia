# Mesh and Surface Elements Guide

## 1. What Are Surface Elements?

Netgen/NGSolve meshes can contain two kinds of elements:

- **Volume elements** (体積要素) -- tetrahedra, hexahedra, etc. that fill the interior of a 3-D domain.
- **Surface elements** (表面要素) -- triangles or quadrilaterals that cover the **boundary** of that domain.

Surface elements represent the outer skin of a mesh. They are stored alongside volume elements inside `.vol` files and are the primary data that the Netgen GUI renders.

---

## 2. Why Are They Important?

### Netgen GUI display behaviour

The Netgen GUI renders **surface elements**, not volume elements directly:

- **Surface elements present** -- the mesh boundary is displayed normally.
- **Volume elements only (no surface elements)** -- nothing is displayed, or a clipping plane is required to see anything.

### Viewer selection guide

| Mesh type | Netgen GUI | ParaView | PyVista | webgui |
|-----------|-----------|----------|---------|--------|
| **Surface elements present** | **Recommended** | Overkill | Overkill | Geometry only |
| **Volume elements only** | Cannot display | **Recommended** | **Recommended** | Needs GridFunction |
| **Geometry check** | **Best** | Requires meshing | Requires meshing | OCC direct |
| **Field visualisation** | Not supported | **Best quality** | Fast | Interactive |

---

## 3. FAQ: Are They Always Required?

### Q: Is it a problem that surface elements are needed?

**A: In practice, no.** In every standard mesh-generation workflow surface elements are created automatically:

| Workflow | Surface elements | Reason |
|----------|-----------------|--------|
| **Netgen direct** (`geo.GenerateMesh()`) | Auto | Boundary mesh generated automatically |
| **NGSolve `Mesh()`** | Auto | STEP/OCC import recognises boundaries |
| **Cubit -> `export_netgen()`** | Auto | Cubit sidesets are converted to boundary elements |
| **GMSH -> NGSolve** | Auto | `.msh` files include boundary elements |

**In short, normal mesh generation requires no extra steps.**

### Q: What if I only have volume elements?

There is almost never a reason to intentionally strip surface elements, but if you end up with a volume-only mesh the options are covered in sections 5 and 6 below.

### Q: Do Cubit meshes include surface elements?

**Yes.** When you define a **sideset** in Cubit and export via `export_netgen()`, the sideset surfaces become surface elements:

```python
import cubit
import cubit_mesh_export
from ngsolve import Mesh
from netgen.gui import StartGUI

# Cubit mesh generation
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'model.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Define sidesets (these become surface elements)
cubit.cmd("sideset 1 surface all")
cubit.cmd("sideset 1 name 'boundary'")

# Export to Netgen (surface elements included)
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# Verify
print(f"Volume elements:  {mesh.ngmesh.ne}")
print(f"Surface elements: {mesh.ngmesh.nse}")  # should be > 0

# Display in Netgen GUI
StartGUI()
mesh.ngmesh.Draw()
```

**Key point**: defining a sideset in Cubit is what creates the corresponding surface elements on export.

---

## 4. Auto-Generation in Standard Workflows

### NGSolve sample meshes

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

**Conclusion**: all NGSolve sample `.vol` files display correctly in the Netgen GUI.

### Recommended standard workflow

```
CAD (STEP) -> Netgen / Cubit -> Mesh generation -> .vol file
                                                      |
                                              Surface elements
                                            generated automatically
                                                      |
                                                Netgen GUI OK
```

---

## 5. Netgen GUI: Limitations & Workarounds

### When Netgen GUI is the right tool

- The mesh contains surface elements.
- You need to inspect geometry or mesh quality.
- You want a lightweight, fast viewer.
- You are following an integrated shape -> mesh -> review workflow.

### When to use ParaView / PyVista instead

- The mesh contains only volume elements.
- You need to visualise internal structure via slicing or clipping.
- You need to visualise field data (B, H, etc.).
- You need publication-quality figures.

### Workaround 1: ParaView slice / clip (recommended for volume-only meshes)

```python
# 1. Export to VTS
import radia as rad
rad.FldUnits('m')
rad.FldVTS(magnet, 'field.vts', ...)

# 2. Open in ParaView
paraview field.vts

# 3. Filters > Slice
#    - Origin: [0, 0, 0.05]
#    - Normal: [0, 0, 1]
#    - Apply

# 4. Filters > Clip
#    - Clip Type: Plane
#    - Normal: [0, 0, 1]
#    - Apply
```

### Workaround 2: PyVista slice / clip

```python
import pyvista as pv

# Load VTS
grid = pv.read('field.vts')

# Create slice at z=0.05 m
slice_z = grid.slice(normal='z', origin=[0, 0, 0.05])
slice_z.plot(scalars='B_magnitude', cmap='coolwarm')

# Or clip half of the domain
clipped = grid.clip(normal='z', origin=[0, 0, 0])
clipped.plot(scalars='B_magnitude', cmap='viridis')
```

### Workaround 3: Regenerate mesh with surface elements via NGSolve

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

### Workaround 4: Add surface elements to an existing volume-only mesh

```python
from ngsolve import *

# Load volume-only mesh
mesh = Mesh('volume_only.vol')

# NGSolve recognises boundaries automatically
# (the original mesh must still contain boundary information)

# Re-export -- surface elements will be included
mesh.ngmesh.Save('with_surface.vol')
```

---

## 6. Troubleshooting

### How to check whether surface elements are present

```bash
python utils/check_vol_surface_elements.py mesh.vol
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

### Netgen GUI shows nothing -- possible causes

| Cause | Diagnosis | Fix |
|-------|-----------|-----|
| **No surface elements** | `mesh.nse == 0` | Regenerate mesh or use ParaView/PyVista (see Section 5) |
| **Mesh too small / too large for viewport** | Elements exist but view is empty | Mouse-wheel zoom, or *View > Center* |
| **Corrupt `.vol` file** | Load raises an error or counts are unexpected | Regenerate the mesh from geometry |

To regenerate a simple test mesh and confirm the GUI works:

```python
from netgen.occ import Box, Pnt, OCCGeometry

geo = OCCGeometry(Box(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
mesh = geo.GenerateMesh(maxh=0.2)
mesh.Save('test.vol')
```

### Quick-reference checklist

- [ ] Run `check_vol_surface_elements.py` to confirm surface element count.
- [ ] Identify the mesh source (Netgen / Cubit / GMSH).
- [ ] If Cubit, verify that sidesets are defined before export.
- [ ] Test with an NGSolve sample mesh to rule out environment issues.
- [ ] If none of the above helps, switch to ParaView or PyVista.

### Practical summary

| Situation | Surface elements | Action |
|-----------|-----------------|--------|
| **Netgen-generated mesh** | Auto-generated | Nothing to do |
| **Cubit -> Netgen** | Auto-converted | Define sidesets |
| **GMSH -> NGSolve** | Auto-converted | Nothing to do |
| **NGSolve samples** | All included | Nothing to do |
| **Volume-only (rare)** | None | Use ParaView / PyVista |

---

**Created**: 2026-02-12
**Applies to**: Netgen GUI, NGSolve `.vol` files, NGSolve mesh workflows
