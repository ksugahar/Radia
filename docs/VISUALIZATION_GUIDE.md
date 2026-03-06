# Visualization Guide: Radia-NGSolve Viewer Selection and Workflows

## Overall Flow

```
┌─────────────────────────────────────────────────────────────────┐
│                    Radia-NGSolve Workflow                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  [1] CAD Modeling                                               │
│      Coreform Cubit / FreeCAD / STEP files                      │
│                          ↓                                       │
│  [2] Mesh Generation                                            │
│      Netgen / Cubit → .vol file (surface elements auto-generated)│
│      GMSH → .msh file                                           │
│                          ↓                                       │
│  [3] Visualization (by purpose)                                  │
│      ├─ .msh inspection: GMSH GUI (native, fastest)             │
│      ├─ Shape inspection: Netgen GUI (lightweight, accurate)    │
│      ├─ Development inspection: PyVista (rapid)                 │
│      ├─ Publication figures: ParaView (high quality)            │
│      └─ Integrated inspection: webgui (browser)                 │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**Important**: Jupyter notebooks (.ipynb) do not work well with Claude Code, so regular Python scripts (.py) are recommended. Both PyVista and webgui can be fully used from .py scripts.

### Jupyter Notebook vs Python Script

**Reasons to prefer Python Scripts (.py)**:

| Aspect | Python Script (.py) | Jupyter Notebook (.ipynb) |
|------|-------------------|-------------------------|
| **Claude Code editing** | ✅ Easy (plain text) | ❌ Difficult (JSON structure) |
| **Version control** | ✅ Clear Git diffs | ❌ Hard to read JSON diffs |
| **File size** | ✅ Small | ❌ Bloated with execution results |
| **Debugging** | ✅ Standard debugger available | ❌ Cell-based only |
| **Automation** | ✅ Easy to run in CI/CD | ❌ Additional setup required |
| **PyVista** | ✅ Interactive window | ⚠️ Inline display (limited) |
| **webgui** | ✅ Auto-launches browser | ⚠️ In-notebook display |

**Conclusion**: Unless there is a specific reason, use regular Python scripts (.py).

```python
# demo.py - Regular Python script
import pyvista as pv

# PyVista: Opens an interactive window (no Jupyter needed)
grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude', cmap='coolwarm')  # ← Opens a window
```

```python
# webgui_demo.py - Regular Python script
from ngsolve.webgui import Draw

# NGSolve webgui: Automatically opens in browser (no Jupyter needed)
Draw(B_gf, mesh, 'B_field')  # ← Opens a browser tab
```

---

## Viewer Comparison Table

| Viewer | Quality | Speed | Purpose | Installation |
|----------|------|------|------|------------|
| **GMSH GUI** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | **Best viewer for .msh (top priority)** | Bundled with GMSH |
| **PyVista** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Python integration (default) | `pip install pyvista` |
| **ParaView** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | **Publication figures (highest quality)** | Standalone application |
| **NGSolve webgui** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Radia field integration | `pip install ngsolve` |
| **Netgen GUI** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | .vol format / shape and mesh inspection | Bundled with NGSolve |

### Viewer Selection Flowchart

```
What is the purpose?
  │
  ├─ .msh mesh inspection and debugging (during development)
  │    → GMSH GUI ⭐⭐⭐⭐⭐ (top priority)
  │
  ├─ Python integration and automation
  │    → PyVista ⭐⭐⭐⭐ (recommended)
  │
  ├─ Publication figure creation
  │    → ParaView ⭐⭐⭐⭐⭐ (highest quality)
  │
  ├─ Radia field exploration
  │    → NGSolve webgui ⭐⭐⭐⭐ (Radia integration)
  │
  └─ Shape and mesh quality check
       → Netgen GUI ⭐⭐⭐ (lightweight, accurate)
```

### Recommendations by Use Case

| Use Case | Recommended Viewer | Reason |
|-------------|--------------|------|
| **.msh mesh inspection** | **GMSH GUI** | Native, fastest, no conversion needed |
| **Debugging** | **GMSH GUI** | Lightweight, instant inspection |
| **Field data (in .msh)** | **GMSH GUI** | View integration |
| **Development inspection** | PyVista | Rapid, script integration |
| **Parameter studies** | PyVista | Batch processing, automation |
| **Mesh quality check** | Netgen GUI | Lightweight, dedicated features |
| **Field distribution exploration** | NGSolve webgui | Interactive |
| **Publication figures** | ParaView | High quality, vector output |
| **Presentation materials** | ParaView | High resolution, beautiful |
| **Animation videos** | ParaView | Keyframe functionality |
| **Jupyter analysis** | PyVista + webgui | Notebook integration |
| **CI/CD regression testing** | PyVista (headless) | Automatable |

---

## Detailed Viewer Evaluations

### 1. GMSH GUI

**Rating**: ⭐⭐⭐⭐⭐ (5/5) - Best for `.msh` files

#### Why GMSH GUI is Optimal for .msh

| Advantage | Details |
|------|------|
| **Native integration** | Designed specifically for `.msh` format, no conversion needed |
| **Lightweight and fast** | Very fast startup and loading |
| **Field data support** | Directly displays solution data within `.msh` |
| **Python API integration** | Automatable via `gmsh.view`, `gmsh.plugin` |
| **Post-processing features** | Built-in isosurfaces, vectors, streamlines |
| **No additional installation** | Bundled with GMSH, no dependencies |

#### Differences from Other Viewers

```
GMSH GUI:     .msh → Direct display ✅
PyVista:      .msh → NGSolve → Conversion → Display ⚠️
ParaView:     .msh → Load → Display ✅ (but lacks GMSH-specific features)
```

#### Usage Examples

```bash
# Simplest and fastest
gmsh geometry.msh

# Or via Python API
python
>>> import gmsh
>>> gmsh.initialize()
>>> gmsh.open('geometry.msh')
>>> gmsh.fltk.run()
```

#### Field Data Integration

GMSH `.msh` files can contain not only meshes but also **field data (Views)**:

```python
import gmsh

gmsh.initialize()
gmsh.open('geometry.msh')

# Add field data (e.g., magnetic flux density)
view_tag = gmsh.view.add("B_field")
gmsh.view.addListData(view_tag, "ST", num_elements, data_list)

# Save (mesh + field data)
gmsh.write('geometry_with_field.msh')

# Open in GMSH GUI → Field is displayed automatically
gmsh.fltk.run()  # Launch GUI
```

**Operations in GMSH GUI**:
- Tools → Visibility: Toggle field display ON/OFF
- Tools → Options → View: Colormap, isosurface settings
- Plugins: Streamlines, Isosurface, Cut, etc.

#### Automation via Python API

**Visualization control from scripts**:

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# Load mesh
gmsh.open('coil_surface.msh')

# Add View (e.g., current density)
view = gmsh.view.add("Current Density")

# Set data (ST = Scalar Triangle)
# Set scalar values for each triangle element
num_triangles = 100  # Example
data = []
for i in range(num_triangles):
    # Triangle vertices (x1,y1,z1, x2,y2,z2, x3,y3,z3)
    # + scalar value
    data.extend([0, 0, 0,  1, 0, 0,  0, 1, 0,  1.5e6])  # Example

gmsh.view.addListData(view, "ST", num_triangles, data)

# Colormap settings
gmsh.view.option.setNumber(view, "ColormapNumber", 2)  # Jet colormap

# Launch GUI
gmsh.fltk.run()
gmsh.finalize()
```

#### Post-Processing Features

**Built-in GMSH GUI features**:

| Feature | Description | Menu |
|------|------|---------|
| **Isosurface** | Isosurface display | Plugins → Isosurface |
| **Streamlines** | Streamline display | Plugins → Streamlines |
| **CutPlane** | Cut plane display | Plugins → CutPlane |
| **Skin** | Show outer surface only | Plugins → Skin |
| **Smooth** | Data smoothing | Plugins → Smooth |

**Usage example**:
```
1. Open .msh in GMSH GUI
2. Plugins → Isosurface
   - View: Select B_field
   - Value: 0.5 (isosurface value)
   - Run
3. A new View is generated (isosurface only)
```

#### Performance Comparison

**Startup time comparison** (Windows):

| Viewer | Startup Time | Memory Usage |
|----------|---------|------------|
| **GMSH GUI** | **<1 sec** | **~50MB** |
| ParaView | ~5 sec | ~200MB |
| PyVista | ~2 sec | ~100MB |

**Large mesh loading** (1 million elements):

| Viewer | Load Time |
|----------|---------|
| **GMSH GUI** | **5 sec** |
| ParaView | 8 sec |
| PyVista | 10 sec (including conversion) |

#### Keyboard Shortcuts

| Key | Function |
|------|------|
| `0` | Toggle mesh display ON/OFF |
| `1-9` | Toggle View 1-9 display |
| `Shift+a` | Toggle axis display ON/OFF |
| `e` | Toggle element edge display ON/OFF |
| `v` | Show View panel |
| `t` | Show Tools panel |

#### Mouse Controls

| Action | Function |
|------|------|
| Left drag | Rotate |
| Middle drag | Pan |
| Scroll wheel | Zoom |
| Double-click | Select object |

#### Script Example: Mesh Quality Check

```python
import gmsh

gmsh.initialize()
gmsh.open('coil_surface.msh')

# Display mesh statistics
gmsh.plugin.setNumber("MeshQuality", "Measure", 1)  # 1=SICN
gmsh.plugin.run("MeshQuality")

# Launch GUI (quality is displayed as a colormap)
gmsh.fltk.run()
gmsh.finalize()
```

#### Example: Mesh Generation → Immediate Visualization

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.model.add("coil")

# Generate coil geometry (using methods described earlier)
# ... (geometry definition)

# Generate surface mesh
gmsh.model.mesh.generate(2)

# Display directly in GUI without saving to file
gmsh.fltk.run()  # ← That's all!

gmsh.finalize()
```

**Advantage**: No file I/O needed, inspect immediately after generation

#### Example: Field Data Visualization

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.open('coil_surface.msh')

# Generate current density data (example)
elements_2d = gmsh.model.mesh.getElements(2)
num_triangles = len(elements_2d[1][0])

# Calculate current density for each triangle (e.g., 1.5e6 A/m^2)
view = gmsh.view.add("Current Density [A/m^2]")
data_list = []

for i in range(num_triangles):
    # Get triangle vertices
    elem_nodes = gmsh.model.mesh.getElement(elements_2d[0][0],
                                             elements_2d[1][0][i])
    coords = []
    for node in elem_nodes[1]:
        coord = gmsh.model.mesh.getNode(node)[0]
        coords.extend(coord)

    # Add scalar value (current density)
    coords.append(1.5e6)
    data_list.extend(coords)

gmsh.view.addListData(view, "ST", num_triangles, data_list)

# Colormap settings
gmsh.view.option.setNumber(view, "ColormapNumber", 2)  # Jet
gmsh.view.option.setNumber(view, "RangeType", 2)  # Custom range
gmsh.view.option.setNumber(view, "CustomMin", 0)
gmsh.view.option.setNumber(view, "CustomMax", 2e6)

# Launch GUI
gmsh.fltk.run()
gmsh.finalize()
```

---

### 2. PyVista

**Rating**: ⭐⭐⭐⭐ (4/5) - Recommended as development default

**Strengths**:
- ✅ Python native - Easy script integration
- ✅ Full Jupyter Notebook/Lab support
- ✅ Rapid visualization (fast development iterations)
- ✅ Pythonic access to all VTK features
- ✅ Easy animation and GIF export
- ✅ Interactive widgets (sliders, checkboxes)
- ✅ Headless execution possible (CI/CD pipelines)
- ✅ Fast rendering
- ✅ Sufficient display quality

**Weaknesses**:
- ❌ Fine-tuning for publication quality is inferior to ParaView
- ❌ Vector graphics output is limited
- ❌ Slow for large-scale data (>10M cells)

**Usage example with Radia-NGSolve integration**:

```python
import pyvista as pv
import radia as rad

# Radia field to VTS output
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
rad.FldVTS(magnet, 'field.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)

# Visualize with PyVista
grid = pv.read('field.vts')
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='B_magnitude', cmap='coolwarm', opacity=0.8)
plotter.add_arrows(grid.points, grid['B_field'], mag=0.01, color='black')
plotter.show()

# Jupyter integration
grid.plot(scalars='B_magnitude', cmap='coolwarm', jupyter_backend='static')
```

**Loading .msh files (via NGSolve)**:

```python
from ngsolve import Mesh
import pyvista as pv
import numpy as np

# Load .msh file (via NGSolve)
mesh = Mesh('geometry.msh')

# Convert to PyVista format
points = []
cells = []
for el in mesh.Elements3D():
    vertices = [mesh.vertices[v.nr].point for v in el.vertices]
    points.extend(vertices)
    cell_indices = list(range(len(cells)*4, len(cells)*4 + 4))
    cells.append([4] + cell_indices)  # [n_points, v0, v1, v2, v3]

points_array = np.array(points)
cells_array = np.hstack(cells)

# Create PyVista mesh
grid = pv.UnstructuredGrid(cells_array, np.array([10]*len(cells)), points_array)

# Display
plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True, color='lightblue')
plotter.show()
```

**Slice display**:

```python
# Slice at the Z=0 plane
slice_z = grid.slice(normal='z', origin=[0, 0, 0])
slice_z.plot(show_edges=True)
```

**Best practices**:
- Visualization inspection during development
- Parameter study automation
- Batch processing (comparing multiple cases)
- Regression test visualization in CI/CD

---

### 3. ParaView

**Rating**: ⭐⭐⭐⭐⭐ (5/5) - Highest publication quality

**Strengths**:
- ✅ Highest quality rendering
- ✅ Vector graphics output (SVG, PDF)
- ✅ High-resolution raster images (300+ DPI)
- ✅ Complex filter chains (Glyph, Contour, StreamTracer)
- ✅ Animation and keyframes
- ✅ Large-scale data support (distributed parallel)
- ✅ Full camera and lighting control

**Weaknesses**:
- ❌ Requires GUI operation (scripting is possible but complex)
- ❌ Not suited for rapid iteration
- ❌ Limited Jupyter integration (pvpython runs as a separate process)

**Usage**:

```bash
# Open in ParaView
paraview geometry.msh
```

**Operation steps**:

1. **File → Open** → `geometry.msh`
2. Click **Apply**
3. Apply filters:
   - **Filters → Slice**: Cross-section display
   - **Filters → Clip**: Clipping display
   - **Filters → Glyph**: Vector display
4. **File → Save Screenshot**: High-resolution export

**Usage example with Radia-NGSolve integration**:

```bash
# 1. Radia field → VTS output (Python)
python generate_field_vts.py

# 2. Open in ParaView
paraview field.vts

# 3. ParaView GUI operations:
#    - Glyph filter: Vector arrows
#    - Contour: Isosurfaces
#    - Slice: Cross-sections
#    - Camera: Angle adjustment
#    - Lighting: Lighting settings

# 4. High-resolution export:
#    File > Save Screenshot
#    - Resolution: 3000x2000 (300 DPI at 10x6.67 inch)
#    - Format: PNG (raster) or SVG (vector)
```

**ParaView script automation** (pvpython):

```python
# publication_figure.py
from paraview.simple import *

# Load VTS
reader = XMLStructuredGridReader(FileName=['field.vts'])

# Add glyph
glyph = Glyph(Input=reader, GlyphType='Arrow')
glyph.ScaleFactor = 0.01
glyph.GlyphMode = 'Every Nth Point'
glyph.Stride = 2

# Render
Show(glyph)
view = GetActiveView()
view.CameraPosition = [0.3, 0.2, 0.5]
view.CameraFocalPoint = [0, 0, 0]

# Save high-res image
SaveScreenshot('figure.png', view, ImageResolution=[3000, 2000])
```

**Best practices**:
- Publication submission figures
- Presentation materials
- High-resolution posters
- Animation videos (MP4, AVI)

---

### 4. NGSolve webgui

**Rating**: ⭐⭐⭐⭐ (4/5) - Interactive exploration

**Strengths**:
- ✅ Full NGSolve integration
- ✅ Simultaneous mesh + field display
- ✅ WebGL - Runs in browser
- ✅ Jupyter integration (within the same notebook)
- ✅ Works with Radia CoefficientFunction
- ✅ Real-time updates (parameter changes)
- ✅ Usable from .py scripts
- ✅ Accurately displays OCC shapes (no approximation)

**Weaknesses**:
- ❌ Cannot directly load VTS files (NGSolve GridFunction only)
- ❌ No advanced filter features
- ❌ Difficult to fine-tune for publication quality
- ❌ Limited export formats

**Usage example with Radia-NGSolve integration**:

```python
from ngsolve import *
from ngsolve.webgui import Draw
import radia as rad

# Cubit mesh → NGSolve
mesh = Mesh('model.msh')

# Radia magnet
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# Radia CoefficientFunction (integrated into radia module since v2.5.0)
B_cf = rad.RadiaField(magnet, 'b')

# Project to GridFunction
fes = HDiv(mesh, order=2)
B_gf = GridFunction(fes)
B_gf.Set(B_cf)

# Interactive display (in Jupyter)
Draw(B_gf, mesh, name='B_field', vectors={'grid_size': 10})

# Also display the mesh
Draw(mesh)
```

**Accurate display of OCC shapes**:

```python
from ngsolve.webgui import Draw
from netgen.occ import Box, Pnt

# OCC shape (accurate)
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
Draw(occ_magnet, name='Magnet')  # Shape is completely accurate

# Field (GridFunction)
Draw(B_gf, mesh, 'B_field')  # Interactive in browser
```

**Best practices**:
- Quick inspection of field distributions
- Mesh quality check
- Real-time feedback during parameter optimization
- Education and demonstration

---

### 5. Netgen GUI

**Rating**: ⭐⭐⭐ (3/5) - Dedicated to .vol format, optimal for shape and mesh inspection

**Strengths**:
- ✅ Netgen/NGSolve native (Tcl/Tk GUI)
- ✅ Displays shapes (OCC) **accurately** (no approximation)
- ✅ Mesh quality visualization (aspect ratio, angles, etc.)
- ✅ Lightweight and fast startup (no browser needed)
- ✅ STL/STEP/IGES loading
- ✅ Integrated workflow (shape inspection → mesh generation → quality check)
- ✅ Usable from regular Python scripts (.py)
- ✅ Accurately displays exterior via surface elements

**Weaknesses**:
- ❌ Cannot directly load `.msh` files
- ❌ Limited field data display
- ❌ No publication-quality rendering
- ⚠️ Older GUI (Tcl/Tk) - but lightweight and stable
- ⚠️ Surface elements required (usually not an issue)

**Usage example**:

```python
from netgen.occ import OCCGeometry, Box, Pnt
from netgen.gui import StartGUI

# Convert Radia magnet to OCC shape
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
geo = OCCGeometry(occ_magnet)

# Inspect in Netgen GUI (opens a native window)
StartGUI()
geo.Draw()  # Displays shape accurately

# Generate mesh
mesh = geo.GenerateMesh(maxh=0.005)
mesh.Draw()  # Inspect mesh quality
```

**.msh → .vol conversion (via NGSolve)**:

```python
from ngsolve import Mesh

mesh = Mesh('geometry.msh')
mesh.ngmesh.Save('geometry.vol')  # Save in .vol format

# Open in Netgen GUI
# python utils/netgen_vol_viewer.py geometry.vol
```

**Choosing between netgen.gui and webgui**:

| Purpose | netgen.gui | ngsolve.webgui |
|------|-----------|---------------|
| Shape inspection | **✅ Recommended** | ⚠️ Browser overhead |
| Mesh quality | **✅ Recommended** | ❌ Limited |
| Field visualization | ❌ Not available | **✅ Recommended** |
| Lightweight and fast | **✅ Native GUI** | ⚠️ Requires browser |

**Windows file association**:

Setting up .vol files to open in Netgen GUI on double-click:

```bash
# Automatic setup (administrator privileges)
cd S:\Radia\01_GitHub\utils
setup_vol_file_association.bat
```

Details: [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md)

**Best practices**:
- Shape inspection before mesh generation (**more accurate than webgui**)
- Mesh quality check (aspect ratio, angles)
- Boundary condition label verification
- Shape verification after CAD import
- Cubit → Netgen workflow verification

---

## Geometry Accuracy: VTS vs OCC

### Problem: Shape Approximation in VTS Export

#### Behavior of rad.FldVTS()

```python
# Radia analysis object (completely accurate)
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# VTS export (computes field values at grid points)
rad.FldVTS(magnet, 'field.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)
```

**Information contained in VTS files**:
- ✅ **Field values**: B, H, A, Phi at grid points (accurate)
- ❌ **Shape information**: Lost (approximated by grid outline)

**Result**:
- When opening `field.vts` in ParaView, the field distribution is accurate, but the magnet shape is approximated by the grid
- A rectangular magnet is displayed with grid-like boundaries (not a perfect rectangle)

### Shape Accuracy by Viewer

| Method | Shape Accuracy | Field Accuracy | Effort | Quality |
|------|---------|--------------|------|-----|
| **PyVista + VTS** | ❌ Approximate | ✅ Accurate | Low | Good |
| **ParaView + VTS** | ❌ Approximate | ✅ Accurate | Medium | Excellent |
| **webgui + OCC** | ✅ Accurate | ✅ Accurate | Medium | Good |
| **ParaView + STL/STEP + VTS** | ✅ Accurate | ✅ Accurate | High | Excellent |

### Solution: ParaView STL+VTS Overlay

To accurately display both shape and field:

```python
# 1. Export shape as STL/STEP (accurate)
from netgen.occ import Box, Pnt, OCCGeometry

box = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
geo = OCCGeometry(box)
mesh = geo.GenerateMesh(maxh=0.002)
mesh.Export('magnet_shape.stl', 'STL Format')

# 2. Export field as VTS
rad.FldVTS(magnet, 'field.vts', ...)

# 3. Overlay both in ParaView for high-quality figures
# - magnet_shape.stl: Display shape semi-transparently
# - field.vts: Display field with colormap
```

### Conclusion

**The user's concern that "shapes are not necessarily preserved when going through ParaView" is correct**:
- VTS format uses structured grids, so shapes are approximated by the grid
- Analytical rectangles and cylinders are represented by grid boundaries

**Solutions**:
1. **Shape not needed**: PyVista/ParaView + VTS (current implementation is sufficient)
2. **Shape important (development)**: NGSolve webgui + OCC shapes
3. **Shape important (publication)**: ParaView + STL/STEP + VTS overlay

---

## Recommended Workflows by Use Case

### Mesh Generation Patterns

#### Pattern A: Direct Netgen Generation (Recommended)

```python
from netgen.occ import Box, Pnt, OCCGeometry
from netgen.gui import StartGUI

# Create shape
box = Box(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05))
geo = OCCGeometry(box)

# Generate mesh (surface elements automatic)
mesh = geo.GenerateMesh(maxh=0.01)

print(f"Volume elements:  {mesh.ne}")
print(f"Surface elements: {mesh.nse}")  # > 0 (auto-generated)

# Save
mesh.Save('magnet.vol')  # Includes surface elements

# Inspect in Netgen GUI
StartGUI()
mesh.Draw()  # ← Displays without issues
```

**Result**: Surface elements are automatically included ✅

#### Pattern B: Cubit → Netgen Conversion

```python
import cubit
import cubit_mesh_export
from ngsolve import Mesh

# Generate mesh in Cubit
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'motor_rotor.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Define sidesets (these become surface elements)
cubit.cmd("sideset 1 surface all")
cubit.cmd("sideset 1 name 'boundary'")

# Convert to Netgen (preserving surface elements)
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# Save
mesh.ngmesh.Save('motor_rotor.vol')  # Includes surface elements
```

**Result**: Cubit sidesets are converted to surface elements ✅

For details on mesh generation and surface elements, see [MESH_GUIDE.md](MESH_GUIDE.md).

#### Pattern C: GMSH → NGSolve

```python
from ngsolve import Mesh

# Load GMSH mesh (NGSolve auto-converts)
mesh = Mesh('geometry.msh')

# Surface elements are automatically recognized
print(f"Surface elements: {mesh.ngmesh.nse}")  # > 0

# Save as .vol
mesh.ngmesh.Save('geometry.vol')
```

**Result**: GMSH boundary elements are converted to Netgen surface elements ✅

### Visualization Workflows (by Purpose)

#### Flow 1: .msh Mesh Inspection (Top Priority)

```
GMSH mesh generate → geometry.msh → **GMSH GUI**
                                        ↓
                                  Instant inspection (no conversion needed)
```

#### Flow 2: Shape and Mesh Inspection

```
Netgen / Cubit → .vol → **Netgen GUI** (lightweight, accurate shapes)
```

```python
from netgen.meshing import Mesh
from netgen.gui import StartGUI

mesh = Mesh()
mesh.Load('magnet.vol')

StartGUI()
mesh.Draw()
```

#### Flow 3: Development and Debugging

```
Cubit → Netgen mesh → Radia solve → VTS export → **PyVista**
                                                     ↓
                                              Rapid inspection and correction
```

```python
import radia as rad
import pyvista as pv

# Radia magnet
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# VTS export
rad.FldVTS(magnet, 'field.vts', ...)

# PyVista visualization (rapid)
grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude', cmap='coolwarm')
```

#### Flow 4: Publication Submission

```
Cubit → Netgen mesh → Radia solve → VTS export → **ParaView**
                                                     ↓
                                         High-resolution PNG/SVG output
```

When both shape and field accuracy are important:
```
Radia → STL/STEP export ───────────┐
     ↓                              ├→ ParaView overlay
     └→ FldVTS() → VTS ────────────┘
                  ↓
            Accurate shape + Accurate field
```

#### Flow 5: Interactive Exploration

```
Cubit → Netgen mesh → Radia CF → NGSolve GridFunction → **webgui**
                                                            ↓
                                                  Instant inspection in browser
```

### Recommended Flowchart (Overall)

```
Mesh generation
    ↓
What is the file format?
    ├─ .msh → GMSH GUI (fastest inspection)
    └─ .vol → Choose by purpose
                ├─ Shape inspection → Netgen GUI (lightest weight)
                ├─ Development inspection → PyVista (rapid)
                ├─ Publication figures → ParaView (high quality)
                └─ Integrated exploration → webgui (accurate)
```

### Frequently Asked Questions

#### Q: Is the lack of surface elements a problem?

**A: In practice, it is not an issue.**

Reasons:
- Netgen-generated meshes: Surface elements are automatically included
- Cubit-converted meshes: Surface elements are generated via sideset definitions
- NGSolve samples: All include surface elements

In rare cases with volume elements only: Use cut-plane display in ParaView/PyVista

#### Q: I want to open files by double-clicking on Windows

**A: This is possible.**

```bash
# Automatic setup (administrator privileges)
cd S:\Radia\01_GitHub\utils
setup_vol_file_association.bat
```

Details: [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md)

#### Q: What about the shape approximation issue in ParaView?

**A: There are two solutions.**

1. **Use Netgen GUI for shape inspection** (accurate, lightweight)
2. **Use STL+VTS overlay for publication figures** (accurate, high quality)

For details, see the "Geometry Accuracy: VTS vs OCC" section above.

---

## Implementation Status

### Phase 1: Basic Visualization (Completed)
- ✅ `rad.FldVTS()` - VTS export
- ✅ Basic plotting with PyVista
- ✅ Manual visualization in ParaView

### Phase 2: Enhanced NGSolve Integration
- ✅ `rad.RadiaField` CoefficientFunction (B, H, A, phi, M)
- ✅ `as_voxel_cf()` for VoxelCoefficient generation
- ✅ Coordinate transform support (origin, u, v, w axes)
- ⏳ webgui integration script collection

### Phase 3: Advanced Visualization (Planned)
- ⏳ ParaView automation scripts (pvpython)
- ⏳ PyVista animation generation
- ⏳ STL/STEP export automation
- ⏳ ParaView overlay automation scripts

### Phase 4: Integrated Viewer (Future)
- ⏳ Custom PyVista interface
- ⏳ Parameter control via Jupyter Widgets
- ⏳ Web application UI (Dash/Streamlit)

---

## References

### Installation Commands

```bash
# GMSH (mesh generation + viewer)
pip install gmsh  # Includes both Python library and GUI

# PyVista (recommended for development)
pip install pyvista

# ParaView (for publications)
# Download from https://www.paraview.org/download/

# NGSolve webgui (field exploration)
pip install ngsolve  # webgui included
```

### Reference Implementation

**EMPY_Field** (`S:\NGSolve\EMPY\EMPY_Field`):
- Example implementation of OCC conversion for Radia analysis objects
- Represents accurate shapes using OCC

### Related Documents

- [MESH_GUIDE.md](MESH_GUIDE.md) - Details on mesh generation and surface elements
- [MESH_GUIDE.md](MESH_GUIDE.md) - GMSH mesh generation workflow
- [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md) - .vol file association setup

---

**Created**: 2026-02-12
**Updated**: 2026-02-22
**Subject**: Visualization guide for the Radia-NGSolve integration framework (viewer selection, shape accuracy, workflows)
