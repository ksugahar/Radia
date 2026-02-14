# Radia-NGSolve Visualization Examples

This directory contains visualization examples for Radia-NGSolve framework.

## Quick Start

### 推奨: 通常のPythonスクリプトワークフロー

**注意**: Jupyterノートブック（.ipynb）はClaude Codeとの相性が悪いため、通常のPythonスクリプト（.py）を推奨します。

### PyVista (Default Viewer)

```bash
# Install PyVista
pip install pyvista

# Run demo - インタラクティブウィンドウが開く
python demo_pyvista_basic.py
```

**Use case:** Daily development, quick field checks, parameter studies

### NGSolve webgui (Interactive Explorer - 推奨)

```bash
# Requires radia_ngsolve.pyd
# ブラウザで自動的に開く
python demo_ngsolve_webgui.py
```

**Use case:** Interactive field exploration, browser-based visualization, education

**利点**:
- ✅ ブラウザベース（Jupyter不要）
- ✅ 通常のPythonスクリプト（.py）から使用可能
- ✅ Claude Codeと相性良好

### ParaView (Publication Quality)

```bash
# Generate VTS file and open in ParaView manually
python generate_paraview_figure.py

# OR: Use pvpython for automated figure generation
pvpython generate_paraview_figure.py
```

**Use case:** Journal papers, presentations, high-resolution figures

### Jupyter Notebook (オプション - 非推奨)

```bash
# Jupyter統合が必要な場合のみ
jupyter notebook jupyter_visualization_demo.ipynb
```

**注意**: Claude Codeでの編集・保守が困難なため、特別な理由がない限り通常の.pyスクリプトを使用してください。

---

## Files

| File | Description | Requirements |
|------|-------------|--------------|
| `demo_pyvista_basic.py` | PyVista basic visualization | pyvista |
| `demo_ngsolve_webgui.py` | NGSolve webgui demo | radia_ngsolve.pyd |
| `generate_paraview_figure.py` | ParaView figure automation | ParaView (optional) |
| `jupyter_visualization_demo.ipynb` | Jupyter notebook demo | pyvista, ipywidgets |
| `README.md` | This file | - |

---

## Viewer Comparison

### PyVista (Recommended Default)

**Pros:**
- ✅ Python-native (easy scripting)
- ✅ Jupyter integration
- ✅ Fast iteration
- ✅ Headless execution (CI/CD)

**Cons:**
- ❌ Publication quality not as good as ParaView
- ❌ Limited vector graphics export

**Best for:**
- Development and debugging
- Parameter studies
- Batch processing
- CI/CD visualization

### ParaView (Publication Quality)

**Pros:**
- ✅ Highest quality rendering
- ✅ Vector graphics (SVG, PDF)
- ✅ High-resolution raster (300+ DPI)
- ✅ Advanced filters (Contour, StreamTracer)

**Cons:**
- ❌ GUI required (scripting is complex)
- ❌ Slower iteration
- ❌ Limited Jupyter integration

**Best for:**
- Journal paper figures
- Presentation slides
- High-resolution posters
- Animation videos

### NGSolve webgui (Interactive)

**Pros:**
- ✅ NGSolve native integration
- ✅ Mesh + field simultaneous display
- ✅ WebGL (browser-based)
- ✅ Jupyter inline display

**Cons:**
- ❌ Cannot load VTS files directly
- ❌ Limited advanced filters
- ❌ Export options limited

**Best for:**
- Interactive field exploration
- Mesh quality checks
- Real-time parameter feedback
- Education and demos

### Netgen GUI (Mesh Checker)

**Pros:**
- ✅ Netgen/NGSolve native
- ✅ Direct OCC shape display
- ✅ Mesh quality visualization
- ✅ Lightweight and fast

**Cons:**
- ❌ Limited field visualization
- ❌ Not suitable for publication
- ❌ Difficult to script

**Best for:**
- Mesh generation verification
- Geometry checking
- Boundary label verification

---

## Workflow Recommendations

### Workflow 1: Development (Default)

```
Cubit → Netgen mesh → Radia solve → VTS export → **PyVista**
                                                     ↓
                                              Quick verification
```

**Tools:** PyVista + Jupyter

### Workflow 2: Publication

```
Cubit → Netgen mesh → Radia solve → VTS export → **ParaView**
                                                     ↓
                                         High-res PNG/SVG export
```

**Tools:** ParaView + pvpython (optional automation)

### Workflow 3: Interactive Exploration

```
Cubit → Netgen mesh → Radia CF → NGSolve GF → **webgui**
                                                  ↓
                                      Browser-based exploration
```

**Tools:** NGSolve webgui + Jupyter

---

## Installation

### PyVista

```bash
pip install pyvista
```

### ParaView

Download from: https://www.paraview.org/download/

**Note:** ParaView includes `pvpython` for scripting.

### NGSolve webgui

Already included in NGSolve installation:

```bash
pip install ngsolve
```

### radia_ngsolve (Optional)

Build from source:

```powershell
# In S:\Radia\01_GitHub
powershell.exe -ExecutionPolicy Bypass -File BuildMSVC.ps1
```

---

## Usage Examples

### PyVista: Basic Plot

```python
import pyvista as pv
import radia as rad

rad.FldUnits('m')

# Create magnet and export field
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
rad.FldVTS(magnet, 'field.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)

# Visualize
grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude', cmap='coolwarm')
```

### NGSolve webgui: Interactive

```python
from ngsolve import *
from ngsolve.webgui import Draw
from radia_ngsolve import RadiaField
import radia as rad

rad.FldUnits('m')

# Create mesh and Radia field
mesh = Mesh('model.msh')
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
B_cf = RadiaField(magnet, 'b')

# Project and display
fes = HDiv(mesh, order=2)
B_gf = GridFunction(fes)
B_gf.Set(B_cf)
Draw(B_gf, mesh, 'B_field', vectors={'grid_size': 10})
```

### ParaView: Manual Workflow

```bash
# 1. Generate VTS
python -c "
import radia as rad
rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
rad.FldVTS(magnet, 'field.vts', [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15], 41, 41, 27, 1, 0, 1.0)
"

# 2. Open in ParaView
paraview field.vts

# 3. Apply filters (Glyph, Contour, etc.) via GUI

# 4. Save Screenshot (3000x2000, PNG or SVG)
```

---

## Troubleshooting

### PyVista: "No module named 'pyvista'"

```bash
pip install pyvista
```

### NGSolve webgui: "Cannot import radia_ngsolve"

Build `radia_ngsolve.pyd`:

```powershell
# S:\Radia\01_GitHub
powershell.exe -ExecutionPolicy Bypass -File BuildMSVC.ps1
```

### ParaView: "pvpython not found"

Add ParaView to PATH or use full path:

```bash
# Windows
"C:\Program Files\ParaView\bin\pvpython.exe" generate_paraview_figure.py

# Linux
/usr/bin/pvpython generate_paraview_figure.py
```

### Jupyter: "PyVista plot not showing"

Set Jupyter backend:

```python
import pyvista as pv
pv.set_jupyter_backend('static')  # or 'trame' for interactive
```

---

## References

- [PyVista Documentation](https://docs.pyvista.org/)
- [ParaView Documentation](https://www.paraview.org/documentation/)
- [NGSolve webgui](https://docu.ngsolve.org/latest/i-tutorials/wgs/wgs.html)
- [Radia Documentation](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)

---

**Last Updated:** 2026-02-12
