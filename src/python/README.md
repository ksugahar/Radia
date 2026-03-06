# Radia Python Utilities

Python utility modules for Radia visualization, data export, and NGSolve integration.

## Files

### rad_ngsolve.cpp

**C++ CoefficientFunction integration for NGSolve (recommended).**

High-performance NGSolve CoefficientFunction wrappers for Radia magnetostatics fields.

**Features:**
- Exact Radia field evaluation in NGSolve
- Automatic unit conversion (NGSolve meters ↔ Radia millimeters)
- Thread-safe Python GIL handling
- Compatible with both 2D and 3D meshes
- Three field types: B-field, H-field, A-field

**Build:**
```bash
cd S:\radia\01_GitHub
.\Build_NGSolve.ps1
```

**Usage:**
```python
from ngsolve import *
import radia as rad
import rad_ngsolve

# Create Radia magnet (mm)
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.Solve(magnet, 0.0001, 10000)

# Create mesh (m)
mesh = Mesh(...)

# CoefficientFunction (exact)
B_cf = rad_ngsolve.RadBfield(magnet)
B = B_cf(mesh(0, 0, 0.02))  # Auto converts m->mm
```

**See also:** `examples/ngsolve_integration/` for complete examples.

### rad_ngsolve_py.py

**Pure Python implementation (alternative).**

Python-only NGSolve wrappers without requiring C++ compilation. Same interface as C++ version.

**Note:** C++ version is recommended for production (better performance, automatic unit conversion).

### radia_coil_builder.py

**Modern fluent interface for constructing complex coil geometries.**

Elegant object-oriented design for multi-segment coil paths with automatic state tracking.

**Features:**
- Fluent method chaining
- Automatic position/orientation tracking
- Type-safe with abstract base classes
- Automatic cross-section transformation
- Direct conversion to Radia objects

**Usage:**
```python
from radia_coil_builder import CoilBuilder

coil = (CoilBuilder(current=1000)
	.set_start([0, 0, 0])
	.set_cross_section(width=20, height=20)
	.add_straight(length=100)
	.add_arc(radius=50, arc_angle=180, tilt=90)
	.add_straight(length=100)
	.add_arc(radius=50, arc_angle=180, tilt=90)
	.to_radia())
```

### radia_vtk_export.py

**VTK export utilities for ParaView visualization.**

Export Radia geometry to VTK Legacy format.

**Usage:**
```python
from radia_vtk_export import exportGeometryToVTK

mag = rad.ObjRecMag([0,0,0], [30,30,10], [0,0,1])
exportGeometryToVTK(mag, 'my_magnet')
```

### radia_pyvista_viewer.py

**Interactive 3D viewer using PyVista.**

Real-time interactive visualization of Radia objects.

**Requirements:**
```bash
pip install pyvista
```

**Usage:**
```python
from radia_pyvista_viewer import view_radia_object

mag = rad.ObjRecMag([0,0,0], [30,30,10], [0,0,1])
view_radia_object(mag)
```

## Visualization Workflow

### Option 1: ParaView (Static Export)
Best for publication-quality figures, batch processing.
```python
from radia_vtk_export import exportGeometryToVTK
exportGeometryToVTK(my_object, 'output')
```

### Option 2: PyVista (Interactive)
Best for quick inspection, interactive exploration.
```python
from radia_pyvista_viewer import view_radia_object
view_radia_object(my_object)
```

### Option 3: NGSolve Integration
For coupled FEM simulations.
```python
import rad_ngsolve
B_cf = rad_ngsolve.RadBfield(magnet)
```

## References

- **NGSolve Integration**: `examples/ngsolve_integration/`
- **Coil Builder Examples**: `examples/complex_coil_geometry/`
- **ParaView**: https://www.paraview.org/
- **PyVista**: https://docs.pyvista.org/

Last Updated: 2025-10-31
