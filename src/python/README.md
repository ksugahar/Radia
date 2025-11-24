# Radia Python Utilities

Python utility modules for Radia visualization, data export, and NGSolve integration.

## Files

### Mesh Import

#### nastran_mesh_import.py

**Unified Nastran mesh import for Radia.**

Imports Nastran mesh files (.bdf, .nas, .dat) and creates Radia geometry.

**Supported Element Types:**
- **CHEXA** - 8-node hexahedron
- **CPENTA** - 6-node wedge/prism
- **CPYRAM** - 5-node pyramid
- **CTETRA** - 4-node tetrahedron
- **CTRIA3** - 3-node triangle (surface mesh, grouped by material ID)
- **GRID/GRID*** - Node definitions (fixed-width format)

**Usage:**
```python
from nastran_mesh_import import create_radia_from_nastran
import radia as rad

rad.FldUnits('mm')

# Import mesh and create Radia objects
mag_obj = create_radia_from_nastran(
    'York.bdf',
    material={'magnetization': [0, 0, 0]},
    units='mm',
    combine=True
)

# Apply material
mat = rad.MatSatIsoFrm([20000, 2], [0.1, 2], [0.1, 2])
rad.MatApl(mag_obj, mat)
```

**CTRIA3 Surface Meshes:**
Surface triangles are automatically grouped by material ID (property ID).
Each material ID creates one polyhedron from all its triangles.

```python
from nastran_mesh_import import import_nastran_mesh

mesh_data = import_nastran_mesh('sphere.bdf', units='mm')
tria_groups = mesh_data['tria_groups']
# Format: {material_id: {'faces': [[n1,n2,n3], ...], 'node_ids': set(...)}}
```

**Migration from nastran_reader.py (removed 2025-11-23):**
- `nastran_reader.py` has been removed and unified into `nastran_mesh_import.py`
- Use `create_radia_from_nastran()` for direct Radia object creation
- Face topologies (TETRA_FACES, WEDGE_FACES, etc.) are in `netgen_mesh_import.py`

#### netgen_mesh_import.py

**Import Netgen/NGSolve meshes to Radia.**

Provides face topology constants and mesh conversion utilities.

**Face Topology Constants:**
```python
from netgen_mesh_import import (
    TETRA_FACES,    # Tetrahedron
    HEX_FACES,      # Hexahedron
    WEDGE_FACES,    # Wedge/Prism
    PYRAMID_FACES   # Pyramid
)
```

### NGSolve Integration

#### radia_ngsolve.cpp

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
import radia_ngsolve

# Create Radia magnet (mm)
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.Solve(magnet, 0.0001, 10000)

# Create mesh (m)
mesh = Mesh(...)

# CoefficientFunction (exact)
B_cf = radia_ngsolve.RadBfield(magnet)
B = B_cf(mesh(0, 0, 0.02))  # Auto converts m->mm
```

**See also:** `examples/ngsolve_integration/` for complete examples.

### radia_ngsolve_py.py

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
import radia_ngsolve
B_cf = radia_ngsolve.RadBfield(magnet)
```

## References

- **NGSolve Integration**: `examples/ngsolve_integration/`
- **Coil Builder Examples**: `examples/complex_coil_geometry/`
- **ParaView**: https://www.paraview.org/
- **PyVista**: https://docs.pyvista.org/

Last Updated: 2025-10-31
