# Radia Examples

Comprehensive collection of Radia examples demonstrating magnetic field computation, material properties, solver methods, and integration with NGSolve FEM.

**Total Examples:** 41 Python scripts across 9 directories

---

## Quick Start

```bash
# Navigate to any example directory
cd examples/simple_problems

# Run any example
python arc_current_with_magnet.py

# View geometry in ParaView
paraview arc_current_with_magnet.vtk
```

---

## Directory Overview

### 1. [simple_problems/](simple_problems/) - Basic Radia Examples
**6 scripts** | **Level:** Beginner

Basic Radia functionality including magnets, coils, materials, and field calculations.

**Key Examples:**
- `arc_current_with_magnet.py` - Arc current element with rectangular magnet
- `arc_current_dual_magnets.py` - Multiple magnets with arc current
- `chamfered_pole_piece.py` - Complex extrusion with chamfer
- `cubic_polyhedron_magnet.py` - Polyhedron magnet creation
- `compare_magpylib.py` - Validation against magpylib library
- `hmatrix_update_magnetization.py` - H-matrix magnetization updates

**Topics Covered:**
- Object creation (`ObjRecMag`, `ObjArcCur`, `ObjPolyhdr`)
- Material properties (`MatLin`, `MatPM`)
- Field computation (`Fld`)
- VTK export for visualization

**Best for:** Learning basic Radia API and concepts

---

### 2. [background_fields/](background_fields/) - Background Field Examples
**4 scripts** | **Level:** Intermediate

Using external background fields with magnetizable materials.

**Key Examples:**
- `quadrupole_analytical.py` - Analytical quadrupole field with magnetizable sphere
- `sphere_in_quadrupole.py` - Python callback background field (rad.ObjBckgCF)
- `permeability_comparison.py` - Material permeability analysis
- `sphere_nastran_analysis.py` - Nastran mesh with background field

**Topics Covered:**
- Background field specification (`ObjBckgCF`)
- Python callback functions for custom fields
- Magnetizable materials in external fields
- Nastran mesh import for complex geometries
- Field/material interaction analysis

**Best for:** Coupling Radia with external field sources or FEM

---

### 3. [electromagnet/](electromagnet/) - Electromagnet Simulation
**3 scripts + data** | **Level:** Advanced

Complete electromagnet with racetrack coil and magnetic yoke.

**Key Files:**
- `magnet.py` - Complete electromagnet simulation
- `racetrack_coil_model.py` - Racetrack coil geometry
- `yoke_model.py` - Nastran mesh import for magnetic yoke
- `York.bdf` - Magnetic yoke mesh (Nastran format, 55KB)

**Topics Covered:**
- Racetrack coil geometry (`ObjRaceTrk`)
- Nastran mesh import (hexahedra + pentahedra)
- Magnetostatic solver for nonlinear materials
- Field distribution calculation (3D grid)
- Graceful degradation (coil-only mode if mesh missing)

**Output:**
- `electromagnet.vtk` - Coil + yoke geometry
- `field_distribution.vtk` - 3D magnetic field vectors

**Best for:** Real-world electromagnet design and analysis

---

### 4. [complex_coil_geometry/](complex_coil_geometry/) - CoilBuilder Examples
**3 scripts** | **Level:** Intermediate

Multi-segment coils using the modern CoilBuilder API.

**Key Examples:**
- `coil_model.py` - 8-segment beam steering coil module
- `visualize_coils.py` - Coil visualization and field verification
- `field_map.py` - 3D magnetic field distribution

**Topics Covered:**
- CoilBuilder fluent API (method chaining)
- Straight and arc segments with automatic state tracking
- Cross-section transformations with tilt
- Field map calculation on structured grids
- Modular coil design (reusable coil modules)

**CoilBuilder Features:**
- 75% less code than manual tracking
- Automatic position/orientation updates
- Type-safe with abstract base classes
- Direct conversion to Radia objects

**Best for:** Building complex multi-segment coil geometries

---

### 5. [NGSolve_Integration/](NGSolve_Integration/) - rad_ngsolve Examples
**9 scripts** | **Level:** Intermediate to Advanced

Coupling Radia magnetic fields with NGSolve finite element analysis.

**Key Examples:**
- `demo_field_types.py` - All field types (b, h, a, m)
- `visualize_field.py` - Field visualization and comparison
- `export_radia_geometry.py` - Geometry export to VTK
- `test_batch_evaluation.py` - Batch field evaluation performance
- `verify_curl_A_equals_B.py` - Verify ∇×A = B mathematically

**Topics Covered:**
- NGSolve CoefficientFunction (`rad_ngsolve.RadiaField`)
- Unit conversion (meters ↔ millimeters)
- Field types: B (flux density), H (field), A (vector potential), M (magnetization)
- GridFunction.Set() for field initialization
- Mesh convergence studies
- Performance optimization

**Unit Convention:**
- NGSolve: meters (m)
- Radia: millimeters (mm)
- Automatic conversion: coordinates × 1000 (m → mm)

**Best for:** Coupling Radia with FEM for multiphysics simulations

---

### 6. [H-matrix/](H-matrix/) - H-matrix Benchmarks
**5 scripts** | **Level:** Advanced

Performance benchmarks for H-matrix accelerated computations.

**Key Examples:**
- `benchmark_solver.py` - H-matrix vs direct solver comparison
- `verify_field_accuracy.py` - H-matrix field accuracy verification
- `benchmark_parallel_construction.py` - OpenMP parallelization benchmarks
- `benchmark_field_evaluation.py` - Single vs batch field evaluation
- `run_all_benchmarks.py` - Automated benchmark suite

**Topics Covered:**
- H-matrix construction and accuracy
- Low-rank approximation (ACA algorithm)
- OpenMP parallel construction
- Field evaluation acceleration
- Accuracy vs performance trade-offs

**Performance Results:**
- H-matrix: O(N log N) complexity
- Direct: O(N²) complexity
- Accuracy: <1% error with ε=1e-6
- Speedup: 10-100x for N>1000 elements

**Best for:** Large-scale magnetostatics problems (1000+ elements)

---

### 7. [solver_time_evaluation/](solver_time_evaluation/) - Solver Benchmarks
**4 scripts** | **Level:** Intermediate

Solver performance analysis and scaling studies.

**Key Examples:**
- `benchmark_linear_material.py` - Linear material solver scaling (O(N^1.15))
- `benchmark_lu_vs_gs.py` - LU decomposition vs Gauss-Seidel comparison
- `benchmark_matrix_construction.py` - Matrix assembly timing (O(N^1.44))
- `benchmark_solver_scaling.py` - General solver scaling analysis

**Topics Covered:**
- Solver complexity analysis
- LU decomposition (exact, O(N^2.5))
- Gauss-Seidel iteration (approximate, O(N^0.5))
- Matrix construction timing
- Performance vs problem size

**Key Findings:**
- Gauss-Seidel preferred for N>27 elements
- LU decomposition better for small problems (N<27)
- Matrix construction: O(N^1.44) scaling

**Best for:** Understanding solver performance characteristics

---

### 8. [solver_benchmarks/](solver_benchmarks/) - Additional Benchmarks
**2 scripts** | **Level:** Intermediate

Additional solver method comparisons and performance tests.

**Key Examples:**
- `benchmark_hmatrix_field.py` - H-matrix field evaluation benchmark
- `benchmark_solver_methods.py` - Comparison of solver methods

**Topics Covered:**
- H-matrix field evaluation
- Solver method selection
- Performance comparison

**Best for:** Choosing the right solver for your problem

---

### 9. [smco_magnet_array/](smco_magnet_array/) - SmCo Magnet Array
**1 script** | **Level:** Intermediate

Samarium-cobalt permanent magnet array simulation.

**Key Example:**
- `smco_array.py` - SmCo magnet array with field calculation

**Topics Covered:**
- Permanent magnet materials (SmCo)
- Magnet array construction
- Field uniformity analysis
- Material properties (Br, Hc)

**Best for:** Permanent magnet array design

---

## Common Patterns

### Path Setup

All examples use consistent path setup:

```python
import sys
import os

# Add Radia module paths
project_root = os.path.join(os.path.dirname(__file__), '..', '..')
sys.path.insert(0, os.path.join(project_root, 'build', 'Release'))
sys.path.insert(0, os.path.join(project_root, 'dist'))
sys.path.insert(0, os.path.join(project_root, 'src', 'python'))

import radia as rad
```

### VTK Export

All simple_problems and demonstration scripts include VTK export:

```python
# VTK Export - Export geometry with same filename as script
try:
    from radia_vtk_export import exportGeometryToVTK

    script_name = os.path.splitext(os.path.basename(__file__))[0]
    vtk_filename = f"{script_name}.vtk"
    vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

    exportGeometryToVTK(g, vtk_path)
    print(f"\n[VTK] Exported: {vtk_filename}")
    print(f"      View with: paraview {vtk_filename}")
except ImportError:
    print("\n[VTK] Warning: radia_vtk_export not available")
```

### Material API

All examples use the new Material API:

```python
# Isotropic linear material (μr = 1000)
mat = rad.MatLin(999)  # chi = μr - 1

# Anisotropic linear material
mat = rad.MatLin([0.06, 0.17], [0, 0, 1])  # [chi_par, chi_perp], easy_axis

# Permanent magnet (NdFeB)
mat = rad.MatPM(1.2, 900000, [0, 0, 1])  # Br, Hc, magnetization_direction

# Saturating material (Steel37)
mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])

# Apply material to object
rad.MatApl(obj, mat)
```

---

## Running Examples

### Prerequisites

**Required:**
- Python 3.12
- Radia (build with `Build.ps1`)
- NumPy

**Optional:**
- NGSolve (for NGSolve_Integration/ examples)
- PyVista (for interactive 3D visualization)
- ParaView (for viewing VTK files)

### Build Radia First

```bash
cd <project_root>
powershell.exe -ExecutionPolicy Bypass -File Build.ps1
```

This creates:
- `dist/radia.pyd` - Core Radia module
- `build/Release/radia.pyd` - Alternative location

### Run Any Example

```bash
cd examples/simple_problems
python arc_current_with_magnet.py
```

### View Output in ParaView

Most examples generate VTK files:

```bash
# Open geometry
paraview arc_current_with_magnet.vtk

# Open field distribution
paraview field_distribution.vtk
```

---

## Example Selection Guide

| Use Case | Recommended Examples |
|----------|---------------------|
| **Learn Radia basics** | `simple_problems/` |
| **Permanent magnets** | `simple_problems/`, `smco_magnet_array/` |
| **Electromagnets** | `electromagnet/`, `complex_coil_geometry/` |
| **Complex coils** | `complex_coil_geometry/` |
| **External fields** | `background_fields/` |
| **FEM coupling** | `NGSolve_Integration/` |
| **Performance analysis** | `solver_time_evaluation/`, `H-matrix/` |
| **Large problems (N>1000)** | `H-matrix/` |
| **Material properties** | `background_fields/`, `solver_time_evaluation/` |

---

## Coordinate System

**All examples use millimeters (mm)** unless otherwise noted.

- X: Horizontal (left-right)
- Y: Horizontal (front-back)
- Z: Vertical (up-down)

**Exception:** NGSolve integration examples use meters (m) for NGSolve meshes, with automatic conversion to mm for Radia.

---

## Visualization

### ParaView (Recommended)

Free, open-source 3D visualization:

1. Open `.vtk` file in ParaView
2. Apply filters:
   - **Glyph** - Vector field arrows
   - **StreamTracer** - Field lines
   - **Contour** - Constant field surfaces
   - **Slice** - Cutting planes

Download: https://www.paraview.org/

### PyVista (Interactive)

Python-based interactive 3D viewer:

```python
from radia_pyvista_viewer import view_radia_object
view_radia_object(mag)
```

Install: `pip install pyvista`

---

## Documentation

Each subdirectory contains a comprehensive README.md with:
- Example descriptions
- API documentation
- Usage instructions
- Troubleshooting tips
- References

**Main Documentation:**
- [README.md](../README.md) - Project overview
- [README_BUILD.md](../README_BUILD.md) - Build instructions
- [docs/OPENMP_PERFORMANCE_REPORT.md](../docs/OPENMP_PERFORMANCE_REPORT.md) - OpenMP benchmarks
- [RAD_NGSOLVE_BUILD_SUCCESS.md](../RAD_NGSOLVE_BUILD_SUCCESS.md) - NGSolve integration

---

## Troubleshooting

### ModuleNotFoundError: No module named 'radia'

**Solution:** Build Radia first:
```bash
cd <project_root>
powershell.exe -ExecutionPolicy Bypass -File Build.ps1
```

### ImportError: DLL load failed

**Cause:** Missing Visual C++ Redistributable

**Solution:** Install Visual C++ 2022 Redistributable:
- Download from Microsoft website
- Both x86 and x64 versions may be needed

### No module named 'rad_ngsolve'

**Solution:** Build NGSolve integration:
```bash
cd <project_root>
powershell.exe -ExecutionPolicy Bypass -File Build_NGSolve.ps1
```

Requires NGSolve installed.

### VTK export not available

**Cause:** `radia_vtk_export.py` not in path

**Solution:** Examples automatically add `src/python/` to path. Verify file exists:
```bash
ls src/python/radia_vtk_export.py
```

---

## Contributing

When adding new examples:

1. **Follow naming conventions:**
   - Use descriptive names: `arc_current_with_magnet.py`
   - Avoid generic names: `test.py`, `example.py`

2. **Include VTK export:**
   - Use consistent VTK export pattern (see above)
   - Export geometry to `{script_name}.vtk`

3. **Add docstrings:**
   - Module-level docstring describing the example
   - Function docstrings for reusable functions

4. **Update README.md:**
   - Add example to subdirectory README.md
   - Include description and key topics

5. **Use new Material API:**
   - `MatLin()` for linear materials
   - `MatPM()` for permanent magnets
   - `MatSatIsoFrm()` for saturating materials

---

## References

- **Original Radia:** https://github.com/ochubar/Radia
- **ESRF Radia Documentation:** https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia
- **NGSolve:** https://ngsolve.org/
- **ParaView:** https://www.paraview.org/

---

**Last Updated:** 2025-11-12
**Total Examples:** 41 Python scripts
**Total Directories:** 9
**Documentation:** 100% coverage (all directories have README.md)
