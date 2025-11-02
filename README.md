# Radia - Python Edition with OpenMP

3D Magnetostatics Computer Code - Optimized for Python 3.12 with OpenMP Parallelization

## Overview

This is a modernized version of Radia focusing on Python integration with performance optimizations:

- **Python 3.12 only** - Streamlined for modern Python
- **OpenMP 2.0 parallelization** - 2.7x speedup on 8-core systems
- **NGSolve integration** - C++ CoefficientFunction for FEM coupling
- **CMake build system** - Modern, cross-platform build
- **Tab indentation** - Consistent code style throughout
- **PyVista viewer** - Modern 3D visualization alternative

## Key Features

- ✓ OpenMP parallel field computation
- ✓ NGSolve C++ CoefficientFunction integration
- ✓ VTK export functionality (`rad.ObjDrwVTK`)
- ✓ PyVista-based 3D viewer (replaces OpenGL viewer)
- ✓ All hexahedron tests passing
- ✓ Comprehensive test suite and benchmarks
- ✓ Removed legacy components (Igor, Mathematica, GLUT, MPI)

## Quick Start

### Installation from PyPI

```bash
pip install radia-ngsolve
```

**Note**: The PyPI package includes pre-built binaries for Windows Python 3.12.

### Build from Source

```bash
# Windows (PowerShell)

# 1. Build radia.pyd (core module)
.\Build.ps1

# 2. Build rad_ngsolve.pyd (optional, for NGSolve integration)
.\build_ngsolve.ps1

# Outputs:
# - dist/radia.pyd
# - build/Release/rad_ngsolve.pyd
```

See [README_BUILD.md](README_BUILD.md) for detailed build instructions.

### Basic Usage

```python
import radia as rad

# Create a rectangular magnet
mag = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])

# Calculate field
field = rad.Fld(mag, 'b', [0,0,20])
print(f"Field: {field} T")
```

### NGSolve Integration

```python
# IMPORTANT: Import ngsolve first
import ngsolve
from ngsolve import Mesh, H1, GridFunction

import radia as rad
import rad_ngsolve

# Create Radia magnet
magnet = rad.ObjRecMag([0,0,0], [20,20,20], [0,0,1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

# Create NGSolve CoefficientFunction
B_field = rad_ngsolve.RadBfield(magnet)

# Use in FEM analysis
gf = GridFunction(fes)
gf.Set(B_field)
```

See [examples/Radia_to_NGSolve_CoefficientFunction/](examples/Radia_to_NGSolve_CoefficientFunction/) for complete examples.

## Documentation

### Core
- [README_BUILD.md](README_BUILD.md) - Build instructions
- [docs/OPENMP_PERFORMANCE_REPORT.md](docs/OPENMP_PERFORMANCE_REPORT.md) - OpenMP benchmarks
- [docs/DIRECTORY_STRUCTURE.md](docs/DIRECTORY_STRUCTURE.md) - Project structure

### NGSolve Integration
- [RAD_NGSOLVE_BUILD_SUCCESS.md](RAD_NGSOLVE_BUILD_SUCCESS.md) - Complete rad_ngsolve documentation
- [examples/Radia_to_NGSolve_CoefficientFunction/README.md](examples/Radia_to_NGSolve_CoefficientFunction/README.md) - NGSolve examples overview
- [examples/Radia_to_NGSolve_CoefficientFunction/EXAMPLES_GUIDE.md](examples/Radia_to_NGSolve_CoefficientFunction/EXAMPLES_GUIDE.md) - Detailed usage guide
- [tests/test_rad_ngsolve.py](tests/test_rad_ngsolve.py) - Integration tests

### Visualization
- [docs/PYVISTA_VIEWER.md](examples/2024_02_03_振分電磁石/PYVISTA_VIEWER.md) - PyVista viewer guide

### Development
- [docs/TAB_CONVERSION_REPORT.md](docs/TAB_CONVERSION_REPORT.md) - Code style conversion
- [docs/CLAUDE.md](docs/CLAUDE.md) - Development notes

## Performance

OpenMP parallelization results (8-core system):

| Threads | Time (sec) | Speedup |
|---------|-----------|---------|
| 1       | 11.67     | 1.00x   |
| 2       | 6.18      | 1.89x   |
| 4       | 3.57      | 3.27x   |
| 8       | 4.33      | 2.70x   |

See [docs/OPENMP_PERFORMANCE_REPORT.md](docs/OPENMP_PERFORMANCE_REPORT.md) for details.

## Examples

Practical examples are available in the `examples/` directory:


### NGSolve Integration
- `examples/Radia_to_NGSolve_CoefficientFunction/` - **Radia → NGSolve: Use Radia fields in FEM**
  - `demo_field_types.py` - All field types demonstration
  - `visualize_field.py` - Field visualization and comparison
  - `export_radia_geometry.py` - Export geometry to VTK

- `examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField/` - **NGSolve → Radia: Background fields**
  - `test_sphere_in_quadrupole.py` - Magnetizable sphere in quadrupole field
  - Uses Python callbacks to define arbitrary background fields

### Magnetostatics
- `examples/simple_problems/` - Basic magnet configurations
- `examples/electromagnet/` - Electromagnet designs
- `examples/complex_coil_geometry/` - Advanced coil geometries

### Legacy Examples
- `examples/2024_02_03_振分電磁石/` - Septum magnet simulation
- `examples/2024_03_02_Rdaiaの6面隊が動作しない調査(要素)/` - Hexahedron tests

## Testing

```bash
# Quick basic test
python tests/test_simple.py

# Comprehensive test suite
python tests/test_radia.py

# Advanced features test
python tests/test_advanced.py

# NGSolve integration test
python tests/test_rad_ngsolve.py

# OpenMP performance test
python tests/test_parallel_performance.py

# Or use pytest to run all tests
pytest tests/

# Run specific test suite
pytest tests/test_rad_ngsolve.py -v
```

See [tests/README.md](tests/README.md) for detailed testing documentation.

## Visualization

### PyVista Viewer (Recommended)

```python
import sys
sys.path.insert(0, 'S:/radia/01_GitHub/src/python')

from radia_pyvista_viewer import view_radia_object

# View Radia object
view_radia_object(mag)
```

### VTK Export

```python
import sys
sys.path.insert(0, 'S:/radia/01_GitHub/src/python')

from radia_vtk_export import export_geometry_to_vtk

# Export to VTK file for Paraview
export_geometry_to_vtk(mag, 'geometry.vtk')
```

## System Requirements

### Build Requirements
- Visual Studio 2022 (MSVC 19.44 or later)
- CMake 3.15 or later
- Python 3.12
- OpenMP 2.0

### Runtime Requirements
- Python 3.12
- NumPy
- NGSolve (optional, for FEM coupling via rad_ngsolve)
- PyVista (optional, for 3D visualization)

## Changes from Original Radia

### Removed Components
- Igor Pro integration
- Mathematica/Wolfram Language integration
- GLUT OpenGL viewer
- MPI support
- C client
- Python 2.7, 3.6, 3.7, 3.8 support

### Added Features
- OpenMP parallelization (2.7x speedup)
- NGSolve C++ CoefficientFunction integration
- PyVista viewer support
- Modern CMake build system
- Comprehensive test suite (including NGSolve tests)
- Performance benchmarks
- Updated documentation

## License

This project is licensed under **LGPL-2.1** (to match NGSolve licensing).

The original RADIA code is licensed under a **BSD-style license** by the European Synchrotron Radiation Facility (ESRF), Copyright © 1997-2018.

Both licenses are included in the `LICENSE` file. The BSD-style license is compatible with and subsumed by the LGPL-2.1 license used for this derivative work.

See:
- `LICENSE` - Complete license text (LGPL-2.1 + Original RADIA BSD-style)
- `COPYRIGHT.txt` - Original RADIA copyright notice

## Credits

**Original Radia**: Pascal Elleaume, Oleg Chubar, and others at ESRF

**This Fork**:
- OpenMP parallelization
- NGSolve C++ integration (rad_ngsolve)
- Python 3.12 optimization
- Build system modernization
- PyVista integration
- Documentation updates

## Related Tools

### Coreform Cubit Mesh Export
[**Coreform_Cubit_Mesh_Export**](https://github.com/ksugahar/Coreform_Cubit_Mesh_Export) - Python library for exporting Coreform Cubit meshes to multiple formats

This tool perfectly complements Radia_NGSolve by providing high-quality mesh generation:

- **Nastran Format Export** - Compatible with Radia's `nastran_reader.py` module
- **Multiple Format Support** - Gmsh, MEG, VTK, and Nastran exports
- **2D/3D Meshing** - Supports complex 3D geometries
- **Second-Order Elements** - High-accuracy mesh generation

**Workflow Example:**
1. Create complex geometry in Coreform Cubit
2. Export to Nastran format using `Coreform_Cubit_Mesh_Export`
3. Import into Radia using `nastran_reader.py`
4. Couple with NGSolve for FEM analysis

See [examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField/](examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField/) for Nastran mesh usage examples.

## Links

- Original Radia: https://github.com/ochubar/Radia
- ESRF Radia Page: https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia
- Coreform Cubit: https://coreform.com/products/coreform-cubit/

---

**Version**: 1.0 (OpenMP + NGSolve Edition)
**Last Updated**: 2025-11-02
