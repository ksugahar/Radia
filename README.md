# Radia - Integrated Electromagnetic Field Environment

3D Magnetostatics and Conductor Analysis - Python 3.12 with OpenMP Parallelization

## Overview

This is more than a simple MMM (Magnetic Moment Method) solver. **Radia provides an integrated electromagnetic field environment** that combines multiple field sources and analysis methods:

### Integrated Field Sources

| Source Type | Method | Frequency Range | Use Case |
|-------------|--------|-----------------|----------|
| **Permanent Magnets** | Analytical (surface charge) | DC | NdFeB, SmCo magnets |
| **Soft Iron** | MMM with demagnetization | DC | Yokes, pole pieces |
| **Coil Currents** | Biot-Savart analytical | DC | Electromagnets, coils |
| **Conductor Analysis** | FastImp (SIBC/pFFT) | AC (kHz-GHz) | Impedance extraction |

### Key Capabilities

1. **Analytical Field Solutions**
   - Permanent magnets: Surface charge method with solid angle integration
   - Current loops/coils: Biot-Savart law implementation
   - Arc currents: Analytical formulas for arc segments

2. **MMM Demagnetization Solver**
   - Self-consistent solution for soft magnetic materials
   - Supports nonlinear B-H curves (`MatSatIsoTab`)
   - H-matrix acceleration for large problems (HACApK)

3. **NGSolve Integration (Reduced Potential Formulation)**
   - Radia fields as `CoefficientFunction` source terms
   - Enables: `curl(curl(A)) = mu_0 * J` with Radia providing J
   - Magnetizable regions computed via FEM with Radia background field

4. **High-Frequency Conductor Analysis (FastImp)**
   - Surface integral equation (EFIE) with pFFT acceleration
   - MQS/EMQS/Full-wave formulations
   - Impedance extraction for interconnects and coils

### Frequency-Dependent Workflow

```
DC Analysis:
  Magnets + Coils ─→ Radia Analytical ─→ NGSolve (reduced potential)
                                         └─→ Current distribution (FEM)

AC Analysis (kHz-GHz):
  Conductors ─→ FastImp (SIBC/pFFT) ─→ Impedance, skin effect
```

**For Radia usage and API documentation**, please refer to the [original Radia repository](https://github.com/ochubar/Radia) and [ESRF Radia documentation](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia).

This repository adds NGSolve integration, FastImp conductor analysis, and performance improvements while maintaining compatibility with the original Radia API.

## Key Features

### Magnetostatics
- ✓ **H-matrix acceleration** - Fast hierarchical matrix solver using [HACApK](https://github.com/Post-Peta-Crest/ppOpenHPC/tree/MATH/HACApK)
- ✓ OpenMP parallel field computation
- ✓ Analytical Biot-Savart for coils and arc currents
- ✓ Surface charge method for permanent magnets
- ✓ MMM demagnetization for soft iron

### Conductor Analysis (FastImp)
- ✓ **pFFT-accelerated BEM** - O(N log N) matrix-vector products
- ✓ Surface discretization: rectangular blocks, loops, wires, spirals
- ✓ DC/MQS/EMQS/Full-wave formulations
- ✓ Impedance extraction with port definition

### Integration & Visualization
- ✓ NGSolve C++ CoefficientFunction integration
- ✓ VTK export functionality (`rad.ObjDrwVTK`)
- ✓ PyVista-based 3D viewer (replaces OpenGL viewer)
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

# 2. Build radia_ngsolve.pyd (optional, for NGSolve integration)
.\Build_NGSolve.ps1

# Outputs:
# - dist/radia.pyd
# - build/Release/radia_ngsolve.pyd
```

See [BUILD.md](BUILD.md) for detailed build instructions.

### Basic Usage

```python
import radia as rad

# Create a hexahedral magnet using ObjHexahedron (8 vertices, faces auto-generated)
vertices = [[-5,-5,-5], [5,-5,-5], [5,5,-5], [-5,5,-5],
            [-5,-5,5], [5,-5,5], [5,5,5], [-5,5,5]]
mag = rad.ObjHexahedron(vertices, [0, 0, 1])  # magnetization in A/m

# Calculate field at external point
field = rad.Fld(mag, 'b', [0, 0, 20])
print(f"Field: {field} T")
```

**Note:** This example uses the standard Radia API. For complete Radia API documentation and examples, see:
- [Original Radia repository](https://github.com/ochubar/Radia)
- [ESRF Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)

### NGSolve Integration

The `radia_ngsolve` module provides a C++ CoefficientFunction interface for using Radia magnetic fields in NGSolve FEM analysis.

**Requirements:**
- NGSolve must be installed separately: `pip install ngsolve`
- **IMPORTANT**: NGSolve must be imported **before** `radia_ngsolve`
- Current build is linked against **NGSolve 6.2.2406** (Windows Python 3.12)

**Installation:**

```bash
# 1. Install NGSolve first
pip install ngsolve

# 2. Install Radia (includes radia_ngsolve.pyd)
pip install radia
```

**Function Specification:**

```python
radia_ngsolve.RadiaField(radia_obj, field_type='b')
```

**Parameters:**
- `radia_obj` (int): Radia object handle returned by `rad.ObjHexahedron()`, `rad.ObjTetrahedron()`, `rad.ObjThckPgn()`, etc.
- `field_type` (str, optional): Field type to compute. Default: `'b'`
  - `'b'`: Magnetic flux density [Tesla]
  - `'h'`: Magnetic field [A/m]
  - `'a'`: Vector potential [T·m]
  - `'m'`: Magnetization [A/m]

**Returns:**
- NGSolve CoefficientFunction (3D vector field)

**Unit Conversion:**
- NGSolve uses **meters**, Radia uses **millimeters**
- Automatic conversion: coordinates are multiplied by 1000 (m → mm)
- Field values remain in SI units (no conversion needed)

**Example:**

```python
# IMPORTANT: Import ngsolve FIRST (required to load NGSolve DLLs)
import ngsolve
from ngsolve import Mesh, H1, GridFunction, HDiv

import radia as rad
from radia import radia_ngsolve  # or: import radia_ngsolve

# Create Radia hexahedral magnet using ObjHexahedron (faces auto-generated)
vertices = [[-10,-10,-10], [10,-10,-10], [10,10,-10], [-10,10,-10],
            [-10,-10,10], [10,-10,10], [10,10,10], [-10,10,10]]  # mm units
Mr = 1.2 / (4 * 3.14159 * 1e-7)  # Br=1.2T -> Mr in A/m
magnet = rad.ObjHexahedron(vertices, [0, 0, Mr])

# Create NGSolve CoefficientFunction for different field types
B_field = radia_ngsolve.RadiaField(magnet, 'b')  # Flux density [T]
H_field = radia_ngsolve.RadiaField(magnet, 'h')  # Magnetic field [A/m]
A_field = radia_ngsolve.RadiaField(magnet, 'a')  # Vector potential [T·m]
M_field = radia_ngsolve.RadiaField(magnet, 'm')  # Magnetization [A/m]

# Use in FEM analysis (NGSolve mesh in meters)
gf = GridFunction(fes)
gf.Set(B_field)  # Automatically converts mesh coordinates m → mm
```

See [examples/ngsolve_integration/](examples/ngsolve_integration/) for complete examples.

## Documentation

### Build & Setup
- [BUILD.md](BUILD.md) - Build instructions (Windows, macOS, Linux)

### User Guides
- [docs/API_REFERENCE.md](docs/API_REFERENCE.md) - Complete Python API reference
- [docs/SOLVER_METHODS.md](docs/SOLVER_METHODS.md) - Solver methods (LU, BiCGSTAB, H-matrix)
- [docs/HMATRIX_USER_GUIDE.md](docs/HMATRIX_USER_GUIDE.md) - H-matrix acceleration guide

### Conductor Analysis (FastImp)
- [docs/FASTIMP_SIBC_INTEGRATION.md](docs/FASTIMP_SIBC_INTEGRATION.md) - FastImp integration guide
- Conductor APIs: `CndRecBlock`, `CndLoop`, `CndWire`, `CndSpiral`
- Frequency range: DC to GHz (MQS/EMQS/Full-wave)

### NGSolve Integration
- [docs/NGSOLVE_USAGE_GUIDE.md](docs/NGSOLVE_USAGE_GUIDE.md) - How to use Radia with NGSolve
- [docs/NGSOLVE_INTEGRATION.md](docs/NGSOLVE_INTEGRATION.md) - Integration overview

### Examples
- [examples/ngsolve_integration/](examples/ngsolve_integration/) - NGSolve integration examples
- [examples/magpylib_integration/](examples/magpylib_integration/) - magpylib background field examples
- [examples/simple_problems/](examples/simple_problems/) - Basic magnet configurations
- [tests/README.md](tests/README.md) - Test suite documentation

## Performance

OpenMP parallelization results (8-core system):

| Threads | Time (sec) | Speedup |
|---------|-----------|---------|
| 1       | 11.67     | 1.00x   |
| 2       | 6.18      | 1.89x   |
| 4       | 3.57      | 3.27x   |
| 8       | 4.33      | 2.70x   |

See [docs/SOLVER_METHODS.md](docs/SOLVER_METHODS.md) for solver performance details.

## Examples

Practical examples are available in the `examples/` directory:


### NGSolve Integration
- `examples/ngsolve_integration/` - **Radia → NGSolve: Use Radia fields in FEM**
  - `demo_field_types.py` - All field types demonstration
  - `visualize_field.py` - Field visualization and comparison
  - `export_radia_geometry.py` - Export geometry to VTK

- `examples/background_fields/` - **NGSolve → Radia: Background fields**
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
python tests/test_radia_ngsolve.py

# OpenMP performance test
python tests/test_parallel_performance.py

# Or use pytest to run all tests
pytest tests/

# Run specific test suite
pytest tests/test_radia_ngsolve.py -v
```

See [tests/README.md](tests/README.md) for detailed testing documentation.

## Visualization

### PyVista Viewer (Recommended)

```python
import sys
import os

# Add src/python to path (adjust path as needed for your installation)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src', 'python'))

from radia_pyvista_viewer import view_radia_object

# View Radia object
view_radia_object(mag)
```

### VTK Export

```python
import sys
import os

# Add src/python to path (adjust path as needed for your installation)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src', 'python'))

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
- NGSolve (optional, for FEM coupling via radia_ngsolve)
- PyVista (optional, for 3D visualization)

## Linear Material with External Field (ObjBckgCF)

Radia supports soft magnetic materials (linear materials) in external fields using `ObjBckgCF`.

### Working Example

```python
import radia as rad
import numpy as np

rad.UtiDelAll()
rad.FldUnits('m')  # Use SI units

# Create soft iron cube (10cm side) using ObjHexahedron (faces auto-generated)
vertices = [[-0.05,-0.05,-0.05], [0.05,-0.05,-0.05], [0.05,0.05,-0.05], [-0.05,0.05,-0.05],
            [-0.05,-0.05,0.05], [0.05,-0.05,0.05], [0.05,0.05,0.05], [-0.05,0.05,0.05]]
cube = rad.ObjHexahedron(vertices, [0, 0, 0])

# Apply isotropic linear material (mu_r = 100)
chi = 99.0  # chi = mu_r - 1
mat = rad.MatLin(chi)  # IMPORTANT: Use single argument for isotropic
rad.MatApl(cube, mat)

# Create external field (1 Tesla in z-direction)
def uniform_B(pos):
    return [0, 0, 1.0]  # Returns B in Tesla
bg = rad.ObjBckgCF(uniform_B)

# Create system and solve
system = rad.ObjCnt([cube, bg])
result = rad.Solve(system, 0.0001, 1000)

# Result: M ~ 2.3 MA/m (correct for cube with mu_r=100 in 1T field)
M = rad.ObjM(cube)
print(f"Magnetization: {M[1]} A/m")
```

### Important Notes

1. **Isotropic materials**: Always use `MatLin(chi)` (single argument)
   - Do NOT use `MatLin(chi, [0, 0, 1e-10])` - this creates a poorly defined anisotropic tensor

2. **Callback returns B in Tesla**: The callback function for `ObjBckgCF` should return magnetic flux density B in Tesla

3. **Units**: When using `rad.FldUnits('m')`, all coordinates are in meters

4. **Tetrahedral elements**: Also works with `ObjTetrahedron` tetrahedral meshes

### Supported Element Types

| Element | API | External Field Support |
|---------|-----|----------------------|
| Hexahedron | `ObjHexahedron` | Yes |
| Tetrahedron | `ObjTetrahedron` | Yes |
| Extruded Polygon | `ObjThckPgn` | Yes |

## Tetrahedral Mesh Import from NGSolve/Netgen

Radia can import tetrahedral meshes from NGSolve/Netgen for complex geometries.

**Features:**
- Direct Netgen/OCC mesh import to Radia MMM
- Independent mesh optimization for FEM vs MMM
- Nastran CTETRA face connectivity standard

**Usage:** See `examples/cube_uniform_field/linear/` for examples.

**Benchmark** (mu_r=100, B0=1T, 100mm cube):

| Mesh | Elements | M_z (MA/m) | Error |
|------|----------|------------|-------|
| 1x1x1 hex | 1 | 2.3171 | 0.00% |
| 2x2x2 hex | 8 | 2.4180 | 4.36% |
| 5x5x5 hex | 125 | 2.6652 | 15.02% |

Analytical: M = 2.3171 MA/m (N = 1/3 for cube)


## Changes from Original Radia

### Removed Components
- Igor Pro integration
- Mathematica/Wolfram Language integration
- GLUT OpenGL viewer
- MPI support
- C client
- Python 2.7, 3.6, 3.7, 3.8 support

### Added Features
- **Integrated electromagnetic environment** - Multiple field sources (magnets, coils, conductors)
- **FastImp conductor analysis** - pFFT-accelerated BEM for impedance extraction
- **H-matrix acceleration** - Hierarchical matrix solver using [HACApK library](https://github.com/Post-Peta-Crest/ppOpenHPC/tree/MATH/HACApK) (9-50x speedup for large problems)
- OpenMP parallelization (2.7x speedup)
- NGSolve C++ CoefficientFunction integration (reduced potential formulation support)
- PyVista viewer support
- Modern CMake build system
- Comprehensive test suite (including NGSolve tests)
- Performance benchmarks
- Updated documentation

## License

This project contains multiple components with different open-source licenses:

### RADIA Core
Licensed under a **BSD-style license** by the European Synchrotron Radiation Facility (ESRF).
- Copyright © 1997-2018, European Synchrotron Radiation Facility
- Permits redistribution and modification with attribution

### HACApK_LH-Cimplm (H-matrix Library)
Licensed under the **MIT License**.
- Copyright © 2015, Akihiro Ida and Takeshi Iwashita
- Location: `src/ext/HACApK_LH-Cimplm/`
- ppOpen-HPC project

Both licenses are permissive open-source licenses that allow free use, modification,
and redistribution.

See:
- `LICENSE` - Complete license text for both components
- `COPYRIGHT.txt` - Original RADIA BSD-style license
- `src/ext/HACApK_LH-Cimplm/LICENSE` - HACApK MIT license

## Credits

**Original Radia**: Pascal Elleaume, Oleg Chubar, and others at ESRF

**This Fork**:
- OpenMP parallelization
- NGSolve C++ integration (radia_ngsolve)
- Python 3.12 optimization
- Build system modernization
- PyVista integration
- Documentation updates

## Related Tools

### Coreform Cubit Mesh Export
[**Coreform_Cubit_Mesh_Export**](https://github.com/ksugahar/Coreform_Cubit_Mesh_Export) - Python library for exporting Coreform Cubit meshes to multiple formats

This tool perfectly complements Radia_NGSolve by providing high-quality mesh generation:

- **Nastran Format Export** - Compatible with Radia's `nastran_mesh_import.py` module
- **Multiple Format Support** - Gmsh, MEG, VTK, and Nastran exports
- **2D/3D Meshing** - Supports complex 3D geometries
- **Second-Order Elements** - High-accuracy mesh generation

**Workflow Example:**
1. Create complex geometry in Coreform Cubit
2. Export to Nastran format using `Coreform_Cubit_Mesh_Export`
3. Import into Radia using `nastran_mesh_import.py`
4. Couple with NGSolve for FEM analysis

See [examples/background_fields/](examples/background_fields/) for Nastran mesh usage examples.

## Links

- Original Radia: https://github.com/ochubar/Radia
- ESRF Radia Page: https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia
- Coreform Cubit: https://coreform.com/products/coreform-cubit/

---

**Version**: 1.4.4 (Integrated EM Environment)
**Last Updated**: 2026-01-08
