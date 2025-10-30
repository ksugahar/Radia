# Radia - Python Edition with OpenMP

3D Magnetostatics Computer Code - Optimized for Python 3.12 with OpenMP Parallelization

## Overview

This is a modernized version of Radia focusing on Python integration with performance optimizations:

- **Python 3.12 only** - Streamlined for modern Python
- **OpenMP 2.0 parallelization** - 2.7x speedup on 8-core systems
- **CMake build system** - Modern, cross-platform build
- **Tab indentation** - Consistent code style throughout
- **PyVista viewer** - Modern 3D visualization alternative

## Key Features

- ✓ OpenMP parallel field computation
- ✓ VTK export functionality (`rad.ObjDrwVTK`)
- ✓ PyVista-based 3D viewer (replaces OpenGL viewer)
- ✓ All hexahedron tests passing
- ✓ Comprehensive test suite and benchmarks
- ✓ Removed legacy components (Igor, Mathematica, GLUT, MPI)

## Quick Start

### Build

```bash
# Windows (PowerShell)
.\build.ps1

# The module will be available in dist/radia.pyd
```

See [README_BUILD.md](README_BUILD.md) for detailed build instructions.

### Usage

```python
import radia as rad

# Create a rectangular magnet
mag = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])

# Calculate field
field = rad.Fld(mag, 'b', [0,0,20])
print(f"Field: {field} T")
```

## Documentation

- [README_BUILD.md](README_BUILD.md) - Build instructions
- [docs/OPENMP_PERFORMANCE_REPORT.md](docs/OPENMP_PERFORMANCE_REPORT.md) - OpenMP benchmarks
- [docs/PYVISTA_VIEWER.md](examples/2024_02_03_振分電磁石/PYVISTA_VIEWER.md) - PyVista viewer guide
- [docs/DIRECTORY_STRUCTURE.md](docs/DIRECTORY_STRUCTURE.md) - Project structure
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

- `examples/2019_11_29_Radia_練習/` - Converted from Wolfram Language
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

# OpenMP performance test
python tests/test_parallel_performance.py

# Or use pytest to run all tests
pytest tests/
```

See [tests/README.md](tests/README.md) for detailed testing documentation.

## Visualization

### PyVista Viewer (Recommended)

```python
from examples.2024_02_03_振分電磁石.radia_pyvista_viewer import view_radia_object

# View Radia object
view_radia_object(mag)
```

### VTK Export

```python
from examples.2024_02_03_振分電磁石.my_module import exportGeometryToVTK

# Export to VTK file for Paraview
exportGeometryToVTK(mag, 'geometry')
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
- OpenMP parallelization
- PyVista viewer support
- Modern CMake build system
- Comprehensive test suite
- Performance benchmarks
- Updated documentation

## License

Original Radia license applies. See source files for details.

## Credits

**Original Radia**: Pascal Elleaume, Oleg Chubar, and others at ESRF

**This Fork**:
- OpenMP parallelization
- Python 3.12 optimization
- Build system modernization
- PyVista integration
- Documentation updates

🤖 Modernization performed with [Claude Code](https://claude.com/claude-code)

## Links

- Original Radia: https://github.com/ochubar/Radia
- ESRF Radia Page: https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia

---

**Version**: 4.32 (OpenMP Edition)
**Last Updated**: 2025-10-29
