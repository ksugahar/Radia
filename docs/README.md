# Radia Documentation

**Version:** 1.3.13
**Date:** 2025-12-11

This folder contains the official documentation for Radia.

## Documentation Organization

### Essential Documentation

#### Getting Started
- [SUPPORTED_ELEMENTS.md](SUPPORTED_ELEMENTS.md) - **Supported magnetic element types** (tet, hex, rectangular)
- [SOLVER_METHODS.md](SOLVER_METHODS.md) - Solver methods (LU, BiCGSTAB)
- [API_REFERENCE.md](API_REFERENCE.md) - Complete Python API reference

#### NGSolve Integration
- [NGSOLVE_USAGE_GUIDE.md](NGSOLVE_USAGE_GUIDE.md) - How to use Radia with NGSolve
- [NGSOLVE_INTEGRATION.md](NGSOLVE_INTEGRATION.md) - Integration overview

#### Mesh Import
- [MESH_MSC_API_DESIGN.md](MESH_MSC_API_DESIGN.md) - Mesh import API design
- [MMM_MSC_IMPLEMENTATION.md](MMM_MSC_IMPLEMENTATION.md) - MMM/MSC implementation details

### Technical Documentation

#### Implementation Details
- [MATERIAL_API_IMPLEMENTATION.md](MATERIAL_API_IMPLEMENTATION.md) - Material property implementation
- [CF_BACKGROUND_FIELD_IMPLEMENTATION.md](CF_BACKGROUND_FIELD_IMPLEMENTATION.md) - Background field implementation
- [NGSOLVE_CF_BACKGROUND_FIELD_DESIGN.md](NGSOLVE_CF_BACKGROUND_FIELD_DESIGN.md) - NGSolve CoefficientFunction design
- [API_EXTENSIONS.md](API_EXTENSIONS.md) - Extended features and new APIs

#### H-Matrix Evaluation (Research)
- [HMATRIX_EVALUATION.md](HMATRIX_EVALUATION.md) - **H-matrix evaluation results** (NOT recommended for typical use)

**Note:** H-matrix acceleration was evaluated and found to provide no benefit for typical Radia use cases (single compact objects). The following docs are archived for reference:
- [HMATRIX_USER_GUIDE.md](HMATRIX_USER_GUIDE.md) - (Archived) User guide
- [HMATRIX_SERIALIZATION.md](HMATRIX_SERIALIZATION.md) - (Archived) Disk cache feature

## Quick Start

### 1. Supported Element Types

Radia currently supports:

| Element Type | API | Status |
|--------------|-----|--------|
| Axis-aligned rectangular | `rad.ObjRecMag()` + `rad.ObjDivMag()` | **Supported** |
| Tetrahedron (4-face MSC) | `rad.ObjPolyhdr()` with `TETRA_FACES` | **Supported** |
| Hexahedron (6-face MSC) | `rad.ObjPolyhdr()` with `HEX_FACES` | **Supported** |
| General polyhedra (>6 faces) | - | **Not Supported** |

See [SUPPORTED_ELEMENTS.md](SUPPORTED_ELEMENTS.md) for details.

### 2. Solver Methods

| Method | Name | Best For |
|--------|------|----------|
| LU Direct | `'lu'` or `0` | Small problems (N < 500) |
| BiCGSTAB | `'bicgstab'` or `1` | General purpose (default) |

See [SOLVER_METHODS.md](SOLVER_METHODS.md) for details.

### 3. Basic Example

```python
import radia as rad
rad.FldUnits('m')

# Create soft iron cube
cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
rad.ObjDivMag(cube, [5, 5, 5])  # 125 elements

# Apply material
mat = rad.MatLin(999.0)  # mu_r = 1000
rad.MatApl(cube, mat)

# External field
ext = rad.ObjBckg([0, 0, 1.0])  # 1 T
grp = rad.ObjCnt([cube, ext])

# Solve
rad.Solve(grp, 0.0001, 1000)

# Get field
B = rad.Fld(grp, 'b', [0, 0, 0])
print(f"B at center: {B}")
```

## Additional Resources

### Examples
Working code examples are in [../examples/](../examples/):
- [simple_problems/](../examples/simple_problems/) - Basic examples
- [solver_benchmarks/](../examples/solver_benchmarks/) - Performance benchmarks
- [ngsolve_integration/](../examples/ngsolve_integration/) - NGSolve integration examples
- [cube_uniform_field/](../examples/cube_uniform_field/) - Cube in uniform field benchmarks

### Tests
Automated tests are in [../tests/](../tests/).

## See Also

- [Main README](../README.md) - Project overview and installation
- [CHANGELOG](../CHANGELOG.md) - Version history

---

**Last Updated:** 2025-12-11
**Project:** Radia Magnetic Field Computation
