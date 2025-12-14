# Radia Documentation

**Version:** 1.3.14
**Date:** 2025-12-15

This folder contains the official documentation for Radia.

## Documentation Organization

### User Guides (Start Here)

| Document | Description |
|----------|-------------|
| [MSC_QUICK_START.md](MSC_QUICK_START.md) | **Quick Start Guide** for mesh-based computation (ObjThckPgn, ObjPolyhdr, Netgen) |
| [SUPPORTED_ELEMENTS.md](SUPPORTED_ELEMENTS.md) | Supported magnetic element types (tet, hex, rectangular) |
| [SOLVER_METHODS.md](SOLVER_METHODS.md) | Solver methods (LU, BiCGSTAB) |
| [API_REFERENCE.md](API_REFERENCE.md) | Complete Python API reference |

### NGSolve Integration

| Document | Description |
|----------|-------------|
| [NGSOLVE_USAGE_GUIDE.md](NGSOLVE_USAGE_GUIDE.md) | How to use Radia with NGSolve |
| [NGSOLVE_INTEGRATION.md](NGSOLVE_INTEGRATION.md) | Integration overview |

### Implementation Details (Advanced)

| Document | Description |
|----------|-------------|
| [MMM_MSC_IMPLEMENTATION.md](MMM_MSC_IMPLEMENTATION.md) | MMM/MSC implementation details |
| [MESH_MSC_API_DESIGN.md](MESH_MSC_API_DESIGN.md) | Mesh import API design |
| [MATERIAL_API_IMPLEMENTATION.md](MATERIAL_API_IMPLEMENTATION.md) | Material property implementation |
| [CF_BACKGROUND_FIELD_IMPLEMENTATION.md](CF_BACKGROUND_FIELD_IMPLEMENTATION.md) | Background field implementation |
| [NGSOLVE_CF_BACKGROUND_FIELD_DESIGN.md](NGSOLVE_CF_BACKGROUND_FIELD_DESIGN.md) | NGSolve CoefficientFunction design |
| [API_EXTENSIONS.md](API_EXTENSIONS.md) | Extended features and new APIs |

### Research Documentation

| Document | Description |
|----------|-------------|
| [HMATRIX_EVALUATION.md](HMATRIX_EVALUATION.md) | H-matrix evaluation (NOT recommended for typical use) |

---

## Quick Start

### 1. Supported Element Types

| Element Type | API | Status |
|--------------|-----|--------|
| Axis-aligned rectangular | `rad.ObjRecMag()` + `rad.ObjDivMag()` | **Supported** |
| Hexahedron (6-face MSC) | `rad.ObjThckPgn()` or `rad.ObjPolyhdr()` | **Supported** |
| Tetrahedron (4-face MSC) | `rad.ObjPolyhdr()` with `TETRA_FACES` | **Supported** |

### 2. Solver Methods

| Method | Name | Best For |
|--------|------|----------|
| LU Direct | `rad.Solve(grp, tol, max_iter, 0)` | Small problems (N < 1000) |
| BiCGSTAB | `rad.Solve(grp, tol, max_iter, 1)` | General purpose (default) |

### 3. Basic Example

```python
import radia as rad
import numpy as np

rad.FldUnits('m')

# Create soft iron cube
cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
rad.ObjDivMag(cube, [5, 5, 5])  # 125 elements

# Apply material
mat = rad.MatLin(999.0)  # mu_r = 1000
rad.MatApl(cube, mat)

# External field
MU_0 = 4 * np.pi * 1e-7
ext = rad.ObjBckg([0, 0, MU_0 * 50000])  # 50,000 A/m
grp = rad.ObjCnt([cube, ext])

# Solve
rad.Solve(grp, 0.001, 1000, 1)

# Get field
B = rad.Fld(grp, 'b', [0, 0, 0])
print(f"B at center: {B}")
```

---

## Additional Resources

### Examples
Working code examples are in [../examples/](../examples/):
- [cube_uniform_field/](../examples/cube_uniform_field/) - **Recommended**: Cube in uniform field benchmarks
- [simple_problems/](../examples/simple_problems/) - Basic examples
- [ngsolve_integration/](../examples/ngsolve_integration/) - NGSolve integration examples

### Tests
Automated tests are in [../tests/](../tests/).

## See Also

- [Main README](../README.md) - Project overview and installation
- [CHANGELOG](../CHANGELOG.md) - Version history

---

**Last Updated:** 2025-12-15
**Project:** Radia Magnetic Field Computation
