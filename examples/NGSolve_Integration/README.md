# NGSolve Integration Examples

Integration of Radia (boundary element method) with NGSolve (finite element method) for magnetostatic field computation.

## Directory Structure

```
ngsolve_integration/
├── mesh_import/          # Mesh import functionality (Netgen/Cubit)
├── field_evaluation/     # Field evaluation tests
├── performance/          # Performance benchmarks
├── verification/         # Physics verification tests
├── demos/                # Simple demonstration scripts
├── utils/                # Visualization and export utilities
└── documentation/        # Technical documentation
```

## Quick Start

### 1. Mesh Import (Hexahedral/Tetrahedral)

Test mesh import from external tools:

```bash
cd mesh_import
python test_hex_mesh_import.py
```

**Features:**
- Nastran (.bdf) hexahedral mesh import
- Netgen tetrahedral mesh import
- Built-in primitive comparison

**→ See: [mesh_import/README.md](mesh_import/README.md)**

### 2. Field Evaluation

Test radia_ngsolve.RadiaField functionality:

```bash
cd field_evaluation
python test_batch_evaluation.py
python test_gridfunction_simple.py
```

**Features:**
- Batch field evaluation at multiple points
- GridFunction.Set() projection
- Coordinate transformation tests

**→ See: [field_evaluation/README.md](field_evaluation/README.md)**

### 3. Verification

Verify physics relationships (curl(A) = B):

```bash
cd verification
python verify_curl_A_equals_B.py
```

**Expected:** curl(A) = B within <5% error

**→ See: [verification/README.md](verification/README.md)**

### 4. Demos

Simple examples for learning:

```bash
cd demos
python demo_batch_evaluation.py
python demo_field_types.py
```

**→ See: [demos/README.md](demos/README.md)**

## Features

### radia_ngsolve Module

The `radia_ngsolve` Python module provides:

```python
from radia_ngsolve import RadiaField

# Create CoefficientFunction from Radia object
B_cf = RadiaField(radia_obj, 'b')  # Magnetic field
A_cf = RadiaField(radia_obj, 'a')  # Vector potential
H_cf = RadiaField(radia_obj, 'h')  # Magnetic field intensity

# Use in NGSolve
from ngsolve import *
gf = GridFunction(HDiv(mesh, order=2))
gf.Set(B_cf)  # Project field to GridFunction
```

**Supported field types:**
- `'b'` - Magnetic field (T)
- `'h'` - Magnetic field intensity (A/m)
- `'a'` - Vector potential (T·m)
- `'m'` - Magnetization (T)

### Mesh Import

Import external meshes to Radia:

```python
# Hexahedral mesh (Nastran format)
from nastran_mesh_import import create_radia_from_nastran
cube = create_radia_from_nastran('cube.bdf', units='m')

# Tetrahedral mesh (NGSolve)
from netgen_mesh_import import netgen_mesh_to_radia
cube = netgen_mesh_to_radia(ngsolve_mesh, units='m')
```

**→ See: [mesh_import/README.md](mesh_import/README.md)**

## Performance

### H-Matrix Acceleration

For large problems (N > 200 elements):

```python
import radia as rad

# Enable H-matrix field evaluation
rad.SetHMatrixFieldEval(1, eps=1e-6)

# Batch evaluation with H-matrix
H_values = rad.FldBatch(obj, 'h', points, use_hmatrix=1)
```

**Performance:**
- O(N log N) complexity vs O(N²) for dense solver
- 10-100x speedup for large problems
- <1% accuracy loss with eps=1e-6

**→ See: [performance/README.md](performance/README.md)**

### GridFunction Performance

```python
# Efficient field projection
B_cf = RadiaField(magnet, 'b')
gf.Set(B_cf)  # Optimized batch evaluation
```

**→ See: [documentation/NGSOLVE_SET_VS_INTERPOLATE.md](documentation/NGSOLVE_SET_VS_INTERPOLATE.md)**

## Best Practices

### Units

**Always use meters for NGSolve integration:**

```python
import radia as rad
rad.FldUnits('m')  # REQUIRED for NGSolve integration
```

NGSolve uses SI units (meters), so Radia must match.

### Finite Element Spaces

**Correct spaces for electromagnetic fields:**

```python
from ngsolve import *

# Vector potential (A) → HCurl
A_space = HCurl(mesh, order=2)
A_gf = GridFunction(A_space)
A_gf.Set(RadiaField(magnet, 'a'))

# Magnetic field (B) → HDiv
B_space = HDiv(mesh, order=2)
B_gf = GridFunction(B_space)
B_gf.Set(RadiaField(magnet, 'b'))
```

**Why:**
- HCurl: Ensures tangential continuity (correct for A)
- HDiv: Ensures normal continuity (correct for B)

### Mesh Resolution

**Field evaluation accuracy depends on mesh size:**

| Distance from magnet | Required mesh size | Expected error |
|---------------------|-------------------|----------------|
| <1 mesh cell | N/A | >10% (avoid) |
| >1 mesh cell | h < 0.015m | <1% |
| >5 mesh cells | h < 0.03m | <0.5% |

**Rule:** Evaluate GridFunction at distances > 1 mesh cell from magnet surfaces.

**→ See: [verification/README.md](verification/README.md)**

## Troubleshooting

### Large errors (>10%)

**Check:**
1. Units: `rad.FldUnits('m')` set?
2. Mesh size: h < 0.015m for 0.1m magnet?
3. Evaluation points: >1 mesh cell from boundaries?
4. FE space: HCurl for A, HDiv for B?

### ModuleNotFoundError: radia_ngsolve

**Cause:** Module not built or not in path.

**Solution:**
```bash
# Build radia_ngsolve module
cd S:/Radia/01_GitHub
cmake --build build --config Release --target radia_ngsolve

# Add to path
import sys
sys.path.insert(0, 'S:/Radia/01_GitHub/build/Release')
```

### GridFunction.Set() hangs

**Cause:** Very fine mesh with many DOFs.

**Solution:**
- Reduce mesh resolution (increase `maxh`)
- Use H-matrix acceleration for large Radia objects
- Check memory usage

## Documentation

Detailed technical documentation:

**→ See: [documentation/INDEX.md](documentation/INDEX.md)**

**Key documents:**
- GridFunction projection: [NGSOLVE_SET_VS_INTERPOLATE.md](documentation/NGSOLVE_SET_VS_INTERPOLATE.md)
- H-matrix analysis: [HMATRIX_ANALYSIS.md](documentation/HMATRIX_ANALYSIS.md)
- Troubleshooting: [HMATRIX_FIELD_EVALUATION_ISSUE.md](documentation/HMATRIX_FIELD_EVALUATION_ISSUE.md)

## Future Directions

Planned additions to `examples/ngsolve_integration/`:

1. **h_formulation/** - H-formulation comparison
   - Compare NGSolve H-formulation solver with Radia
   - Benchmark accuracy and performance
   - Hybrid solver workflows

2. **magnetization_import/** - Import NGSolve-computed magnetization
   - Read magnetization from NGSolve GridFunction
   - Apply to Radia geometry
   - Coupled Radia-NGSolve simulations

## Contributing

When adding new examples:
1. Choose appropriate subdirectory
2. Add README.md if creating new category
3. Update this main README.md
4. Follow existing code style (relative paths, error handling)

## Related

- `src/python/radia_ngsolve.cpp` - C++ pybind11 implementation
- `src/python/netgen_mesh_import.py` - Tetrahedral mesh importer
- `src/python/nastran_mesh_import.py` - Hexahedral mesh importer
- `tests/` - Unit tests for integration features

---

**Author**: Radia Development Team
**Last Updated**: 2025-11-22
