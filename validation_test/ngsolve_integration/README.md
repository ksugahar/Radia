# NGSolve Integration Examples

Integration of Radia HDiv-VIM with NGSolve finite-element workflows for magnetostatic field computation.

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
- Cubit -> Netgen direct hex mesh import
- Netgen tetrahedral mesh import
- Built-in primitive comparison

**→ See: [mesh_import/README.md](mesh_import/README.md)**

### 2. Field Evaluation

Test rad.RadiaField functionality:

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

### RadiaField (integrated into main radia module since v2.5.0)

`RadiaField` is accessed as `rad.RadiaField()` from the main radia module:

```python
import radia as rad

# Create CoefficientFunction from Radia object
B_cf = rad.RadiaField(radia_obj, 'b')  # Magnetic field
A_cf = rad.RadiaField(radia_obj, 'a')  # Vector potential
H_cf = rad.RadiaField(radia_obj, 'h')  # Magnetic field intensity

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
# Hexahedral mesh (via Cubit -> Netgen direct export)
import cubit
import cubit_mesh_export
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import geometry 'model.step'")
cubit.cmd("mesh volume all")
ngmesh = cubit_mesh_export.extract_curved_mesh()
mesh = Mesh(ngmesh)
cube = netgen_mesh_to_radia(mesh, material={'magnetization': [0,0,0]}, units='m')

# Tetrahedral mesh (NGSolve)
from netgen_mesh_import import netgen_mesh_to_radia, extract_elements, compute_element_centroid
cube = netgen_mesh_to_radia(ngsolve_mesh, units='m')

# Custom processing with extract_elements
elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']  # Correctly extracted coordinates
    centroid = compute_element_centroid(vertices)
```

**NOTE**: Nastran BDF format is REMOVED. Use Cubit -> Netgen direct export for hex meshes.

**CRITICAL POLICY - NGSolve Mesh Access**:

| Rule | Description |
|------|-------------|
| **ALWAYS** | Use functions from `netgen_mesh_import.py` |
| **NEVER** | Directly access `mesh.ngmesh.Points()`, `mesh.vertices[]`, or `el.vertices[].nr` |
| **NO EXCEPTIONS** | Applies to all scripts including examples, tests, and debugging code |

**Why?** NGSolve has TWO different indexing schemes:
- `mesh.ngmesh.Points()[i]` is **1-indexed** (valid: 1 to nv)
- `mesh.vertices[i]` is **0-indexed** (valid: 0 to nv-1)
- `el.vertices[i].nr` returns **0-indexed** value (for use with `mesh.vertices[]` only)

Mixing these causes off-by-one errors that are difficult to debug.

**→ See: [mesh_import/README.md](mesh_import/README.md)**

## Performance

### H-Matrix Acceleration

For large problems (N > 200 elements):

```python
import radia as rad

# Batch field evaluation (auto-detected from input shape)
H_values = rad.Fld(obj, 'h', points)  # points shape (N, 3) -> batch
```

**Performance:**
- O(N log N) complexity vs O(N²) for dense solver
- 10-100x speedup for large problems
- <1% accuracy loss with eps=1e-6

**→ See: [performance/README.md](performance/README.md)**

### GridFunction Performance

```python
# Efficient field projection
B_cf = rad.RadiaField(magnet, 'b')
gf.Set(B_cf)  # Optimized batch evaluation
```

**→ See: [documentation/NGSOLVE_SET_VS_INTERPOLATE.md](documentation/NGSOLVE_SET_VS_INTERPOLATE.md)**

## Best Practices

### Units

**Radia always uses meters, which is compatible with NGSolve SI units.**

### Finite Element Spaces

**Correct spaces for electromagnetic fields:**

```python
from ngsolve import *

# Vector potential (A) → HCurl
A_space = HCurl(mesh, order=2)
A_gf = GridFunction(A_space)
A_gf.Set(rad.RadiaField(magnet, 'a'))

# Magnetic field (B) → HDiv
B_space = HDiv(mesh, order=2)
B_gf = GridFunction(B_space)
B_gf.Set(rad.RadiaField(magnet, 'b'))
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
1. Units: all dimensions in meters?
2. Mesh size: h < 0.015m for 0.1m magnet?
3. Evaluation points: >1 mesh cell from boundaries?
4. FE space: HCurl for A, HDiv for B?

### RadiaField not found

**Cause:** Radia module not built with NGSolve support, or outdated version.

**Solution:** Since v2.5.0, `RadiaField` is integrated into the main `radia` module:
```python
import radia as rad
B_cf = rad.RadiaField(magnet, 'b')  # No separate module needed
```

If using an older version, rebuild:
```bash
cd S:/Radia/01_GitHub
.\Build.ps1
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

Planned additions to `validation_test/ngsolve_integration/`:

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

- `src/radia/radia_pybind.cpp` - C++ pybind11 implementation (includes RadiaField)
- `src/radia/netgen_mesh_import.py` - Mesh importer (tet/hex)
- `src/cubit_plugin/` - Cubit mesh export plugin (C++ + pybind11)
- `tests/` - Unit tests for integration features

---

**Author**: Radia Development Team
**Last Updated**: 2026-01-16
**Mesh Workflow**: Cubit -> Netgen (direct) -> Radia
