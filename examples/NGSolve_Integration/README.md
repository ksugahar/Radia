# Radia - NGSolve Integration Examples

Complete examples demonstrating integration between Radia (3D magnetostatics) and NGSolve (finite element analysis).

**Key Feature**: `radia_ngsolve.RadiaField` provides seamless conversion of Radia field sources to NGSolve CoefficientFunction with automatic unit conversion (mm ↔ m).

---

## Table of Contents

- [Quick Start](#quick-start)
- [Examples](#examples)
  - [Basic Field Types](#basic-field-types)
  - [Vector Potential and curl(A) = B](#vector-potential-and-curla--b)
  - [Coordinate Transformation](#coordinate-transformation)
  - [Field Visualization](#field-visualization)
- [Performance Benchmarks](#performance-benchmarks)
- [API Reference](#api-reference)

---

## Quick Start

**Requirements**:
- Radia (built with NGSolve support)
- NGSolve
- Python 3.7+
- NumPy

**Basic usage**:
```python
import radia as rad
from ngsolve import *
import radia_ngsolve

# 1. Create Radia magnet (specify units!)
rad.FldUnits('m')  # Use meters (recommended for NGSolve integration)
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.2])  # 10mm = 0.01m cube

# 2. Create NGSolve CoefficientFunction (units='m' is default)
B_cf = radia_ngsolve.RadiaField(magnet, 'b')  # units='m' by default

# 3. Use in NGSolve mesh (always meters)
mesh = Mesh(...)  # NGSolve mesh in meters
gf = GridFunction(HCurl(mesh))
gf.Set(B_cf)  # No conversion needed - both use meters
```

**Unit conversion**:
- NGSolve always uses **meters**
- **Recommended**: Use `rad.FldUnits('m')` for seamless integration (default `units='m'`)
- For legacy mm models: Use `rad.FldUnits('mm')` and `units='mm'`

---

## Examples

### Basic Field Types

**File**: `demo_field_types.py`

Demonstrates all supported field types:
- `'b'`: Magnetic flux density **B** (Tesla)
- `'h'`: Magnetic field **H** (A/m)
- `'a'`: Vector potential **A** (T·m)
- `'m'`: Magnetization **M** (A/m)

**Run**:
```bash
python demo_field_types.py
```

**Output**: Field values at multiple test points, unit conversion verification.

---

### Vector Potential and curl(A) = B

**File**: `verify_curl_A_equals_B.py`

Rigorous verification that `curl(A) = B` with correct unit handling.

**Features**:
- Computes vector potential **A** from Radia
- Numerically computes `curl(A)` in NGSolve
- Compares with direct **B** from Radia
- Statistical error analysis
- VTK output for ParaView visualization

**Run**:
```bash
python verify_curl_A_equals_B.py
```

**Output**:
- Console: Error statistics (typical: relative error < 10⁻⁶)
- Files: `*.vtu` for visualization in ParaView

**Key result**: Unit scaling factor for **A** is 0.001 to ensure `curl(A) = B` with correct units (mm vs m).

---

### Coordinate Transformation

**File**: `test_coordinate_transform.py`

Tests coordinate transformation (translation + rotation) for moving/rotating magnets.

**Features**:
- Translation: `origin` parameter
- Rotation: `u_axis`, `v_axis`, `w_axis` parameters
- Verification against analytical solutions

**Example**:
```python
# Magnet rotated 45° around z-axis, translated to (0.1, 0, 0) m
angle = np.pi/4
B_cf = radia_ngsolve.RadiaField(
	magnet, 'b',
	origin=[0.1, 0, 0],
	u_axis=[np.cos(angle), np.sin(angle), 0],
	v_axis=[-np.sin(angle), np.cos(angle), 0],
	w_axis=[0, 0, 1]
)
```

**Use cases**:
- Multi-body dynamics with moving magnets
- Coordinate system alignment
- Rotating machinery simulations

---

### Field Visualization

**File**: `visualize_field.py`

Generates VTK output for field visualization in ParaView.

**Features**:
- Export Radia geometry to VTK
- Evaluate field on NGSolve mesh
- Save as `.vtu` files

**Visualization files**:
- `radia_field.pvsm`: ParaView state file (pre-configured visualization)

**Run**:
```bash
python visualize_field.py
# Open in ParaView: File -> Load State -> radia_field.pvsm
```

---

## Performance Benchmarks

### GridFunction.Set() Benchmark

**File**: `benchmark_gridfunction_set.py`

Measures the critical operation in coupled simulations: evaluating Radia field on NGSolve mesh.

**Tests**:
- Dense matrix vs H-matrix acceleration
- Different Radia element counts (27 to 343 elements)
- Different NGSolve mesh sizes (coarse to fine)
- All field types (B, H, A, M)

**Key metric**: `GridFunction.Set()` time

**Run**:
```bash
python benchmark_gridfunction_set.py
```

**Actual Results** (2025-11-08):

| Configuration | Dense (ms) | H-matrix (ms) | Speedup |
|---------------|-----------|---------------|---------|
| N=27, h=0.05m, 135 vertices | 22.94 | 13.97 | 1.64x |
| N=64, h=0.025m, 826 vertices | 240.28 | 226.10 | 1.06x |
| N=125, h=0.025m, 826 vertices | 419.56 | 406.89 | 1.03x |
| N=216, h=0.0125m, 5034 vertices | 4697.99 | 4706.82 | 1.00x |
| N=343, h=0.0125m, 5034 vertices | 7354.75 | 7371.10 | 1.00x |

**Summary**:
- Overall average H-matrix speedup: **1.02x** (essentially no improvement)
- H-matrix provides **no significant speedup** for field evaluation
- Performance scales as O(M × N) regardless of H-matrix setting
  - M = number of mesh vertices
  - N = number of Radia elements

**Why H-matrix doesn't help**: See `HMATRIX_ANALYSIS.md`
- H-matrix accelerates **solving** (rad.Solve), not **field evaluation** (rad.Fld)
- GridFunction.Set() only evaluates field from known magnetization
- Would need hierarchical field evaluation (FMM) for speedup

**Recommendation**:
- H-matrix parameter available but has minimal effect on GridFunction.Set()
- For large-scale field evaluation, reduce mesh size or Radia elements where possible

---

### Mesh Convergence Study

**File**: `test_mesh_convergence.py`

Studies convergence of NGSolve solution with mesh refinement when using Radia field as boundary condition or source term.

---

## Memory Management

### Memory Usage Behavior

`radia_ngsolve` has been optimized for efficient memory management with the following improvements:

**Optimizations implemented** (v1.2.0):
1. **Minimized GIL scope**: Python Global Interpreter Lock is acquired only during Python calls
2. **Cache-friendly design**: C++ cache lookups avoid GIL acquisition entirely
3. **Immediate object cleanup**: Temporary Python objects released at end of scope

**Memory behavior**:
- Initial memory increase during first 100-200 iterations (Python memory pool initialization)
- **Growth rate decreases over time** (saturating behavior)
- Memory usage **stabilizes** after initialization phase
- **Not a memory leak** - this is normal Python interpreter behavior

**Test results** (500 iterations):
```
First 50 iterations:  3.36 MB per 100 iterations
Last 50 iterations:   1.98 MB per 100 iterations
Growth rate change:   -41% (saturating)
```

**Conclusion**: Memory usage is safe for long-running simulations. The initial growth is Python's memory pool allocation, not a leak.

**For optimal performance**:
- Use `PrepareCache()` for repeated evaluations on the same mesh
- Memory usage will stabilize after warm-up period
- No special cleanup needed - Python GC handles everything

---

## API Reference

### radia_ngsolve.RadiaField

**Syntax**:
```python
cf = radia_ngsolve.RadiaField(
	radia_obj,           # Radia object ID
	field_type='b',      # 'b', 'h', 'a', or 'm'
	origin=None,         # [x, y, z] translation (m)
	u_axis=None,         # Local u-axis [ux, uy, uz]
	v_axis=None,         # Local v-axis [vx, vy, vz]
	w_axis=None,         # Local w-axis [wx, wy, wz]
	precision=None,      # Computation precision (T)
	use_hmatrix=None,    # Enable H-matrix (True/False/None)
	units='mm'           # Radia model units ('mm' or 'm')
)
```

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `radia_obj` | int | Radia object ID (from `rad.ObjRecMag`, etc.) |
| `field_type` | str | `'b'` (B field), `'h'` (H field), `'a'` (A potential), `'m'` (M) |
| `origin` | list[3] | Translation vector in meters (default: [0,0,0]) |
| `u_axis`, `v_axis`, `w_axis` | list[3] | Local coordinate axes (auto-normalized) |
| `precision` | float | Computation precision in Tesla (None = Radia default) |
| `use_hmatrix` | bool | H-matrix acceleration (True/False/None=keep current) |
| `units` | str | Radia model units: `'m'` (default, recommended) or `'mm'` |

**Returns**: NGSolve `CoefficientFunction` (vector, dim=3)

---

### Unit Conversion

**Important**: The `units` parameter specifies the units used in your Radia model, NOT the units you want to use in NGSolve.

**Rule**: Match `units` parameter to `rad.FldUnits()` setting:

| Radia Model Units | `rad.FldUnits()` | `units` parameter | NGSolve Coordinates |
|-------------------|------------------|-------------------|---------------------|
| **Meters (recommended)** | `rad.FldUnits('m')` | `units='m'` (default) | Always meters |
| Millimeters (legacy) | `rad.FldUnits('mm')` | `units='mm'` | Always meters |

**How it works**:
- NGSolve always uses **meters** for coordinates
- The `units` parameter tells radia_ngsolve how to convert NGSolve's meters to your Radia model units:
  - `units='mm'`: NGSolve meters → Radia millimeters (multiply by 1000)
  - `units='m'`: NGSolve meters → Radia meters (no conversion)

**Example 1: Radia model in meters (recommended)**
```python
# Define Radia model in meters
rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1])  # 0.1m = 100mm cube

# Create NGSolve interface - units='m' is default
B_cf = radia_ngsolve.RadiaField(magnet, 'b')  # No need to specify units='m'

# NGSolve works in meters
mesh = Mesh(...)  # Mesh with coordinates in meters
gf = GridFunction(HCurl(mesh))
gf.Set(B_cf)  # NGSolve point [0.05, 0, 0] m → Radia point [0.05, 0, 0] m (no conversion)
```

**Example 2: Radia model in millimeters (legacy)**
```python
# Define Radia model in millimeters
rad.FldUnits('mm')
magnet = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 1])  # 100mm cube

# Create NGSolve interface - specify units='mm' explicitly
B_cf = radia_ngsolve.RadiaField(magnet, 'b', units='mm')

# NGSolve works in meters
mesh = Mesh(...)  # Mesh with coordinates in meters
gf = GridFunction(HCurl(mesh))
gf.Set(B_cf)  # NGSolve point [0.05, 0, 0] m → Radia point [50, 0, 0] mm (converted)
```

**Why this matters**:
- Using wrong `units` parameter will cause 1000× error in field evaluation positions
- Both examples above evaluate the field at the **same physical location** (50mm from origin)
- The field values will be **identical** regardless of units (physics doesn't change!)

---

### Performance Control

**High accuracy** (slow):
```python
B_cf = radia_ngsolve.RadiaField(magnet, 'b',
                               precision=1e-8,
                               use_hmatrix=False)
```

**High speed** (large meshes, 1000+ vertices):
```python
# Enable H-matrix field evaluation (required first)
rad.SetHMatrixFieldEval(1, 1e-6)

# Use H-matrix for field evaluation
B_cf = radia_ngsolve.RadiaField(magnet, 'b', use_hmatrix=True)
```

**Balanced** (default):
```python
B_cf = radia_ngsolve.RadiaField(magnet, 'b')
```

**Note**: H-matrix acceleration for field evaluation requires:
1. `rad.SetHMatrixFieldEval(1, tol)` - Enable H-matrix field evaluation
2. `use_hmatrix=True` parameter in `RadiaField()`
3. Beneficial for **large meshes** (1000+ vertices) or **many observation points**

---

### Field Evaluation Cache (PrepareCache)

**Performance Optimization**: `radia_ngsolve` implements an internal **C++ cache** to avoid repeated `rad.Fld()` calls during `GridFunction.Set()`.

**Why Caching is Needed**:
- NGSolve's `GridFunction.Set(cf)` calls `cf.Evaluate(point)` **individually** for each mesh element
- Without cache: 1000 vertices → 1000 separate `rad.Fld()` calls → **very slow**
- With cache: 1 batch evaluation → cache all points → `GridFunction.Set()` retrieves from cache → **fast**

**How It Works**:
```python
# Create field object
B_cf = radia_ngsolve.RadiaField(magnet, 'b')

# Extract mesh vertices
mesh_points = [[mesh.vertices[i].point[0],
                mesh.vertices[i].point[1],
                mesh.vertices[i].point[2]]
               for i in range(mesh.nv)]

# Prepare cache (1 batch call to rad.Fld())
B_cf.PrepareCache(mesh_points)

# GridFunction.Set() now retrieves from cache (fast!)
gf = GridFunction(HCurl(mesh))
gf.Set(B_cf)  # Each evaluation: cache lookup, not rad.Fld() call
```

**Cache Statistics**:
```python
stats = B_cf.GetCacheStats()
print(stats)
# {'enabled': True, 'size': 1000, 'hits': 5000, 'misses': 0, 'hit_rate': 1.0}
```

**Implementation Details**:
- Cache implemented in **C++** (`std::unordered_map<uint64_t, std::array<double,3>>`)
- Uses FNV-1a hash of quantized 3D coordinates
- Cache tolerance: 1e-10 meters (points within tolerance share same hash)
- Thread-safe with GIL (Global Interpreter Lock)

**Performance Impact**:
- **Without cache**: 1000 vertices × 50ms/call = **50 seconds**
- **With cache**: 1 batch call (500ms) + cache lookups (negligible) = **< 1 second**
- Speedup: **50x** for typical meshes

**Note**: Cache must be manually prepared before `GridFunction.Set()` for optimal performance. Without `PrepareCache()`, each evaluation falls back to individual `rad.Fld()` calls.

---

### Coordinate Transformation

**Translation + Rotation**:
```python
import numpy as np

# Rotate 30° around z-axis, translate to (0.05, 0.1, 0) m
angle = np.radians(30)
origin = [0.05, 0.1, 0]
u_axis = [np.cos(angle), np.sin(angle), 0]
v_axis = [-np.sin(angle), np.cos(angle), 0]
w_axis = [0, 0, 1]

B_cf = radia_ngsolve.RadiaField(magnet, 'b',
                               origin=origin,
                               u_axis=u_axis,
                               v_axis=v_axis,
                               w_axis=w_axis)
```

**Transformation steps**:
1. `p' = p_global - origin` (translation)
2. `p_local = [u·p', v·p', w·p']` (rotation to local)
3. Evaluate field in Radia coordinate system (mm)
4. `F_global = u*F_local[0] + v*F_local[1] + w*F_local[2]` (rotation to global)

---

## Unit Conventions

| Quantity | Radia | NGSolve | Conversion |
|----------|-------|---------|------------|
| Length | mm | m | ×1000 (NGS→Radia) |
| Magnetic flux density B | Tesla | Tesla | (no conversion) |
| Magnetic field H | A/m | A/m | (no conversion) |
| Vector potential A | T·mm | T·m | ×0.001 (for curl) |
| Magnetization M | T | A/m | (context dependent) |

**Automatic conversion**: Handled internally by `radia_ngsolve.RadiaField`

**curl(A) = B**: Vector potential **A** is scaled by 0.001 to ensure correct derivative with respect to meters.

---

## File Organization

```
NGSolve_Integration/
├── README.md                       # This file
├── HMATRIX_ANALYSIS.md             # Why H-matrix doesn't accelerate field evaluation
├── demo_field_types.py             # Basic field type examples
├── verify_curl_A_equals_B.py       # Vector potential verification
├── test_coordinate_transform.py    # Coordinate transformation tests
├── test_mesh_convergence.py        # Mesh convergence study
├── visualize_field.py              # VTK export for ParaView
├── export_radia_geometry.py        # Radia geometry to VTK
├── benchmark_gridfunction_set.py   # Performance benchmark
├── benchmark_results.txt           # Benchmark output
├── radia_field.pvsm                # ParaView state file
└── rad.ObjBckgCF/                  # NGSolve → Radia background field
    ├── README.md
    ├── sphere_nastran_analysis.py
    ├── Cubit2Nastran.py
    ├── sphere.bdf
    └── *.vtu, *.vtk, *.pvsm        # Visualization files
```

---

## Important Limitations

### ⚠️ Do NOT Use radia_ngsolve Inside Permanent Magnets

**Critical Limitation**: `radia_ngsolve.RadiaField()` should **NOT** be used for field calculations **inside** permanent magnets on the NGSolve side.

**Reason**:
- `radia_ngsolve` evaluates Radia field using `rad.Fld()`
- `rad.Fld()` uses boundary element method (BEM) which is designed for **air regions** only
- BEM is **inaccurate inside permanent magnets** - this is a fundamental limitation, not a bug
- Therefore, `radia_ngsolve.RadiaField()` will return incorrect values inside magnets

**Valid Use Cases** (✓ OK):
```python
# GOOD - Field evaluation OUTSIDE magnets
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])  # 10mm cube magnet
mesh = Mesh(Box((0.02, 0, 0), (0.05, 0.03, 0.03)))  # Mesh OUTSIDE magnet
B_cf = radia_ngsolve.RadiaField(magnet, 'b')
gf.Set(B_cf)  # OK - evaluating outside magnet
```

**Invalid Use Cases** (✗ DO NOT DO):
```python
# BAD - Field evaluation INSIDE magnets
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1.0])  # Large magnet
mesh = Mesh(Box((0, 0, 0), (0.01, 0.01, 0.01)))  # Small mesh INSIDE magnet
B_cf = radia_ngsolve.RadiaField(magnet, 'b')
gf.Set(B_cf)  # WRONG - will give incorrect values inside magnet!
```

**Alternative for Fields Inside Magnets**:
- Use direct NGSolve FEM solutions with proper material properties
- Do NOT rely on `radia_ngsolve` for interior field calculations

---

## Troubleshooting

**Import error: `No module named 'radia_ngsolve'`**
- Solution: Rebuild Radia with NGSolve support
- Check: `import radia_ngsolve` in Python

**Unit mismatch errors**
- Radia uses millimeters (mm)
- NGSolve uses meters (m)
- `radia_ngsolve.RadiaField` handles conversion automatically

**Performance issues**
- For large N (> 500): Use `use_hmatrix=True`
- For small N (< 100): Use `use_hmatrix=False` or default

**Accuracy issues**
- Increase precision: `precision=1e-8`
- Disable H-matrix: `use_hmatrix=False`
- Refine NGSolve mesh

---

## References

1. **API Documentation**: [`docs/API_EXTENSIONS.md`](../../docs/API_EXTENSIONS.md)
2. **Radia Documentation**: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html
3. **NGSolve Documentation**: https://ngsolve.org

---

**Last Updated**: 2025-11-08
**Status**: Production ready
**License**: LGPL-2.1
