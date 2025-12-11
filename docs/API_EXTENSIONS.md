# Radia API Extensions

This document describes custom extensions to the original ESRF Radia API.

**Original Documentation**: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html

**Date**: 2025-12-11
**Version**: 1.3.13
**Status**: Active development

> **Note (v1.3.13):** The solver architecture was significantly updated in v1.3.13:
> - **Recommended API**: Use `rad.Solve()` for all new code
> - **Deprecated API**: `RlxPre()`, `RlxMan()`, `SetRelaxSubInterval()` are legacy
> - **Solver Methods**: LU (Method 0) and BiCGSTAB (Method 1) replaced original Gauss-Seidel
> - **H-matrix**: Under research, NOT available in v1.3.13
>
> See [SOLVER_METHODS.md](SOLVER_METHODS.md) for current solver documentation.

---

## Table of Contents

- [Background Field Sources](#background-field-sources)
  - [ObjBckg](#objbckg)
  - [ObjBckgCF](#objbckgcf)
- [Solver Extensions](#solver-extensions)
  - [SetRelaxSubInterval](#setrelaxsubinterval)
  - [RlxPre](#rlxpre)
  - [RlxMan - Method 5 Support](#rlxman---method-5-support)
- [Performance Features](#performance-features)
  - [SolverHMatrixDisable/Enable](#solverhmatrixdisableenable)
- [NGSolve Integration](#ngsolve-integration)
  - [radia_ngsolve.RadiaField](#radia_ngsolveradiafield)

---

## Background Field Sources

### ObjBckg

**Purpose**: Create a source of uniform background magnetic field.

**Syntax**:
```python
field_src = rad.ObjBckg([Bx, By, Bz])
```

**Parameters**:
- `[Bx, By, Bz]`: Uniform magnetic field vector in Tesla

**Returns**:
- Integer key identifying the field source object

**Usage Example**:
```python
import radia as rad

# Create uniform background field of 0.5 T in z-direction
bg_field = rad.ObjBckg([0, 0, 0.5])

# Use in field computation
grp = rad.ObjCnt([magnet, bg_field])
B = rad.Fld(grp, 'b', [0, 0, 0])
```

**Notes**:
- Can be combined with other objects using `rad.ObjCnt()`
- Useful for applying external field to magnetostatic problems
- Field is uniform throughout entire space

---

### ObjBckgCF

**Purpose**: Create a source of arbitrary spatially-varying background field using a Python callback function.

**Syntax**:
```python
field_src = rad.ObjBckgCF(callback_function)
```

**Parameters**:
- `callback_function`: Python function with signature `f([x,y,z]) -> [Bx,By,Bz]`
  - Input: Position `[x, y, z]` in millimeters
  - Output: Field `[Bx, By, Bz]` in Tesla

**Returns**:
- Integer key identifying the field source object

**Usage Example**:
```python
import radia as rad
import math

# Define spatially-varying field (e.g., gradient field)
def gradient_field(pos):
	x, y, z = pos
	gradient = 0.01  # T/mm
	Bx = gradient * x
	By = gradient * y
	Bz = 0.5  # Constant z-component
	return [Bx, By, Bz]

# Create background field source
bg_field = rad.ObjBckgCF(gradient_field)

# Use in computation
grp = rad.ObjCnt([magnet, bg_field])
B = rad.Fld(grp, 'b', [10, 20, 0])  # Will call gradient_field([10, 20, 0])
```

**NGSolve Integration**:

For use with NGSolve's CoefficientFunction:

```python
import ngsolve as ngs
import radia_ngsolve

# Create NGSolve coefficient function (in meters)
def ngsolve_field_wrapper(x, y, z):
	# Convert NGSolve coordinates (m) to Radia coordinates (mm)
	pos_mm = [x*1000, y*1000, z*1000]
	B = gradient_field(pos_mm)
	return B

# Create background field
bg_cf = radia_ngsolve.BackgroundFieldCF(ngsolve_field_wrapper)
```

**Notes**:
- Callback is invoked for every field evaluation point
- Performance depends on callback function efficiency
- Coordinate system: Radia uses millimeters, NGSolve uses meters
- Can implement any analytical or numerical field model

---

## Solver Extensions

### SetRelaxSubInterval

**Purpose**: Configure element grouping for LU decomposition solver (Method 5).

**Syntax**:
```python
rad.SetRelaxSubInterval(intrc, start_idx, end_idx, relax_together=1)
```

**Parameters**:
- `intrc`: Interaction matrix key (from `rad.RlxPre()`)
- `start_idx`: Starting element index (0-based, inclusive)
- `end_idx`: Ending element index (0-based, inclusive)
- `relax_together`: Integer flag
  - `1` = RelaxTogether (use LU decomposition for this group) - **default**
  - `0` = RelaxApart (use Gauss-Seidel for this group)

**Returns**:
- `1` on success

**Usage Example**:

**Basic - Enable LU for all elements**:
```python
import radia as rad

# Create geometry with N elements
elements = []
for i in range(N):
	elem = rad.ObjRecMag([x, y, z], [dx, dy, dz], [0, 0, 0.1])
	rad.MatApl(elem, material)
	elements.append(elem)

grp = rad.ObjCnt(elements)

# Build interaction matrix
intrc = rad.RlxPre(grp, grp)

# Enable LU decomposition for all elements
rad.SetRelaxSubInterval(intrc, 0, N-1, 1)

# Solve using Method 5 (LU decomposition)
rad.RlxMan(intrc, 5, 100, 1.0)
```

**Advanced - Mixed solver strategies**:
```python
# Group 1 (elements 0-49): Use LU decomposition
rad.SetRelaxSubInterval(intrc, 0, 49, 1)

# Group 2 (elements 50-99): Use Gauss-Seidel
rad.SetRelaxSubInterval(intrc, 50, 99, 0)

# Solve with Method 5
# - Group 1 solved with direct LU inversion
# - Group 2 solved with iterative Gauss-Seidel
rad.RlxMan(intrc, 5, 100, 1.0)
```

**Performance Considerations**:
- **LU decomposition (together=1)**: O(N³) per solve
  - Best for: Small groups, direct solve required, few iterations
  - Memory: O(N²) for matrix storage

- **Gauss-Seidel (together=0)**: O(N²) per iteration
  - Best for: Large problems, many iterations, iterative convergence
  - Memory: O(N²) for interaction matrix only

**Notes**:
- Must be called **after** `rad.RlxPre()` and **before** solver methods
- Can be called multiple times to define different groups
- Overlapping intervals: behavior is undefined (avoid)
- See [benchmark results](../examples/solver_benchmarks/README.md) for scaling analysis

---

### RlxPre

**Purpose**: Build interaction matrix for relaxation solver.

**Syntax**:
```python
intrc = rad.RlxPre(obj, srcobj=0)
```

**Parameters**:
- `obj`: Main object key (magnetizable material)
- `srcobj`: Additional external field source object key (optional, default=0)
  - If `srcobj=0` or omitted, only `obj` contributes to field
  - If `srcobj` specified, it acts as additional field source

**Returns**:
- Integer key identifying the interaction matrix

**Usage Example**:

**Basic - Self-consistent problem**:
```python
# Single object self-field
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0.1])
rad.MatApl(magnet, material)

intrc = rad.RlxPre(magnet, magnet)  # Self-interaction
```

**With external field source**:
```python
# Magnet in external coil field
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,0.1])
rad.MatApl(magnet, material)

coil = rad.ObjRaceTrk([0,0,20], [30,20], [5,5], 1000, 3)  # Current coil

# Build interaction matrix with external source
intrc = rad.RlxPre(magnet, coil)
```

**Notes**:
- Matrix construction time: O(N²) where N = number of elements
- Matrix is reusable for multiple `RlxMan()` or `RlxAuto()` calls
- Use `rad.RlxUpdSrc()` to update external source without rebuilding matrix

**Related Functions**:
- `rad.SetRelaxSubInterval()` - Configure LU decomposition groups
- `rad.RlxMan()` - Manual relaxation with interaction matrix
- `rad.RlxAuto()` - Automatic relaxation with interaction matrix
- `rad.RlxUpdSrc()` - Update external source data

---

### RlxMan - Method 5 Support

**Purpose**: Execute manual relaxation procedure with extended solver methods.

**Original ESRF API**: Methods 0-4
**Extended API**: Methods 0-5 (adds LU decomposition)

**Syntax**:
```python
rad.RlxMan(intrc, method, iter_num, relax_param)
```

**Parameters**:
- `intrc`: Interaction matrix key (from `rad.RlxPre()`)
- `method`: Solver method number (0-5)
- `iter_num`: Number of iterations to perform
- `relax_param`: Relaxation parameter value (typically 1.0)

**Solver Methods**:

| Method | Name | Scaling | Description |
|--------|------|---------|-------------|
| 0 | Simple iteration | O(N²) | Basic fixed-point iteration |
| 1 | Over-relaxation | O(N²) | Accelerated fixed-point |
| 2 | Under-relaxation | O(N²) | Damped iteration |
| 3 | Gauss-Seidel | O(N²) | Sequential element updates |
| 4 | Gauss-Seidel | O(N²) | (Same as Method 3) |
| **5** | **LU decomposition** | **O(N³)** | **Direct matrix inversion** ⭐ NEW |

**Method 5 - LU Decomposition (Extended)**:

**When to use**:
- Direct solve required (no iteration uncertainty)
- Small problem size (N < 200)
- Few solver calls (matrix decomposition cost amortized)

**Requirements**:
- Must call `rad.SetRelaxSubInterval(intrc, start, end, 1)` first
- Otherwise behaves like Method 3 (Gauss-Seidel)

**Usage Example**:

**Method 5 - LU decomposition**:
```python
import radia as rad
from time import perf_counter

# Build geometry
elements = []
for i in range(N):
	elem = rad.ObjRecMag([x, y, z], [dx, dy, dz], [0, 0, 0.1])
	rad.MatApl(elem, material)
	elements.append(elem)

grp = rad.ObjCnt(elements)

# Measure matrix construction time (O(N²))
t0 = perf_counter()
intrc = rad.RlxPre(grp, grp)
t_matrix = perf_counter() - t0
print(f"Matrix construction: {t_matrix*1000:.2f} ms")

# Enable LU decomposition for all elements
rad.SetRelaxSubInterval(intrc, 0, N-1, 1)

# Measure solver time (O(N³) for LU)
t0 = perf_counter()
rad.RlxMan(intrc, 5, 100, 1.0)  # Method 5
t_solver = perf_counter() - t0
print(f"LU solver: {t_solver*1000:.2f} ms")
```

**Method 4 - Gauss-Seidel (comparison)**:
```python
# Same matrix construction
intrc = rad.RlxPre(grp, grp)

# Gauss-Seidel solver (O(N²) per iteration)
t0 = perf_counter()
rad.RlxMan(intrc, 4, 100, 1.0)  # Method 4
t_solver = perf_counter() - t0
print(f"Gauss-Seidel: {t_solver*1000:.2f} ms")
```

**Performance Comparison**:

See detailed benchmark results in [`examples/solver_benchmarks/README.md`](../examples/solver_benchmarks/README.md).

**Summary**:
- **Small N (< 64)**: LU and GS comparable
- **Medium N (64-343)**: LU becomes 100-1000× slower per iteration
- **Large N (> 343)**: LU becomes 1000-10000× slower per iteration
- **Recommendation**: Use Method 4 (Gauss-Seidel) for most problems

**Notes**:
- Method 5 without `SetRelaxSubInterval()` falls back to Method 3
- LU decomposition is performed once, then reused for all iterations
- For nonlinear materials, LU decomposition is re-computed each iteration
- Matrix construction time O(N²) is identical for all methods

**Related Documentation**:
- Original ESRF Methods 0-4: https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html
- Benchmark analysis: [`examples/solver_benchmarks/README.md`](../examples/solver_benchmarks/README.md)

---

## Performance Features

### SolverHMatrixDisable/Enable (DEPRECATED)

> **NOT AVAILABLE in v1.3.13**
>
> H-matrix acceleration was evaluated and found to provide **NO benefit for typical Radia use cases**
> (single compact objects). These APIs are not available in the current release.
>
> See [HMATRIX_EVALUATION.md](HMATRIX_EVALUATION.md) for evaluation results.

~~**Purpose**: Control H-matrix (hierarchical matrix) acceleration for field computation.~~

~~**Syntax**:~~
```python
# NOT AVAILABLE in v1.3.13
# rad.SolverHMatrixDisable()
# rad.SolverHMatrixEnable()
```

---

### SolverTetraMethod

**Purpose**: Set the field computation method for tetrahedral mesh elements.

**Syntax**:
```python
rad.SolverTetraMethod(method)
```

**Parameters**:
- `method`: Integer selecting the computation method
  - `0`: Original Radia polygon-based method (default)
  - `1`: Analytical method using magnetic charge formulation

**Returns**: None

**Background**:

Tetrahedral elements in Radia can compute their magnetic field contribution using two different algorithms:

| Method | Description | Best For |
|--------|-------------|----------|
| 0 | Original Radia polygon method | General purpose, standard accuracy |
| 1 | Analytical method | High-permeability materials (mu_r > 100) |

**Method 0 (Original)**:
- Uses Radia's standard polygon-based field calculation
- Treats each tetrahedral face as a magnetic charge sheet
- Well-tested and reliable for most applications

**Method 1 (Analytical)**:
- Uses analytical formulas for magnetic field from surface charges
- May provide better accuracy for high-permeability soft magnetic materials
- Mathematically equivalent to Method 0 but with different numerical implementation

**Usage Example**:

```python
import radia as rad
from netgen_mesh_import import netgen_mesh_to_radia

# Use analytical method for high-permeability material
rad.SolverTetraMethod(1)

# Import tetrahedral mesh with soft iron (mu_r = 4000)
mesh_data = netgen_mesh_to_radia(
    'soft_iron.vol',
    material={'magnetization': [0, 0, 0]},  # Zero initial M
    units='m'
)

# Apply linear material
mat = rad.MatLin(3999)  # ksi = mu_r - 1
rad.MatApl(mesh_data, mat)

# Solve
rad.Solve(mesh_data, 0.0001, 10000)

# Compute field
B = rad.Fld(mesh_data, 'b', [0.1, 0.1, 0.1])
```

**Migration from Environment Variable**:

Previously, the tetrahedral method was controlled via the `RADIA_TETRA_METHOD` environment variable. This has been replaced by the `SolverTetraMethod()` API for better discoverability and control.

```python
# Old way (deprecated):
import os
os.environ['RADIA_TETRA_METHOD'] = '1'  # No longer works

# New way (recommended):
import radia as rad
rad.SolverTetraMethod(1)
```

**Important Notes**:
- Setting is global and persists for the session
- Must be called **before** `rad.Solve()` or `rad.RlxPre()`
- Both methods should produce identical results (within numerical precision)
- Method 0 and Method 1 have been verified to produce 0.00% difference

---

## Version History

### v1.0.8 (2025-11-27)
- Added `SolverTetraMethod()` API for tetrahedral mesh computation method control
- Fixed ANALYTICAL method (Method 1) to match standard polygon method
- Deprecated `RADIA_TETRA_METHOD` environment variable

### v1.0.7 (2025-11-08)
- Added `SetRelaxSubInterval()` for LU decomposition control
- Extended `RlxMan()` to support Method 5 (LU decomposition)
- Updated documentation with benchmark results

### v1.0.6 (2025-11-02)
- Fixed vector potential coordinate conversion
- Added `ObjBckgCF()` for arbitrary background fields
- Improved NGSolve integration

### v1.3.13 (2025-12-11)
- **Major solver refactoring**: Replaced Implicit SS (Gauss-Seidel) with BiCGSTAB
- Added `rad.Solve()` as recommended API
- Deprecated `RlxPre()`, `RlxMan()`, `SetRelaxSubInterval()` (still available for legacy)
- H-matrix APIs removed (under research)
- Added OpenBLAS and OpenMP optimizations

### v1.0.5 (2025-10-30)
- Added `SolverHMatrixDisable/Enable()` controls (REMOVED in v1.3.13)
- Performance improvements for H-matrix (REMOVED in v1.3.13)

---

## NGSolve Integration

### radia_ngsolve.RadiaField

**Purpose**: Create NGSolve CoefficientFunction for Radia magnetic field with full control over computation accuracy and performance.

**Module**: `radia_ngsolve` (C++ extension, included in `radia` package)

> **CRITICAL: Import Order**
>
> You **MUST** import modules in this exact order:
> 1. `import radia` - Load radia first
> 2. `import ngsolve` - Load NGSolve **BEFORE** radia_ngsolve
> 3. `from radia import radia_ngsolve` - Load radia_ngsolve last
>
> Incorrect order causes `ImportError: DLL load failed`.

**Syntax**:
```python
# Correct import order
import radia as rad
import ngsolve
from ngsolve import *
from radia import radia_ngsolve

cf = radia_ngsolve.RadiaField(
    radia_obj,
    field_type='b',
    origin=None,
    u_axis=None,
    v_axis=None,
    w_axis=None,
    precision=None,
    use_hmatrix=None
)
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `radia_obj` | int | *required* | Radia object ID |
| `field_type` | str | `'b'` | Field type: `'b'` (flux density), `'h'` (magnetic field), `'a'` (vector potential), `'m'` (magnetization) |
| `origin` | list[3] | `None` | Translation vector [x, y, z] in meters |
| `u_axis` | list[3] | `None` | Local u-axis (auto-normalized) |
| `v_axis` | list[3] | `None` | Local v-axis (auto-normalized) |
| `w_axis` | list[3] | `None` | Local w-axis (auto-normalized) |
| `precision` | float | `None` | Computation precision in Tesla (None = Radia default) |
| `use_hmatrix` | bool | `None` | Enable H-matrix acceleration (None = keep current setting) |

**Returns**:
- NGSolve `CoefficientFunction` (vector-valued, dim=3)

---

#### Field Types

| Type | Description | Units | Notes |
|------|-------------|-------|-------|
| `'b'` | Magnetic flux density **B** | Tesla (T) | Default, most common |
| `'h'` | Magnetic field **H** | A/m | For material analysis |
| `'a'` | Vector potential **A** | T·m | Use `curl(A)` to get **B** |
| `'m'` | Magnetization **M** | A/m | Inside magnetic materials |

---

#### Coordinate Transformation

If `origin`, `u_axis`, `v_axis`, or `w_axis` are specified, coordinate transformation is applied:

1. **Translation**: `p' = p_global - origin`
2. **Rotation to local**: `p_local = [u·p', v·p', w·p']`
3. **Field evaluation**: Radia computes field in local coordinates (mm)
4. **Rotation to global**: `F_global = u*F_local[0] + v*F_local[1] + w*F_local[2]`

**Use cases**:
- Moving/rotating magnets in FEM simulation
- Coordinate system alignment between NGSolve and Radia
- Multi-body dynamics coupling

---

#### Performance Control

**precision** parameter:

Controls accuracy vs speed trade-off. Smaller values = more accurate but slower.

```python
# High accuracy (slow)
B_cf = radia_ngsolve.RadiaField(magnet, 'b', precision=1e-8)

# Normal accuracy (default)
B_cf = radia_ngsolve.RadiaField(magnet, 'b')  # Uses Radia default

# Fast evaluation (less accurate)
B_cf = radia_ngsolve.RadiaField(magnet, 'b', precision=1e-4)
```

**Internally calls**: `rad.FldCmpPrc()` with specified precision for B, H, A, M fields.

---

**use_hmatrix** parameter:

Controls H-matrix (hierarchical matrix) acceleration.

```python
# Disable H-matrix (accurate for small N)
B_cf = radia_ngsolve.RadiaField(magnet, 'b', use_hmatrix=False)

# Enable H-matrix (fast for large N > 1000)
B_cf = radia_ngsolve.RadiaField(magnet, 'b', use_hmatrix=True)

# Keep current setting (default)
B_cf = radia_ngsolve.RadiaField(magnet, 'b')  # use_hmatrix=None
```

**When to enable H-matrix**:
- ✅ Large number of elements (N > 1000)
- ✅ Many field evaluation points
- ✅ Memory constraints
- ❌ Small problems (N < 100, overhead not worth it)
- ❌ Highest accuracy required

**Internally calls**: `rad.SolverHMatrixEnable()` or `rad.SolverHMatrixDisable()`.

---

#### Usage Examples

**Basic usage**:
```python
import radia as rad
from ngsolve import *
from netgen.occ import *
import radia_ngsolve

# Create Radia magnet
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1.2])

# Create NGSolve CoefficientFunction
B_cf = radia_ngsolve.RadiaField(magnet, 'b')

# Use in NGSolve mesh
mesh = Mesh(...)
gf = GridFunction(HCurl(mesh))
gf.Set(B_cf)
```

**High accuracy with coordinate transformation**:
```python
# Magnet rotated 45 degrees around z-axis, translated to (0.1, 0, 0) m
import numpy as np

angle = np.pi/4
origin = [0.1, 0, 0]  # meters
u_axis = [np.cos(angle), np.sin(angle), 0]
v_axis = [-np.sin(angle), np.cos(angle), 0]
w_axis = [0, 0, 1]

B_cf = radia_ngsolve.RadiaField(
    magnet, 'b',
    origin=origin,
    u_axis=u_axis,
    v_axis=v_axis,
    w_axis=w_axis,
    precision=1e-8,
    use_hmatrix=False
)
```

**Fast evaluation for large system**:
```python
# Create large magnet array
elements = []
for i in range(20):
    for j in range(20):
        for k in range(20):
            elem = rad.ObjRecMag([i*10, j*10, k*10], [8,8,8], [0,0,1])
            elements.append(elem)

magnet_array = rad.ObjCnt(elements)  # 8000 elements

# Use H-matrix for fast evaluation
B_cf = radia_ngsolve.RadiaField(magnet_array, 'b', use_hmatrix=True)

# Evaluate on large mesh
mesh = Mesh(...)  # Large mesh
gf = GridFunction(HCurl(mesh))
gf.Set(B_cf)  # Fast thanks to H-matrix
```

**Vector potential with curl**:
```python
# Get vector potential A
A_cf = radia_ngsolve.RadiaField(magnet, 'a')

# Compute B = curl(A)
B_from_curl = Curl(A_cf)

# Compare with direct B
B_direct = radia_ngsolve.RadiaField(magnet, 'b')

# Should match (within numerical precision)
```

---

#### Performance Comparison

See [`examples/ngsolve_integration/README.md`](../examples/ngsolve_integration/README.md) for detailed benchmarks.

**Typical speedup with H-matrix** (N > 500):
- Field evaluation: 2-10× faster
- Negligible accuracy loss (relative error < 10⁻¹⁰)

---

#### Notes

- **Unit conversion**: NGSolve uses meters, Radia uses millimeters → automatic conversion (×1000)
- **Vector potential scaling**: `A` is scaled by 0.001 so that `curl(A) = B` with correct units
- **Thread safety**: Settings (`precision`, `use_hmatrix`) are applied at object creation
- **Global side effects**: `precision` and `use_hmatrix` modify global Radia settings

---

#### Related Functions

- `rad.Fld()` - Direct Radia field computation
- `rad.FldCmpPrc()` - Set global computation precision
- `rad.SolverHMatrixEnable()`/`Disable()` - Control H-matrix globally

---

#### Version History

- **v0.08 (2025-11-08)**: Added `precision` and `use_hmatrix` parameters
- **v0.07 (2025-11-02)**: Added coordinate transformation support
- **v0.06 (2025-10-30)**: Initial release with basic field types

---

## References

1. **Original ESRF Radia Documentation**
   https://www.esrf.fr/home/Accelerators/instrumentation--equipment/Software/Radia/Documentation/ReferenceGuide.html

2. **Benchmark Results**
   [`examples/solver_benchmarks/README.md`](../examples/solver_benchmarks/README.md)

3. **NGSolve Integration Guide**
   [`docs/NGSOLVE_USAGE_GUIDE.md`](NGSOLVE_USAGE_GUIDE.md)

4. **Source Code**
   - Python bindings: `src/python/radpy.cpp`
   - Solver implementation: `src/core/radapl2.cpp`
   - Interaction matrix: `src/core/radintrc.cpp`

---

**Last Updated**: 2025-12-11
**Maintained By**: Radia Development Team
**License**: LGPL-2.1 (modifications), BSD-style (original RADIA from ESRF)
