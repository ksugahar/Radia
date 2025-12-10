# NGSolve Integration for Radia

## Overview

This document describes the NGSolve integration module (`radia_ngsolve`) that provides CoefficientFunction wrappers for Radia magnetic field calculations.

## Implementation Summary

### Files Created/Modified

#### New Files
1. **`src/radia/radia_ngsolve.cpp`**
   - NGSolve CoefficientFunction wrappers for Radia
   - Implements `RadiaField` unified interface (replaces legacy RadBfield/RadHfield/RadAfield)
   - Supports coordinate transformation (origin, u/v/w axes)
   - Uses pybind11 for Python bindings

2. **`examples/ngsolve_integration/test_radia_ngsolve.py`**
   - Comprehensive test script
   - Demonstrates all features of the radia_ngsolve module

3. **`examples/ngsolve_integration/README.md`**
   - User documentation
   - API reference
   - Usage examples

#### Modified Files
1. **`src/lib/radpy_pyapi.cpp`**
   - Python API bindings for Radia
   - Updated to work with NGSolve integration

2. **`CMakeLists.txt`**
   - Added optional NGSolve module build target
   - Detects NGSolve installation automatically

## Architecture

### CoefficientFunction Classes

The implementation follows the NGSolve CoefficientFunction pattern:

```cpp
namespace ngfem {
	class RadiaBFieldCF : public CoefficientFunction {
	    // 3D vector field (dimension = 3)
	    virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                         FlatVector<> result) const override;
	};
}
```

### Key Design Decisions

1. **Template-free implementation**: Unlike EMPY_Field which uses templates, we use direct class implementations since Radia objects are referenced by integer indices.

2. **Three separate classes**:
   - `RadiaBFieldCF` - B-field (magnetic flux density)
   - `RadiaHFieldCF` - H-field (magnetic field intensity)
   - `RadiaAFieldCF` - A-field (vector potential)

3. **Direct RadFld calls**: Each evaluation directly calls Radia's `RadFld` function with appropriate field type.

4. **Error handling**: Returns zero field on error to avoid exceptions during FEM integration.

## Function Signatures

### RadiaField (Unified Interface - Recommended)

```python
RadiaField(radia_obj: int, field_type: str = 'b',
           origin: list = None, u_axis: list = None,
           v_axis: list = None, w_axis: list = None,
           precision: float = None, units: str = 'm') -> CoefficientFunction
```

Creates a field coefficient function from a Radia object with optional coordinate transformation.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `radia_obj` | int | (required) | Radia object index |
| `field_type` | str | `'b'` | Field type: `'b'` (flux density), `'h'` (H-field), `'a'` (vector potential), `'m'` (magnetization) |
| `origin` | list[3] | `None` | Origin of local coordinate system `[x, y, z]` |
| `u_axis` | list[3] | `None` | U-axis direction vector (will be normalized) |
| `v_axis` | list[3] | `None` | V-axis direction vector (will be normalized) |
| `w_axis` | list[3] | `None` | W-axis direction vector (will be normalized) |
| `precision` | float | `None` | Radia field calculation precision |
| `units` | str | `'m'` | Length unit: `'m'` (meters) or `'mm'` (millimeters) |

**Returns:**
- NGSolve CoefficientFunction representing 3D vector field

**Coordinate Transformation:**

When `origin` and axis vectors are specified, the field is evaluated in a local coordinate system:

1. **Translation**: Point is translated by subtracting `origin`
2. **Rotation**: Translated point is projected onto local axes (u, v, w)
3. **Field Transform**: Computed field is transformed back to global coordinates

Example:
```python
# Evaluate field in a rotated coordinate system
# Origin at (1, 0, 0), local z-axis aligned with global x-axis
B_local = radia_ngsolve.RadiaField(magnet, 'b',
                                    origin=[1.0, 0.0, 0.0],
                                    u_axis=[0, 1, 0],  # local x -> global y
                                    v_axis=[0, 0, 1],  # local y -> global z
                                    w_axis=[1, 0, 0])  # local z -> global x
```

### Cache Methods

RadiaField supports caching for improved performance in repeated evaluations:

```python
# Prepare cache with evaluation points
cf = radia_ngsolve.RadiaField(magnet, 'b')
cf.PrepareCache(points)  # points: list of [x, y, z] coordinates

# Get cache statistics
stats = cf.GetCacheStats()  # Returns dict with 'size', 'hits', 'misses'

# Clear cache
cf.ClearCache()
```

### Legacy Functions (Deprecated)

The following functions are deprecated. Use `RadiaField` instead:

### RadBfield

```python
RadBfield(radia_obj: int, field_comp: str = 'b') -> CoefficientFunction
```

Creates a B-field coefficient function from a Radia object.

**Parameters:**
- `radia_obj`: Radia object index (from rad.ObjRecMag, rad.ObjCnt, etc.)
- `field_comp`: Field component identifier ('b', 'h', or 'a')

**Returns:**
- NGSolve CoefficientFunction representing 3D vector field

### RadHfield

```python
RadHfield(radia_obj: int) -> CoefficientFunction
```

Creates an H-field coefficient function from a Radia object.

### RadAfield

```python
RadAfield(radia_obj: int) -> CoefficientFunction
```

Creates a vector potential coefficient function from a Radia object.

## Build System Integration

### CMake Configuration

The NGSolve module is optional and built only if NGSolve is detected:

```cmake
find_package(NGSolve CONFIG QUIET)

if(NGSolve_FOUND)
	add_library(radia_ngsolve MODULE
	    ${RADIA_LIB_SOURCES}
	    ${RADIA_NGSOLVE_SOURCES}
	)
	target_link_libraries(radia_ngsolve PRIVATE
	    Python3::Python
	    NGSolve::ngsolve
	)
endif()
```

### Build Requirements

- NGSolve installed in Python 3.12 environment
- pybind11 (included with NGSolve)
- Radia library sources

## Usage Examples

### Basic Field Evaluation

```python
import radia as rad
from ngsolve import *
import radia_ngsolve

# Create Radia geometry
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1000])
rad.Solve(magnet, 0.0001, 10000)

# Create NGSolve mesh
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))

# Create coefficient function
B = radia_ngsolve.RadBfield(magnet)

# Use in NGSolve
B_integral = Integrate(B, mesh)
Draw(B, mesh, "B_field")
```

### FEM Assembly

```python
# Define finite element space
fes = HCurl(mesh, order=2)
u, v = fes.TnT()

# Bilinear form
a = BilinearForm(fes)
a += curl(u) * curl(v) * dx

# Linear form with Radia field as source
f = LinearForm(fes)
f += B * v * dx

# Solve
a.Assemble()
f.Assemble()
gfu = GridFunction(fes)
gfu.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec
```

## Important Limitations

### ⚠️ Do NOT Use Inside Permanent Magnets

**Critical Limitation**: `radia_ngsolve` should **NOT** be used for field calculations **inside** permanent magnets.

**Reason**:
- `radia_ngsolve` evaluates Radia field using `rad.Fld()`
- `rad.Fld()` uses Magnetic Moment Method (MMM) which is designed for **air regions** only
- MMM is **inaccurate inside permanent magnets** - this is a fundamental limitation of Radia, not a bug
- Therefore, using `radia_ngsolve.RadiaField()` inside magnets will return incorrect values

**Valid Use Cases** (✓ OK):
- Field evaluation **outside** permanent magnets (air regions)
- Field evaluation in non-magnetic materials
- Background field sources for FEM problems in air regions

**Invalid Use Cases** (✗ DO NOT DO):
- Field evaluation **inside** permanent magnets
- NGSolve mesh overlapping with Radia magnet geometry

**Alternative for Interior Fields**:
- Use direct NGSolve FEM solutions with proper material properties
- Do NOT rely on `radia_ngsolve` for calculations inside magnetic materials

---

## Coordinate Systems and Units

### Unit System Configuration

**CRITICAL**: NGSolve typically uses meters (m), but Radia defaults to millimeters (mm).

**Recommended approach**:
```python
import radia as rad
rad.FldUnits('m')  # Set to meters for NGSolve integration
```

### Radia Default Units
- **Length**: millimeters (mm) - can be changed with `rad.FldUnits('m')`
- **Magnetization**: kA/m
- **B-field**: Tesla (T)
- **H-field**: A/m

### NGSolve
- **Length**: typically meters (m)
- **Fields**: matches Radia units

**Important**: Always call `rad.FldUnits('m')` at the start of your script when using NGSolve integration.

## Performance Considerations

### Field Evaluation Cost

Each integration point evaluation calls `RadFld`:
- For a mesh with N elements and Q quadrature points per element
- Total field evaluations: N × Q
- Consider Radia precision settings to balance accuracy/speed

### Optimization Tips

1. **Mesh refinement**: Use adaptive mesh refinement near field sources
2. **Precision control**: Set appropriate Radia precision with `rad.FldCmpCrt`
3. **Presolving**: Always call `rad.Solve` before creating coefficient functions
4. **Caching**: For static problems, Radia internally caches field calculations

## Testing

### Test Script

Run the comprehensive test:

```bash
python examples/ngsolve_integration/test_radia_ngsolve.py
```

### Test Coverage

The test script verifies:
- ✓ Module import
- ✓ CoefficientFunction creation
- ✓ Field evaluation at integration points
- ✓ Integration over mesh
- ✓ Field arithmetic operations
- ✓ Visualization (if available)

## Comparison with EMPY_Field

| Feature | EMPY_Field | radia_ngsolve |
|---------|-----------|-------------|
| Physics engine | Custom analytical | Radia (MMM) |
| Field sources | Coils, magnets | Any Radia object |
| Template usage | Yes | No |
| Object management | C++ objects | Integer indices |
| Field types | B, H, A, Ω | B, H, A |

## Implemented Features (v1.3+)

Features implemented in the current version:

1. **Field caching**: Cache field values at integration points (`PrepareCache`, `ClearCache`, `GetCacheStats`)
2. **Batch evaluation**: Evaluate multiple points efficiently
3. **Complex geometries**: Support for Radia groups and transformations via coordinate system parameters
4. **Unit conversion**: Automatic mm<->m conversion via `units` parameter
5. **Coordinate transformation**: Local coordinate system support (origin + u/v/w axes)

## Future Enhancements

Potential improvements:

1. **Field derivatives**: Implement gradient computation
2. **Automatic precision**: Adaptive precision based on mesh size

## References

- [NGSolve CoefficientFunction Documentation](https://docu.ngsolve.org/latest/i-tutorials/unit-2.1-coefficient/coefficientfunction.html)
- [Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
- [NGSolve Integration Examples](../examples/ngsolve_integration/)

## Authors

- Implementation: Claude Code (October 2025)
- Based on patterns from EMPY_Field NGSolve integration
- Radia core: O. Chubar, P. Elleaume, et al.

## License

Part of the Radia project. See main LICENSE file for details.
