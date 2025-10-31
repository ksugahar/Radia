# NGSolve Integration for Radia

## Overview

This document describes the NGSolve integration module (`rad_ngsolve`) that provides CoefficientFunction wrappers for Radia magnetic field calculations.

## Implementation Summary

### Files Created/Modified

#### New Files
1. **`src/python/rad_ngsolve.cpp`**
   - NGSolve CoefficientFunction wrappers for Radia
   - Implements `RadBfield`, `RadHfield`, and `RadAfield` classes
   - Uses pybind11 for Python bindings

2. **`examples/ngsolve_integration/test_rad_ngsolve.py`**
   - Comprehensive test script
   - Demonstrates all features of the rad_ngsolve module

3. **`examples/ngsolve_integration/README.md`**
   - User documentation
   - API reference
   - Usage examples

#### Modified Files
1. **`src/python/radpy.cpp` → `src/python/rad_py.cpp`**
   - Renamed for consistency with naming convention
   - Updated file header comment

2. **`CMakeLists.txt`**
   - Updated reference from `radpy.cpp` to `rad_py.cpp`
   - Added optional NGSolve module build target
   - Detects NGSolve installation automatically

3. **Documentation files**
   - Updated references to `radpy.cpp` in all documentation
   - Files: CHANGELOG.md, CODE_QUALITY_REPORT.md, docs/DIRECTORY_STRUCTURE.md, etc.

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
    add_library(rad_ngsolve MODULE
        ${RADIA_LIB_SOURCES}
        ${RADIA_NGSOLVE_SOURCES}
    )
    target_link_libraries(rad_ngsolve PRIVATE
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
import rad_ngsolve

# Create Radia geometry
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1000])
rad.Solve(magnet, 0.0001, 10000)

# Create NGSolve mesh
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))

# Create coefficient function
B = rad_ngsolve.RadBfield(magnet)

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

## Coordinate Systems and Units

### Radia
- **Length**: millimeters (mm)
- **Magnetization**: kA/m
- **B-field**: Tesla (T)
- **H-field**: A/m

### NGSolve
- **Length**: user-defined (typically meters or millimeters)
- **Fields**: matches Radia units

**Important**: Ensure coordinate system consistency between Radia geometry and NGSolve mesh.

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
python examples/ngsolve_integration/test_rad_ngsolve.py
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

| Feature | EMPY_Field | rad_ngsolve |
|---------|-----------|-------------|
| Physics engine | Custom analytical | Radia (BEM) |
| Field sources | Coils, magnets | Any Radia object |
| Template usage | Yes | No |
| Object management | C++ objects | Integer indices |
| Field types | B, H, A, Ω | B, H, A |

## Future Enhancements

Potential improvements:

1. **Field caching**: Cache field values at integration points
2. **Batch evaluation**: Evaluate multiple points in single Radia call
3. **Field derivatives**: Implement gradient computation
4. **Complex geometries**: Support for Radia groups and transformations
5. **Unit conversion**: Automatic mm↔m conversion

## References

- [NGSolve CoefficientFunction Documentation](https://docu.ngsolve.org/latest/i-tutorials/unit-2.1-coefficient/coefficientfunction.html)
- [Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
- [EMPY_Field Implementation](../02_EMPY_Field/)

## Authors

- Implementation: Claude Code (October 2025)
- Based on patterns from EMPY_Field NGSolve integration
- Radia core: O. Chubar, P. Elleaume, et al.

## License

Part of the Radia project. See main LICENSE file for details.
