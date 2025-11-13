# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and refactoring policies for the Radia project when working with Claude Code.

## Source File Organization

### Refactoring Policy: C++ Source File Names Should Reflect API Functionality

**Goal**: Improve code maintainability by organizing C++ source files to clearly correspond to Python API functional groups, making it easier to locate implementation code for specific APIs.

**Rationale**:
- Makes it easy to find which file contains the implementation for a specific API function
- Reduces time spent searching through cryptically-named files (radapl1.cpp, radapl2.cpp, etc.)
- Improves onboarding experience for new developers
- Facilitates targeted code maintenance and debugging

**Important**:
- **Python API names remain unchanged** (backward compatibility maintained)
- **C++ function names remain unchanged** (no breaking changes to internal interfaces)
- **Only C++ source file names are refactored** (organizational improvement only)

### Current File Structure Issues

**Current files** (cryptic naming):
- `radapl1.cpp` - What does "apl1" mean? Contains object creation functions
- `radapl2.cpp` - What does "apl2" mean? Contains material functions
- `radapl3.cpp` - What does "apl3" mean? Contains transformation functions
- `radapl4.cpp` - What does "apl4" mean? Contains field computation functions
- `radapl7.cpp` - What does "apl7" mean? Contains multi-generation extrusion

**Problem**: Developer must memorize or search to know which file contains which API implementation.

### Target File Organization

**Proposed file structure** (functional naming):

| Current File | New File Name | API Category | Contains |
|--------------|---------------|--------------|----------|
| `radapl1.cpp` | `radobj_creation.cpp` | Object Creation | `ObjRecMag`, `ObjCylMag`, `ObjPolyhdr`, `ObjThckPgn`, etc. |
| `radapl2.cpp` | `radmaterial.cpp` | Materials | `MatLin`, `MatStd`, `MatSatIsoFrm`, `MatApl`, etc. |
| `radapl3.cpp` | `radtransform.cpp` | Transformations | `TrfTrsl`, `TrfRot`, `TrfPlSym`, `TrfMlt`, etc. |
| `radapl4.cpp` | `radfield_compute.cpp` | Field Computation | `Fld`, `FldLst`, `FldInt`, `FldEnr`, `FldFrc`, etc. |
| `radapl7.cpp` | `radobj_subdivision.cpp` | Subdivision | `ObjDivMag`, `ObjCutMag`, multi-generation functions |
| `radentry.cpp` | ✓ Keep (clear) | API Entry Points | C API entry points (RadXxx functions) |
| `radinter.cpp` | ✓ Keep (clear) | Interface Layer | Wrapper functions between C API and C++ implementation |

### File Naming Convention

**Pattern**: `rad<functional_area>.cpp`

**Functional areas**:
- `obj_creation` - Object creation APIs (`ObjRecMag`, `ObjCylMag`, etc.)
- `obj_subdivision` - Object subdivision APIs (`ObjDivMag`, `ObjCutMag`, etc.)
- `material` - Material definition and application (`MatLin`, `MatApl`, etc.)
- `transform` - Geometric transformations (`TrfTrsl`, `TrfRot`, etc.)
- `field_compute` - Field computation (`Fld`, `FldLst`, `FldEnr`, etc.)
- `solver` - Relaxation and solving (`RlxPre`, `RlxAuto`, `Solve`, etc.)
- `container` - Container operations (`ObjCnt`, `ObjAddToCnt`, etc.)
- `query` - Object property queries (`ObjM`, `ObjGeoVol`, etc.)

### Implementation Strategy

When refactoring file names:

1. **Create mapping table**: Document current file → new file mapping
2. **Update CMakeLists.txt**: Change file references in build system
3. **Update #include statements**: Change includes to use new file names
4. **Verify build**: Ensure project compiles after renaming
5. **Update documentation**: Reflect new file organization in comments
6. **Test thoroughly**: Run full test suite to ensure no regressions

### API to File Mapping Examples

**Object Creation APIs** → `radobj_creation.cpp`:
- `rad.ObjRecMag()` → Implementation in `radobj_creation.cpp`
- `rad.ObjCylMag()` → Implementation in `radobj_creation.cpp`
- `rad.ObjPolyhdr()` → Implementation in `radobj_creation.cpp`

**Material APIs** → `radmaterial.cpp`:
- `rad.MatLin()` → Implementation in `radmaterial.cpp`
- `rad.MatStd()` → Implementation in `radmaterial.cpp`
- `rad.MatApl()` → Implementation in `radmaterial.cpp`

**Field Computation APIs** → `radfield_compute.cpp`:
- `rad.Fld()` → Implementation in `radfield_compute.cpp`
- `rad.FldLst()` → Implementation in `radfield_compute.cpp`
- `rad.FldEnr()` → Implementation in `radfield_compute.cpp`

### Benefits

1. **Faster development**: "Where is `ObjRecMag` implemented?" → "Look in `radobj_creation.cpp`"
2. **Easier maintenance**: All material-related code in one clearly-named file
3. **Better organization**: Functional grouping matches API documentation structure
4. **Reduced cognitive load**: No need to memorize arbitrary file numbering scheme

## Memory Management

### Exception Safety

All functions that allocate memory with `new` must follow this pattern for exception safety:

```cpp
Type* ptr = nullptr;
try {
	ptr = new Type(...);
	Handle h(ptr);
	ptr = nullptr;  // Ownership transferred to handle
	...
}
catch(...) {
	if(ptr) delete ptr;  // Cleanup if exception before ownership transfer
	Initialize();
	return 0;
}
```

**Key Points**:
- Initialize raw pointers to `nullptr` before `try` block
- Set to `nullptr` immediately after ownership transfer to handle (radThg, etc.)
- Clean up in `catch(...)` block if pointer is still non-null
- This prevents memory leaks when exceptions occur before ownership transfer

### RAII (Resource Acquisition Is Initialization)

Prefer RAII containers over manual memory management:

```cpp
// Good - RAII with std::vector
std::vector<radTPolygon> polygons;

// Avoid - Manual memory management
radTPolygon* polygons = new radTPolygon[n];  // Requires manual delete[]
```

## Material Specification (MatLin)

### Overview

`rad.MatLin()` defines **anisotropic linear magnetic materials** with two susceptibility values (parallel and perpendicular to the easy magnetization axis).

**Usage:** Magnetic materials that respond to external fields (NOT for permanent magnets)

### Two Forms

#### Form 1: Magnitude Only (Easy Axis from Object)
```python
mat = rad.MatLin([ksi_par, ksi_perp], mr)
```
- `mr`: Magnitude of remanent magnetization (scalar)
- Easy axis direction taken from the magnetization vector specified at object creation

```python
# Example: Easy axis from object's initial magnetization direction
obj = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])  # M in z-direction
mat = rad.MatLin([0.06, 0.17], 100)  # mr magnitude = 100, axis = z from object
rad.MatApl(obj, mat)
```

#### Form 2: Explicit Vector (Easy Axis Specified)
```python
mat = rad.MatLin([ksi_par, ksi_perp], [mrx, mry, mrz])
```
- `[mrx, mry, mrz]`: Remanent magnetization vector (defines easy axis explicitly)
- Easy axis direction independent of object's magnetization

```python
# Example: Easy axis explicitly defined
mat = rad.MatLin([0.06, 0.17], [0, 0, 100])  # Easy axis in z-direction, |Mr| = 100
rad.MatApl(obj, mat)
```

### Isotropic Linear Materials

For **isotropic materials** (same susceptibility in all directions):
```python
# Use equal susceptibilities
mat = rad.MatLin([ksi, ksi], [0, 0, 100])  # Mr still needed as reference
```

### Important Notes

1. **Permanent Magnets**: Do NOT use MatLin for permanent magnets (fixed magnetization)
   - Permanent magnets: `rad.ObjRecMag([x,y,z], [dx,dy,dz], [mx,my,mz])` only
   - No material needed for permanent magnets

2. **Remanent Magnetization (Mr)**:
   - Defines the easy magnetization axis (parallel direction)
   - Must be non-zero (either as magnitude or vector)
   - For isotropic materials, direction is arbitrary but still required

3. **Two Susceptibilities**:
   - ksi_par: Susceptibility parallel to easy axis
   - ksi_perp: Susceptibility perpendicular to easy axis
   - For isotropic: use `[ksi, ksi]`

## H-Matrix Field Evaluation Testing

### Test Requirements

H-matrix field evaluation tests **must** meet these minimum requirements:

```python
# Enable H-matrix explicitly
rad.SetHMatrixFieldEval(1, 1e-6)  # enable=1, eps=1e-6

# Minimum problem sizes
num_elements = 100      # At least 100 magnetic elements
num_points = 100        # At least 100 observation points

# Use H-matrix with explicit flag
H_field = rad.FldBatch(container, 'h', obs_points, 1)  # use_hmatrix=1
```

### Why These Requirements?

- **100+ elements**: Ensures H-matrix construction with both low-rank and full blocks
- **100+ observation points**: Provides meaningful accuracy statistics
- **Explicit enable**: User controls when to use H-matrix (no automatic thresholds)

### Test Structure

```python
# Test 1: Direct calculation (baseline)
H_direct = rad.FldBatch(container, 'h', obs_points, 0)  # use_hmatrix=0

# Test 2: H-matrix calculation
rad.SetHMatrixFieldEval(1, 1e-6)
H_hmat = rad.FldBatch(container, 'h', obs_points, 1)  # use_hmatrix=1

# Accuracy check
rel_error = np.linalg.norm(H_hmat - H_direct) / np.linalg.norm(H_direct)
assert rel_error < 0.01, f"H-matrix error {rel_error*100:.2f}% > 1%"
```

### Example Test

See `test_hmatrix_large.py` for a complete example with:
- 10×10 grid = 100 magnets
- 10×10 grid = 100 observation points
- Accuracy verification < 1%
- Performance comparison

---

**Last Updated**: 2025-11-09
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation

## Example Scripts VTK Output Policy

### Requirement: Consistent VTK File Naming

**Goal**: All example scripts in the `examples/` folder should export VTK files with the same base name as the Python script for easy identification and visualization.

**Policy**:
- Every example script (`*.py`) in `examples/` folder **must** generate a corresponding VTK file (`*.vtk`)
- VTK filename must match the Python script basename
- Example: `arc_current_with_magnet.py` → `arc_current_with_magnet.vtk`

**Rationale**:
- Easy identification: Users can immediately find the VTK file for each example
- Consistent workflow: Run example → View corresponding VTK in ParaView
- Prevents filename confusion and file accumulation
- Facilitates automated visualization pipelines

### Implementation Template

Add this code at the end of each example script:

```python
# VTK Export - Export geometry with same filename as script
try:
	from radia_vtk_export import exportGeometryToVTK
	import os

	# Get script basename without extension
	script_name = os.path.splitext(os.path.basename(__file__))[0]
	vtk_filename = f"{script_name}.vtk"
	vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

	# Export geometry
	exportGeometryToVTK(g, vtk_path)  # g = Radia object
	print(f"\n[VTK] Exported: {vtk_filename}")
	print(f"      View with: paraview {vtk_filename}")
except ImportError:
	print("\n[VTK] Warning: radia_vtk_export not available (VTK export skipped)")
except Exception as e:
	print(f"\n[VTK] Warning: Export failed: {e}")
```

### Required Changes

**Folders to update**:
1. `examples/simple_problems/` - Basic Radia examples
2. `examples/background_fields/` - Background field examples
3. `examples/electromagnet/` - Electromagnet simulations
4. `examples/complex_coil_geometry/` - Complex coil models
5. `examples/smco_magnet_array/` - Magnet array examples

**Exceptions**:
- Benchmark scripts in `examples/solver_benchmarks/` (optional)
- Test scripts in `examples/NGSolve_Integration/test_*.py` (optional)
- Analysis scripts in `examples/solver_benchmarks/` (optional)

### File Naming Examples

| Python Script | VTK Output | Status |
|--------------|------------|---------|
| `arc_current_with_magnet.py` | `arc_current_with_magnet.vtk` | ✓ Required |
| `quadrupole_analytical.py` | `quadrupole_analytical.vtk` | ✓ Required |
| `sphere_in_quadrupole.py` | `sphere_in_quadrupole.vtk` | ✓ Required |
| `magnet.py` | `magnet.vtk` | ✓ Required |
| `benchmark_solver.py` | N/A | Optional |

---

**Last Updated**: 2025-11-11

## Mesh File Preservation Policy

### Requirement: Preserve Mesh and Journal Files

**Goal**: Preserve mesh files and their source files used in examples to enable regeneration and modifications.

**Policy**:
- **NEVER DELETE** mesh files (`.bdf`, `.nas`, `.msh`, `.vtk`) from `examples/` directories
- **NEVER DELETE** Cubit journal files (`.jou`, `.journal`) if present
- **NEVER DELETE** mesh generation scripts (`cubit.py`, `generate_mesh.py`, etc.) if present
- These files are critical for:
  - Reproducing results
  - Understanding mesh generation process
  - Modifying geometries for new problems
  - Educational purposes

**Protected Files**:
- `examples/electromagnet/York.bdf` - Magnetic yoke mesh (Nastran format)
- `*.jou`, `*.journal` - Cubit journal files (mesh generation scripts)
- `cubit.py` - Python scripts that call Cubit API
- Any other mesh files (`.nas`, `.msh`, `.mesh`, `.vtk`) in examples/

**Rationale**:
- Mesh files are difficult to recreate without original CAD or mesh generation tools
- Journal files document the exact steps used to create meshes
- Enables users to modify geometries for their own problems
- Some meshes require expensive commercial software (Cubit, Ansys, etc.)

**If Mesh Files Are Missing**:
- Examples should gracefully degrade (e.g., run in "coil only" mode)
- Provide clear error messages indicating which files are optional
- Document in README where to obtain or how to generate missing mesh files

**Example - Graceful Degradation Pattern**:
```python
nas_file = os.path.join(os.path.dirname(__file__), 'York.bdf')

if os.path.exists(nas_file):
    yoke = create_yoke_from_nastran(nas_file)
    print(f"[OK] Loaded yoke mesh from {nas_file}")
else:
    print(f"[INFO] Mesh file not found: {nas_file}")
    print(f"[INFO] Continuing without yoke (coil-only mode)")
    yoke = None
```

---

**Last Updated**: 2025-11-12

## Windows Console Encoding (cp932) Compatibility

### Requirement: ASCII-Only Characters in Python Scripts

**Goal**: Ensure all Python scripts run correctly on Windows systems with cp932 (Japanese) console encoding without UnicodeEncodeError.

**Policy**:
- **NEVER use Unicode mathematical symbols** in print statements or string literals
- **ALWAYS use ASCII equivalents** for mathematical and special symbols
- This applies to all Python scripts, especially:
  - Benchmark scripts
  - Example scripts
  - Test scripts
  - Documentation strings

**Forbidden Unicode Characters and Their ASCII Replacements**:

| Unicode | Symbol | ASCII Replacement | Example |
|---------|--------|-------------------|---------|
| `\u00b2` | ² | `^2` | `N²` → `N^2` |
| `\u00b3` | ³ | `^3` | `N³` → `N^3` |
| `\u03b1` | α | `alpha` | `α = 2.0` → `alpha = 2.0` |
| `\u2192` | → | `->` | `A → B` → `A -> B` |
| `\u2248` | ≈ | `~=` | `x ≈ 2` → `x ~= 2` |
| `\u2264` | ≤ | `<=` | `N ≤ 100` → `N <= 100` |
| `\u2265` | ≥ | `>=` | `N ≥ 250` → `N >= 250` |
| `\u00d7` | × | `x` | `3×3` → `3x3` |
| `\u00b1` | ± | `+/-` | `±0.5` → `+/- 0.5` |
| `\u221e` | ∞ | `inf` | `→ ∞` → `-> inf` |

**Rationale**:
- Windows console (cmd.exe) defaults to cp932 encoding in Japanese environments
- cp932 cannot encode mathematical Unicode symbols
- This causes `UnicodeEncodeError` when print() tries to output these characters
- ASCII equivalents are universally readable and display correctly on all systems

**Prevention Pattern**:
```python
# BAD - Uses Unicode symbols
print("  Complexity: O(N²)")
print("  Range: 100 ≤ N < 250")
print("  Arrow: α → 2.0")

# GOOD - Uses ASCII equivalents
print("  Complexity: O(N^2)")
print("  Range: 100 <= N < 250")
print("  Arrow: alpha -> 2.0")
```

**Testing**:
- Before committing Python scripts, verify they run on Windows cp932 console
- Test command: `python script.py` (should complete without UnicodeEncodeError)
- If errors occur, replace Unicode symbols with ASCII equivalents

**Common Locations to Check**:
1. Print statements in benchmark scripts
2. Error messages and warnings
3. Summary and report outputs
4. Docstrings with mathematical notation
5. Comments with special symbols (less critical but still recommended)

---

**Last Updated**: 2025-11-13
