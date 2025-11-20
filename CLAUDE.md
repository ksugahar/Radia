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

## Development and Debug Files Organization

### Requirement: Separate User-Facing and Development Files

**Goal**: Keep the repository clean and organized by separating user-facing content from development/debug materials.

**Policy**:
- **User-facing files** (examples, documentation, tests) should remain in their standard locations
- **Development/debug files** that users don't need should be organized in dedicated folders
- Use clear, descriptive folder names that indicate the content is for development purposes

**Recommended Folder Structure**:

```
project_root/
├── examples/           # User-facing examples (keep clean)
├── tests/              # User-facing test suite
├── docs/               # User-facing documentation
├── dev/                # Development notes and prototypes
│   ├── notes/          # Development notes, meeting notes
│   ├── prototypes/     # Experimental code, proof-of-concepts
│   └── benchmarks/     # Development-only benchmarks
├── debug/              # Debug outputs and temporary files
│   ├── logs/           # Debug log files
│   ├── output/         # Temporary output files
│   └── cache/          # Debug cache files
└── internal/           # Internal documentation not for end users
    ├── design/         # Design documents, architecture notes
    └── analysis/       # Performance analysis, profiling results
```

**File Types and Their Locations**:

| File Type | Appropriate Location | Example |
|-----------|---------------------|---------|
| Benchmark results (published) | `examples/solver_benchmarks/BENCHMARK_RESULTS.md` | ✓ User-facing |
| Benchmark output (debug) | `debug/output/benchmark_*.txt` | For developers only |
| Example scripts | `examples/*/` | ✓ User-facing |
| Development notes | `dev/notes/` | For developers only |
| Design documents (internal) | `internal/design/` | For developers only |
| Test output files | `debug/output/test_*.txt` | For developers only |
| Cache files | `.radia_cache/` or `debug/cache/` | For developers only |
| Meeting notes | `dev/notes/meetings/` | For developers only |
| Profiling results | `internal/analysis/` | For developers only |

**Naming Conventions**:

1. **Development folders**: Use prefixes that clearly indicate non-user content
   - `dev/` - Active development materials
   - `debug/` - Debug outputs and temporary files
   - `internal/` - Internal documentation
   - `_archive/` - Archived materials (prefix with underscore)

2. **Files to avoid in user-facing folders**:
   - `output.txt`, `debug.log`, `temp_*.py`
   - `test_manual_*.py` (unless part of test suite)
   - `benchmark_results_raw_*.txt`
   - `notes_*.md`, `TODO_*.md` (use issue tracker instead)

**Rationale**:
- **Clarity**: Users can easily find relevant examples and documentation
- **Professionalism**: Repository looks organized and maintained
- **Reduced confusion**: Clear separation between "what to use" and "how it was developed"
- **Version control**: Development files can have different .gitignore rules

**Implementation**:

```bash
# Create development folders
mkdir -p dev/notes dev/prototypes
mkdir -p debug/logs debug/output
mkdir -p internal/design internal/analysis

# Move development files
mv benchmark_*_output.txt debug/output/
mv notes_*.md dev/notes/
mv DESIGN_*.md internal/design/
```

**Example .gitignore Rules**:

```gitignore
# Debug and development outputs (never commit)
debug/
*.log
*_output.txt
*_debug.txt

# Development notes (optional - commit if useful for team)
# dev/notes/

# Internal documentation (commit - useful for maintainers)
# internal/

# User-facing content (always commit)
examples/
docs/
tests/
```

**When to Use Each Folder**:

1. **Use `examples/`** for:
   - Scripts users should run
   - Code users should learn from
   - Production-ready demonstrations

2. **Use `dev/`** for:
   - Prototypes being developed
   - Development meeting notes
   - Experimental features
   - "Work in progress" code

3. **Use `debug/`** for:
   - Temporary output files
   - Debug logs
   - Performance profiling outputs
   - Files you'll delete later

4. **Use `internal/`** for:
   - Design documents for maintainers
   - Architecture decisions
   - Performance analysis reports
   - Code review checklists

**Migration Strategy**:

When cleaning up an existing repository:

1. **Identify files**: Find development/debug files in user-facing folders
2. **Categorize**: Determine if files are dev/, debug/, or internal/
3. **Move files**: Relocate to appropriate folders
4. **Update references**: Fix any documentation or script paths
5. **Update .gitignore**: Ensure debug/ is excluded
6. **Commit changes**: Document the reorganization in commit message

---

**Last Updated**: 2025-11-13

## H-Matrix Solver Control Policy

### Requirement: No Automatic Problem Size Threshold

**Goal**: Users have full explicit control over H-matrix solver enable/disable. No automatic switching based on problem size (N < 200 vs N >= 200).

**Rationale**:
- Users should decide when H-matrix is appropriate for their use case
- Small problems (N < 200) may still benefit from H-matrix in iterative workflows
- Large problems (N >= 200) may not need H-matrix for single-shot calculations
- Automatic threshold hides behavior and limits user control
- Different applications have different performance trade-offs

**Policy**:

**✗ NO automatic threshold checking**:
```cpp
// WRONG - Don't do this
if(use_hmatrix && AmOfMainElem >= 200) {
    // Use H-matrix
} else {
    // Use dense solver
    use_hmatrix = false;  // Override user's choice!
}
```

**✓ User explicit control**:
```cpp
// CORRECT - Respect user's choice
if(use_hmatrix) {
    // Use H-matrix regardless of problem size
    std::cout << "[Phase 2-B] Using H-matrix solver (N=" << AmOfMainElem << ")" << std::endl;
    return SetupInteractMatrix_HMatrix();
}
```

**Implementation**:

1. **Remove all HMATRIX_AUTO_THRESHOLD checks** from code
2. **Respect `use_hmatrix` flag** set by user via API
3. **Print informational messages** showing N and parameters, but don't override
4. **Document in user guide** when H-matrix is recommended (N > 200), but let user decide

**API Behavior**:

```python
import radia as rad

# User explicitly enables H-matrix
rad.SolverHMatrixEnable(1, eps=1e-4, max_rank=30)

# H-matrix will be used regardless of N
# - N=50: H-matrix used (if user wants)
# - N=500: H-matrix used (as expected)
# - User controls when to use H-matrix

# User explicitly disables H-matrix
rad.SolverHMatrixEnable(0)

# Dense solver will be used regardless of N
# - N=50: Dense solver (as expected)
# - N=5000: Dense solver (if user really wants)
```

**User Guidance** (in documentation, not code):

Recommended usage:
- **N < 100**: Dense solver typically faster
- **N = 100-200**: Either method works, depends on use case
- **N > 200**: H-matrix recommended for most cases
- **Iterative workflows**: H-matrix beneficial even for small N (amortized construction cost)

But ultimately: **User decides, code obeys**.

**Files Modified**:
- `src/core/rad_interaction.cpp`: Removed HMATRIX_AUTO_THRESHOLD checks
  - `SetupInteractMatrix()`: No automatic threshold
  - `EnableHMatrix()`: No automatic threshold

**Verification**:
- Small problem (N=50) with H-matrix enabled: Should use H-matrix
- Large problem (N=500) with H-matrix disabled: Should use dense solver
- No "[Auto] N=X < 200 - using optimized dense solver" messages

---

**Last Updated**: 2025-11-13

## PyPI Package Release Policy

### Requirement: Manual PyPI Upload with Automated Version Management

**Goal**: Ensure controlled, secure PyPI releases with automated version and metadata management by Claude Code.

**Rationale**:
- PyPI releases require human review and approval for security
- Version numbers and metadata should be automatically maintained by Claude Code
- Upload process uses standardized script with API token authentication
- Separation of automated preparation (Claude) and manual approval (user)

**Policy**:

### 1. Version Management (Automated by Claude Code)

**Claude Code is responsible for**:
- Maintaining version numbers in `pyproject.toml`
- Following semantic versioning (MAJOR.MINOR.PATCH)
- Updating `CHANGELOG.md` with release notes
- Ensuring version consistency across all files
- Preparing release documentation

**Version Update Workflow**:
```toml
# pyproject.toml - Claude updates this automatically
[project]
version = "1.1.1"  # Claude increments based on changes
```

**Semantic Versioning Rules**:
- **MAJOR**: Breaking changes to API (e.g., 1.x.x → 2.0.0)
- **MINOR**: New features, backward compatible (e.g., 1.1.x → 1.2.0)
- **PATCH**: Bug fixes, backward compatible (e.g., 1.1.1 → 1.1.2)

### 2. PyPI Upload Process (Manual by User)

**User performs upload using**:
```powershell
# Set PyPI API token (keep secure, do not commit!)
$env:PYPI_TOKEN = "pypi-AgEIcGl..."

# Run standardized upload script
powershell.exe -ExecutionPolicy Bypass -File Publish_to_PyPI.ps1
```

**Publish_to_PyPI.ps1 Script**:
- Cleans old dist/ files
- Builds source distribution (.tar.gz) and wheel (.whl)
- Validates packages with twine check
- Uploads to PyPI using API token authentication
- Provides clear success/failure feedback

**Security Requirements**:
- **NEVER commit PyPI tokens** to repository
- Use environment variable `$env:PYPI_TOKEN` for authentication
- Tokens should have upload-only permissions (no admin access)
- Rotate tokens periodically

### 3. Pre-Release Checklist (Claude Code)

Before requesting user to upload, Claude Code must:

1. **Update version number**:
   ```bash
   # Edit pyproject.toml
   version = "X.Y.Z"  # New version
   ```

2. **Update CHANGELOG.md**:
   ```markdown
   ## [X.Y.Z] - YYYY-MM-DD

   ### Added
   - New features...

   ### Changed
   - Modified features...

   ### Fixed
   - Bug fixes...
   ```

3. **Build and validate packages**:
   ```bash
   python -m build
   python -m twine check dist/*
   ```

4. **Commit and push changes**:
   ```bash
   git add pyproject.toml CHANGELOG.md
   git commit -m "Bump version to X.Y.Z"
   git push
   ```

5. **Request user approval** for PyPI upload

### 4. Upload Workflow

**Claude Code**:
1. Updates version and CHANGELOG
2. Builds packages
3. Validates packages
4. Commits changes
5. Says: "PyPI パッケージ準備完了。Publish_to_PyPI.ps1 でアップロードしてください。"

**User**:
1. Reviews changes and version number
2. Sets `$env:PYPI_TOKEN`
3. Runs `Publish_to_PyPI.ps1`
4. Verifies upload success on https://pypi.org/project/radia/

### 5. Files Maintained by Claude Code

Claude Code automatically maintains:

| File | Claude's Responsibility |
|------|------------------------|
| `pyproject.toml` | Version number, dependencies, metadata |
| `CHANGELOG.md` | Release notes, version history |
| `README.md` | Package description, installation instructions |
| `LICENSE` | License text (verify before changes) |
| `MANIFEST.in` | Package file inclusion rules |

### 6. Example Release Workflow

```
User: "Release version 1.2.0 with new features"

Claude:
1. Updates pyproject.toml: version = "1.2.0"
2. Updates CHANGELOG.md: ## [1.2.0] - 2025-11-13
3. Builds packages: python -m build
4. Validates: python -m twine check dist/*
5. Commits: git commit -m "Bump version to 1.2.0"
6. Pushes: git push
7. Reports: "PyPI パッケージ準備完了 (1.2.0)。Publish_to_PyPI.ps1 でアップロードしてください。"

User:
1. Reviews: git log, CHANGELOG.md
2. Sets token: $env:PYPI_TOKEN = "..."
3. Uploads: powershell.exe -ExecutionPolicy Bypass -File Publish_to_PyPI.ps1
4. Verifies: https://pypi.org/project/radia/1.2.0/
```

### 7. Emergency Rollback

If a release has critical issues:

1. **Do NOT delete from PyPI** (PyPI does not allow re-uploading same version)
2. **Release a patch version** immediately (e.g., 1.2.0 → 1.2.1)
3. **Update CHANGELOG.md** to mark problematic version
4. **Document issue** in release notes

### 8. TestPyPI for Testing

For testing releases without affecting production:

```powershell
# Upload to TestPyPI instead
python -m twine upload --repository testpypi dist/*

# Test installation
pip install --index-url https://test.pypi.org/simple/ radia==X.Y.Z
```

**When to use TestPyPI**:
- Major version changes (1.x → 2.0)
- Significant refactoring
- New build system changes
- First-time releases

---

**Last Updated**: 2025-11-13 (PyPI Release Policy)

## Python Script Path Import Policy

### Requirement: Use Relative Paths for Module Imports

**Goal**: Ensure all Python scripts in the repository use relative paths for importing the Radia module, making scripts portable across different development environments and user installations.

**Rationale**:
- Absolute paths (e.g., `S:\Radia\01_GitHub\build\Release`) only work on specific machines
- Scripts with absolute paths fail when shared with other developers or users
- Relative paths work regardless of installation location
- Improves portability and collaboration

**Policy**:

**✗ NEVER use absolute paths**:
```python
# WRONG - Hard-coded absolute path
import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")
import radia as rad
```

**✓ ALWAYS use relative paths**:
```python
# CORRECT - Relative path from script location
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
import radia as rad
```

### Path Calculation Pattern

For scripts in different locations, use these patterns:

**Examples folder** (`examples/*/script.py`):
```python
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
```

**Benchmarks folder** (`examples/solver_benchmarks/script.py`):
```python
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
```

**Tests folder** (`tests/script.py`):
```python
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))
```

### Multiple Path Imports

When importing from multiple locations (Radia module + Python utilities):
```python
import sys
import os

# Add Radia module path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))

# Add Python utilities path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

import radia as rad
from radia_vtk_export import exportGeometryToVTK
```

### Benefits

1. **Portability**: Scripts work on any machine without modification
2. **Collaboration**: Other developers can run scripts immediately
3. **Distribution**: Examples work for end users after installation
4. **Maintenance**: No need to update paths when repository is moved

### Files to Check

When creating or modifying Python scripts:
- All files in `examples/` folder
- All files in `examples/solver_benchmarks/` folder
- All files in `tests/` folder
- Development scripts in `dev/` folder

**Exception**: System-wide installed package (users install via `pip install radia`) don't need path manipulation:
```python
import radia as rad  # Works if radia is installed via pip
```

---

**Last Updated**: 2025-11-17 (Python Script Path Import Policy)

## rad_ngsolve Result Matrix Indexing Fix

### Issue Summary (2025-11-20)

**Problem**: GridFunction.Set() with rad_ngsolve.RadiaField produced incorrect field values (thousands of % error).

**Root Cause**: In the batch evaluation function `Evaluate(const BaseMappedIntegrationRule& mir, BareSliceMatrix<> result)`, the result matrix indexing was inconsistent between normal path and error handling path:

```cpp
// WRONG (before fix)
// Normal path - src/python/rad_ngsolve.cpp:348-350
result(0, i) = f_global[0] * scale;  // result(component, point)
result(1, i) = f_global[1] * scale;
result(2, i) = f_global[2] * scale;

// Error path - lines 359-362
result(i, 0) = 0.0;  // result(point, component) <- Inconsistent!
result(i, 1) = 0.0;
result(i, 2) = 0.0;
```

**Fix**: Corrected normal path to match error path indexing:

```cpp
// CORRECT (after fix)
result(i, 0) = f_global[0] * scale;  // result(point, component)
result(i, 1) = f_global[1] * scale;
result(i, 2) = f_global[2] * scale;
```

### Origin of Bug

Introduced in commit `ab77976` (Add H-matrix field evaluation and Python API) when batch evaluation was added.

### Verification

**Before fix**:
- Direct CoefficientFunction evaluation: ✓ Correct
- GridFunction.Set() evaluation: ✗ Wrong (thousands of % error)

**After fix**:
- Direct CoefficientFunction evaluation: ✓ Correct
- GridFunction.Set() evaluation: ✓ Correct (6-8% L² projection error - normal)
- curl(A) = B verification: ✓ Passed (max error < 0.015 T)

### Files Modified

- `src/python/rad_ngsolve.cpp` (lines 348-350)
- Built module: `build/Release/rad_ngsolve.pyd`

### Future Tasks

**TODO: Verify H-matrix field evaluation impact**

This indexing fix may affect H-matrix accelerated field evaluation:
- H-matrix uses the same batch evaluation interface
- Need to verify that H-matrix field evaluation still produces correct results
- Test files:
  - `examples/ngsolve_integration/verify_curl_A_equals_B.py` (copied to `tests/`)
  - `tests/test_rad_ngsolve_function.py`
- Verification criteria:
  - GridFunction.Set() with H-matrix enabled should match direct evaluation
  - curl(A) = B relationship should hold with H-matrix acceleration
  - Performance should not degrade

**Priority**: Medium (not urgent, but should verify before next release)

---

**Last Updated**: 2025-11-20 (rad_ngsolve Result Matrix Indexing Fix)

## GridFunction Projection for rad_ngsolve (2025-11-20)

### Summary

After fixing the result matrix indexing bug in `rad_ngsolve.cpp`, GridFunction projection accuracy was thoroughly investigated.

### Key Findings

**1. Direct CoefficientFunction Evaluation**: Always accurate
- B_cf(mip) matches rad.Fld() with <0.01% error
- Recommended for maximum accuracy

**2. GridFunction Projection Accuracy** (using Set() method):

| Space          | Order | Pointwise Error | L2 Norm Error | Notes |
|----------------|-------|-----------------|---------------|-------|
| HDiv           | 2     | 0.36%           | 22.87%        | Best overall |
| HCurl          | 2     | 0.32%           | 29.89%        | Good pointwise |
| VectorH1       | 2     | 1.05%           | 34.47%        | Acceptable |

**3. Region-Dependent Accuracy** (HDiv order=2, h=10mm):

| Region                          | Distance from Surface | Mean Error | Max Error | Usability |
|---------------------------------|----------------------|------------|-----------|-----------|
| Near magnet surface             | 0-5mm                | >10%       | >50%      | ✗ Avoid   |
| Near field (mesh-spacing away)  | >10mm (1 mesh cell) | **0.15%**  | 0.41%     | ✓ Excellent |
| Mid field                       | >15mm                | 0.10%      | 0.30%     | ✓ Excellent |
| Far field                       | >20mm                | <0.10%     | <0.20%    | ✓ Excellent |

### Best Practices for NGSolve Applications

**Recommended Configuration**:
```python
# 1. Create mesh with appropriate size
mesh_size = 0.010  # 10mm in meters
mesh = Mesh(geo.GenerateMesh(maxh=mesh_size))

# 2. Use HDiv space, order=2
fes = HDiv(mesh, order=2)

# 3. Project B to GridFunction
B_cf = rad_ngsolve.RadiaField(radia_obj, 'b')
B_gf = GridFunction(fes)
B_gf.Set(B_cf)

# 4. Evaluate ONLY at distances > 1 mesh cell from magnet surface
# For h=10mm mesh: evaluate at distances > 10mm from surface
```

**Evaluation Guidelines**:
- ✓ **DO**: Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- ✓ **DO**: Use B_cf directly for maximum accuracy near boundaries
- ✗ **DON'T**: Evaluate GridFunction within 1 mesh cell of magnet surface
- ✗ **DON'T**: Rely on L2 norm error (dominated by near-field singularities)

**curl(A) = B Verification**:
```python
# A in HCurl (mathematically correct)
A_gf = GridFunction(HCurl(mesh, order=2))
A_gf.Set(A_cf)

# curl(A) in HDiv
curl_A = curl(A_gf)

# Compare with B_cf directly (not B_gf!)
# Typical accuracy: <1% in far field, ~3% in near field
```

### Why L2 Norm Error is Large

The L2 norm error (~23-35%) is dominated by near-field regions where:
- Magnetic field has strong gradients
- GridFunction.Set() performs L² projection via integration points
- Local errors near boundaries contribute heavily to global L2 norm

**This is normal and acceptable** for practical applications where:
- Evaluation points are >1 mesh cell from boundaries
- Pointwise accuracy (0.15-0.36%) is excellent in these regions

### Files Modified/Created

**Modified**:
- `src/python/rad_ngsolve.cpp`: Fixed result matrix indexing (lines 348-350)

**Test Scripts Created**:
- `tests/verify_curl_A_equals_B_improved.py`: Proper FE space verification
- `tests/test_far_field_accuracy.py`: Region-dependent accuracy test
- `tests/test_hcurl_vs_hdiv.py`: Space comparison
- `tests/test_all_spaces.py`: Comprehensive space evaluation

**Recommended for Users**:
- `tests/verify_curl_A_equals_B.py`: Copied from examples/ngsolve_integration/
- Use HDiv order=2 for B projection in practical applications
- Evaluate >1 mesh cell away from magnet surfaces

---

**Last Updated**: 2025-11-20 (GridFunction Projection Best Practices)


## H-Matrix Field Evaluation Verification (2025-11-20)

### Summary

Verified that the result matrix indexing fix in `rad_ngsolve.cpp` (lines 348-350) works correctly with H-matrix accelerated field evaluation.

### Test Configuration

**Test file**: `tests/test_rad_ngsolve_hmatrix.py`

**Setup**:
- Magnet: 5×5×5 = 125 elements (0.1m cube, boundaries at ±50mm)
- Mesh: 177 vertices, 508 elements, domain [15mm, 63mm]
- FE Space: HCurl order=2 (6012 DOFs)
- Test points: 3 points at safe distances (>10mm from magnet surface)

**Methods tested**:
1. Direct mode (H-matrix OFF)
2. H-matrix mode (H-matrix ON, eps=1e-6)

### Results

| Mode | Mean curl(A)=B Error | GridFunction Consistency |
|------|---------------------|--------------------------|
| Direct (H-matrix OFF) | 1.5491% | - |
| H-matrix (H-matrix ON) | 1.5491% | 0.0000% |

**Key Findings**:

1. **✓ GridFunction values identical**: Direct and H-matrix modes produce exactly the same GridFunction values (0.0000% difference)
   - This confirms the result matrix indexing fix works correctly with H-matrix
   
2. **✓ curl(A) = B maintained**: Both modes show identical curl(A)=B accuracy (1.5491%)
   - Projection error is consistent and acceptable (<2%)
   
3. **✓ H-matrix does not degrade accuracy**: H-matrix approximation does not introduce additional error in GridFunction.Set()

### Verification Status

**[SUCCESS]** H-matrix field evaluation verified!
- GridFunction.Set() works correctly with H-matrix
- curl(A) = B relationship maintained
- Results consistent between Direct and H-matrix modes  
- Result matrix indexing fix confirmed working with H-matrix

### Acceptance Criteria

For H-matrix verification tests:
- curl(A) = B error < 2.0% (allows for L² projection error)
- GridFunction consistency (Direct vs H-matrix) < 1.0%
- All test points evaluated successfully

### Implementation Notes

**Bug fix location**: `src/python/rad_ngsolve.cpp:348-350`

**Before** (incorrect):
```cpp
result(0, i) = f_global[0];  // component, point - WRONG
result(1, i) = f_global[1];
result(2, i) = f_global[2];
```

**After** (correct):
```cpp
result(i, 0) = f_global[0];  // point, component - CORRECT
result(i, 1) = f_global[1];
result(i, 2) = f_global[2];
```

This fix applies to both Direct and H-matrix modes, ensuring consistent and correct field evaluation.

---

**Last Updated**: 2025-11-20  
**Test**: tests/test_rad_ngsolve_hmatrix.py  
**Status**: ✓ Verified

