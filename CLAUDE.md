# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and policies for the Radia project when working with Claude Code.

---

## Critical Policies

### Green's Function: Laplace Kernel Only (MQS/Darwin)

**POLICY**: Radia uses **Laplace kernel only**: $G(r) = 1/(4\pi r)$. Target regime is MQS (Magneto-Quasi-Static) to Darwin approximation.

**Do NOT**:
- Add Helmholtz kernel ($e^{-jkr}/r$) to any Green's function
- Use wave number $k$ in field calculations (except for skin depth)
- Implement full-wave EFIE or MFIE formulations

Skin depth is computed from frequency for SIBC, but field propagation uses quasi-static approximation.

**Affected Components**: `rad_green_fullwave.h/cpp`, `rad_conductor.cpp` (`GreenFunction()`), `rad_hacapk.cpp`.

### Matrix Storage: Row-Major (C-style)

**POLICY**: All interaction matrices use **row-major [target][source] format**.
- `A[i][j]` stored at `i * stride + j`; represents effect ON target i FROM source j
- All BLAS calls use `CblasRowMajor`
- Python interface returns NumPy C-contiguous (row-major) arrays

**Source Files**: `rad_interaction.cpp`, `rad_relaxation_methods.cpp`, `rad_hacapk.cpp`.

### Binary File Policy

**POLICY**: No binary files (`.pyd`, `.dll`, `.so`, `.lib`, `.exe`) in the git repository.
- Hosted on GitHub Releases (tag: `binaries`)
- Pre-push hook auto-uploads `.pyd` on `git push`
- After cloning, run `./download_binaries.sh` to fetch binaries
- `.png`, `.pdf` allowed in repository; `.msh`, `.vtu`, `.vtk`, `.vol` are gitignored

### File Placement Policy

**POLICY**: Generated output files (`.png`, `.msh`, `.vtu`, `.vtk`, `.vol`, `.vts`) must be placed **next to their corresponding `.py` script**.
- Example outputs belong in `examples/<category>/` alongside their script
- Do NOT place generated files at the repository root
- `.msh` files in `examples/**/gmsh_models/` are tracked (pre-generated mesh definitions)
- Build output goes to `build*/` or `dist/` (both gitignored)

### Unit System Policy

**POLICY**: Radia always uses **meters**. There is no unit conversion in C++. All coordinates are in meters, all current densities in A/m^2.

**`FldUnits` is removed**: Do NOT call `rad.FldUnits()` in any code. Radia always uses meters with no configuration needed.

```python
# CORRECT
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])  # meters, A/m

# WRONG - hard-coded conversion
x_mm = x_m * 1000.0  # DO NOT DO THIS
```

**Radia Units** (always meters, no conversion):
- All coordinates in meters
- B in Tesla, H in A/m, A in T*m
- Current density J in A/m^2
- Physical constants in `rad_constants.h`: `MU_0_OVER_FOUR_PI = 1e-7`, `INV_FOUR_PI = 1/(4*pi)`

### Magnetization Units: A/m (NOT Tesla)

**POLICY**: Radia uses **M in A/m**. Common conversion: `M = Br / mu_0` (e.g., Br=1.2T -> M=954930 A/m).

Do NOT confuse M (A/m) with J (magnetic polarization, Tesla): J = mu_0 * M.

### Windows Console Encoding (cp932)

**POLICY**: NEVER use Unicode mathematical symbols in print statements. Use ASCII equivalents: `^2` not `²`, `->` not `→`, `<=` not `≤`, etc. Windows console defaults to cp932 in Japanese environments.

### Naming Policy: External Project References

**POLICY**: Do NOT use "ELF" or "ELF_MAGIC" in Radia source code, documentation, or comments. Radia is an independent project. Academic citations are allowed.

### No Console Output from C++ Code

**POLICY**: No `printf`/`cout`/`cerr` in C++ code for logging. All user-facing output through Python. Allowed: error messages via `Send.ErrorMessage(...)` and `#ifdef DEBUG_...` guards.

### Field Comparison: Vector Difference

**POLICY**: Compare magnetic fields using **vector difference** `norm(B1 - B2)`, not scalar magnitude difference `abs(|B1| - |B2|)`. Magnetic field is a vector quantity.

### FMM (Fast Multipole Method): Removed (2026-03-06)

**ExaFMM-t was removed from the repository**. Do NOT re-implement FMM acceleration.

**Why FMM failed for Radia**:

1. **Dipole approximation accuracy is poor for MSC elements**: MSC (surface charge) elements have distributed charge on 4-8 faces. A single dipole m=M*V approximates this poorly at intermediate distances (r ~ 2-5 element sizes). The O((a/r)^2) error is unacceptable for engineering accuracy.

2. **FMM Solve (Method 3) was useless**: Compact geometries (C-type magnets, iron yokes) have 87% near-field pairs. Near-field correction memory equals the full dense matrix, eliminating FMM's O(N log N) advantage. HACApK (H-matrix, Method 2) is 10-100x faster because ACA+ compression works on the same near-field blocks.

3. **FMM field evaluation had no benefit over direct**: For typical Radia models (N < 10,000 elements), direct B_genComp with TaskManager parallelization is fast enough. FMM overhead (tree build, M2L translation) exceeds direct computation time for these sizes.

4. **HACApK covers all large-scale needs**: H-matrix acceleration (ACA+) provides O(N log N) memory and O(N log^2 N) MatVec for the interaction matrix, which is the actual bottleneck.

**Lesson**: FMM is effective for point charges/dipoles in unbounded space (N-body). It is NOT effective for BEM/MSC where source distributions are extended (face integrals) and geometries are compact.

---

## Architecture Overview

### Development Strategy: Complement NGSolve

Radia's role is to **complement NGSolve**, not compete with it. Focus on areas where FEM is weak.

```
┌─────────────────────────────────────────────────────────────────┐
│                    Electromagnetic Analysis                      │
├─────────────────────────────────────────────────────────────────┤
│  NGSolve (FEM)              │  Radia (BEM/Integral Methods)     │
│  ───────────────────────────│──────────────────────────────────│
│  OK: Bounded domains        │  OK: Unbounded domains (open BC) │
│  OK: Complex geometry       │  OK: Permanent magnets (no mesh) │
│  OK: Nonlinear materials    │  OK: Thin conductors (PEEC)      │
│  OK: Transient analysis     │  OK: SPICE circuit extraction    │
│  OK: Multi-physics coupling │  OK: Model order reduction (MOR) │
│  WEAK: Open boundary (PML)  │  OK: Natural open boundary       │
│  WEAK: Thin structures      │  OK: Surface impedance (SIBC)    │
│  WEAK: Circuit parameters   │  OK: L, R, C, M extraction       │
└─────────────────────────────────────────────────────────────────┘
```

**Do NOT Implement** (use existing libraries):
- FEM solvers (use NGSolve)
- General sparse solvers (use MKL/MUMPS)
- Full-wave BEM (use ngbem for high frequency)
- CAD geometry kernels (use OpenCASCADE via NGSolve)
- PEEC from scratch (use PAMELA)
- Custom H-matrix algorithms (use HACApK)

**Radia C++ Core** (maintain and enhance):
1. MMM - Magnetic Moment Method for permanent magnets and soft iron
2. MSC - Magnetic Surface Charge for hexahedra/tetrahedra
3. Field computation - B, H, A, Phi in unbounded domains
4. NGSolve integration - RadiaField CoefficientFunction

### Solver Methods: MMM and MSC

| Method | Element | DOF | Description |
|--------|---------|-----|-------------|
| **MMM** | Tetrahedra (4 faces) | 3 (Mx, My, Mz) | Magnetic dipole distributions |
| **MSC** | Hexahedra (6 faces) | 6 (sigma/face) | Surface charge solid angle integration |
| **MSC** | Wedges (5 faces) | 5 (sigma/face) | Transition elements |

**Mixed Element Support**: All solvers (LU, BiCGSTAB, HACApK) support mixed hex+wedge+tet meshes. Variable DOF offset arrays: `m_elemDOF`, `m_elemDOFOffset`, `m_totalDOF`.

**BiCGSTAB Block Jacobi**: Automatically switches to block Jacobi preconditioner when diagonal ratio > 10 or min dominance < 0.1 (distorted elements). Uses LAPACK `dgetrf_`/`dgetri_` for block inversion.

**Interaction Matrix Blocks** (mixed elements):
- **3x3** (tet-tet), **5x5** (wedge-wedge), **6x6** (hex-hex)
- **5x6 / 6x5** (wedge-hex cross), **3x6 / 3x5 / 6x3 / 5x3** (tet-hex/wedge cross)
- Implementation: `SetupInteractMatrix_VariableDOF()`, compile flag `RADIA_MSC_SUPPORT`

### Unified Field Computation Architecture

**POLICY**: All field computation MUST use `rad_field_unified.h/cpp`.

```
┌─────────────────────────────────────────────────────────────────┐
│                    rad_field_unified.h/cpp                       │
│  ─────────────────────────────────────────────────────────────  │
│  ComputeFieldSingle()     - Single point, static field          │
│  ComputeFieldBatch()      - Batch points, TaskManager parallelized │
│  ComputeComplexFieldSingle() - Complex (AC) field               │
│  ComputeComplexFieldBatch()  - Complex batch with TaskManager   │
│  IsPointInsideAnyElement() - Inside/outside classification      │
│  ComputeBFromMagnetization() - Dipole field from M (complex)    │
└─────────────────────────────────────────────────────────────────┘
                              │
           ┌──────────────────┼──────────────────┐
           ▼                  ▼                  ▼
    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
    │ rad.Fld()   │    │ rad.FldVTS()│    │ CplMagFld() │
    │ rad.FldBatch│    │ VTS export  │    │ PEEC+MMM    │
    └─────────────┘    └─────────────┘    └─────────────┘
```

**Users**: `rad.Fld()`, `rad.FldBatch()`, `rad.FldVTS()`, `rad.RadiaField()`, `rad_particle_trajectory`, `CplMagFld()`.

**Key Features**: Inside/outside classification (solid angle method), TaskManager parallelized batch, complex field support (PEEC+MMM AC).

### Field Calculation: Surface Current vs Surface Charge

- **ObjRecMag**: Surface current model (rectangular blocks). 8-corner BufVect formula, efficient and non-cancelling on symmetry axes.
- **ObjHexahedron/ObjTetrahedron**: Surface charge model (general polyhedra). Face-based solid angle integration. A field may be zero on symmetry axes (mathematical cancellation, not a bug).

**rad.Fld() inside materials**: MMM gives dipole approximations inside materials; MSC gives uniform field per element. For validation, compare sigma values or external field points, not internal fields.

### Vector Potential A Field

A field is **implemented** for all element types using face integration (Wilton et al. formula). Formula: `A = (mu_0/4pi) * (M x BufVect)`. Satisfies `B = curl(A)` (verified numerically). Verification script: `examples/ngsolve_integration/verify_curl_A_equals_B/`.

### User-Facing Element APIs

- `rad.ObjRecMag(center, dimensions, magnetization)` -- Rectangular magnets (optimized formulas)
- `rad.ObjHexahedron(vertices, magnetization)` -- Arbitrary hexahedra (8 vertices)
- `rad.ObjTetrahedron(vertices, magnetization)` -- Tetrahedra (4 vertices)
- `rad.ObjWedge(vertices, magnetization)` -- Wedges (6 vertices)
- Mesh import functions (`netgen_mesh_to_radia`, `gmsh_to_radia`) for complex geometries

### EIEM2 Evaluation Point Convention

**POLICY**: The MSC interaction matrix evaluation point for face `i` is:
```cpp
EvalPt = 0.5 * (FaceCenter[i] + ElementCenter)
```
Do NOT change this. This matches ELF's EIEM2 convention exactly.

**MSC Source Files**: `rad_polyhedron.cpp` (element dispatch), `rad_poly_analytical.cpp` (triangle/quad integration), `rad_interaction.cpp` (interaction matrix, `PrecomputeHexaGeometry()`).

See `docs/MSC_QUICK_START.md` for quick start guide.

---

## API Guardrails

### Common Mistakes Checklist

**1. ObjBckg Requires Callable (CRITICAL)**
```python
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])  # CORRECT
bkg = rad.ObjBckg([0, 0, 0.1])             # WRONG - not a callable
```

**2. UtiDelAll() Cleanup**: Every script must call `rad.UtiDelAll()` before exiting.

**3. Relative Path Imports**:
```python
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))  # CORRECT
sys.path.insert(0, r'S:\Radia\01_GitHub\src\radia')  # WRONG - machine-specific
```

**4. MatLin Usage**: For isotropic materials, ALWAYS use single-argument form `MatLin(mu_r)`. MatLin is for soft magnetic materials only -- permanent magnets specify magnetization directly in `ObjHexahedron(vertices, [Mx, My, Mz])`.

**5. Docstring Units**: Use "in constructor length units", not "in mm".

**6. State Mutation**: Computation methods must NOT leave object state inconsistent on exception.

### Background Field API

```python
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])      # Uniform 0.1T in z
bkg = rad.ObjBckg(quadrupole_field_function)    # Spatially varying
container = rad.ObjCnt([mag_obj, bkg])
rad.Solve(container, 0.0001, 1000, 1)
```
Legacy `ObjBckg([Bx, By, Bz])` array form is NOT supported. Callback receives `[x, y, z]` in current units and returns `[Bx, By, Bz]` in Tesla.

### Memory Management

```cpp
// Exception-safe pattern
Type* ptr = nullptr;
try {
    ptr = new Type(...);
    Handle h(ptr);
    ptr = nullptr;  // Ownership transferred
} catch(...) {
    if(ptr) delete ptr;
    Initialize();
    return 0;
}
```
Prefer RAII containers (`std::vector`) over manual `new`/`delete`.

### Deprecated Relaxation API

| Deprecated | Replacement |
|------------|-------------|
| `RlxPre()`, `RlxMan()`, `RlxAuto()` | `rad.Solve(obj, prec, maxiter, method)` |
| `RlxUpdSrc()`, `SetRelaxSubInterval()` | `rad.Solve()` |

---

## Build & Release

### Build: MSVC + Intel MKL

**POLICY**: Use **MSVC** compiler with **Intel MKL**. Intel oneAPI compiler (icx-cl) is NOT compatible with NGSolve linking.

```powershell
powershell.exe -ExecutionPolicy Bypass -File "Build.ps1"
powershell.exe -ExecutionPolicy Bypass -File "Build.ps1" -Rebuild  # Clean rebuild
```

**Required Software**: Visual Studio 2022 (MSVC), Intel oneAPI Base Toolkit (MKL only, NOT the compiler).

### BLAS/LAPACK: Intel MKL Only

**POLICY**: OpenBLAS is NOT supported. MKL provides optimized BLAS/LAPACK. MKL internally uses Intel OpenMP (`libiomp5md.dll`) for its own threading, but Radia no longer links it directly.

**Required MKL DLLs** (loaded at runtime via pip dependency): `mkl_rt.*.dll`, `mkl_core.*.dll`, `mkl_intel_thread.*.dll`, `mkl_def.*.dll`, `mkl_avx2.*.dll`, `mkl_vml_*.dll`, `libiomp5md.dll` (MKL dependency), `libmmd.dll`, `svml_dispmd.dll`.

### Parallelization: NGSolve TaskManager

**POLICY**: Use **NGSolve TaskManager** for thread-level parallelization, NOT raw OpenMP parallel regions.

NGSolve's TaskManager provides work-stealing task-based parallelism that integrates with MKL and avoids nested OpenMP issues. All new parallel code in Radia should use TaskManager.

```cpp
// CORRECT: NGSolve TaskManager
#include <ngstd.hpp>
TaskManager::CreateJob([&](const TaskInfo& ti) {
    // work-stealing parallel loop
    for (size_t i = ti.task_nr; i < n; i += ti.ntasks) {
        // compute...
    }
});

// AVOID: Raw OpenMP parallel for (legacy code only)
#pragma omp parallel for
for (int i = 0; i < n; i++) { ... }
```

**When to use TaskManager**:
- Field computation loops (ComputeFieldBatch)
- Interaction matrix assembly
- Any embarrassingly parallel loop

**When OpenMP is acceptable**:
- MKL internal threading (controlled by `mkl_set_num_threads`)
- Legacy code not yet migrated

### PyPI Release Workflow

1. Claude Code bumps version in `pyproject.toml` and `src/radia/__init__.py`, updates `CHANGELOG.md`
2. **CRITICAL**: Run `Build.ps1` BEFORE `python -m build`
3. **MANDATORY**: Verify wheel contains correct `.pyd` before upload:
   ```python
   import zipfile
   whl = zipfile.ZipFile('dist/radia-X.Y.Z-py3-none-any.whl')
   for info in whl.infolist():
       if info.filename == 'radia/radia.pyd':
           print(f'Wheel .pyd: {info.file_size} bytes')  # Must match build output
   ```
4. **User manually executes** `Publish_to_PyPI.ps1` (Claude Code does NOT upload)
5. NEVER commit PyPI tokens to repository

**Build Path**: Only `build-msvc/` is supported. `build/Release/` is REMOVED.

**Testing Checklist** (before PyPI upload):
1. `python -m build`
2. Inspect wheel contents (must show `radia/radia.pyd`, NOT `python/radia.pyd`)
3. `pip install dist/radia-x.y.z-py3-none-any.whl`
4. `python -c "import radia; print(radia.__version__)"`

### MKL DLL Policy: Do NOT Bundle

**POLICY**: PyPI packages MUST NOT bundle Intel MKL DLLs. `pyproject.toml` declares `mkl>=2024.2.0` as dependency; pip installs MKL DLLs to `{sys.prefix}/Library/bin/`. `__init__.py` adds the path via `os.add_dll_directory()`.

**Do NOT**: Copy MKL/Intel OpenMP DLLs into `src/radia/` or include `*.dll` in `package_data`.

### Package Structure

```
src/radia/
  __init__.py           # DLL path setup + re-export from C++ module
  _radia_pybind.pyd     # Main C++ extension (includes RadiaField CoefficientFunction)
  cln_core.pyd          # CLN transient solver
  mmm_core.pyd          # MMM solver
  peec_matrices.pyd     # PEEC matrix assembly
  *.py                  # Python utility modules
  # NO .dll files
```

**Always use `Build.ps1`** for building. Never use manual cmake commands -- the script handles CMake configure + build + `.pyd` copy to `src/radia/`.

---

## Mesh & NGSolve Integration

### NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2601** or later. Version 6.2.2406~6.2.2501 had a regression in Periodic Boundary Conditions (`Identify()` lost during `Glue()`), fixed in 6.2.2601+.

Reference: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

```bash
pip install radia  # NGSolve is a required dependency (>=6.2.2601)
```

### NGSolve Recommended Configuration

```python
fes = HDiv(mesh, order=2)  # Best accuracy
B_gf = GridFunction(fes)
B_gf.Set(rad.RadiaField(radia_obj, 'b'))  # C++ CoefficientFunction in _radia_pybind.pyd
```

- Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- Use CoefficientFunction directly for maximum accuracy near boundaries
- Avoid GridFunction evaluation within 1 mesh cell of magnet surface

### NGSolve Magnetization → Radia Open Boundary Field Evaluation

NGSolve FEM solves M(x) inside bounded domains but struggles with open boundary (PML needed). Radia provides natural open boundary evaluation using **exact analytical formulas** (NOT dipole approximation).

```
NGSolve FEM Solve → M per element → netgen_mesh_to_radia() → Radia objects → rad.Fld()
```

**POLICY**: Do NOT use dipole approximation (m=M*V) for NGSolve → Radia pipeline. Register elements as proper Radia ObjHexahedron/ObjTetrahedron with solved magnetization. Radia's surface charge/surface current analytical formulas are exact for constant M per element, with no approximation error at any distance.

**Use cases**:
- External field from FEM-solved nonlinear iron core (no PML needed)
- Stray field evaluation at large distances (exact, not approximate)
- Particle trajectory through FEM-solved magnet assembly
- NGSolve CoefficientFunction for coupling back into FEM

**Workflow**:
```python
import radia as rad
from ngsolve import *
from radia.netgen_mesh_import import netgen_mesh_to_radia

rad.UtiDelAll()

# 1. NGSolve solves nonlinear problem → M per element
# (user's FEM solve code here)

# 2. Convert mesh to Radia objects with per-element magnetization
def material_from_ngsolve(el_idx):
    M = get_element_magnetization(gf_M, mesh, el_idx)  # user function
    return {'magnetization': M.tolist()}

container = netgen_mesh_to_radia(mesh, material=material_from_ngsolve, units='m')
# No Solve() needed - M is already known from NGSolve

# 3. Evaluate field at arbitrary external points (exact analytical formulas)
B = rad.Fld(container, 'b', [0, 0, 0.1])          # single point
B_batch = rad.FldBatch(container, obs_points, 0)   # batch
rad.FldVTS(container, 'field.vts', ...)             # VTK export
```

**Why Radia objects, not dipoles**:
- Surface charge model: exact for constant M, zero approximation error
- Near-field: no distance limitation (dipoles fail at r < 2*element_size)
- `netgen_mesh_to_radia()` already supports per-element material via callable

### NGSolve Mesh Access Policy

**POLICY**: All mesh access MUST use functions from `src/radia/netgen_mesh_import.py`. NEVER directly access `mesh.ngmesh.Points()` or `el.vertices[].nr` -- NGSolve has two indexing schemes (0-indexed vs 1-indexed) that cause off-by-one errors.

```python
# CORRECT
from netgen_mesh_import import netgen_mesh_to_radia, extract_elements
radia_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0,0,0]}, units='m')

# WRONG - index confusion
pt = mesh.ngmesh.Points()[v.nr]  # Off-by-one!
```

### Mesh Generation Policy

**2-Tier workflow by element order**:

| Element Order | Workflow | Tools |
|---------------|----------|-------|
| **1st + 2nd order** | GMSH direct | GMSH, or Cubit -> GMSH export |
| **3rd order+** | Netgen + Cubit via STEP | `mesh.Curve(order)` with OCC geometry (NGSolve FEM only) |

**Radia supports 1st order only** (Tet4, Hex8, Wedge6). 2nd order planned.

### Mesh Import Paths

```
Path A: GMSH direct (no NGSolve needed)
  CAD -> GMSH -> .msh -> gmsh_mesh_import.py -> Radia

Path B: GMSH via NGSolve
  CAD -> GMSH -> .msh -> NGSolve Mesh() -> netgen_mesh_import.py -> Radia

Path C: Cubit for complex hex
  CAD -> Cubit -> GMSH export -> Path A or B
```

**Key import functions**:

| Module | Function | Purpose |
|--------|----------|---------|
| `gmsh_mesh_import` | `gmsh_to_radia(file, mu_r, ...)` | GMSH .msh -> Radia (no NGSolve) |
| `netgen_mesh_import` | `netgen_mesh_to_radia(mesh, ...)` | NGSolve mesh -> Radia |
| `netgen_mesh_import` | `create_hex_mesh_grid(...)` | Structured hex (no external tool) |

### Coreform Cubit Mesh Export

For complex hexahedral meshes, use the **Coreform Cubit Mesh Export** tool.

**Repository**: https://github.com/ksugahar/Coreform_Cubit_Mesh_Export

```python
# Cubit -> GMSH -> Radia (recommended for complex hex)
from gmsh_mesh_import import gmsh_to_radia
core = gmsh_to_radia('cubit_exported.msh', mu_r=1000)
```

Cubit workflow for journal files: define blocks before export, use `cubit_mesh_export` utilities.

### PEEC Conductor Mesh

PEEC conductors use **surface mesh only** (SIBC handles skin effect). Use GMSH `generate(2)` for surface meshes. Supported: Tri3, Quad4 (1st order), Tri6, Quad8/9 (2nd order).

### Nastran Format: REMOVED

Nastran BDF support is **REMOVED**. Use Cubit -> Netgen direct export. Cubit can read legacy `.bdf` files if needed.

### Mesh Operations: Dropped APIs

`ObjDivMag`, `ObjDivMagPln`, `ObjCutMag` are NOT supported. All mesh operations use external tools (GMSH, Cubit).

### Mesh File Preservation

**NEVER DELETE** mesh files (`.bdf`, `.nas`, `.msh`, `.vtk`), Cubit journal files (`.jou`), or mesh generation scripts. These are difficult to recreate.

### Available Mesh Access Functions

From `src/radia/netgen_mesh_import.py`:
- `netgen_mesh_to_radia()` -- Convert entire mesh to Radia (recommended)
- `extract_elements()` -- Extract element data for custom processing
- `compute_element_centroid()` -- Centroid from vertex list
- `create_radia_tetrahedron()` / `create_radia_hexahedron()` -- Single elements
- `create_hex_mesh_grid()` -- Structured hex grid (no external tool)
- Constants: `TETRA_FACES`, `HEX_FACES`, `WEDGE_FACES`, `PYRAMID_FACES` (1-indexed face topology)

---

## H-Matrix Acceleration (HACApK)

### Policy: Use HACApK Only

**POLICY**: Do NOT implement custom H-matrix algorithms. Use the HACApK library at `src/ext/HACApK_LH-Cimplm/` (MIT license).

**Solver Methods**:

| Method | Name | Use Case |
|--------|------|----------|
| 0 | LU | Small problems (N < 500), guaranteed convergence |
| 1 | BiCGSTAB | General purpose, medium problems |
| 2 | HACApK | Large problems (N > 1000), O(N log N) memory |

**ソルバー選択ガイドライン**:
- **小規模 (N<500)**: LU推奨 (確実な収束)
- **中規模 (500<N<2000)**: BiCGSTAB推奨 (最速)
- **大規模 (N>2000)**: HACApK推奨 (メモリ効率)

### HACApK Parameters

```python
rad.SetHACApKParams(eps, leaf_size, eta)
```

| Parameter | C++ Default | Recommended | Description |
|-----------|-------------|-------------|-------------|
| eps | 1e-4 | 1e-4 | ACA tolerance (1e-6 to 1e-2) |
| leaf_size | 32 | **10** | Minimum cluster size (elements). 10 for MSC 6DOF hex (60 DOF/leaf), ELF-compatible |
| eta | 2.0 | 2.0 | Admissibility parameter |

**leaf_size=10 の根拠**: C++ デフォルトは 32 だが、MSC 6DOF ヘキサヘドロンでは leaf_size=10 → 実際のリーフは最大 11 要素 × 6DOF = **66 DOF/leaf** となる（二分木分割の端数）。C++ デフォルト 32 では 32×6=192 DOF/leaf で leaf が大きすぎ、ACA+ の低ランク近似の恩恵が小さくなる。leaf_size=10 (66 DOF) が圧縮効率と per-leaf 計算コストのバランスが良い。全ベンチマークで `rad.SetHACApKParams(1e-4, 10, 2.0)` を使用。

See `docs/HMATRIX_EVALUATION.md` for full evaluation report.

### Under-Relaxation for Nonlinear Problems

```python
rad.SetRelaxParam(0.3)  # 30% damping (0.0 = full step)
rad.Solve(container, 0.0001, 1000, 1)
rad.SetRelaxParam(0.0)  # Reset to full step
```

---

## IMA (Image Method of Analysis)

### IMA Sign Selection Policy

| Field vs Mirror Plane | IMA Sign |
|----------------------|----------|
| Field **parallel** to mirror | **+** (symmetric) |
| Field **perpendicular** to mirror | **-** (antisymmetric) |

```python
# Z-field, X-Z quarter model
rad.Solve(container, 0.0001, 100, 0, image='+x-z')  # Bz parallel to X-mirror, perp to Z-mirror

# X-field, X-Z quarter model
rad.Solve(container, 0.0001, 100, 0, image='-x+z')
```

### IMA Boundary Element Limitation

IMA produces incorrect results for **boundary elements** (faces ON symmetry plane) when observation points are **also on the symmetry plane** (~0.5x magnitude). Off-plane observation points work correctly (fixed 2026-02-04).

**Workaround**: Use explicit element duplication for models with boundary elements.

**When IMA is Safe**:
1. Non-boundary elements only (offset from symmetry planes)
2. Observation points off-plane
3. MMM (tetrahedra) -- no limitation

---

## PEEC & Conductor Solver

### Architecture Overview

**Approach**: PEEC (Partial Element Equivalent Circuit) with SIBC (Surface Impedance Boundary Condition) and ESIM (Effective Surface Impedance Method).

**Target**: Induction heating (1-500 kHz), WPT (6.78/13.56 MHz), power electronics (DC-1 MHz).

See `docs/` for detailed PEEC documentation.

### Filament-Panel Architecture (FastImp Style)

```
Surface Mesh -> Face -> Panel (Star: charge)  -> P matrix
             -> Edge -> Filament (Loop: current) -> L, R matrices
```

Loop-Star basis transformation is NOT needed -- filaments and panels are inherently separate in PEEC.

### PEEC System Equation

```
[R + jwL + Zs    jwM_LS  ] [I_filament]   [V]
[jwM_LS^T        P/(jw)  ] [Q_panel   ] = [0]
```

### Node-Segment Topology API

```python
from peec_matrices import PyPEECBuilder
from peec_topology import PEECCircuitSolver

builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n2)
topo = builder.build_topology()

solver = PEECCircuitSolver(topo)
Z = solver.compute_port_impedance(freq=1e6)
```

### Multi-Filament (nwinc/nhinc)

Use `nwinc`/`nhinc` parameters to subdivide conductor cross-sections for skin/proximity effect:
```python
builder.add_connected_segment(n1, n2, 3e-3, 3e-3, sigma=5.8e7, nwinc=3, nhinc=3)
```

Guidelines: DC=1x1, moderate skin (d/delta~2-5)=3x3, strong skin (d/delta>5)=5x5+.

### FastHenry .inp Parser

```python
from fasthenry_parser import FastHenryParser
parser = FastHenryParser()
parser.parse_file('inductor.inp')
result = parser.solve()
```

Supports: `.Units`, `N`/`E` definitions, `.external`, `.freq`, `.default`, `.equiv`, `.magnetic` blocks, line continuation `+`.

### Coupled PEEC + MMM

```python
from peec_coupled import CoupledPEECSolver
solver = CoupledPEECSolver(topology_dict, magnetic_objects=[core_id])
solver.compute_coupling_matrix()  # N_seg Radia Solve calls
Z = solver.compute_port_impedance(freq)
Z_sweep = solver.frequency_sweep(freqs)
L_total = solver.get_effective_inductance()  # L_air + Delta_L
```

For linear materials, `Delta_L` is frequency-independent (computed once).

**Physics**: `Z_eff(f) = diag(R + Zs(f)) + jw * (L_air + Delta_L)` where Delta_L comes from Biot-Savart -> `rad.ObjBckg()` + `rad.Solve()` -> vector potential A -> mutual inductance.

**FastHenry .magnetic Block** for coupled simulations:
```
.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  mu_r=1000
.endmagnetic
```

### SIBC Implementations

| Conductor Type | Method | Library |
|---------------|--------|---------|
| Circular | Bessel I0/I1 | `scipy.special.iv` |
| Rectangular (d << w) | Dowell formula | C++ rad_peec_surface_impedance.cpp |
| Nonlinear magnetic | ESIM cell problem | `esim_cell_problem.py` |

**Bessel**: Use `scipy.special.iv` (modified Bessel), NOT `jv` (regular Bessel). MKL does not provide Bessel functions.

### ESIM (Effective Surface Impedance Method)

ESIM solves 1D cell problem for H-dependent surface impedance: `d/dz[(1/mu(z)) * dH/dz] = jw*sigma*H`.

Supports complex permeability: `mu = mu' - j*mu"` for magnetic hysteresis/grain eddy current losses.

Use for: induction heating workpieces, nonlinear iron cores, lossy ferrite at high frequency.

Reference: `src/radia/esim_cell_problem.py`, `src/radia/esim_coupled_solver.py`.

### Deleted Legacy PEEC APIs

The following C++ APIs are **REMOVED**: `CndLoop`, `CndRecBlock`, `CndLoopFromHelix`, `CplMagCreate`, `CplMagSolve`, `CplMagSetFrequency`, `CndHexahedron`, `CndWire`, `CndSpiral`, `MatSIBC`. Use `PEECBuilder` and `CoupledPEECSolver` instead.

### PRIMA Model Order Reduction

**POLICY**: Use PRIMA (not CLN/Cauer) terminology. Both use Lanczos tridiagonalization; PRIMA (1998, IEEE TCAD) is the standard reference.

Key classes: `SPICEExtractionConfig`, `PRIMASchurExtractor`, `LoopStarMagneticCoupled` in `lanczos_reduction.py`.

### ngbem Integration

Radia PEEC works alongside NGSolve ngbem:

| Range | Solver |
|-------|--------|
| DC - 1 MHz | Radia PEEC + SIBC |
| DC - 1 MHz | ngbem (Weggler EFIE, low-freq stable) |
| 1 MHz - GHz | ngbem (Helmholtz) |

Radia PEEC unique features: direct circuit extraction (L, R, C), native SPICE netlist, Lanczos MOR, MMM coupling.

### Integration Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  NGSolve Geometry & Mesh                                        │
└──────────────────────────┬──────────────────────────────────────┘
                           │
           ┌───────────────┴───────────────┐
           ▼                               ▼
┌─────────────────────┐         ┌─────────────────────┐
│  Radia PEEC         │         │  ngbem              │
│  - Loop-Star        │         │  - EFIE/MFIE        │
│  - SIBC/ESIM        │  <--->  │  - H-matrix         │
│  - Lanczos MOR      │ coupling│  - Helmholtz/Laplace│
│  - SPICE output     │         │  - Low-freq Weggler │
└─────────────────────┘         └─────────────────────┘
           │                               │
           └───────────────┬───────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  Unified Solution (NGSolve GridFunction)                        │
└─────────────────────────────────────────────────────────────────┘
```

### PEEC Source Files

**C++ Core**:
- `src/core/rad_peec_matrices.h/cpp` -- PEECSegment, PEECPort, PEECMatrices, MutualInductance
- `src/core/rad_peec_surface_impedance.cpp` -- Dowell formula
- `src/lib/rad_peec_matrices_api.cpp` -- pybind11 bindings

**Python**:
- `src/radia/peec_topology.py` -- PEECCircuitSolver (MNA nodal admittance)
- `src/radia/peec_coupled.py` -- CoupledPEECSolver
- `src/radia/fasthenry_parser.py` -- FastHenry .inp parser
- `src/radia/esim_cell_problem.py` -- ESIM cell problem solver
- `src/radia/lanczos_reduction.py` -- PRIMA model order reduction

---

## Material Specification

### MatLin - Linear Materials

```python
mat = rad.MatLin(mu_r)                       # Isotropic (preferred)
mat = rad.MatLin([mu_r_par, mu_r_perp], [ex, ey, ez])  # Anisotropic
```

For isotropic materials, ALWAYS use single-argument form. MatLin is for soft magnetic materials only.

### MatSatIsoTab - Nonlinear (B-H Curve)

```python
BH_DATA = [[0.0, 0.0], [100.0, 0.1], [1000.0, 1.2], [50000.0, 2.0]]
mat = rad.MatSatIsoTab(BH_DATA)  # [[H(A/m), B(T)], ...]
```

### Permanent Magnets

For fixed magnetization PM, specify directly -- no `Solve()` needed:
```python
pm = rad.ObjHexahedron(vertices, [0, 0, 954930])  # M in A/m
B = rad.Fld(pm, 'b', [0, 0, 0.1])
```

Call `Solve()` only when soft iron is present alongside permanent magnets.

PM material classes (`MatMagFixed`, `MatMagLinear`, `MatMagCurve`) are available but currently all behave as fixed magnetization. Full demagnetization is planned.

See `docs/ELF_CONVENTIONS.md` for detailed unit system documentation.

### Permanent Magnet + Soft Iron Interaction

When combining PM with soft iron, use `Solve()`:
```python
pm = rad.ObjHexahedron(pm_vertices, [0, 0, 954930])  # Fixed PM
iron = rad.ObjHexahedron(iron_vertices, [0, 0, 0])    # Zero initial M
mat_iron = rad.MatLin(1000)
rad.MatApl(iron, mat_iron)
assembly = rad.ObjCnt([pm, iron])
result = rad.Solve(assembly, 0.0001, 1000, 0)  # LU solver
B = rad.Fld(assembly, 'b', [0, 0, 0.1])
```

---

## File & Naming Conventions

### Python Script Path Import

```python
# From examples/: use ../../src/radia
# From tests/: use ../src/radia
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
```

Import from `src/radia` package (not build directories).

### Script Naming Convention

Use **snake_case** with functional prefixes:

| Prefix | Purpose | Example |
|--------|---------|---------|
| `demo_` | Educational demonstration | `demo_batch_evaluation.py` |
| `benchmark_` | Performance measurement | `benchmark_solver_scaling.py` |
| `verify_` | Correctness verification | `verify_curl_A_equals_B.py` |
| `compare_` | Method comparison | `compare_radia_ngsolve.py` |
| (none) | Physical model name | `sphere_in_quadrupole.py` |

### VTK Export

All example scripts should export VTK files with the same basename as the script. Use `rad.FldVTS()` for field data export.

### Benchmark Policy

ベンチマーク実行ルール:
1. 1ケース毎に実行（並列実行しない、メモリ測定の正確性のため）
2. メモリ使用量を記録（`psutil` を使用、`tracemalloc` はC++メモリを追跡しない）
3. 結果JSONに含める: `n_elements`, `solver_method`, `t_solve`, `peak_memory_mb`, `iterations`

```python
import psutil, os
def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return mem.peak_wset / (1024 * 1024) if hasattr(mem, 'peak_wset') else mem.rss / (1024 * 1024)
```

---

## Visualization Policy

### Tool Selection

| Purpose | Tool |
|---------|------|
| Quick/interactive | **PyVista** (default) |
| Publication figures | **ParaView** |
| Geometry | Netgen OCC + NGSolve Draw() |
| Field data | `rad.FldVTS()` -> PyVista/ParaView |

**Do NOT** implement custom visualization in Radia C++ code.

**Removed APIs**: `rad.ObjDrwVTK()`, `exportGeometryToVTK()`, `radia_pyvista_viewer.py`.

### VTS Field Export

```python
rad.FldVTS(magnet, 'field_output.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)
```

### Visualization Stack

```
┌─────────────────────────────────────────────────────────────────┐
│                    Visualization Framework                       │
├─────────────────────────────────────────────────────────────────┤
│  Geometry (CAD)        │  Field Data          │  Interactive    │
│  ──────────────────────│─────────────────────│────────────────│
│  Netgen OCC shapes     │  rad.FldVTS()       │  NGSolve Draw() │
│  STEP import (Cubit)   │  PyVista meshes     │  Netgen GUI     │
│  ObjRecMag -> OCC      │  ParaView VTS/VTU   │  webgui         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Universal Relaxation Network (URN)

All URN examples, data, and scripts in `examples/Universal_Relaxation_Network/`.

**Policy**:
- Synthetic data MUST be clearly marked as synthetic
- Real-world datasets MUST include license and citation info
- All paper results reproducible from scripts in this directory

---

**Last Updated**: 2026-03-01
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
