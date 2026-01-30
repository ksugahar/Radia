# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and refactoring policies for the Radia project when working with Claude Code.

## Green's Function Policy: Laplace Kernel Only (MQS/Darwin)

### Helmholtz Kernel Removed (2026-01-09)

**CRITICAL**: Radia uses **Laplace kernel only** for all Green's function computations.

**Policy**:
- **Use Laplace kernel**: $G(r) = 1/(4\pi r)$ for all integral equation formulations
- **Helmholtz kernel REMOVED**: $G(r) = e^{-jkr}/(4\pi r)$ is NOT supported
- **Frequency regime**: MQS (Magneto-Quasi-Static) to Darwin approximation

**Rationale**:
1. **Target Applications**: MagLev, WPT, Induction Heating - all operate in quasi-static regime
2. **Validity**: Valid when wavelength >> problem size ($kL << 1$)
3. **Performance**: Laplace kernel enables efficient FMM/H-matrix acceleration
4. **Simplicity**: Single kernel reduces code complexity and potential bugs

**Affected Components**:
- `rad_green_fullwave.h/cpp` - All Green's functions use Laplace kernel
- `rad_conductor.cpp` - `GreenFunction()` returns $1/(4\pi r)$ for all formulations
- `rad_hacapk.cpp` - HACApK H-matrix uses Laplace kernel
- `rad_exafmm.h/cpp` - FMM uses Laplace kernel ($1/r^3$ for dipoles)

**Do NOT**:
- Add Helmholtz kernel ($e^{-jkr}/r$) to any Green's function
- Use wave number $k$ in field calculations (except for skin depth)
- Implement full-wave EFIE or MFIE formulations

**Skin Effect Handling**:
Skin depth is computed from frequency for SIBC, but field propagation uses quasi-static approximation:
```cpp
// Skin depth calculation (OK)
double delta = std::sqrt(2.0 / (omega * mu * sigma));

// Green's function (Laplace only)
double G = 1.0 / (4.0 * M_PI * r);  // NOT exp(-jkr) / (4*pi*r)
```

---

## Development Strategy: Complement NGSolve (2026-01-16)

### Radia Focuses on What NGSolve Cannot Do Well

**CRITICAL**: Radia's role is to **complement NGSolve**, not compete with it. Focus development on areas where NGSolve (FEM) is weak.

**Strategic Positioning**:

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

**Radia Core Competencies** (where NGSolve is weak):

| Capability | Why NGSolve Struggles | Radia's Approach |
|------------|----------------------|------------------|
| **Open boundaries** | Requires PML/ABC, adds DOFs | Natural with BEM |
| **Permanent magnets** | Needs volume mesh | Analytical (ObjRecMag) |
| **Thin conductors** | Mesh aspect ratio issues | PEEC (surface only) |
| **Circuit extraction** | Post-processing needed | Direct L,R,C,M output |
| **SPICE export** | Not supported | Verilog-A generation |
| **Model order reduction** | Manual implementation | PRIMA/Lanczos built-in |

### No Reinventing the Wheel

**Policy**: Use established libraries, do NOT implement from scratch.

| Component | Decision | Library |
|-----------|----------|---------|
| **PEEC Solver** | External | PAMELA |
| **H-matrix/ACA** | External | HACApK (integrated) |
| **BLAS/LAPACK** | External | Intel MKL |
| **FEM** | External | NGSolve |
| **MOR** | Python | scipy + numpy |

**Radia C++ Core** (maintain and enhance):
1. **MMM** - Magnetic Moment Method for permanent magnets and soft iron
2. **MSC** - Magnetic Surface Charge for hexahedra/tetrahedra
3. **Field computation** - B, H, A, Phi in unbounded domains
4. **NGSolve integration** - RadiaField CoefficientFunction

**Do NOT Implement**:
- FEM solvers (use NGSolve)
- General sparse solvers (use MKL/MUMPS)
- Full-wave BEM (use ngbem for high frequency)
- CAD geometry kernels (use OpenCASCADE via NGSolve)
- PEEC from scratch (use PAMELA)

**Rationale**:
1. NGSolve already excels at FEM - don't duplicate
2. Focus resources on unique value: BEM + circuit extraction
3. Integration > reinvention

---

## Build Policy: MSVC + Intel MKL

### Compiler Requirement (2025-12-27)

**CRITICAL**: Use **MSVC (Visual Studio C++ compiler)** with **Intel MKL** for building Radia.

**Policy**:
- **Use MSVC** for building both `radia.pyd` and `radia_ngsolve.pyd`
- **Use Intel MKL** for BLAS/LAPACK operations (mkl_rt.lib)
- **radia_ngsolve REQUIRES MSVC** due to ABI compatibility with MSVC-compiled NGSolve libraries
- Intel oneAPI compiler (icx-cl) is NOT compatible with NGSolve linking

**Build Command**:
```powershell
powershell.exe -ExecutionPolicy Bypass -File "BuildMSVC.ps1"
# Or for clean rebuild:
powershell.exe -ExecutionPolicy Bypass -File "BuildMSVC.ps1" -Rebuild
```

**Why MSVC instead of Intel Compiler**:
1. **NGSolve Compatibility**: NGSolve is compiled with MSVC; Intel compiler produces LLVM bitcode objects that are incompatible with MSVC-compiled libraries
2. **radia_ngsolve Linking**: The `add_ngsolve_python_module` CMake function requires MSVC ABI compatibility
3. **Intel MKL Still Used**: mkl_rt.lib works with both MSVC and Intel compilers, providing fast BLAS/LAPACK

**Rationale**:
1. **NGSolve Integration**: Required for RadiaField CoefficientFunction support
2. **MKL Performance**: Intel MKL provides optimized BLAS/LAPACK regardless of compiler
3. **Build Simplicity**: Single build produces both radia.pyd and radia_ngsolve.pyd

**Required Software**:
- Visual Studio 2022 (MSVC compiler)
- Intel oneAPI Base Toolkit (for Intel MKL only, NOT the compiler)

---

## BLAS/LAPACK Policy: Intel MKL Only

### OpenBLAS Dropped (2025-12-28)

**Policy**: Radia uses **Intel MKL only** for BLAS/LAPACK operations. **OpenBLAS is NOT supported**.

**Rationale**:
1. **Performance**: Intel MKL is faster than OpenBLAS on Intel CPUs
2. **OpenMP Integration**: Intel MKL uses Intel OpenMP (libiomp5md.dll) which integrates well with parallel Radia code
3. **Simplicity**: Single BLAS library reduces complexity and potential conflicts

**Removed Files** (2025-12-28):
- `src/ext/openblas/` - OpenBLAS headers, libraries, and binaries
- `src/radia/libopenblas.dll` - OpenBLAS runtime DLL

**Required Intel MKL DLLs** (auto-copied by BuildMSVC.ps1):
- `mkl_rt.*.dll` - MKL SDL (Single Dynamic Library) runtime
- `mkl_core.*.dll` - MKL core
- `mkl_intel_thread.*.dll` - MKL threading
- `mkl_def.*.dll`, `mkl_avx2.*.dll` - CPU kernels
- `mkl_vml_def.*.dll`, `mkl_vml_avx2.*.dll` - Vector math library
- `libiomp5md.dll` - Intel OpenMP runtime
- `libmmd.dll`, `svml_dispmd.dll` - Intel compiler runtime

**Note**: The DLL patterns use wildcards (e.g., `mkl_rt.*.dll`) for version-agnostic compatibility with future Intel oneAPI releases.

---

## OpenMP Policy: Intel OpenMP Only

### Use Intel OpenMP Instead of MSVC OpenMP (2025-12-28)

**Policy**: Radia uses **Intel OpenMP (libiomp5md.dll)** only. **MSVC OpenMP (vcomp140.dll) is NOT used**.

**Rationale**:
1. **Intel MKL Compatibility**: Intel MKL uses Intel OpenMP internally; mixing with MSVC OpenMP causes conflicts
2. **Performance**: Intel OpenMP provides better threading performance on Intel CPUs
3. **Consistency**: Single OpenMP runtime avoids DLL conflicts and undefined behavior

**CMake Configuration**:
- Use standard `/openmp` compiler flag for MSVC
- Link against `libiomp5md.lib` directly instead of CMake's `OpenMP::OpenMP_CXX`
- CMake detects Intel OpenMP and links it instead of MSVC's vcomp140.dll

**Required DLLs** (auto-copied by BuildMSVC.ps1):
- `libiomp5md.dll` - Intel OpenMP runtime (from Intel oneAPI compiler directory)

**Verification**:
After building, verify only Intel OpenMP is loaded:
```python
import psutil, os
process = psutil.Process(os.getpid())
for dll in process.memory_maps():
    if 'omp' in dll.path.lower():
        print(dll.path)
# Should show: libiomp5md.dll
# Should NOT show: vcomp140.dll
```

**Note**: If both `libiomp5md.dll` and `vcomp140.dll` are loaded, there is a configuration error that will cause OpenMP parallelization to malfunction.

---

## Radia Solver Methods: MMM and MSC

Radia supports two solver methods:

### MMM (Magnetic Moment Method) - Tetrahedra
- Used for **tetrahedral elements** (4 faces)
- **3 DOF per element**: Magnetization vector (Mx, My, Mz)
- Represents magnetic objects as distributions of magnetic dipoles

### MSC (Magnetic Surface Charge) - Hexahedra
- Used for **hexahedral elements** (6 faces)
- **6 DOF per element**: Surface charge density (sigma) per face
- Computes field from surface charges using solid angle integration
- Use `ObjHexahedron()` for hexahedral elements

**Note**: Radia does NOT use BEM (Boundary Element Method). The MSC method uses surface charges but differs from classical BEM.

### Mixed Hex/Tetra Element Support (2025-12-26)

Radia supports **mixed meshes** containing both hexahedral (6DOF) and tetrahedral (3DOF) elements.

**Solver Compatibility**:

| Solver | Mixed Elements | Notes |
|--------|----------------|-------|
| LU (Method 0) | **Supported** | Dense LU with variable DOF blocks |
| BiCGSTAB (Method 1) | **Supported** | Iterative solver with variable DOF |
| HACApK (Method 2) | NOT Supported | HACApK cluster tree requires uniform DOF |

**Note**: HACApK's cluster tree algorithm assumes fixed DOF per element. For mixed meshes, use LU (small problems) or BiCGSTAB (large problems).

**Interaction Matrix Blocks**:
- **3x3 block** (tetra-tetra): Standard demagnetization tensor
- **6x6 block** (hex-hex): Surface charge interaction (MSC)
- **3x6 block** (tetra from hex): H-field at tetra center from hex face charges
- **6x3 block** (hex from tetra): Normal dot N_matrix at hex eval points

**Usage**:
```python
import radia as rad

rad.FldUnits('m')

# Create mixed container with hex and tetra elements
hex_vertices = [[0,0,0], [0.1,0,0], [0.1,0.1,0], [0,0.1,0],
                [0,0,0.1], [0.1,0,0.1], [0.1,0.1,0.1], [0,0.1,0.1]]
hex_obj = rad.ObjHexahedron(hex_vertices, [0, 0, 0])   # 6DOF MSC
tetra_obj = rad.ObjTetrahedron(tetra_vertices, [0, 0, 0])  # 3DOF MMM

container = rad.ObjCnt([hex_obj, tetra_obj])
mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(container, mat)

# Solve with LU or BiCGSTAB (NOT HACApK)
rad.Solve(container, 0.001, 100, 0)  # Method 0 = LU
```

**Implementation Details**:
- `RADIA_MSC_SUPPORT` compile flag enables mixed element support
- Variable DOF offset arrays: `m_elemDOF`, `m_elemDOFOffset`, `m_totalDOF`
- Block interaction computation in `SetupInteractMatrix_VariableDOF()`

---

## Field Calculation Methods: Surface Current vs Surface Charge

### ObjRecMag - Surface Current Model (Rectangular Blocks Only)

`ObjRecMag` uses the **surface current approximation** for field calculations:

- **Applicable to**: Rectangular blocks (parallelepipeds) only
- **B/H field**: 8-corner analytical formula with arctangent and logarithm integrals
- **A field (vector potential)**: Uses unified BufVect formula: `A = (1/4π) * M × BufVect`
- **Phi field (scalar potential)**: Uses `Phi = (1/4π) * M · BufVect`

**Key advantage**: The 8-corner BufVect formula is computationally efficient and does NOT cancel on symmetry axes.

```python
import radia as rad
rad.FldUnits('m')

# Create rectangular permanent magnet
rec_mag = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 954930])

# All field types work correctly
B = rad.Fld(rec_mag, 'b', [0.05, 0, 0])  # Magnetic field
A = rad.Fld(rec_mag, 'a', [0, 0, 0.05])  # Vector potential (non-zero on axis)
```

### ObjHexahedron/ObjTetrahedron - Surface Charge Model (General Polyhedra)

`ObjHexahedron` and `ObjTetrahedron` use **surface charge integration**:

- **Applicable to**: Arbitrary hexahedra (6 quadrilateral faces) and tetrahedra (4 triangular faces)
- **B/H field**: Face-based solid angle integration: `H = (1/4π) * Σ σ_i * Ω_i`
- **A field (vector potential)**: Face-based integration: `A = (1/4π) * Σ (M × n_i) * I_i`

**Limitation**: On symmetry axes, face-based A integration gives A=0 due to symmetric cancellation. This is mathematically correct for the face-based formula but differs from ObjRecMag.

```python
import radia as rad
rad.FldUnits('m')

# Create hexahedral magnet (arbitrary shape)
vertices = [
    [-0.02, -0.02, -0.03], [0.02, -0.02, -0.03],
    [0.02, 0.02, -0.03], [-0.02, 0.02, -0.03],
    [-0.02, -0.02, 0.03], [0.02, -0.02, 0.03],
    [0.02, 0.02, 0.03], [-0.02, 0.02, 0.03],
]
hex_mag = rad.ObjHexahedron(vertices, [0, 0, 954930])

# B/H fields work correctly everywhere
B = rad.Fld(hex_mag, 'b', [0.05, 0, 0])

# A field: May be zero on symmetry axes (mathematical cancellation)
A = rad.Fld(hex_mag, 'a', [0, 0, 0.05])  # Could be ~0 on z-axis
```

### API Summary (2025-12-31)

**User-facing APIs** for creating magnetic elements:
- `ObjRecMag(center, dimensions, magnetization)` - Rectangular magnets (optimized formulas)
- `ObjHexahedron(vertices, magnetization)` - Arbitrary hexahedra (8 vertices, auto-generates faces)
- `ObjTetrahedron(vertices, magnetization)` - Tetrahedra (4 vertices, auto-generates faces)
- Mesh import functions (`netgen_mesh_to_radia`, `create_radia_from_nastran`) for complex geometries

---

## Background Field Policy (2025-12-31)

### Background Field API

**Policy**: Use `ObjBckg(callback)` for all background field applications.

**API**:
- `rad.ObjBckg(callback)` - Background field via Python callback function
  - Callback receives `[x, y, z]` in current units and returns `[Bx, By, Bz]` in Tesla

**Uniform Background Field**:
```python
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])  # 0.1 T in z-direction
```

**NOT Supported** (Do NOT use):
- Solve-time background field specification
- Legacy `ObjBckg([Bx, By, Bz])` array form - use `lambda p: [Bx, By, Bz]` instead

**Usage**:

```python
import radia as rad

rad.FldUnits('m')

# Create magnetic object
mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
mat = rad.MatLin(999)  # mu_r = 1000
rad.MatApl(mag_obj, mat)

# Uniform background field
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])  # 0.1 T in z-direction

# Quadrupole background field
def quadrupole_field(point):
    x, y, z = point
    G = 10.0  # T/m gradient
    return [G * y, G * x, 0]

bkg = rad.ObjBckg(quadrupole_field)

# Add background to container and solve
container = rad.ObjCnt([mag_obj, bkg])
rad.Solve(container, 0.0001, 1000, 1)
```

**Rationale**:
- Single unified API for all background field types
- Uniform fields expressed as `lambda p: [Bx, By, Bz]`
- Consistent design: all fields are callback-based

---

## Memory Management

### Exception Safety

All functions that allocate memory with `new` must follow this pattern:

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
- Set to `nullptr` immediately after ownership transfer
- Clean up in `catch(...)` block if pointer is still non-null

### RAII (Resource Acquisition Is Initialization)

Prefer RAII containers over manual memory management:

```cpp
// Good - RAII with std::vector
std::vector<radTPolygon> polygons;

// Avoid - Manual memory management
radTPolygon* polygons = new radTPolygon[n];  // Requires manual delete[]
```

---

## Unit System Policy

### Always Use Meters (SI Units)

**Policy**:
- **All examples** in `examples/` folder MUST use `rad.FldUnits('m')`
- **NGSolve integration** ALWAYS requires `rad.FldUnits('m')`

**Rationale**:
- Radia default: millimeters (mm)
- NGSolve default: meters (m)
- Without `rad.FldUnits('m')`, coordinates are off by 1000x

**Correct workflow**:
```python
import radia as rad
rad.FldUnits('m')  # REQUIRED for NGSolve integration

# Hexahedral magnet using ObjHexahedron (8 vertices, magnetization in A/m)
vertices = [[-0.02,-0.02,-0.03], [0.02,-0.02,-0.03], [0.02,0.02,-0.03], [-0.02,0.02,-0.03],
            [-0.02,-0.02,0.03], [0.02,-0.02,0.03], [0.02,0.02,0.03], [-0.02,0.02,0.03]]
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])  # meters, A/m
```

---

## Radia Field Computation Limitations

### rad.Fld() Accuracy Inside Magnets

**Important Limitation**: `rad.Fld()` does **NOT** accurately compute field values **inside** permanent magnets.

**Rationale**:
- Radia MMM is designed for field calculation in **air regions** (outside magnetic materials)
- Inside magnets, `rad.Fld()` returns inaccurate values (known limitation, not a bug)

**Testing Strategy**:
- X **Avoid**: Direct comparison of `rad.Fld()` inside magnets
- OK **Use**: Large magnet with small mesh region (field approximately uniform)

### Vector Potential A Field Implementation (2025-12-31)

**Status**: Vector potential A is now **IMPLEMENTED** for all ObjHexahedron/ObjTetrahedron elements.

**Implementation Details**:
- Uses **face integration** (not dipole approximation) for accurate results
- Formula: `A = (mu_0/4pi) * (M x BufVect)` with mm-to-m conversion factor
- Matches the analytical formula used in radTRecMag for rectangular blocks
- Extended to arbitrary triangular and quadrilateral faces using the Wilton et al. formula

**Usage**:
```python
import radia as rad
rad.FldUnits('m')

# Create hexahedral magnet
vertices = [[-0.02,-0.02,-0.03], [0.02,-0.02,-0.03], ...]  # 8 vertices
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])

# Get vector potential A at a point
A = rad.Fld(magnet, 'a', [0.03, 0.02, 0.05])  # Returns [Ax, Ay, Az] in T*m
```

**Maxwell Equation Consistency**:
Vector potential A satisfies `B = curl(A)` (verified numerically).

**Verification Script**:
- `examples/ngsolve_integration/verify_curl_A_equals_B/` - Verifies curl(A) = B

### Radia Internal Unit System: SI Meters (2025-12-31 Refactoring)

**POLICY**: Radia now uses **meters (m)** as the internal base length unit, matching ELF.

**Key Design Principles**:
1. **FldUnits controls ONLY input geometry** - Converts user geometry coordinates to meters at input time
2. **One-time conversion only** - After initial conversion, all internal processing uses meters
3. **B, H, and evaluation points are NEVER scaled** - Fixed to SI units (Tesla, A/m, meters)
4. **FldUnits functionality is LIMITED** - Only geometry input uses unit conversion

**How it Works**:
1. `rad.FldUnits('m')` - Uses meters directly (default, recommended)
2. `rad.FldUnits('mm')` - Converts user geometry coordinates from mm to m at input time
3. **Evaluation points always in meters** - `rad.Fld()` point coordinates are always in meters
4. **Field values are always in SI**: B in Tesla, H in A/m, A in T*m

**Physical Constants Used**:

| Constant | Value | Usage |
|----------|-------|-------|
| `MU_0_OVER_FOUR_PI` | `1.0e-7` H/m | Vector potential A, Biot-Savart |
| `INV_FOUR_PI` | `1/(4*pi)` | Scalar potential Phi, solid angle |
| `MU_0` | `4*pi*1e-7` H/m | B-H relations |

**Field Calculation Formulas (SI Units)**:

| Field | Formula | Constant |
|-------|---------|----------|
| B field | `B = -mu_0 * grad(Omega)` | Uses solid angle |
| H field | `H = -grad(Phi)` | Uses `INV_FOUR_PI` |
| Phi (scalar) | `Phi = (1/4pi) * M . BufVect` | `INV_FOUR_PI` |
| **A field** | `A = (mu_0/4pi) * M x BufVect` | `MU_0_OVER_FOUR_PI = 1e-7` |

**Comparison with ELF (Now Matching)**:

| Aspect | Radia (v1.4.3+) | ELF |
|--------|-----------------|-----|
| Internal units | **m (SI)** | m (SI) |
| A field constant | `mu_0/(4*pi) = 1e-7` | `mu_0/(4*pi) = 1e-7` |
| Unit conversion bugs | **None** | None |
| Maxwell equations | B = curl(A), H = -grad(Phi) | Same |

**Policy for C++ Developers**:
1. **All coordinates in meters** - no mm assumptions in calculations
2. **Use SI constants from rad_constants.h**: `MU_0_OVER_FOUR_PI`, `INV_FOUR_PI`
3. **Test with curl(A) = B** for vector potential implementations
4. **No hardcoded conversion factors** - all units are consistent

**Source Files with Field Constants**:
- `rad_constants.h`: Central location for physical constants
- `rad_rectangular_block.cpp`: ObjRecMag A/Phi/B field computation
- `rad_poly_analytical.cpp`: ObjHexahedron/ObjTetrahedron A/Phi field computation
- `rad_arc_current.cpp`: Biot-Savart field from arc currents

### Unified Field Computation Architecture (2026-01-16)

**POLICY**: All field computation MUST use `rad_field_unified.h/cpp` as the central module.

**Architecture**:
```
┌─────────────────────────────────────────────────────────────────┐
│                    rad_field_unified.h/cpp                       │
│  ─────────────────────────────────────────────────────────────  │
│  ComputeFieldSingle()     - Single point, static field          │
│  ComputeFieldBatch()      - Batch points, OpenMP parallelized   │
│  ComputeComplexFieldSingle() - Complex (AC) field               │
│  ComputeComplexFieldBatch()  - Complex batch with OpenMP        │
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
           │                  │                  │
           ▼                  ▼                  ▼
    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
    │ radentry.cpp│    │ radentry.cpp│    │ rad_peec_   │
    │ Python API  │    │ VTS output  │    │ mmm_coupled │
    └─────────────┘    └─────────────┘    └─────────────┘
```

**Users of rad_field_unified**:
| Component | Function Used | Purpose |
|-----------|--------------|---------|
| `rad.Fld()` | `ComputeFieldSingle` | Python API single-point field |
| `rad.FldBatch()` | `ComputeFieldBatch` | Python API batch field |
| `rad.FldVTS()` | `ComputeFieldBatch` | VTS export grid field |
| `radia_ngsolve` | `ComputeFieldSingle` | RadiaField CoefficientFunction |
| `rad_particle_trajectory` | `ComputeFieldForTrajectory` | Beam tracking |
| `CplMagFld()` | `ComputeComplexFieldSingle` | PEEC+MMM coupled field |

**Key Features**:
1. **Inside/Outside Classification**: Uses solid angle method for accurate determination
2. **OpenMP Parallelization**: Batch computations parallelized with `#pragma omp parallel for`
3. **Complex Field Support**: For PEEC+MMM AC analysis with complex magnetization
4. **FMM Acceleration**: Optional dipole approximation for large problems

**Inside/Outside Handling**:
```cpp
// Automatic inside/outside handling
RadFieldUnified::ComputeConfig config;
config.check_inside = true;           // Enable inside check
config.return_internal_field = true;  // Return M*mu0 for inside points

RadFieldUnified::FieldResult result =
    RadFieldUnified::ComputeFieldSingle(g3dPtr, point, FIELD_B, config);

if (result.status == RadFieldUnified::STATUS_INSIDE) {
    // Point is inside element result.element_id
    // result.Bx/By/Bz contains internal field (mu0 * M)
}
```

**Complex Field Computation (PEEC+MMM)**:
```cpp
// For coupled PEEC+MMM with complex magnetization
std::complex<double>* M_complex = ...;  // From coupled solver solution
int n_elements = ...;

RadFieldUnified::ComplexFieldResult result =
    RadFieldUnified::ComputeComplexFieldSingle(
        g3dPtr, point, M_complex, n_elements, config);

// result.Bx, By, Bz are std::complex<double>
```

---

## Material Specification

### MatLin - Linear Materials

`rad.MatLin()` defines **linear magnetic materials** (soft magnetic materials, NOT permanent magnets).

**IMPORTANT**: MatLin is for **linear materials only**. For permanent magnets, use `ObjHexahedron()` or `ObjRecMag()` with magnetization vector.

### API Forms (Industry Standard)

```python
# Form 1: Isotropic linear material
mat = rad.MatLin(mu_r)  # Relative permeability (mu_r >= 1)

# Form 2: Anisotropic linear material with easy axis
mat = rad.MatLin([mu_r_par, mu_r_perp], [ex, ey, ez])
```

**Parameters**:
- **mu_r**: Relative permeability (industry standard, mu_r >= 1)
- **[mu_r_par, mu_r_perp]**: Parallel and perpendicular relative permeabilities
- **[ex, ey, ez]**: Easy axis direction vector (does NOT need normalization)

**Important Notes**:
1. **Linear materials ONLY**: MatLin is for soft magnetic materials (iron, steel, mu-metal, etc.)
2. **Permanent magnets**: Do NOT use MatLin - define magnetization directly in `ObjHexahedron(vertices, [Mx,My,Mz])`
3. **Isotropic materials**: **ALWAYS prefer single-argument form `MatLin(mu_r)`** for isotropic materials.
4. **Easy axis**: For anisotropic materials, the easy axis vector must have significant magnitude (e.g., `[0, 0, 1]`)

**Example**:
```python
import radia as rad
rad.FldUnits('m')

# Soft iron cube (isotropic, mu_r=4000)
iron_vertices = [[0,0,0], [0.1,0,0], [0.1,0.1,0], [0,0.1,0],
                 [0,0,0.1], [0.1,0,0.1], [0.1,0.1,0.1], [0,0.1,0.1]]
cube = rad.ObjHexahedron(iron_vertices, [0, 0, 0])  # Zero magnetization
mat = rad.MatLin(4000)  # mu_r = 4000
rad.MatApl(cube, mat)

# Anisotropic material with easy axis in z-direction
iron_vertices2 = [[0.2,0,0], [0.3,0,0], [0.3,0.1,0], [0.2,0.1,0],
                  [0.2,0,0.1], [0.3,0,0.1], [0.3,0.1,0.1], [0.2,0.1,0.1]]
cube2 = rad.ObjHexahedron(iron_vertices2, [0, 0, 0])
mat2 = rad.MatLin([5001, 101], [0, 0, 1])  # Easy axis along z
rad.MatApl(cube2, mat2)
```

### MatSatIsoTab - Nonlinear Materials (B-H Curve)

`rad.MatSatIsoTab()` defines **nonlinear isotropic magnetic materials** using a B-H curve.

**API (Industry Standard)**:

```python
# B-H curve: [[H1, B1], [H2, B2], ...] where H is in A/m and B is in Tesla
mat = rad.MatSatIsoTab(BH_data)
```

**Example**:
```python
# Steel B-H curve data
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [500.0, 0.8],
    [1000.0, 1.2],
    [5000.0, 1.7],
    [50000.0, 2.0],
]

mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(steel_obj, mat)
```

**Note**: Radia internally converts B-H to M-H using: M = B/mu_0 - H

### Permanent Magnet Materials

Radia provides several APIs for permanent magnet materials:

#### Method 1: ObjHexahedron with Magnetization (Recommended for Fixed PM)

For permanent magnets where demagnetization is negligible, specify magnetization directly in `ObjHexahedron`:

```python
import radia as rad
rad.FldUnits('m')

# Define hexahedral vertices (8 corners)
vertices = [
    [-0.05, -0.05, -0.05],
    [0.05, -0.05, -0.05],
    [0.05, 0.05, -0.05],
    [-0.05, 0.05, -0.05],
    [-0.05, -0.05, 0.05],
    [0.05, -0.05, 0.05],
    [0.05, 0.05, 0.05],
    [-0.05, 0.05, 0.05],
]

# NdFeB magnet: Br = 1.2 T = 954930 A/m (= Br / mu_0)
Mr = 954930  # A/m
pm = rad.ObjHexahedron(vertices, [0, 0, Mr])

# Compute field (NO Solve needed for fixed PM)
B = rad.Fld(pm, 'b', [0, 0, 0.1])  # Field at z=0.1m
```

**Key Points**:
- Use `ObjHexahedron` for arbitrary hexahedral shapes (8 vertices)
- Magnetization is specified in A/m: Mr = Br / mu_0
- **No `Solve()` required** for fixed magnetization permanent magnets
- Only call `Solve()` when soft iron is present in the model

#### Method 2: PM Material Classes (Future Demagnetization Support)

Three material classes are available for future demagnetization implementation:

```python
# MatMagFixed - Fixed magnetization (no demagnetization)
# Magnetization [Mx, My, Mz] in A/m
mat = rad.MatMagFixed([0, 0, 954930])

# MatMagLinear - Linear demagnetization (Br/Hc model)
# Br [T], Hc [A/m], easy axis [ex, ey, ez]
mat = rad.MatMagLinear(1.2, 955000, [0, 0, 1])

# MatMagCurve - User-defined B-H demagnetization curve
# [[H1,B1], [H2,B2], ...], easy axis [ex, ey, ez]
BH_curve = [[0.0, 1.2], [-500000, 0.6], [-955000, 0.0]]
mat = rad.MatMagCurve(BH_curve, [0, 0, 1])
```

**Note**: Currently all three material classes behave as fixed magnetization. Full demagnetization implementation is planned for future versions.

#### Permanent Magnet + Soft Iron Interaction

When combining permanent magnets with soft iron, use `Solve()` to calculate the induced magnetization in the iron:

```python
import radia as rad
rad.FldUnits('m')

# PM magnet (with fixed magnetization)
pm_vertices = [...]  # 8 vertex coordinates
pm = rad.ObjHexahedron(pm_vertices, [0, 0, 954930])

# Soft iron yoke (zero initial magnetization)
iron_vertices = [...]  # 8 vertex coordinates
iron = rad.ObjHexahedron(iron_vertices, [0, 0, 0])
mat_iron = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(iron, mat_iron)

# Create assembly and solve
assembly = rad.ObjCnt([pm, iron])
result = rad.Solve(assembly, 0.0001, 1000, 0)  # LU solver

# Now compute fields
B = rad.Fld(assembly, 'b', [0, 0, 0.1])
```

---

## Windows Console Encoding (cp932) Compatibility

**Policy**: **NEVER use Unicode mathematical symbols** in print statements.

**Forbidden Unicode → ASCII Replacements**:

| Unicode | Symbol | ASCII | Example |
|---------|--------|-------|---------|
| `\u00b2` | ² | `^2` | `N²` → `N^2` |
| `\u00b3` | ³ | `^3` | `N³` → `N^3` |
| `\u2192` | → | `->` | `A → B` → `A -> B` |
| `\u2248` | ≈ | `~=` | `x ≈ 2` → `x ~= 2` |
| `\u2264` | ≤ | `<=` | `N ≤ 100` → `N <= 100` |
| `\u2265` | ≥ | `>=` | `N ≥ 250` → `N >= 250` |

**Rationale**: Windows console (cmd.exe) defaults to cp932 encoding in Japanese environments, causing `UnicodeEncodeError` for Unicode symbols.

---

## Python Script Path Import Policy

**Policy**: Import from `src/radia` package directory (not build directories).

```python
# ✓ CORRECT - Import from src/radia package
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad

# ✗ WRONG - Import from build directory (may be outdated)
sys.path.insert(0, r"S:\Radia\01_GitHub\build-msvc")
import radia as rad
```

**Path patterns**:
- Examples folder: `'../../src/radia'`
- Tests folder: `'../src/radia'`

**Rationale**:
- `src/radia/` contains the latest .pyd files copied by BuildMSVC.ps1
- Build directories may contain outdated or version-tagged .pyd files
- `build/Release/` is REMOVED - do not use

---

## File Organization Policies

### Mesh File Preservation

**Policy**:
- **NEVER DELETE** mesh files (`.bdf`, `.nas`, `.msh`, `.vtk`)
- **NEVER DELETE** Cubit journal files (`.jou`, `.journal`)
- **NEVER DELETE** mesh generation scripts

**Rationale**: Mesh files are difficult to recreate without original CAD or mesh generation tools.

### Cubit Mesh Generation

**Policy**: Use `cubit_mesh_export` utilities for Cubit-based workflows.

**Requirements**:
- Journal files (`.jou`) MUST define blocks before export
- Use `cubit_mesh_export` for NASTRAN format conversion

**Correct workflow**:
```python
# In Cubit journal file (.jou):
# 1. Create geometry
# 2. Generate mesh
# 3. Define blocks (REQUIRED)
block 1 volume 1
block 1 element type hex8

# 4. Export using cubit_mesh_export
export nastran "geometry.bdf" dimension 3 overwrite

# In Python script:
from cubit_mesh_export import export_cubit_mesh
export_cubit_mesh('geometry.bdf', blocks={'1': {'material': 'NdFeB', 'Mr': [0, 0, 1.2]}})
```

**Common mistakes**:
```python
# WRONG - No block definition in .jou file
export nastran "geometry.bdf"  # Missing block assignment!

# WRONG - Block defined after export
export nastran "geometry.bdf"
block 1 volume 1  # Too late!
```

### VTK Export Policy

All example scripts should export VTK files with the same basename as the script.

---

## H-Matrix Acceleration Policy (HACApK ACA+)

### Policy: Do NOT Implement Custom H-Matrix - Use HACApK Only

**Date**: 2025-12-19
**Status**: Development Policy

**CRITICAL POLICY**:
1. **Do NOT implement custom H-matrix algorithms** (ACA, ACA+, or any low-rank approximation)
2. **If H-matrix acceleration is needed**, use the HACApK library at `src/ext/HACApK_LH-Cimplm/`
3. **Custom implementations have been removed** - rad_hmatrix*.cpp/h files were deleted (2025-12-18)

**Rationale**:
- Custom ACA implementations are prone to bugs and difficult to validate
- HACApK is a proven, MIT-licensed library with established correctness
- Benchmarks showed NO speedup for typical Radia use cases (single compact objects)

**Reference Implementation (for future use)**:
- **Library**: HACApK (Hierarchical Approximation with ACA+ for Krylov methods)
- **Source**: `src/ext/HACApK_LH-Cimplm/` (C implementation)
- **Original**: ppOpen-HPC project (MIT License)
- **Key Files**:
  - `cHACApK_base.h` - Core H-matrix structures and cluster tree
  - `cHACApK_base.c` - ACA+ implementation
  - `mpi_stub.h` - MPI stub for single-process execution (created for Radia)

**Why HACApK ACA+** (if needed in future):
1. **Proven Algorithm**: ACA+ is a well-established low-rank approximation method
2. **MIT License**: Compatible with Radia's license
3. **C Implementation**: Direct integration possible (no Fortran dependency)
4. **MPI Stub Available**: Can run without MPI dependency using provided stub

**Current Status (2025-12-23)**:

HACApK H-matrix solver is **implemented and available** as Method 2.

**When to Use HACApK**:
- Large-scale problems (N > 1,000 elements, >6,000 DOF)
- Memory-constrained environments (H-matrix uses O(N log N) memory vs O(N^2) for dense)
- NOT recommended for small problems where LU is faster

**Current Solver Methods** (v1.3.15):

| Method | Name | Description | Use Case |
|--------|------|-------------|----------|
| 0 | LU | Dense LU decomposition with LAPACK dgesv | Small problems (N < 500), guaranteed convergence |
| 1 | BiCGSTAB | Iterative BiCGSTAB with Jacobi preconditioner | General purpose, medium problems |
| 2 | HACApK | BiCGSTAB with H-matrix acceleration (ACA+) | Large problems (N > 1000), memory efficiency |

**Note**: Original Radia used Implicit SS (Gauss-Seidel) which had slow convergence for nonlinear materials. This was replaced with BiCGSTAB in v1.3.13. HACApK was added in v1.3.15.

### Deprecated Relaxation API (2026-01-14)

The following legacy relaxation functions are **DEPRECATED** and will emit warnings:

| Function | Status | Replacement |
|----------|--------|-------------|
| `RlxPre()` | Deprecated | Use `Solve()` - handles matrix construction internally |
| `RlxMan()` | Deprecated | Use `Solve(obj, prec, maxiter, method)` |
| `RlxAuto()` | Deprecated | Use `Solve(obj, prec, maxiter, method)` |
| `RlxUpdSrc()` | Deprecated | Use `Solve()` for re-solving |
| `SetRelaxSubInterval()` | Deprecated | Use `Solve(obj, prec, maxiter, 0)` for LU |

**Migration Example**:

```python
# OLD (deprecated)
intrc = rad.RlxPre(container, container)
rad.RlxMan(intrc, 5, 1, 1.0)  # Method 5 = LU

# NEW (recommended)
rad.Solve(container, 0.0001, 1000, 0)  # Method 0 = LU
```

### HACApK Parameters

Configure HACApK parameters using `rad.SetHACApKParams(eps, leaf_size, eta)` or `rad.SetHMatrixEpsilon(eps)`:

| Function | Description |
|----------|-------------|
| `SetHACApKParams(eps, leaf_size, eta)` | Set all H-matrix parameters |
| `SetHMatrixEpsilon(eps)` | Set only ACA tolerance (ELF-compatible: `magic.set_hmatrix_epsilon`) |
| `GetHACApKStats()` | Get H-matrix statistics after solve |

| Parameter | Description | Default | Range |
|-----------|-------------|---------|-------|
| **eps** | ACA tolerance (lower = more accurate, higher = faster) | 1e-4 | 1e-6 to 1e-2 |
| **leaf_size** | Minimum cluster size for leaf nodes | 10 | 5 to 50 |
| **eta** | Admissibility parameter (higher = more low-rank blocks) | 2.0 | 1.0 to 3.0 |

**Parameter Rationale**:

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `eps` | 1e-4 | Balance between accuracy and compression. ELF-compatible default. |
| `leaf_size` | 10 | Minimum cluster size. Smaller values allow deeper tree but increase overhead. 10 is ELF-compatible and provides good balance for typical problems. Values < 5 rarely improve performance; values > 50 reduce H-matrix benefit. |
| `eta` | 2.0 | Standard admissibility: clusters are "well-separated" when `dist(c1,c2) >= eta * max(diam(c1), diam(c2))`. eta=2.0 is conservative and ELF-compatible. Lower values (1.0-1.5) allow more aggressive compression but may reduce accuracy. |

**ACA Tolerance (eps) Guidelines**:

| eps | Accuracy | Speed | Use Case |
|-----|----------|-------|----------|
| 1e-6 | Highest | Slowest | Validation, reference solutions |
| 1e-4 | Good | Good | **Default, ELF-compatible** |
| 1e-3 | Moderate | Fast | Quick iterations, exploration |
| 1e-2 | Lower | Fastest | Initial guess, previews |

**Recommended Usage**:
```python
import radia as rad

# Method 0: Dense LU (direct solver, default)
rad.Solve(container, 0.0001, 1000, 0)

# Method 1: BiCGSTAB (iterative solver, better for large problems)
rad.Solve(container, 0.0001, 1000, 1)

# Method 2: HACApK (H-matrix accelerated BiCGSTAB)
rad.SetHACApKParams(1e-4, 10, 2.0)  # eps=1e-4 (ELF-compatible)
rad.Solve(container, 0.0001, 1000, 2)

# ELF-compatible: Set only epsilon (other params unchanged)
rad.SetHMatrixEpsilon(1e-3)  # magic.set_hmatrix_epsilon(1e-3)
rad.Solve(container, 0.0001, 1000, 2)

# Under-relaxation for difficult nonlinear problems
rad.SetRelaxParam(0.3)  # 30% damping (0.0 = full step, 0.0-1.0)
rad.Solve(container, 0.0001, 1000, 1)
rad.SetRelaxParam(0.0)  # Reset to full step
```

**Benchmark Script Usage**:
```bash
# Run HACApK benchmark with default eps=1e-4
python benchmark_hexahedron_msc.py --hacapk 5 10 15

# Run with custom ACA tolerance
python benchmark_hexahedron_msc.py --hacapk --eps 1e-3 5 10 15
```

**Documentation**:
- `docs/HMATRIX_EVALUATION.md` - Full evaluation report
- `examples/cube_uniform_field/nonlinear/` - Benchmark scripts and results

### Performance Comparison with ELF (2025-12-30)

**All solvers now match or exceed ELF performance** (Hex N=10, 1000 elements):

| Solver | ELF | Radia | Radia/ELF |
|--------|-----|-------|-----------|
| LU | 14.0s | 14.1s | 1.01x (同等) |
| BiCGSTAB | 5.1s | **3.2s** | **0.62x (Radia高速)** |
| HACApK | 3.0s | **3.4s** | **1.13x (ほぼ同等)** |

**最適化履歴**:
- 2025-12-29: HACApK 10.2s (3.4x遅い) - 線形ソルブボトルネック特定
- 2025-12-30: UpdateDiagonal fast method有効化 → 3.4s (**3.0x高速化**)

**ボトルネック修正**:
1. ✅ `UpdateDiagonal`: fast method有効化 (O(N) vs O(N^2))
2. ✅ `HACApK_update_diagonal_wrapper`: OpenMP並列化
3. ✅ ACA+ fill: 既にOpenMP並列化済み

**ソルバー選択ガイドライン**:
- **小規模 (N<500)**: LU推奨 (確実な収束)
- **中規模 (500<N<2000)**: BiCGSTAB推奨 (最速)
- **大規模 (N>2000)**: HACApK推奨 (メモリ効率、O(N log N))

### Known Performance Issue: Tetrahedron HACApK (2025-12-30)

**問題**: TetrahedronメッシュでHACApKがELFより4-7倍遅い

| メッシュ | ELF (s) | Radia (s) | Radia/ELF |
|---------|---------|-----------|-----------|
| maxh=0.10 (4994要素) | 19.5 | 80.3 | **4.1x遅い** |
| maxh=0.15 (2211要素) | 3.8 | 25.6 | **6.7x遅い** |

**ボトルネック分析** (maxh=0.10):

| 項目 | ELF (s) | Radia (s) | Radia/ELF |
|------|---------|-----------|-----------|
| H-matrix構築 | 6.1 | 68.9 | **11.3x遅い** |
| 線形ソルブ | 13.4 | 11.2 | 0.84x (Radia高速) |

**原因**:
- H-matrix構築 (ACA+ fill) がTetrahedronで特に遅い
- ELFはFortran + OpenMPで高度に最適化
- Radiaの三角形積分計算がTetraで重い

**TODO (次の改善項目)**:
1. [ ] TetrahedronのACA+ fill並列化の最適化
2. [ ] 三角形積分のキャッシュ効率改善
3. [ ] `radTInteraction::B_comp_tetrahedron_MSC()` のプロファイリング
4. [ ] ELFの実装との比較調査

**Note**: HexahedronではRadiaはELFと同等以上の性能を達成済み。

---

## NGSolve Integration Best Practices

### NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2405** only.

**Version Constraint**: `ngsolve==6.2.2405`

**Reason**: NGSolve 6.2.2406+ has a regression bug in Periodic Boundary Conditions.
The `Identify()` information is lost during mesh generation with `Glue()`.

**Reference**: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

**Installation**:
```bash
pip install radia[ngsolve]  # Installs with correct NGSolve version constraint
# OR
pip install ngsolve==6.2.2405
```

### Recommended Configuration

```python
fes = HDiv(mesh, order=2)  # Best accuracy
B_gf = GridFunction(fes)
B_gf.Set(radia_ngsolve.RadiaField(radia_obj, 'b'))
```

**Evaluation guidelines**:
- Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- Use CoefficientFunction directly for maximum accuracy near boundaries
- Avoid GridFunction evaluation within 1 mesh cell of magnet surface

---

## PyPI Package Release Policy

### Version Management (Automated by Claude Code)

Claude Code is responsible for:
- Maintaining version numbers in `pyproject.toml` and `src/radia/__init__.py`
- Following semantic versioning (MAJOR.MINOR.PATCH)
- Updating `CHANGELOG.md` with release notes
- Committing and pushing version bump changes

### PyPI Upload (Manual by User)

**IMPORTANT**: Claude Code does NOT execute PyPI upload. The user must manually run the upload script.

**Workflow**:
1. Claude Code prepares the release (version bump, CHANGELOG, commit, push)
2. **CRITICAL: Run `BuildMSVC.ps1` BEFORE `python -m build`** to ensure latest .pyd files
3. Claude Code runs `python -m build` to create wheel package
4. **CRITICAL: Verify wheel contains correct .pyd** (see verification below)
5. **User manually executes** `Publish_to_PyPI.ps1` with their API token

### Wheel Verification (MANDATORY before PyPI upload)

**ALWAYS verify the wheel contains the latest .pyd before upload!**

```python
# Verification script - run BEFORE PyPI upload
import zipfile
import os
from datetime import datetime

whl_path = 'dist/radia-X.Y.Z-py3-none-any.whl'  # Update version
whl = zipfile.ZipFile(whl_path)

# Get wheel .pyd info
for info in whl.infolist():
    if info.filename == 'radia/radia.pyd':
        print(f'Wheel radia.pyd: {info.file_size} bytes, {info.date_time}')
        break

# Compare with build-msvc .pyd
msvc_pyd = 'build-msvc/radia.cp312-win_amd64.pyd'
if os.path.exists(msvc_pyd):
    stat = os.stat(msvc_pyd)
    print(f'Build radia.pyd: {stat.st_size} bytes, {datetime.fromtimestamp(stat.st_mtime)}')

# SIZES MUST MATCH! If wheel is smaller, it's using old .pyd
```

**Common mistake**: `python -m build` uses sdist which may contain old .pyd.
Solution: Always run `BuildMSVC.ps1` immediately before `python -m build`.

**Build Path Priority** (setup.py):
1. `build-msvc/` (MSVC + Intel MKL - ONLY supported path)

**Note**: `build/Release/` is REMOVED. Only `build-msvc/` is supported.

```powershell
# Set PyPI API token (keep secure!)
$env:PYPI_TOKEN = "pypi-AgEIcGl..."

# Run upload script
powershell.exe -ExecutionPolicy Bypass -File Publish_to_PyPI.ps1
```

**Security**: NEVER commit PyPI tokens to repository. NEVER ask user for their token.

### MKL DLL Bundling Policy (2025-12-27)

**Policy**: PyPI packages MUST include Intel MKL runtime DLLs for user convenience.

**Required DLLs** (Windows):
- `mkl_core.2.dll` - MKL core library
- `mkl_intel_thread.2.dll` - MKL Intel OpenMP threading
- `mkl_def.2.dll` - MKL definitions
- `mkl_rt.2.dll` - MKL runtime (optional, for dynamic linking)
- `libiomp5md.dll` - Intel OpenMP runtime

**Source Location** (Intel oneAPI):
```
C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\
C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\
```

**Package Structure**:
```
src/radia/
  __init__.py
  radia.pyd           # C++ extension
  mkl_core.2.dll      # MKL runtime (bundled)
  mkl_intel_thread.2.dll
  mkl_def.2.dll
  libiomp5md.dll      # Intel OpenMP
```

**Rationale**:
- Users should NOT need to install Intel oneAPI to use Radia
- Bundling DLLs ensures consistent behavior across environments
- Reduces support burden from "DLL not found" errors

**Build Integration**:
- `BuildWithIntel.ps1` should copy required DLLs to `src/radia/`
- `setup.py` / `pyproject.toml` should include DLLs in `package_data`

**License Consideration**:
Intel MKL runtime DLLs are redistributable under Intel oneAPI EULA.
Include appropriate license notices in the package.

---

## Nastran Format Policy: REMOVED (2026-01-16)

### Nastran BDF Support Removed from Radia

**CRITICAL**: Nastran BDF format support is **REMOVED** from Radia. Use **Coreform Cubit → Netgen direct export** exclusively.

**Policy**:
- **Nastran import modules REMOVED** from Radia codebase
- **All Nastran workflows REMOVED** - no backwards compatibility
- **Use Coreform Cubit** for all mesh operations (hex/tet)
- **Cubit can read legacy .bdf** files if needed

**Rationale**:
1. **Legacy format**: Nastran BDF is a decades-old format with limitations
2. **Complexity**: Fixed-width fields, multiple continuation styles, error-prone parsing
3. **Better alternatives**: Cubit `export_netgen()` provides direct in-memory mesh transfer
4. **No users yet**: No backwards compatibility concerns

**Removed Files**:
- `src/radia/nastran_mesh_import.py` - REMOVED
- `src/radia/nastran_reader.py` - REMOVED (already deprecated)
- All examples using `import_nastran_mesh()` - REMOVED or refactored

**Correct Workflow** (2026-01-16):
```
Cubit geometry → export_netgen() → Netgen mesh → Radia
                 (direct, no file)
```

**For Legacy .bdf Files**:
```
Legacy .bdf → Cubit (import nastran) → export_netgen() → Netgen → Radia
              (Cubit handles .bdf reading)
```

### Cubit → Netgen Workflow

**Standard workflow** (recommended):
```python
# RECOMMENDED - Use Cubit direct export
import cubit
import cubit_mesh_export
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

# Generate mesh in Cubit, export directly to Netgen
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import geometry 'model.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Direct export to Netgen (no intermediate file)
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# Convert to Radia
mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

### Reading Legacy Nastran Files via Cubit

**If you have existing Nastran .bdf files**, use Coreform Cubit to read them and add boundary condition labels before exporting to Netgen:

```python
# RECOMMENDED - Read legacy Nastran via Cubit, add labels, export to Netgen
import cubit
import cubit_mesh_export
from ngsolve import Mesh

cubit.init(['cubit', '-nojournal', '-batch'])

# Read legacy Nastran mesh into Cubit
cubit.cmd("import nastran 'legacy_mesh.bdf'")

# Add boundary condition labels (sidesets/nodesets) in Cubit
cubit.cmd("sideset 1 surface 1 2 3")
cubit.cmd("sideset 1 name 'conductor'")
cubit.cmd("sideset 2 surface 4 5 6")
cubit.cmd("sideset 2 name 'ferrite'")
cubit.cmd("sideset 3 surface 7 8")
cubit.cmd("sideset 3 name 'shield'")

# Export to Netgen with boundary labels preserved
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# Boundary labels are now accessible as mesh.GetBoundaries()
print(mesh.GetBoundaries())  # ['conductor', 'ferrite', 'shield']
```

**Benefits of Cubit as intermediary**:
1. **Label management**: Add/modify boundary condition names
2. **Mesh repair**: Fix mesh quality issues
3. **Visualization**: Inspect mesh before export
4. **Format flexibility**: Read various legacy formats (Nastran, ANSYS, Abaqus)
5. **CAD format support**: Import STEP, IGES, SAT, and other CAD formats directly
6. **Modern workflow**: No custom parsers needed in Radia
```

### CAD Format Support via Cubit

Coreform Cubit supports direct import of CAD formats, eliminating need for custom importers:

**Supported CAD Formats**:
- STEP (.step, .stp) - ISO standard, recommended
- IGES (.iges, .igs) - Legacy CAD exchange
- Parasolid (.x_t, .x_b) - Siemens/NX native
- ACIS SAT (.sat) - Spatial native
- STL (.stl) - Triangulated surface
- BREP (.brep) - OpenCascade native

**Workflow with CAD files**:
```python
import cubit
import cubit_mesh_export
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

cubit.init(['cubit', '-nojournal', '-batch'])

# Import CAD file (STEP, IGES, etc.)
cubit.cmd("import step 'motor_rotor.step' heal")

# Set up mesh
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("volume all size auto factor 5")
cubit.cmd("mesh volume all")

# Add boundary condition labels
cubit.cmd("sideset 1 surface with z_coord > 0")
cubit.cmd("sideset 1 name 'north_pole'")

# Export directly to Netgen
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# Convert to Radia
mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

---

## Unit System Policy: No Hard-Coded Unit Conversions

### Requirement: Centralized Unit Control via rad.FldUnits() and radia_ngsolve

**Goal**: All unit conversions must be controlled through explicit API calls (`rad.FldUnits()` or `radia_ngsolve` constructor), never through hard-coded conversion factors in user code.

**Policy**:

**✓ ALLOWED - Explicit unit control**:
```python
# Method 1: Set Radia units globally
import radia as rad
rad.FldUnits('m')  # All Radia operations now use meters

# Create hexahedral magnet with ObjHexahedron
vertices = [[-0.05,-0.05,-0.05], [0.05,-0.05,-0.05], [0.05,0.05,-0.05], [-0.05,0.05,-0.05],
            [-0.05,-0.05,0.05], [0.05,-0.05,0.05], [0.05,0.05,0.05], [-0.05,0.05,0.05]]
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])  # 0.1m cube, 1.2T equivalent

# Method 2: Specify units in radia_ngsolve constructor
from radia_ngsolve import RadiaField
B_cf = RadiaField(magnet, 'b', units='m')  # Explicitly use meters
```

**✗ FORBIDDEN - Hard-coded unit conversions**:
```python
# WRONG - Hard-coded mm to m conversion
for pt in obs_points:
    f.write(f'{pt[0]/1000.0} {pt[1]/1000.0} {pt[2]/1000.0}\n')  # ✗ DO NOT DO THIS

# WRONG - Hard-coded scaling factors
x_mm = x_m * 1000.0  # ✗ DO NOT DO THIS
field_m = field_mm / 1000.0  # ✗ DO NOT DO THIS
```

**Rationale**:
- **Single source of truth**: Units controlled by `rad.FldUnits()` only
- **Consistency**: All code uses same unit system set at initialization
- **Maintainability**: Changing units requires one line change, not searching for conversion factors
- **Error prevention**: Hard-coded conversions cause bugs when unit system changes

**Unit Detection**:

Use `rad.FldUnits()` without arguments to get current unit system:

```python
import radia as rad

# Get current units (returns multi-line string)
units_str = rad.FldUnits()
# Parse to detect length unit
if 'Length:  mm' in units_str:
    length_unit = 'mm'
    length_scale = 1.0  # No conversion needed
elif 'Length:  m' in units_str:
    length_unit = 'm'
    length_scale = 0.001  # mm to m
else:
    raise ValueError(f"Unknown length unit in: {units_str}")
```

**No Exceptions**:

All code, including `radia_vtk_export.py`, `radia_ngsolve.cpp`, `nastran_mesh_import.py`, must:
1. Query current unit system via `rad.FldUnits()`
2. Convert based on detected units, not hard-coded assumptions
3. Never assume Radia is using mm or m

This ensures code works regardless of user's `rad.FldUnits()` setting.

**Implementation Pattern**:

```python
# CORRECT - Use rad.FldUnits() to control units
import radia as rad

# Set unit system once at start
rad.FldUnits('mm')  # or 'm' for NGSolve integration

# Create hexahedral magnet with ObjHexahedron (100mm cube)
vertices = [[-50,-50,-50], [50,-50,-50], [50,50,-50], [-50,50,-50],
            [-50,-50,50], [50,-50,50], [50,50,50], [-50,50,50]]
magnet = rad.ObjHexahedron(vertices, [0, 0, 954930])  # 100mm, A/m
field = rad.Fld(magnet, 'b', [50, 50, 50])  # 50mm point

# Export - automatically uses correct units
from radia_vtk_export import exportGeometryToVTK
exportGeometryToVTK(magnet, 'output.vtk')  # Handles units internally
```

**Migration Guidelines**:

When removing hard-coded unit conversions:

1. **Identify conversion factors**: Search for `/1000`, `*1000`, `0.001`, etc.
2. **Determine intended unit system**: mm or m?
3. **Add `rad.FldUnits()` at script start**: Set unit system explicitly
4. **Remove conversion factors**: Use values directly in chosen unit system
5. **Update comments**: Document which unit system is used

**Files to Check**:

When writing or modifying code:
- ✓ Check for hard-coded `*1000`, `/1000`, `*0.001`, `/0.001`
- ✓ Ensure `rad.FldUnits()` is called at script start
- ✓ Verify no manual coordinate scaling
- ✓ Use `radia_vtk_export.py` for VTK export (handles units)

**Example - Corrected Code**:

Before (hard-coded conversions):
```python
# BAD - Hard-coded unit conversion
x_range = np.linspace(-90, 90, 21)  # mm
for pt in obs_points:
    f.write(f'{pt[0]/1000.0} {pt[1]/1000.0} {pt[2]/1000.0}\n')  # Manual mm->m
```

After (unit-aware):
```python
# GOOD - Use rad.FldUnits() and exportGeometryToVTK
rad.FldUnits('m')  # Set to meters
x_range = np.linspace(-0.09, 0.09, 21)  # m (no conversion needed)

# Use radia_vtk_export for geometry
from radia_vtk_export import exportGeometryToVTK
exportGeometryToVTK(magnet, 'output.vtk')  # Automatic unit handling

# For field data, use same unit system
for pt in obs_points:
    f.write(f'{pt[0]} {pt[1]} {pt[2]}\n')  # Already in meters
```

---

## NGSolve Mesh Access Policy

### Centralized Mesh Access via netgen_mesh_import.py

**Policy**:
- **All NGSolve mesh access** MUST use functions from `src/radia/netgen_mesh_import.py`
- **NEVER** directly access `mesh.ngmesh.Points()`, `mesh.vertices[]`, or `el.vertices[].nr` in any script
- **ALWAYS** import mesh handling functions from `netgen_mesh_import.py`
- **NO EXCEPTIONS**: This applies to all scripts including examples, tests, and debugging code

**Enforcement**:
- Direct mesh access is a bug source due to index confusion
- All new code MUST use `extract_elements()` or `netgen_mesh_to_radia()`
- Existing code with direct access MUST be refactored

**Rationale**:

NGSolve has two different indexing schemes that cause off-by-one errors:

| Access Method | Indexing | Notes |
|--------------|----------|-------|
| `mesh.ngmesh.Points()[i]` | **1-indexed** | Index 0 raises error, valid: 1 to nv |
| `mesh.vertices[i]` | **0-indexed** | Valid: 0 to nv-1 |
| `el.vertices[i].nr` | Returns value for **0-indexed** `mesh.vertices[]` | Use with `mesh.vertices[]` only |

**Common Bug Pattern**:
```python
# WRONG - Using 0-indexed .nr with 1-indexed ngmesh.Points()
for v in el.vertices:
    pt = mesh.ngmesh.Points()[v.nr]  # Off-by-one error!

# CORRECT - Use 0-indexed consistently
for v in el.vertices:
    vertex = mesh.vertices[v.nr]
    pt = vertex.point
```

**Correct Usage**:

```python
# Import from centralized module
from netgen_mesh_import import netgen_mesh_to_radia, extract_elements, TETRA_FACES

# Option 1: Direct conversion to Radia (recommended)
radia_obj = netgen_mesh_to_radia(mesh,
                                  material={'magnetization': [0, 0, 0]},
                                  units='m',
                                  material_filter='magnetic')

# Option 2: Extract elements for custom processing
elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']  # Already extracted correctly
    # ...
```

**Module Location**: `src/radia/netgen_mesh_import.py`

**Available Functions**:
- `netgen_mesh_to_radia()`: Convert entire mesh to Radia geometry (recommended)
- `extract_elements()`: Extract element data for custom processing
- `compute_element_centroid()`: Compute centroid from vertex list
- `create_radia_tetrahedron()`: Create single Radia tetrahedron
- `create_radia_hexahedron()`: Create single Radia hexahedron

**Available Constants**:
- `TETRA_FACES`: 1-indexed face topology for tetrahedra
- `HEX_FACES`: 1-indexed face topology for hexahedra
- `WEDGE_FACES`: 1-indexed face topology for wedges
- `PYRAMID_FACES`: 1-indexed face topology for pyramids

---

## Example Script Naming Convention

### Requirement: Consistent snake_case Naming with Functional Prefixes

**Goal**: All example scripts in `examples/` folder must follow a consistent naming convention for easy identification and organization.

**Policy**:

**Naming pattern**: `<prefix>_<description>.py`

- Use **snake_case** (all lowercase with underscores)
- Use **functional prefix** to indicate script purpose
- Use **descriptive names** that explain what the script does

**Standard prefixes**:

| Prefix | Purpose | Example |
|--------|---------|---------|
| `demo_` | Educational demonstration of a feature | `demo_batch_evaluation.py` |
| `example_` | Complete working example | `example_hmatrix_cache_usage.py` |
| `benchmark_` | Performance measurement script | `benchmark_solver_scaling.py` |
| `test_` | Validation/verification script | `test_batch_evaluation.py` |
| `verify_` | Correctness verification | `verify_curl_A_equals_B.py` |
| `compare_` | Comparison between methods | `compare_radia_ngsolve_cube.py` |
| `visualize_` | Visualization script | `visualize_field.py` |
| `run_` | Runner/orchestrator script | `run_all_benchmarks.py` |
| (none) | Descriptive physical model name | `sphere_in_quadrupole.py`, `arc_current_with_magnet.py` |

**Naming rules**:

1. **✓ CORRECT - snake_case**:
   ```
   sphere_in_quadrupole.py
   benchmark_solver_scaling.py
   demo_batch_evaluation.py
   verify_curl_A_equals_B.py
   ```

2. **✗ INCORRECT - CamelCase or PascalCase**:
   ```
   Cubit2Nastran.py        # Should be: cubit_to_nastran.py
   York_cubit_mesh.py      # Should be: york_cubit_mesh.py (already correct case, but York should be lowercase)
   CompareResults.py       # Should be: compare_results.py
   ```

3. **✗ INCORRECT - No prefix for functional scripts**:
   ```
   accuracy.py             # Should be: verify_accuracy.py or benchmark_accuracy.py
   plot.py                 # Should be: visualize_results.py or plot_benchmark_results.py
   ```

**Directory-specific guidelines**:

| Directory | Typical Prefixes | Notes |
|-----------|------------------|-------|
| `examples/simple_problems/` | (none), `demo_` | Physical model names preferred |
| `examples/solver_benchmarks/` | `benchmark_`, `run_`, `plot_`, `verify_` | Performance focus |
| `examples/ngsolve_integration/` | `demo_`, `example_`, `test_`, `verify_` | Educational + validation |
| `examples/background_fields/` | (none), `compare_` | Physical model names |
| `examples/electromagnet/` | `main_`, `visualize_` | Workflow scripts |

**Migration checklist**:

When renaming files:
1. Use `git mv` to preserve history: `git mv OldName.py new_name.py`
2. Update imports in other files
3. Update README.md references
4. Update documentation
5. Commit with clear message: `"Rename OldName.py to new_name.py (naming convention)"`

**Examples of good names**:

```
# Physical models - descriptive, no prefix needed
sphere_in_quadrupole.py           # Clear physics description
arc_current_with_magnet.py        # Clear what it models
cubic_polyhedron_magnet.py        # Clear geometry + physics

# Functional scripts - prefix required
demo_batch_evaluation.py          # Demo of batch feature
benchmark_solver_scaling.py       # Benchmark solver performance
verify_curl_A_equals_B.py         # Verify Maxwell equation
compare_radia_ngsolve.py          # Compare two methods
visualize_field.py                # Visualize field data
run_all_benchmarks.py             # Orchestrator script
```

**Files to rename** (current violations):

1. `background_fields/Cubit2Nastran.py` → `background_fields/cubit_to_nastran.py`
2. `electromagnet/York_cubit_mesh.py` → `electromagnet/york_cubit_mesh.py`

**Rationale**:
- **Consistency**: Easy to scan and find scripts by purpose
- **Clarity**: Prefix immediately indicates script type
- **Python convention**: PEP 8 recommends snake_case for module names
- **Sorting**: Related scripts group together alphabetically

---

## Mesh Generation Tools

### CAD Modeling and Mesh Generation Policy (2026-01-23)

**CRITICAL**: Use **Coreform Cubit for CAD modeling**, then import to **Netgen for meshing**.

**Workflow**:
```
Cubit (CAD modeling) → STEP export → Netgen (mesh import) → NGSolve (Curve for high-order)
```

**Rationale**:
1. **Cubit excels at CAD**: Complex geometry creation, Boolean operations, parametric modeling
2. **Netgen excels at meshing**: High-quality tetrahedral/hexahedral mesh generation
3. **STEP as interchange**: Standard CAD format, preserves geometry accurately
4. **High-order elements**: Use `mesh.Curve(order)` for curved boundaries

**Policy**:
- **Complex CAD**: Create in Cubit, export to STEP, import to Netgen
- **Simple geometry**: Can use Netgen OCC directly (`netgen.occ.Box`, `Sphere`, etc.)
- **Curved elements**: Always call `mesh.Curve(order)` after importing STEP geometry
- **Mesh quality**: Use Cubit's mesh controls for element size grading

**Implementation**:

```python
from netgen.occ import OCCGeometry
from ngsolve import Mesh

# Import STEP file (created in Cubit)
geo = OCCGeometry('model.step')

# Generate mesh
ngmesh = geo.GenerateMesh(maxh=0.05)

# Create NGSolve mesh and curve for high-order accuracy
mesh = Mesh(ngmesh)
mesh.Curve(3)  # 3rd order curved elements for accurate geometry

# Now use with Radia or NGSolve FEM
```

**Note**: The `mesh.Curve(order)` method is essential for accurate representation of curved boundaries when using high-order finite elements.

---

### Tool Selection by Element Type

| Element Type | Tool | Notes |
|--------------|------|-------|
| **Tetrahedral** | **Netgen** | Simple geometry. Uses `netgen.occ.Box` + `OCCGeometry.GenerateMesh()` |
| Tetrahedral | **Coreform Cubit** | Complex geometry. Uses `cubit_mesh_export.export_netgen()` |
| Tetrahedral | **GMSH via NGSolve** | Import .msh files using `ngsolve.Mesh()` |
| **Hexahedral** | **Coreform Cubit** | Required. Netgen cannot generate 3D hex meshes |
| Mixed (hex+tet) | **Coreform Cubit** | Required for mixed element meshes |

**Note**: Nastran BDF format is **REMOVED**. Use Cubit to read legacy .bdf files if needed.

### GMSH Mesh Import via NGSolve

GMSH meshes can be imported via NGSolve's `Mesh()` function:

```python
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

rad.FldUnits('m')

# Import GMSH mesh via NGSolve
# NGSolve automatically reads .msh format
mesh = Mesh('geometry.msh')

# Convert to Radia
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')
```

**Note**: GMSH meshes must use NGSolve-compatible format (Gmsh 2.2 ASCII or later).

### Netgen Tetrahedral Mesh Generation

```python
from netgen.occ import Box, Pnt, OCCGeometry
from netgen.meshing import MeshingParameters

# Create box geometry
p1 = Pnt(-0.5, -0.5, -0.5)
p2 = Pnt(0.5, 0.5, 0.5)
box = Box(p1, p2)

# Generate mesh - MUST wrap in OCCGeometry first
geo = OCCGeometry(box)
mp = MeshingParameters(maxh=0.3)  # Maximum element size
mesh = geo.GenerateMesh(mp)

# Extract nodes and tetrahedra
nodes = [[p[0], p[1], p[2]] for p in mesh.Points()]
tetrahedra = [[v.nr - 1 for v in el.vertices] for el in mesh.Elements3D()]
```

**Important Notes**:
- Use `OCCGeometry(box).GenerateMesh(mp)`, NOT `box.GenerateMesh(mp)`
- Use `v.nr - 1` for vertex indices, NOT `int(v) - 1`

### Coreform Cubit Mesh Export

For complex hexahedral meshes, use the **Coreform Cubit Mesh Export** tool:

**Repository**: https://github.com/ksugahar/Coreform_Cubit_Mesh_Export

**Features**:
- Export Cubit meshes directly to Python (no file I/O needed)
- Export to GMSH format (.msh) for NGSolve import
- Export to Netgen format for visualization and verification
- Supports hexahedral, tetrahedral, and mixed element meshes

**Usage with Radia**:

```python
# Option 1: Direct Python export from Cubit
from coreform_cubit_mesh_export import get_mesh_data
mesh_data = get_mesh_data()  # Get mesh directly from Cubit session

# Option 2: Export to GMSH, then import via NGSolve
# In Cubit: export_gmsh("mesh.msh")
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

mesh = Mesh('mesh.msh')
mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

**Note**: This is the recommended workflow for complex hexahedral geometries that Netgen cannot mesh directly.

### Coreform Cubit Policy for CplMag (PEEC-MMM Coupling)

**Policy (2026-01-11)**: For CplMag (coupled PEEC conductor + MMM magnetic material) simulations, use **Coreform Cubit** for hexahedral mesh generation.

**Rationale**:
1. **High-quality hex meshes**: Cubit produces better quality hexahedral meshes than automatic generators
2. **Complex geometries**: Magnetic cores often have complex shapes (E-cores, pot cores, toroids)
3. **Mesh control**: Cubit allows precise control over element size and distribution
4. **Multi-element MMM**: CplMag requires multiple MMM elements for accurate coupling (NOT single dipole)

**Workflow**:

```python
# Step 1: Generate hex mesh in Cubit (save as .cub or export to Nastran .bdf)

# Step 2: Import mesh using cubit_mesh_export
import sys
sys.path.insert(0, 'S:/CoreformCubit/01_GitHub')
from cubit_mesh_export import get_mesh_data

# Get hex elements from Cubit session
mesh_data = get_mesh_data()
hex_elements = mesh_data['hex_elements']  # List of [8 vertices] per element

# Step 3: Create Radia objects
import radia as rad
rad.FldUnits('m')

sub_cores = []
mat = rad.MatLin(mu_r)
for vertices in hex_elements:
    sub_core = rad.ObjHexahedron(vertices, [0, 0, 0])
    rad.MatApl(sub_core, mat)
    sub_cores.append(sub_core)

# Create container for multi-element core
core_container = rad.ObjCnt(sub_cores)

# Step 4: Create CplMag solver with multi-element core
coil = rad.CndLoop([0, 0, 0], radius, [0, 0, 1], 'r', w, h, sigma, n_radial, n_azimuthal)
solver = rad.CplMagCreate(coil, core_container)
rad.CplMagSetFrequency(solver, freq)
rad.CplMagSetMu(solver, mu_r, mu_r_imag)
result = rad.CplMagSolve(solver)
```

**Key Points**:
- **ObjCnt**: Use container to group multiple hex elements for CplMag
- **Multi-element**: Each element contributes to MMM interaction matrix
- **No single dipole**: Full element-to-element coupling is computed

**Repository**: `S:\CoreformCubit\01_GitHub` contains `cubit_mesh_export` utilities

### Coreform Cubit Policy for PEEC Conductor Mesh

**Policy (2026-01-11)**: For PEEC conductor meshes, use **Coreform Cubit** to generate surface meshes exported to Netgen format.

**Rationale**:
1. **Unified workflow**: Same tool for both MMM (magnetic) and PEEC (conductor) meshes
2. **Quality control**: Cubit provides better mesh quality for complex conductor geometries
3. **Wedge/Prism elements**: Cubit supports wedge elements for thin skin layers (induction heating)
4. **Curved elements**: Future support via NGSolve high-order curving (PR submitted)

**Current Limitation**:
- Radia currently supports **1st order elements only** (linear)
- Curved surface approximation uses piecewise-linear facets
- High-order curving (SetDeformation) is planned for future versions

**Workflow**:

```python
# Step 1: Generate conductor mesh in Cubit
# - Create conductor geometry (coil, wire, trace)
# - Mesh with surface elements (TRI, QUAD)
# - Export to Gmsh format

# Step 2: Import via NGSolve
from ngsolve import Mesh
from netgen.read_gmsh import ReadGmsh

mesh = Mesh(ReadGmsh("conductor_mesh.msh"))

# Step 3: Create PEEC conductor from mesh (future API)
# conductor = rad.CndFromMesh(mesh, sigma=5.8e7)  # Planned API

# Current workaround: Use CndLoop for simple geometries
coil = rad.CndLoop([0, 0, 0], radius, [0, 0, 1], 'r', w, h, sigma, n_radial, n_azimuthal)
```

**Note**: Full mesh-based PEEC conductor creation (`CndFromMesh`) is planned for future implementation.

### Hexahedral Mesh Import Functions (netgen_mesh_import.py)

**Functions** for importing hexahedral meshes into Radia:

| Function | Purpose |
|----------|---------|
| `cubit_hex_to_radia(hex_elements, ...)` | Convert Cubit hex element data to Radia geometry |
| `create_hex_mesh_grid(center, size, divisions, ...)` | Create structured hex mesh without Cubit |

**Usage with CplMag**:

```python
import radia as rad
from netgen_mesh_import import create_hex_mesh_grid, cubit_hex_to_radia

rad.FldUnits('m')

# Create coil conductor
coil = rad.CndLoop([0, 0, 0], 0.05, [0, 0, 1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)

# Create multi-element magnetic core (no Cubit needed)
core = create_hex_mesh_grid(
    center=[0, 0, 0],
    size=[0.03, 0.03, 0.03],  # 30mm cube
    divisions=[3, 3, 3],       # 27 elements
    mu_r=1000                  # mu_r = 1000
)

# Solve coupled system
solver = rad.CplMagCreate(coil, core)
rad.CplMagSetFrequency(solver, 1000)
rad.CplMagSetMu(solver, 1000, 0)
result = rad.CplMagSolve(solver)
```

**Usage with Cubit**:

```python
import radia as rad
from netgen_mesh_import import cubit_hex_to_radia

rad.FldUnits('m')

# Get hex elements from Cubit (via cubit_mesh_export)
# hex_elements = [[[x1,y1,z1], [x2,y2,z2], ..., [x8,y8,z8]], ...]

# Convert to Radia geometry
core = cubit_hex_to_radia(hex_elements, mu_r=1000)
```

---

## Mesh Operations Policy: Netgen with Coreform Cubit Integration

### Mesh APIs Dropped (2026-01-11)

**CRITICAL**: Radia's internal mesh operation APIs are **NOT SUPPORTED**.

**Dropped APIs**:
- `ObjDivMag` - Internal mesh subdivision (REMOVED)
- `ObjDivMagPln` - Plane-based subdivision (REMOVED)
- `ObjCutMag` - Cutting objects by plane (REMOVED from Python API)

**Policy**:
- **All mesh operations** must use **Netgen with Coreform Cubit integration**
- Coreform Cubit provides geometry and high-quality hex meshing
- Netgen/NGSolve provides the mesh import interface to Radia
- **Repository**: `S:\CoreformCubit\01_GitHub` contains `cubit_mesh_export` utilities

**Key Functions** (from `netgen_mesh_import.py`):
| Function | Purpose |
|----------|---------|
| `create_hex_mesh_grid()` | Simple structured hex mesh (no Cubit needed) |
| `cubit_hex_to_radia()` | Import Cubit hex mesh to Radia |
| `netgen_mesh_to_radia()` | Import Netgen/NGSolve mesh to Radia |

**Rationale**:
1. **Cubit produces higher quality meshes**: Professional meshing tool with quality control
2. **Unified workflow**: Same tool for both conductor (PEEC) and magnetic (MMM) meshes
3. **Flexibility**: Cubit supports complex geometries, grading, and boundary layers
4. **NGSolve integration**: Cubit meshes exported to Netgen/NGSolve format via `cubit_mesh_export`
5. **Maintenance**: Reduces Radia C++ codebase complexity

**Workflow**:

```python
import radia as rad
from netgen_mesh_import import create_hex_mesh_grid, cubit_hex_to_radia

rad.FldUnits('m')

# Method 1: Simple structured mesh (no Cubit needed)
core = create_hex_mesh_grid(
    center=[0, 0, 0],
    size=[0.1, 0.1, 0.1],
    divisions=[3, 3, 3],
    mu_r=1000
)

# Method 2: Complex geometry via Cubit + Netgen
# 1. Create geometry and mesh in Cubit
# 2. Export via cubit_mesh_export.export_netgen()
# 3. Import to Radia
core = cubit_hex_to_radia(hex_elements_from_cubit, mu_r=1000)

# Method 3: Tetrahedral mesh via Netgen
# 1. Create geometry in Cubit, export to STEP
# 2. Import STEP to Netgen, generate tet mesh
# 3. Import to Radia
from netgen_mesh_import import netgen_mesh_to_radia
core = netgen_mesh_to_radia(ngsolve_mesh, material={'magnetization': [0,0,0]})
```

**Note**: Legacy examples using `ObjDivMag` or `ObjCutMag` are DEPRECATED and will not run.

---

## MSC (Magnetic Surface Charge) Method

### Overview

Radia uses **MSC (Magnetic Surface Charge)** for all hexahedral elements:

| Element Type | Faces | DOF | Python API | Use Case |
|--------------|-------|-----|------------|----------|
| **Tetrahedron** | 4 triangular | 3 (Mx, My, Mz) | `netgen_mesh_to_radia()` | Complex curved geometry |
| **Hexahedron** | 6 quadrilateral | 6 (sigma per face) | `netgen_mesh_to_radia()` | Permanent magnets, soft iron |

**Policy (2025-12-27, updated 2025-12-31)**:
- **Python API**: Use `ObjHexahedron()` and `ObjTetrahedron()` for individual elements
- **Mesh import**: Use `netgen_mesh_to_radia()` for Netgen meshes
- **Tetrahedron**: 3 DOF (Mx, My, Mz) - MMM method with uniform magnetization
- **Hexahedron**: 6 DOF (sigma per face) - MSC method with surface charges
- 3 DOF hexahedron (MMM) is NOT supported - all hexahedra use 6 DOF MSC

**Python APIs for element creation**:
- `rad.ObjHexahedron(vertices, magnetization)` - 8 vertices, auto-generates faces
- `rad.ObjTetrahedron(vertices, magnetization)` - 4 vertices, auto-generates faces
- `netgen_mesh_import.netgen_mesh_to_radia()` - Batch import from Netgen mesh

### Implementation

**Source files** (`src/core/`):
- `rad_polyhedron.cpp`: Element type detection and dispatch
- `B_comp_tetrahedron_MSC()`: 4-face tetrahedral field computation
- `B_comp_hexahedron_MSC()`: 6-face hexahedral field computation
- `FieldFromChargedTriangle()`: Analytic field from charged triangle (solid angle formula)
- `FieldFromQuadFace()`: Quad face split into 2 triangles
- `rad_poly_analytical.cpp`: `RadFieldFromTriangleFaceGlobal()` for triangle integration

**Key Features**:
- Uses **global coordinates** directly (no local coordinate transformations)
- Computes field using **solid angle integration** formula (van Oosterom & Strackee, 1983)
- Handles **outward normal orientation** automatically

### Quick Start

**Tetrahedral mesh (Netgen)**:

```python
import radia as rad
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')

# Generate tetrahedral mesh
cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube_solid.mat('magnetic')
geo = OCCGeometry(cube_solid)
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Import to Radia (uses MSC automatically)
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')

# Apply material and solve
mat = rad.MatLin(999)  # mu_r = 1000
rad.MatApl(mag_obj, mat)
rad.Solve(mag_obj, 0.0001, 1000, 1)  # Method 1 = BiCGSTAB
```

**Hexahedral mesh (Netgen)**:

```python
import radia as rad
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')

# Generate hexahedral mesh using Netgen
# Note: Netgen generates tetrahedral meshes by default
# For true hexahedral meshes, use structured mesh generation
cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube_solid.mat('magnetic')
geo = OCCGeometry(cube_solid)
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Import to Radia
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')

# Apply material and solve
mat = rad.MatLin(999)  # mu_r = 1000
rad.MatApl(mag_obj, mat)
rad.Solve(mag_obj, 0.0001, 1000, 1)
```

### Benchmark Results (2025-12-19)

**Nonlinear Material (H_ext = 50,000 A/m)**:

| Element Type | N_elem | Solver | Iterations | M_avg_z (A/m) |
|--------------|--------|--------|------------|---------------|
| Hexahedral MSC | 125 | LU | 3 | 173,400 |
| Tetrahedral MSC | ~200 | LU | 2 | ~190,000 |

**Notes**:
- Tetrahedron: 3 DOF (Mx, My, Mz) - uniform magnetization
- Hexahedron: 6 DOF (sigma per face) - surface charge density
- LU solver (Method 0) and BiCGSTAB (Method 1) both work

### TrfMlt REMOVED (2026-01-31)

**CRITICAL**: `TrfMlt`, `TrfPlSym`, `TrfZerPara`, and `TrfZerPerp` have been **REMOVED** from Radia.

**Reason**: The shared-DOF design in TrfMlt was fundamentally incompatible with MSC 6DOF hexahedra. Element-based management (IMA) is the correct approach.

**Replacement**: Use **IMA (Image) Symmetry** for plane symmetry:

```python
import radia as rad

rad.FldUnits('m')

# Build full model geometry
hex_objects = [rad.ObjHexahedron(verts, [0,0,0]) for verts in all_vertices]
for h in hex_objects:
    rad.MatApl(h, rad.MatLin(mu_r))
container = rad.ObjCnt([coil] + hex_objects)

# Enable IMA x-mirror (half model)
intrc = rad.PreRelax(container, container)
n_ima = rad.SetIMASymmetry(intrc, 'x')  # 'x', 'y', or 'z' mirror
rad.BuildIMAMatrix(intrc)

# Solve with reduced DOF
rad.Solve(container, 0.0001, 100, 0)
B = rad.Fld(container, 'b', [0, 0, 0])
```

**Design Principle**: Element-based management (independent DOFs per element) is essential for correct physics. Face-based management (shared DOFs) causes errors.

See [docs/IMA_SYMMETRY_DESIGN.md](docs/IMA_SYMMETRY_DESIGN.md) for implementation details.

### Documentation

- [docs/MSC_QUICK_START.md](docs/MSC_QUICK_START.md): Quick start guide
- [examples/cube_uniform_field/nonlinear/](examples/cube_uniform_field/nonlinear/): Nonlinear benchmark

---

## PyPI Package Structure (v1.3.8+)

### Package Directory: src/radia (NOT src/python)

**Critical Requirement**: The package directory MUST be named `src/radia` for `pip install radia` to work correctly.

**Rationale**:
- `setuptools.find_packages(where="src")` discovers packages by directory name
- If the directory is named `python`, the wheel creates `python/` folder, not `radia/`
- Users would need `import python` instead of `import radia` - completely wrong!

**Correct Structure**:
```
src/
  radia/              # Package directory (import radia works)
    __init__.py       # Re-exports symbols from radia.pyd
    radia.pyd         # Core C++ extension module
    radia_ngsolve.pyd # Optional NGSolve integration
    *.py              # Python utility modules
```

**__init__.py Requirements**:
```python
# Radia Python package
__version__ = "1.3.8"

# Import all symbols from the C++ extension module
try:
    from radia.radia import *
except ImportError:
    from .radia import *
```

**Key Points**:
1. **Directory name = Package name**: `src/radia/` -> `import radia`
2. **C++ module import**: `__init__.py` must re-export symbols from `radia.pyd`
3. **setup.py package_data**: Use `"radia"` not `"python"` in package_data dict
4. **Version sync**: Keep version consistent across `__init__.py`, `setup.py`, `pyproject.toml`

### Build Policy: Always Use BuildRadiaInternal.ps1

**CRITICAL**: Always use the standard build script to ensure .pyd files are correctly copied.

**Policy**:
1. **ALWAYS use `BuildRadiaInternal.ps1`** for building - NEVER use manual cmake commands
2. The script automatically copies .pyd to `src/radia/radia.pyd` after build
3. The script verifies the copy was successful

**Standard Build Command**:
```powershell
# From project root directory:
powershell.exe -ExecutionPolicy Bypass -File BuildRadiaInternal.ps1

# With options:
powershell.exe -ExecutionPolicy Bypass -File BuildRadiaInternal.ps1 -Rebuild  # Clean rebuild
powershell.exe -ExecutionPolicy Bypass -File BuildRadiaInternal.ps1 -Verbose  # Verbose output
```

**Why This Matters**:
- CMake outputs `radia.cp312-win_amd64.pyd` (version-tagged name)
- Python package expects `src/radia/radia.pyd` (simple name)
- Manual cmake builds do NOT copy/rename the .pyd file
- Using old .pyd causes "AttributeError: module has no attribute" errors

**What the Script Does**:
1. Runs CMake configure and build
2. Finds `radia.cp312-win_amd64.pyd` in `build-msvc/`
3. Copies it to `src/radia/radia.pyd` (with rename)
4. Verifies the copy with timestamp check

**NEVER Do This**:
```powershell
# WRONG - Does NOT copy .pyd to src/radia
cmake --build build-msvc --config Release --target radia

# WRONG - Manual copy forgets the rename
copy build-msvc\radia.cp312-win_amd64.pyd src\radia\
```

**If You Must Use Manual Commands** (not recommended):
```powershell
# Build
cmake --build "s:\Radia\01_GitHub\build" --config Release --target radia

# MUST copy with rename
Copy-Item "build\Release\radia.cp312-win_amd64.pyd" "src\radia\radia.pyd" -Force

# Verify
(Get-Item "src\radia\radia.pyd").LastWriteTime
```

**Rationale**:
- `src/radia/` is the package source directory
- `python -m build` uses files from `src/radia/` to create wheels
- Ensures wheel contains the latest built binaries
- Avoids confusion about which .pyd is the "current" version

**Note**: The `.gitignore` excludes `.pyd` files, so these copied files are NOT committed to git.

**Testing Checklist** (before PyPI upload):
1. Build wheel: `python -m build`
2. Inspect wheel: `unzip -l dist/radia-x.y.z-py3-none-any.whl`
3. Verify structure: Contents should show `radia/radia.pyd`, NOT `python/radia.pyd`
4. Test install: `pip install dist/radia-x.y.z-py3-none-any.whl`
5. Test import: `python -c "import radia; print(radia.__version__)"`

---

## Benchmark Policy

### ベンチマーク実行ルール

1. **1ケース毎に実行**: 複数のソルバーやパラメータを比較する場合でも、1ケースずつ実行して結果を確認
2. **シーケンシャル実行のみ**: ベンチマークは並列実行せず、必ず1つずつ順番に実行する（メモリ測定の正確性のため）
3. **メモリ使用量を記録**: ピークメモリ使用量を測定し、結果JSONに含める

### メモリ測定方法

**psutil を使用**（C++拡張のメモリも含めて測定可能）:

```python
import psutil
import os

def get_peak_memory_mb():
    """Get peak memory usage in MB (Windows: peak_wset)"""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    # Windows: peak_wset = peak working set size
    # Linux/Mac: rss = resident set size (ピークは別途追跡必要)
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)  # MB
    else:
        return mem_info.rss / (1024 * 1024)  # MB (fallback)
```

**注意**: `tracemalloc` はPython内部のメモリのみ追跡するため、C++拡張モジュール (radia.pyd) のメモリは測定されない。

### 結果JSONフォーマット

```json
{
  "n_elements": 1000,
  "solver_method": 0,
  "t_solve": 5.2,
  "peak_memory_mb": 450.0,
  "iterations": 6,
  "M_avg_z": 149846.0
}
```

---

## Naming Policy: External Project References

### Do NOT Use "ELF" or "ELF_MAGIC" in Radia Codebase

**Policy (2025-12-16)**:
- **Do NOT use "ELF" or "ELF_MAGIC"** as terminology in Radia source code, documentation, or comments
- The MSC (Magnetic Surface Charge) method in Radia is a **Radia-native implementation**
- References to external projects should be kept to academic citations only

**Rationale**:
- Radia is an independent project with its own implementation
- Avoid confusion with external codebases
- Maintain clear intellectual property boundaries

**Allowed**:
- Academic citations in documentation (e.g., "Reference: Yano & Sugahara, J. Magn. Soc. Jpn., 2023")
- General algorithm descriptions (e.g., "MSC method", "surface charge method")

**Not Allowed**:
- Variable/function names containing "ELF" (e.g., `ELF_MAGIC_compatible`)
- Comments like "following ELF_MAGIC convention"
- Documentation referring to "ELF_MAGIC format"

**Migration**: If existing code contains "ELF" references, replace with Radia-native terminology.

---

## Logging and Debugging Policy

### No Console Output from C++ Code

**Policy (2025-12-21)**:
- **Do NOT add printf/cout/cerr output in C++ code** for logging or debugging
- All user-facing output should be handled through **Python scripts**
- C++ code should be silent unless there is a critical error that requires immediate attention

**Rationale**:
- Python stdout/stderr can be captured and formatted properly
- C++ printf may not appear in Python console depending on buffering
- Keeps separation of concerns: C++ for computation, Python for user interaction

**Allowed**:
- Error messages via Radia's error handling system (`Send.ErrorMessage(...)`)
- Debug output controlled by compile-time flags (`#ifdef DEBUG_...`)

**Not Allowed**:
- Unconditional `printf()`, `fprintf(stderr, ...)`, `std::cout`, `std::cerr`
- Timing output from C++ code (use Python's `time.time()` instead)

**Debugging Approach**:
1. Use Python-side timing: `time.time()` before/after calls
2. Use conditional compile flags: `#ifdef RADIA_DEBUG_TIMING`
3. Write to file only when explicitly enabled by environment variable

---

## Conductor Solver Policy: PEEC + Surface Impedance

**Approach**:
- **PEEC (Partial Element Equivalent Circuit)**: Loop-Star decomposition for coils and conductors
- **SIBC (Surface Impedance Boundary Condition)**: Skin effect via analytical formulas
- **ESIM (Effective Surface Impedance Method)**: Nonlinear/H-dependent surface impedance

**Target Applications** (PEEC + SIBC covers):
- **Induction heating**: 1 kHz - 500 kHz
- **Wireless power transfer (WPT)**: 6.78 MHz, 13.56 MHz
- **Power electronics**: DC - 1 MHz (inverters, converters, transformers)
- **Eddy current analysis**: Low to medium frequency

### NGSolve ngbem Integration

Radia PEEC is designed to work alongside **NGSolve ngbem** for unified electromagnetic analysis.

**ngbem** features:
- Helmholtz kernel (full-wave)
- H-matrix acceleration (HLib/H2Lib)
- EFIE/MFIE/CFIE formulations
- Native NGSolve integration

**Current frequency allocation**:

| Range | Solver | Use Case |
|-------|--------|----------|
| DC - 1 MHz | Radia PEEC + SIBC | Power electronics, WPT, transformers |
| 1 MHz - GHz | ngbem | RF heating, antennas, EMC/shielding |

**Future: ngbem low-frequency support** (requested via NGSolve issue):
- Loop-Star decomposition for MQS stability
- Laplace kernel option for quasi-static problems
- SIBC/ESIM integration

### Radia PEEC Unique Features

Even with ngbem low-frequency support, Radia PEEC provides:

| Feature | Radia PEEC | ngbem |
|---------|------------|-------|
| **Circuit extraction** | Direct (L, R, C) | Post-processing needed |
| **SPICE netlist** | Native output | Conversion needed |
| **Lanczos MOR** | Implemented | Separate implementation |
| **KAN/CFE learning** | Implemented | Not available |
| **MMM coupling** | CplMag | FEM-BEM coupling |
| **Schur complement** | Port extraction | Manual extraction |

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
│  - SPICE output     │         │                     │
└─────────────────────┘         └─────────────────────┘
           │                               │
           └───────────────┬───────────────┘
                           ▼
┌─────────────────────────────────────────────────────────────────┐
│  Unified Solution (NGSolve GridFunction)                        │
└─────────────────────────────────────────────────────────────────┘
```

### PRIMA Model Order Reduction

**Policy**: Use PRIMA (not CLN/Cauer) terminology for model order reduction.

**PRIMA vs CLN Equivalence**:
- Both use **Lanczos tridiagonalization** to produce RL ladder networks
- Mathematically identical: tridiagonal matrix -> series RL ladder
- PRIMA (1998, IEEE TCAD) is the standard academic reference
- CLN is a later repackaging with potential patent ambiguity

**Implementation** (in `lanczos_reduction.py`):
```python
# PRIMASchurExtractor class handles:
# - PRIMA Lanczos with re-orthogonalization (higher accuracy)
# - Per-group ACA tolerance (magnetic, dielectric, conductor)
# - Schur complement for port impedance extraction
# - SPICE netlist generation
```

**Key Classes**:
- `SPICEExtractionConfig`: Configuration with Lanczos order and ACA tolerances
- `PRIMASchurExtractor`: Full SPICE extraction workflow
- `LoopStarMagneticCoupled`: Loop-Star basis transformation

**Configuration Example**:
```python
config = SPICEExtractionConfig(
    n_lanczos_loop=20,          # Lanczos order for conductor loops
    n_lanczos_star=10,          # Lanczos order for capacitive nodes
    aca_tol_magnetic=1e-3,      # ACA tolerance for magnetic elements
    aca_tol_dielectric=1e-4,    # ACA tolerance for dielectric elements
    aca_tol_conductor=1e-5,     # ACA tolerance for conductor/shield
    port_indices=[0, 1],
)
```

**Rationale**:
1. **Low-frequency stability**: Eliminates breakdown at DC (jomega*L -> 0)
2. **MQS validity**: At power frequencies, displacement current is negligible
3. **Circuit extraction**: Direct mapping to RLC ladder elements
4. **Passivity**: PRIMA Lanczos preserves passivity
5. **Re-orthogonalization**: Higher accuracy than plain Lanczos

---

## Complex Permeability and ESIM Policy

### Complex Permeability Support (2026-01-09)

Radia SIBC (Surface Impedance Boundary Condition) supports **complex permeability**:

```
μ = μ' - jμ"
```

Where:
- **μ'** (real part): Energy storage, determines reactive power
- **μ"** (imaginary part): Energy loss from magnetic hysteresis/eddy currents in grains
- **Loss tangent**: tan(δ_m) = μ" / μ'

**Physical Effects**:

| Component | Physical Meaning | Power |
|-----------|------------------|-------|
| μ' | Magnetic energy storage | Reactive: Q' = (ω/2) * μ' * |H|^2 |
| μ" | Magnetic loss | Active: P_mag = (ω/2) * μ" * |H|^2 |

**API Usage**:

```python
from esim_cell_problem import ESIMCellProblemSolver

# Constant complex permeability (ferrite at 50 kHz)
solver = ESIMCellProblemSolver(
    sigma=0.01,           # S/m (ferrite is nearly insulating)
    frequency=50000,      # Hz
    complex_mu=(2000, 200)  # (μ'_r, μ"_r)
)

# H-dependent complex permeability (saturable ferromagnetic)
complex_mu_data = [
    [0, 2000, 200],      # [H (A/m), μ'_r, μ"_r]
    [1000, 1000, 100],
    [10000, 200, 20],
]
solver = ESIMCellProblemSolver(
    sigma=2e6,
    frequency=50000,
    complex_mu=complex_mu_data
)
```

### ESIM (Effective Surface Impedance Method)

**Policy**: ESIM is the recommended method for nonlinear magnetic materials in conductor problems.

**Field-Dependent Surface Impedance**:
```
Zs(H) = Re{Zs(H)} + j·Im{Zs(H)}
```

The surface impedance is field-dependent, not just frequency-dependent.

**1D Cell Problem**:
ESIM solves the 1D cell problem in the depth direction:

```
d/dz[(1/μ(z)) · dH/dz] = jωσ·H
```

where μ(z) = μ(H(z)) for nonlinear materials (e.g., from the B-H curve of a specific material).

**What ESIM Does**:
1. Solves 1D Cell Problem in depth direction for each surface H-field value
2. Computes effective surface impedance Zs(H0) that accounts for:
   - Skin effect (nonuniform current distribution)
   - Magnetic saturation (H-dependent μ)
   - Magnetic losses (complex μ = μ' - jμ")
3. Builds lookup table for fast 3D solver iteration

**When to Use ESIM**:
- Induction heating (workpiece heating analysis)
- Nonlinear iron cores in transformers/motors
- Lossy ferrite components at high frequency
- Any case where μ(H) varies significantly

**ESIM vs Standard SIBC**:

| Aspect | Standard SIBC | ESIM |
|--------|---------------|------|
| μ assumption | Constant | H-dependent |
| Magnetic loss | Optional (μ") | Included via μ"(H) |
| Accuracy in saturation | Poor | Good |
| Computational cost | Low | Higher (1D solve per H0) |

**Reference Implementation**:
- `src/radia/esim_cell_problem.py`: Cell problem solver with complex μ
- `src/radia/esim_coupled_solver.py`: 3D coupled solver using ESIM

**Literature Reference**:
K. Hollaus, V. Hanser, and M. Schobinger, "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics, 2025.

---

## Visualization Policy: NGSolve + VTK (2026-01-16)

### Unified Visualization Framework

**CRITICAL**: Radia uses **NGSolve/Netgen** and **VTK (PyVista/ParaView)** for all visualization needs.

**Policy**:
- **Default visualization**: PyVista (quick, interactive, Jupyter-friendly)
- **Publication-quality**: ParaView (fine-tuned rendering, high-resolution export)
- **Geometry visualization**: Netgen OCC + NGSolve Draw()
- **Field visualization**: VTS export + PyVista (default) / ParaView (publication)
- **Mesh visualization**: NGSolve mesh + Netgen GUI
- **DO NOT** implement custom visualization in Radia C++ code

**Tool Selection**:
| Purpose | Tool | Notes |
|---------|------|-------|
| Quick visualization | **PyVista** | Default, Jupyter integration |
| Interactive exploration | **PyVista** | Python scripting |
| Publication figures | **ParaView** | Fine control over rendering |
| High-resolution export | **ParaView** | Vector graphics (SVG, PDF) |
| Animation | **ParaView** | Keyframe animation |

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
│  ObjCylinder -> OCC    │                     │                 │
└─────────────────────────────────────────────────────────────────┘
```

### Analytical Objects → OCC Shapes (TODO)

**Goal**: Export Radia analytical objects (no mesh) as OCC shapes for unified visualization.

**Planned Implementation**:
```python
from netgen.occ import Box, Cylinder, Sphere
import radia as rad

# Radia analytical object
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# Export as OCC shape for visualization
occ_shape = rad.ExportOCC(magnet)  # Returns netgen.occ shape

# Combine with Cubit-imported geometry
from ngsolve import Mesh
mesh = Mesh(...)  # From Cubit export_netgen

# Unified visualization
from ngsolve.webgui import Draw
Draw(mesh)
Draw(occ_shape)  # Analytical object as CAD
```

**Supported Conversions** (TODO):
| Radia Object | OCC Shape | Notes |
|--------------|-----------|-------|
| ObjRecMag | Box | Rectangular permanent magnet |
| ObjCylMag | Cylinder | Cylindrical permanent magnet |
| ObjSphMag | Sphere | Spherical permanent magnet |
| ObjArcCur | Torus section | Arc current coil |
| ObjRaceTrk | Composite | Racetrack coil |

**Reference**: EMPY_Field implementation at `S:\NGSolve\EMPY\EMPY_Field`

### VTS Field Export

**Policy**: Use `rad.FldVTS()` for field data export.

```python
import radia as rad

rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# Export field grid to VTS
rad.FldVTS(magnet, 'field_output.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)
```

### PyVista Integration (Default)

**Policy**: Use PyVista as the **default** visualization tool.

```python
import pyvista as pv

# Read VTS field data
grid = pv.read('field_output.vts')

# Quick visualization (default workflow)
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='B_magnitude', cmap='coolwarm')
plotter.add_arrows(grid.points, grid['B_field'], mag=0.01)
plotter.show()

# Jupyter notebook integration
grid.plot(scalars='B_magnitude', cmap='coolwarm', jupyter_backend='static')
```

**Why PyVista as default**:
- Python-native, integrates with Jupyter notebooks
- Quick iterative visualization during development
- Scriptable for batch processing
- Good enough quality for most use cases

### ParaView (Publication Quality)

**Policy**: Use ParaView for **publication-quality** figures.

```bash
# Open VTS file in ParaView
paraview field_output.vts
```

**ParaView workflow for publications**:
1. Open VTS file in ParaView
2. Apply filters (Glyph, Contour, Slice, StreamTracer)
3. Adjust rendering (lighting, camera, colormap)
4. Export high-resolution image (PNG, TIFF) or vector graphics (SVG, PDF)

**When to use ParaView**:
- Journal paper figures (fine control over appearance)
- High-resolution exports (>300 DPI)
- Vector graphics export (SVG, PDF for LaTeX)
- Complex visualizations (streamlines, isosurfaces)
- Animation sequences

### NGSolve webgui

**Policy**: Use NGSolve webgui for interactive visualization.

```python
from ngsolve import *
from ngsolve.webgui import Draw

# Mesh from Cubit
mesh = Mesh(...)

# Field from Radia
from radia_ngsolve import RadiaField
B_cf = RadiaField(magnet, 'b')

# GridFunction projection
B_gf = GridFunction(HDiv(mesh, order=2))
B_gf.Set(B_cf)

# Interactive visualization
Draw(B_gf, mesh, name='B_field')
```

### Removed Legacy Visualization

**Removed APIs** (2026-01-09):
- `rad.ObjDrwVTK()` - Use NGSolve Draw() instead
- `exportGeometryToVTK()` - Use OCC export instead
- `radia_pyvista_viewer.py` - Use PyVista directly

---

## Universal Relaxation Network (URN) Policy (2026-01-19)

### URN Examples Directory

**CRITICAL**: All URN-related examples, data, and scripts MUST be placed in:
```
examples/Universal_Relaxation_Network/
```

**Directory Structure**:
```
examples/Universal_Relaxation_Network/
  data/
    synthetic/                    # Synthetic benchmark data
      liion_battery_eis.csv       # Physics-based synthetic Li-ion EIS
      mnzn_ferrite_impedance.csv  # Physics-based synthetic ferrite data
    real_world/                   # Publicly available real datasets
      nasa_battery/               # NASA Li-ion Battery Aging Dataset
      mendeley_eis/               # Mendeley SoC EIS Dataset
  universal_relaxation_network.py # Main URN implementation
  validate_urn_vs_vf.py           # Validation script (URN vs VF comparison)
  demo_spice_timedomain.py        # Time-domain SPICE simulation demo
```

**Policy**:
1. **All paper data here**: Data mentioned in `docs/paper/urn_paper.tex` MUST exist in this directory
2. **Synthetic data labeled**: Synthetic data MUST be clearly marked as synthetic in file headers
3. **Real data with attribution**: Real-world datasets MUST include license and citation info
4. **Reproducibility**: All paper results MUST be reproducible from scripts in this directory

**Do NOT**:
- Place URN examples in `examples/peec_integration/` (legacy location)
- Claim synthetic data as real measurements
- Use proprietary datasets without proper licensing

**Paper-Data Consistency**:
Any data file referenced in `urn_paper.tex` MUST:
1. Exist in `examples/Universal_Relaxation_Network/data/`
2. Have matching parameters (frequency range, impedance values)
3. Include header comments explaining data source

---

## Publication-Quality Figure Generation Policy (2026-01-19)

### Matplotlib Settings for IEEE/Academic Papers

**CRITICAL**: When generating figures for academic papers, use the following matplotlib settings for publication-quality PDF output.

**Required Settings**:

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Font settings: Times New Roman, 10pt at 8cm width
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 10
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 8

# High quality output
rcParams['figure.dpi'] = 300
rcParams['savefig.dpi'] = 300
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.02  # Minimal margins

# PDF font embedding: Type 42 (TrueType) for Acrobat compatibility
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# Tick settings: INWARD on ALL sides
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['xtick.bottom'] = True
rcParams['ytick.left'] = True
rcParams['ytick.right'] = True

# Line widths
rcParams['axes.linewidth'] = 0.5
rcParams['xtick.major.width'] = 0.5
rcParams['ytick.major.width'] = 0.5
rcParams['xtick.minor.width'] = 0.3
rcParams['ytick.minor.width'] = 0.3

# Figure size: 8cm width (standard single-column)
CM_TO_INCH = 1 / 2.54
FIG_WIDTH = 8 * CM_TO_INCH   # 8cm = 3.15 inches
FIG_HEIGHT = 6 * CM_TO_INCH  # Adjustable

fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))
# ... plot code ...
plt.savefig('figure.pdf', format='pdf')
```

**Key Requirements**:

| Setting | Value | Rationale |
|---------|-------|-----------|
| Font | Times New Roman | IEEE standard |
| Font size | 10pt at 8cm | Readable in print |
| Tick direction | Inward | Professional appearance |
| Ticks | All 4 sides | Complete axis frame |
| Margins | Minimal (0.02 in) | Maximize data area |
| Output | PDF | Vector graphics, scalable |
| DPI | 300 | High quality |

**Figure Dimensions**:

| Column Type | Width (cm) | Width (inch) |
|-------------|-----------|--------------|
| Single column | 8.0 | 3.15 |
| Double column | 17.0 | 6.69 |
| Full page | 19.0 | 7.48 |

**Do NOT**:
- Use PNG for paper figures (use PDF)
- Use default matplotlib fonts (use Times New Roman)
- Use outward ticks (use inward)
- Leave large margins (use minimal padding)
- Forget ticks on top/right axes

**Example Script Location**:
- `examples/Universal_Relaxation_Network/generate_paper_figures.py`

---

**Last Updated**: 2026-01-19 (Added Publication Figure Policy)
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation