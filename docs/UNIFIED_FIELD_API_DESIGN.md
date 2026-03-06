# Radia Unified Field API Design

## Overview

This document describes the unified field computation architecture for Radia,
integrating static magnetostatics, AC conductor analysis (RWG-EFIE), and hybrid
methods (PEEC+RWG) under a **single unified API**: `rad.Fld()`.

## Design Philosophy

**Core Principle**: All field computations use `rad.Fld(obj, field_type, point)`.
The function automatically detects the object type and dispatches to the appropriate solver.

```python
# Unified API - rad.Fld() handles everything
B = rad.Fld(obj, 'b', [0, 0, 0.1])

# For AC analysis, set frequency first
rad.CndSetFrequency(obj, 50000)  # 50 kHz
B_complex = rad.Fld(obj, 'b', [0, 0, 0.1])
# Returns [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im] for conductors
```

## Implementation Status (2026-01-08)

### Completed

- **CndFld() API removed**: All conductor field computation now uses `rad.Fld()`
- **Unified RadFld()**: Modified to detect conductor handles (>= 10000) and dispatch
- **IsConductorHandle()**: Helper function to check if handle is a conductor
- **ComputeConductorField()**: Conductor field computation with complex return values

### Return Value Convention

| Object Type | Return Format |
|------------|---------------|
| Magnetic (PM, iron) | `[Bx, By, Bz]` (3 real values) |
| Conductor (AC) | `[Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]` (6 values) |

### Handle Ranges

| Handle Range | Object Type |
|-------------|-------------|
| 1 - 9999 | Magnetic objects (radTg3d) |
| 10000+ | Conductor objects (radTConductor) |
| 20000+ | SIBC materials |

## Previous State (Before Unification)

| API | Target | Status |
|-----|--------|--------|
| `rad.Fld(obj, 'b', pt)` | PM, soft iron | Active |
| `rad.CndFld(cnd, 'b', pt)` | Conductors | **REMOVED** |
| `solver.compute_B_field()` | RWG surfaces | To be integrated |

### Problems Solved

1. **Single API**: Users now use `rad.Fld()` for all field computations
2. **Automatic dispatch**: Object type detection handles routing internally
3. **Consistent interface**: Same API for DC and AC analysis

## Unified Architecture

```
                         User API
    ┌─────────────────────────────────────────────────────────┐
    │                                                         │
    │     rad.Fld(obj, field_type, point)                    │
    │                                                         │
    │     - Detects object type automatically                │
    │     - Returns real for DC, complex for AC              │
    │     - Handles containers with mixed object types       │
    │                                                         │
    └────────────────────────┬────────────────────────────────┘
                             │
                             ▼
    ┌─────────────────────────────────────────────────────────┐
    │              Field Source Dispatcher (C++)              │
    │                                                         │
    │   radTg3d?          -> Static field (Biot-Savart/MSC)  │
    │   radTConductor?    -> AC field (RWG-EFIE)             │
    │   radTConductorGroup? -> Hybrid (PEEC+RWG)             │
    │   ObjCnt?           -> Sum all sources                 │
    │                                                         │
    └────────────────────────┬────────────────────────────────┘
                             │
         ┌───────────────────┼───────────────────┐
         │                   │                   │
         ▼                   ▼                   ▼
┌─────────────────┐ ┌─────────────────┐ ┌─────────────────┐
│ Static Magneto  │ │ RWG-EFIE Solver │ │ Hybrid PEEC+RWG │
│ (Biot-Savart)   │ │ (Surface IE)    │ │ (Coil+Work)     │
│                 │ │                 │ │                 │
│ - ObjRecMag     │ │ - radTConductor │ │ - ConductorGroup│
│ - ObjHexahedron │ │ - RWGMesh       │ │ - Coupled solve │
│ - ObjTetrahedron│ │ - MQS/Full-wave │ │ - Mutual M      │
└─────────────────┘ └─────────────────┘ └─────────────────┘
```

## Implementation Plan

### Phase 1: C++ RWG-EFIE Integration (Priority: High)

**Goal**: Make `rad.CndFld()` use C++ RWG-EFIE solver instead of Python.

**Files to modify**:
- `src/core/rad_rwg_basis.cpp`: Complete `RWGEFIESolver::Solve()` and `ComputeB()`
- `src/core/rad_conductor.cpp`: Use `RWGEFIESolver` in `radTConductorSolver`
- `src/lib/rad_conductor_api.cpp`: Connect `CndFld*` functions to RWG solver

**Key Implementation**:

```cpp
// rad_conductor.cpp
void radTConductorSolver::Solve() {
    // Build RWG mesh from conductor panels
    BuildRWGMesh();

    // Create and configure RWG solver
    rwgSolver_ = std::make_unique<RWGEFIESolver>();
    rwgSolver_->SetMesh(rwgMesh_);
    rwgSolver_->SetConductivity(conductors_[0]->GetConductivity());
    rwgSolver_->SetFrequency(frequency_);

    // Define port from conductor definition
    rwgSolver_->DefinePort(GetPortEdges());
    rwgSolver_->SetVoltageExcitation(GetExcitationVoltage());

    // Solve EFIE
    rwgSolver_->Solve();

    // Extract results
    portImpedance_ = rwgSolver_->GetImpedance();
    totalCurrent_ = rwgSolver_->GetCurrent();
}

void radTConductorSolver::ComputeB(const TVector3d& point,
                                    std::complex<double>& Bx,
                                    std::complex<double>& By,
                                    std::complex<double>& Bz) const {
    if (rwgSolver_) {
        rwgSolver_->ComputeB(point, Bx, By, Bz);
    }
}
```

### Phase 2: Unified Field API (Priority: High)

**Goal**: Single `rad.Fld()` API that handles both static and AC cases.

**Selected Approach: Extend rad.Fld() (Option A)**

```python
# Case 1: Static magnetostatics (existing behavior, unchanged)
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1e6])
B = rad.Fld(magnet, 'b', [0, 0, 0.1])  # Returns [Bx, By, Bz]

# Case 2: AC conductor analysis (NEW)
coil = rad.CndLoop([0, 0, 0], 0.05, [0, 0, 1], 'c', 0.002, 0.002, 5.8e7, 8, 30)
rad.CndSetFrequency(coil, 50000)  # 50 kHz
rad.CndSolve(coil)
B_complex = rad.Fld(coil, 'b', [0, 0, 0.1])  # Returns [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]

# Case 3: Mixed static + AC (NEW)
# Container with both magnet and conductor
container = rad.ObjCnt([magnet])
rad.ObjCntAdd(container, coil)  # Add conductor to same container
B_total = rad.Fld(container, 'b', [0, 0, 0.1])  # Sum of all fields
```

**Return Value Convention**:
- Static objects: `[Bx, By, Bz]` (3 elements, real)
- AC conductors: `[Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]` (6 elements, complex)
- Mixed container: Returns 6 elements if any AC conductor present

**Implementation in C++**:

```cpp
// radentry.cpp - Unified Fld() dispatcher
void Fld(double* pB, int obj, const char* fieldId, double* pP) {
    // Check object type
    radTg3d* pObj = GetObject(obj);

    if (IsConductor(obj)) {
        // AC conductor field
        radTConductor* pCond = GetConductor(obj);
        std::complex<double> Bx, By, Bz;
        pCond->ComputeB(TVector3d(pP[0], pP[1], pP[2]), Bx, By, Bz);

        // Return 6 values: [real, real, real, imag, imag, imag]
        pB[0] = Bx.real(); pB[1] = By.real(); pB[2] = Bz.real();
        pB[3] = Bx.imag(); pB[4] = By.imag(); pB[5] = Bz.imag();
    }
    else if (pObj != nullptr) {
        // Static field (existing code)
        TVector3d B;
        pObj->B_comp(TVector3d(pP[0], pP[1], pP[2]), B);
        pB[0] = B.x; pB[1] = B.y; pB[2] = B.z;
    }
}
```

**Why Option A (unified rad.Fld) instead of separate rad.FldAC**:
- Consistent with Radia design philosophy: single API for field computation
- No confusion about which function to use
- Natural extension: conductors are just another field source
- Container support: mix magnets and conductors in same ObjCnt

### Phase 3: Hybrid PEEC+RWG (Priority: Medium)

**Goal**: Unified coil-workpiece analysis with automatic method selection.

**New Class**: `radTConductorGroup`

```cpp
class radTConductorGroup {
public:
    // Add coil (uses PEEC automatically for wire geometry)
    void AddCoil(std::shared_ptr<radTConductor> coil);

    // Add workpiece (uses RWG automatically for surface geometry)
    void AddWorkpiece(std::shared_ptr<radTConductor> workpiece);

    // Solve coupled system
    void Solve();

    // Get results
    std::complex<double> GetCoilImpedance() const;
    std::complex<double> GetMutualImpedance() const;
    double GetWorkpiecePower() const;

    // Field computation (sum of all sources)
    void ComputeB(const TVector3d& point, ...);

private:
    std::vector<std::shared_ptr<radTConductor>> coils_;
    std::vector<std::shared_ptr<radTConductor>> workpieces_;

    // Coupled impedance matrix
    // [Z_cc  Z_cw] [I_c]   [V_c]
    // [Z_wc  Z_ww] [I_w] = [0  ]
};
```

**Python API**:

```python
# Create group
group = rad.CndGroup()

# Add coil (PEEC)
coil = rad.CndLoop([0, 0, 0.02], 0.04, [0, 0, 1], 'c', 0.002, 0.002, 5.8e7, 8, 30)
rad.CndGroupAddCoil(group, coil)

# Add workpiece (RWG)
work = rad.CndRecBlock([0, 0, 0], [0.1, 0.1, 0.005], 2e6, 8)
rad.CndGroupAddWorkpiece(group, work)

# Set frequency and solve
rad.CndSetFrequency(group, 50000)
rad.CndSolve(group)

# Get results
Z = rad.CndGetImpedance(group)
P = rad.CndGetPower(group)  # Power in workpiece

# Field computation
B = rad.CndFld(group, 'b', [0.05, 0, 0])
```

### Phase 4: Loop-Star Decomposition (Priority: Low)

**Goal**: Low-frequency stable MQS formulation in C++.

Move `LoopStarDecomposition` and `LoopStarEFIESolver` from Python to C++.

```cpp
class LoopStarDecomposition {
public:
    void Build(const RWGMesh& mesh);

    // Transform RWG coefficients to Loop-Star
    void ToLoopStar(const std::vector<std::complex<double>>& rwg,
                    std::vector<std::complex<double>>& loop,
                    std::vector<std::complex<double>>& star);

    // Transform back
    void FromLoopStar(const std::vector<std::complex<double>>& loop,
                      const std::vector<std::complex<double>>& star,
                      std::vector<std::complex<double>>& rwg);

private:
    // Loop functions: one per triangle (div J = 0)
    // Star functions: one per vertex minus one (curl J = 0)
    Eigen::SparseMatrix<double> T_loop_;
    Eigen::SparseMatrix<double> T_star_;
};
```

## Migration Path

### For Existing Python RWG Users

```python
# Old (Python solver)
from radia.rwg_efie_solver import RWGEFIESolver, RWGMesh
mesh = RWGMesh()
mesh.create_wire_mesh(path, wire_radius)
solver = RWGEFIESolver(mesh)
solver.set_frequency(50000)
solver.solve()
B = solver.compute_B_field([0, 0, 0.1])

# New (C++ integrated)
import radia as rad
cond = rad.CndWire(path, 'c', wire_radius, wire_radius, 5.8e7, 8)
rad.CndSetFrequency(cond, 50000)
rad.CndSolve(cond)
B = rad.CndFld(cond, 'b', [0, 0, 0.1])
```

### For Hybrid Coil-Workpiece Analysis

```python
# Old (Python hybrid solver)
from radia.rwg_efie_solver import HybridPEECRWGSolver
solver = HybridPEECRWGSolver(coil_path, wire_radius, workpiece_mesh, ...)
result = solver.solve()

# New (C++ integrated)
import radia as rad
group = rad.CndGroup()
coil = rad.CndSpiral(...)
work = rad.CndRecBlock(...)
rad.CndGroupAddCoil(group, coil)
rad.CndGroupAddWorkpiece(group, work)
rad.CndSolve(group)
```

## API Reference (Proposed)

### Conductor Creation

| Function | Description |
|----------|-------------|
| `CndRecBlock(center, size, sigma, n)` | Rectangular block conductor |
| `CndLoop(center, R, normal, cs, w, h, sigma, na, nl)` | Circular loop |
| `CndSpiral(center, Ri, Ro, pitch, turns, axis, cs, w, h, sigma, na)` | Spiral coil |
| `CndWire(path, cs, w, h, sigma, na)` | Wire along path |
| `CndGroup()` | Container for coupled conductors |

### Analysis

| Function | Description |
|----------|-------------|
| `CndSetFrequency(cnd, freq)` | Set analysis frequency |
| `CndSolve(cnd)` | Solve impedance |
| `CndGetImpedance(cnd)` | Get port impedance |
| `CndGetPower(cnd)` | Get power dissipation |

### Field Computation (Unified API)

| Function | Description |
|----------|-------------|
| `Fld(obj, type, pt)` | **Unified field API** - handles all object types |
| `FldBatch(obj, type, pts)` | Batch field computation at multiple points |

**Return values**:
- Static objects (radTg3d): `[Bx, By, Bz]` (3 real values)
- AC conductors (radTConductor): `[Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]` (6 values)
- Mixed container: Returns 6 values if any AC conductor present

**Deprecated APIs** (kept for backward compatibility):
- `CndFld()` - Use `Fld()` instead
- `CndFldBatch()` - Use `FldBatch()` instead

### Usage Examples

```python
import radia as rad

# 1. Static magnetostatics (existing)
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1e6])
B_static = rad.Fld(magnet, 'b', [0, 0, 0.15])  # [Bx, By, Bz]

# 2. AC conductor
coil = rad.CndLoop([0, 0, 0], 0.05, [0, 0, 1], 'c', 0.002, 0.002, 5.8e7, 8, 30)
rad.CndSetFrequency(coil, 50000)
rad.CndSolve(coil)
B_ac = rad.Fld(coil, 'b', [0, 0, 0.1])  # [Bx_re, By_re, Bz_re, Bx_im, By_im, Bz_im]

# 3. Combined static + AC
container = rad.ObjCnt([magnet])
rad.ObjCntAdd(container, coil)
B_total = rad.Fld(container, 'b', [0, 0, 0.1])  # Sum of all fields

# 4. Coupled coil-workpiece
group = rad.CndGroup()
rad.CndGroupAddCoil(group, coil)
work = rad.CndRecBlock([0, 0, 0], [0.1, 0.1, 0.005], 2e6, 8)
rad.CndGroupAddWorkpiece(group, work)
rad.CndSolve(group)
B_coupled = rad.Fld(group, 'b', [0, 0, 0.1])
P_heat = rad.CndGetPower(group)  # Heating power in workpiece
```

## Implementation Timeline

| Phase | Description | Priority | Status |
|-------|-------------|----------|--------|
| 1 | C++ RWG-EFIE integration | High | Pending |
| 2 | Unified `rad.Fld()` API | High | **Design complete** |
| 3 | Hybrid PEEC+RWG | Medium | Pending |
| 4 | Loop-Star decomposition | Low | Python impl done |

## Next Steps

1. **Modify `radentry.cpp`**: Add conductor type detection in `Fld()` function
2. **Complete `RWGEFIESolver::ComputeB()`** in C++
3. **Add `IsConductor()` helper** to check object type
4. **Update Python bindings** in `radpy_pyapi.cpp`
5. **Deprecate `CndFld()`** with warning message pointing to `Fld()`

## References

1. S.M. Rao et al., "Electromagnetic scattering by surfaces of arbitrary shape," IEEE TAP, 1982
2. A.E. Ruehli, "Equivalent Circuit Models for 3D Multiconductor Systems," IEEE Trans. MTT, 1974
3. J.R. Phillips, J.K. White, "A precorrected-FFT method," IEEE TCAD, 1997
4. G. Vecchi, "Loop-Star Decomposition of Basis Functions in EFIE," IEEE TAP, 1999
