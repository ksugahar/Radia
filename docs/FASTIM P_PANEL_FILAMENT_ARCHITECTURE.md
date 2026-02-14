# FastImp Panel + Filament Architecture

**Date**: 2026-02-13
**Confirmed**: SIBC/ESIM applies to panel elements (surface)

---

## Dual Representation: Filament + Panel

FastImp uses a **two-level discretization** for conductors:

```
Conductor (3D solid)
    │
    ├── Filaments (1D) ← Loop elements, carry current I
    │   └── Centerline segments (線上要素)
    │
    └── Panels (2D) ← Star elements, carry surface charge σ
        └── Surface patches (表面要素)
```

### Why Two Levels?

| Purpose | Element Type | DOF | Formula |
|---------|-------------|-----|---------|
| **Magnetic field** (inductive) | Filament (Loop) | Current I | Biot-Savart: `B = ∫ (I dl × r) / r³` |
| **Electric field** (capacitive) | Panel (Star) | Charge σ | Coulomb: `E = ∫ σ dS / (4πε₀r)` |

---

## Loop-Star Decomposition

### System Equation

```
[Z_LL   Z_LS] [I_L]   [V_L]
[Z_SL   Z_SS] [I_S] = [V_S]

where:
  Z_LL = R + jω*L + Z_s  (Loop-Loop: inductive + resistive + skin effect)
  Z_SS = P / (jω*ε)      (Star-Star: capacitive)
  Z_LS = jω*M_LS         (Loop-Star coupling)
```

### Matrix Components

| Matrix | Dimension | Physical Meaning | Computed From |
|--------|-----------|------------------|---------------|
| **L** | `n_loop × n_loop` | Inductance | Filament-filament Neumann formula |
| **R** | `n_loop × n_loop` | DC resistance | Filament geometry + σ |
| **Z_s** | `n_loop × n_loop` | Surface impedance | **SIBC/ESIM on panels** |
| **P** | `n_star × n_star` | Potential coefficient | Panel-panel integration |
| **M_LS** | `n_loop × n_star` | Loop-Star coupling | Filament-panel integration |

---

## Surface Impedance Boundary Condition (SIBC)

### Where SIBC is Applied

**CRITICAL**: SIBC/ESIM is applied to **PANEL elements**, not filaments.

```cpp
// rad_peec_matrices.h:248
void PEECSolver::SetSurfaceImpedance(
    const std::vector<std::complex<double>>& Zs_diag
);

// Effect on Z_LL matrix:
Z_LL[i][i] = R[i] + jω*L[i][i] + Z_s[i]
                                   ↑
                           Surface Impedance (from SIBC/ESIM)
```

### SIBC Formula

For a conductor panel with:
- Conductivity: `σ` [S/m]
- Permeability: `μ = μ₀μᵣ` [H/m]
- Frequency: `f` [Hz]

**Surface impedance**:
```
Z_s = (1 + j) / (σ * δ)

where δ = sqrt(2 / (ω*μ*σ))  (skin depth)
```

### ESIM (Effective Surface Impedance Method)

For **nonlinear materials** (μ depends on H):
- Solve 1D cell problem in depth direction
- Compute effective Z_s(H₀) at each panel
- Build lookup table for fast 3D iteration

**Reference**: K. Hollaus et al., IEEE Trans. Magnetics, 2025

---

## Current Implementation Status

### Implemented Components

| Component | File | Status |
|-----------|------|--------|
| **Panel structure** | `rad_conductor.h:50` | ✅ Implemented |
| **SIBC API** | `rad_peec_matrices.h:248` | ✅ Implemented |
| **Filament structure** | `rad_peec_matrices.h:52` | ✅ Implemented |

### Code Locations

**1. Panel Definition** (`rad_conductor.h`):
```cpp
struct SurfacePanel {
    TVector3d center;         // Panel center
    TVector3d normal;         // Outward normal
    double area;              // Panel area
    std::vector<TVector3d> vertices;  // 3 or 4 vertices
    enum Type { Triangle, Quadrilateral } type;
};
```

**2. Filament Definition** (`rad_peec_matrices.h`):
```cpp
struct PEECSegment {
    TVector3d center;       // Segment center [m]
    TVector3d direction;    // Unit direction vector
    double length;          // Segment length [m]
    double width;           // Cross-section width [m]
    double height;          // Cross-section height [m]
    double sigma;           // Conductivity [S/m]
};
```

**3. SIBC Interface** (`rad_peec_matrices.cpp:294`):
```cpp
void PEECSolver::SetSurfaceImpedance(
    const std::vector<std::complex<double>>& Zs_diag
) {
    Zs_ = Zs_diag;
    hasSurfaceImpedance_ = !Zs_.empty();
}
```

---

## Problem: GMD Approximation

### Current Implementation (WRONG)

**Location**: `rad_peec_matrices.cpp:207-223`

```cpp
double PEECMatrixBuilder::SelfInductance(const PEECSegment& seg) const {
    // GMD approximation for self-inductance
    // L = (mu_0 / 2*pi) * l * (ln(2*l/GMD) - 1)
    // GMD for rectangular cross-section: GMD ~ 0.2235 * (w + h)

    double l = seg.length;
    double gmd = 0.2235 * (seg.width + seg.height);  // ← WRONG

    return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l *
           (std::log(2.0 * l / gmd) - 1.0);
}
```

**Why it's wrong**:
1. **Converts rectangular to circular**: GMD = 0.2235 * (w + h) ≈ equivalent radius
2. **Introduces error**: ~6% in inductance for square cross-section
3. **Not FastImp approach**: FastImp uses panel integration, not GMD

### Correct Implementation (Neumann Formula)

**Replace with segment integration**:

```cpp
double PEECMatrixBuilder::SelfInductance(const PEECSegment& seg) const {
    // Neumann formula with numerical integration
    // L_ii = (mu_0/4π) * ∫∫ (d_i · d_j) / |r_i - r_j| dl_i dl_j

    // For straight segment: use Grover formula or Gauss quadrature
    // NO GMD approximation
}
```

**Grover Formula** (exact for straight segment):
```
L = (mu_0/2π) * l * [ln(2*l/sqrt(w²+h²)) + 0.25 + (w²+h²)/(12*l²)]
```

---

## Integration Plan

### Phase 1: Remove GMD (Current)

- [x] Locate GMD approximation code
- [x] Verify SIBC applies to panels (confirmed)
- [x] Document panel+filament architecture
- [ ] **Replace `SelfInductance()` with Neumann integration**

### Phase 2: Panel Integration

1. **Add panel support to PEECMatrixBuilder**:
   ```cpp
   class PEECMatrixBuilder {
       void AddPanel(const SurfacePanel& panel);
       void ComputePanelP();  // Star-Star matrix
       void ComputeFilamentPanelM_LS();  // Loop-Star coupling
   };
   ```

2. **Update GMSH workflow**:
   - Centerline mesh → Filaments (already working)
   - Surface mesh → Panels (NEW)

### Phase 3: SIBC/ESIM Integration

1. **Compute Z_s for each panel**:
   ```python
   # For each panel:
   Z_s[i] = compute_sibc(sigma, mu_r, freq, H_magnitude[i])

   # Or with ESIM:
   Z_s[i] = esim_solver.compute_surface_impedance(H0[i])
   ```

2. **Apply to Z_LL matrix**:
   ```cpp
   // In BuildImpedanceMatrix():
   Z_LL[i][i] += Zs_[i];  // Add surface impedance to diagonal
   ```

---

## Workflow Comparison

### Old (GMD Approximation):

```
1D Mesh (edges) → Filaments → GMD → L matrix (6% error)
                                ↑
                           WRONG (approximation)
```

### New (Panel-based FastImp):

```
1D Mesh (centerline) → Filaments → Neumann formula → L matrix
2D Mesh (surface)    → Panels    → SIBC/ESIM      → Z_s
                                                       ↓
                                              Z_LL = R + jωL + Z_s
```

---

## References

1. Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, vol. 24, no. 7, pp. 981-998, 2005
2. G. Vecchi, "Loop-Star Decomposition", IEEE TAP, vol. 47, no. 2, pp. 339-346, 1999
3. K. Hollaus et al., "Effective Surface Impedance in Scalar Potential Formulation", IEEE Trans. Magnetics, 2025
4. F. W. Grover, "Inductance Calculations", Dover Publications, 1946

---

**Status**: Architecture verified
**Next**: Implement Neumann formula to replace GMD approximation
