# FastImp Loop-Star PEEC Integration Plan

**Date**: 2026-02-13
**Task**: Replace GMD approximation with panel-based FastImp Loop-Star formulation

---

## Current Implementation Status

### Existing Components

1. **SurfacePanel Structure** (`rad_conductor.h`):
   ```cpp
   struct SurfacePanel {
       TVector3d center;         // Panel center
       TVector3d normal;         // Outward normal
       double area;              // Panel area
       std::vector<TVector3d> vertices;  // Panel vertices (3 or 4)
       enum Type { Triangle, Quadrilateral } type;
   };
   ```

2. **PEECSegment Structure** (`rad_peec_matrices.h`):
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

3. **GMD Approximation** (`rad_peec_matrices.cpp:207-223`):
   - **PROBLEM**: Uses GMD approximation for self-inductance
   - Formula: `L = (mu_0/2pi) * l * (ln(2*l/GMD) - 1)` where `GMD = 0.2235 * (w + h)`
   - **Result**: 6% error in inductance (too high)

### References in Codebase

| File | Line | Reference |
|------|------|-----------|
| `rad_conductor.h` | 10 | `[1] Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005` |
| `radentry.cpp` | 1907 | `// Conductor Analysis API Implementation (FastImp-based)` |
| `rad_conductor.cpp` | 759 | `// If totalCurrent_ is set (from Solve), use filament Biot-Savart` |
| `rad_conductor.cpp` | 914 | `// For wire conductors, use filament approximation along wire centerline` |

---

## FastImp Approach (from literature)

### Two-Level Discretization

FastImp uses a **filament + panel** approach:

1. **Filaments** (Loop elements):
   - 1D line segments along conductor centerline
   - Carry current I (Loop DOF)
   - Used for magnetic field computation (Biot-Savart)

2. **Panels** (Star elements):
   - 2D surface patches on conductor surface
   - Carry surface charge sigma (Star DOF)
   - Used for electric field computation (capacitance)

### Filament-Panel PEEC System (Loop-Star Transformation NOT Needed)

**CRITICAL**: In the Filament+Panel formulation, the Loop-Star basis transformation (Vecchi 1999) is **NOT needed**.

In MoM/BEM with RWG basis functions, a single set of basis functions must be decomposed
into solenoidal (Loop) and irrotational (Star) parts for low-frequency numerical stability.
However, in the PEEC Filament+Panel formulation, Filament = Loop elements and
Panel = Star elements are **inherently separate from the start**.

```
PEEC System (Filament+Panel):
  [Z_LL   Z_LS] [I_filament]   [V]
  [Z_SL   Z_SS] [Q_panel   ] = [0]

where:
  I_filament = current through filaments (Loop unknowns)
  Q_panel    = charge on panels (Star unknowns)

  Z_LL = R + jω*L + Z_s  (Filament-Filament: inductive + resistive)
  Z_SS = P / (jω)        (Panel-Panel: capacitive, P = potential coefficient)
  Z_LS = jω*M_LS         (Filament-Panel coupling)
  Z_SL = Z_LS^T          (reciprocity)
```

**Why Loop-Star Transformation is NOT Needed**:

| Method | Loop-Star Separation | Reason |
|--------|---------------------|--------|
| **MoM/BEM (RWG)** | **Required** | Single basis set must be decomposed into solenoidal/irrotational parts |
| **PEEC Filament+Panel** | **NOT needed** | Filament=Loop, Panel=Star are inherently separate from the formulation |

**Reference**: FastImp (Zhu et al., 2005) uses this same architecture.

### Matrix Computation (FastImp)

1. **L matrix (inductance)**: Neumann formula with numerical integration over filament pairs
   ```
   L_ij = (mu_0/4pi) * integral_i integral_j (d_i · d_j) / r
   ```
   - NOT GMD approximation
   - Full segment-to-segment integration

2. **P matrix (potential coefficient)**: Panel-to-panel integration
   ```
   P_ij = (1/4pi*eps_0) * integral_i integral_j dS_i dS_j / r
   ```

3. **M_LS matrix** (Loop-Star coupling): Filament-to-panel integration
   ```
   M_LS[i][j] = (mu_0/4pi) * integral_filament_i integral_panel_j (d_i · n_j) / r
   ```

---

## Problems with Current Implementation

### 1. GMD Approximation Error

**Location**: `rad_peec_matrices.cpp:207-223`

```cpp
double PEECMatrixBuilder::SelfInductance(const PEECSegment& seg) const {
    // GMD approximation for self-inductance
    // L = (mu_0 / 2*pi) * l * (ln(2*l/GMD) - 1)
    // GMD for rectangular cross-section: GMD ~ 0.2235 * (w + h)

    double l = seg.length;
    double gmd = 0.2235 * (seg.width + seg.height);  // ← APPROXIMATION

    return (PEEC_MU_0 / (2.0 * RadConst::PI)) * l * (std::log(2.0 * l / gmd) - 1.0);
}
```

**Why it's wrong**:
- GMD converts rectangular cross-section to equivalent circular radius
- Introduces ~6% error in inductance
- FastImp does NOT use this approximation

### 2. Point-Matching Approximation

**Location**: `rad_peec_matrices.cpp:225-243`

```cpp
double PEECMatrixBuilder::MutualInductance(const PEECSegment& seg_i,
                                            const PEECSegment& seg_j) const {
    // Neumann formula approximation (point matching)
    // L_ij = (mu_0 / 4*pi) * (d_i . d_j) * l_i * l_j / r_ij

    // ← POINT MATCHING: uses center-to-center distance only
    double r = distance(seg_i.center, seg_j.center);
    return (PEEC_MU_0 * PEEC_INV_FOUR_PI) * dot * seg_i.length * seg_j.length / r;
}
```

**Why it's wrong**:
- Uses center-to-center distance only (point matching)
- Should use segment-to-segment integration
- Loses accuracy for close or parallel segments

---

## Integration Plan

### Phase 1: Document FastImp Algorithm (CURRENT)

- [x] Locate SurfacePanel implementation
- [x] Locate GMD approximation code
- [x] Find FastImp references in codebase
- [ ] **TODO**: Access FastImp paper for exact formulas
- [ ] **TODO**: Document panel-based inductance formulas

### Phase 2: Implement Panel-Based Inductance

1. **Replace `SelfInductance()` function**:
   - Remove GMD approximation
   - Implement proper Neumann formula with segment integration
   - Handle rectangular cross-section accurately (no circular conversion)

2. **Replace `MutualInductance()` function**:
   - Remove point-matching approximation
   - Implement segment-to-segment integration (Neumann formula)
   - Use Gauss quadrature for accuracy

### Phase 3: Add Panel Support

1. **Extend `PEECMatrixBuilder` to accept panels**:
   ```cpp
   void AddPanel(const SurfacePanel& panel);
   void ComputePanelP();  // Panel potential coefficient matrix
   void ComputeFilamentPanelM_LS();  // Filament-panel coupling
   ```

2. **Update GMSH import workflow**:
   - Generate 1D centerline mesh (filaments) ← ALREADY WORKING
   - Generate 2D surface mesh (panels) ← NEW
   - Link filaments to panels

### Phase 4: Validation

1. **Test on circular coil**:
   - Compare against analytical formula
   - Target: < 1% error (down from current 6%)

2. **Test against FastImp reference**:
   - If FastImp source code available, run same geometry
   - Compare L, R, P matrices

---

## File Modification Plan

| File | Action | Priority |
|------|--------|----------|
| `rad_peec_matrices.cpp` | Remove GMD approximation (lines 207-223) | **HIGH** |
| `rad_peec_matrices.cpp` | Implement segment integration for L matrix | **HIGH** |
| `rad_peec_matrices.h` | Add panel-based API | MEDIUM |
| `demo_peec_from_1d_mesh.py` | Test new implementation | HIGH |
| `generate_1d_coil_mesh.py` | Add surface panel generation | MEDIUM |

---

## Next Steps

1. **Access FastImp paper** (Zhu et al., IEEE TCAD 2005)
   - Semantic Scholar: https://www.semanticscholar.org/paper/e6b9cadae5ac4a036faec2f00a1215799d636533
   - Extract exact formulas for inductance matrix

2. **Check FastImp GitHub repository**
   - Repository: https://github.com/ediloren/FastImp
   - Look for C++ source code implementing Loop-Star decomposition

3. **Implement segment integration**
   - Replace GMD approximation with proper Neumann formula
   - Add numerical integration (Gauss quadrature)

---

## References

1. Z. Zhu, B. Song, and J. White, "Algorithms in FastImp: A Fast and Wideband Impedance Extraction Program for Complicated 3-D Geometries," IEEE Trans. Computer-Aided Design, vol. 24, no. 7, pp. 981-998, July 2005.

2. G. Vecchi, "Loop-Star Decomposition of Basis Functions in the Discretization of the EFIE," IEEE Trans. Antennas and Propagation, vol. 47, no. 2, pp. 339-346, Feb. 1999.

3. A. E. Ruehli, "Equivalent Circuit Models for Three-Dimensional Multiconductor Systems," IEEE Trans. Microwave Theory and Techniques, vol. 22, no. 3, pp. 216-221, Mar. 1974.

---

**Status**: Phase 1 (Documentation) - In Progress
**Next Action**: Access FastImp paper for exact inductance formulas
