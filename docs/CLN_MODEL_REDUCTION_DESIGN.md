# CLN (Cauer Ladder Network) Model Order Reduction Design

## Overview

This document describes the design for CLN-based model order reduction in Radia's coupled PEEC-MMM solver.

## 1. Matrix Symmetry Analysis

### 1.1 PEEC Part (Conductor)

The PEEC formulation uses Loop-Star decomposition for low-frequency stability:

```
[Z_LL   Z_LS] [I_L]   [V_L]
[Z_SL   Z_SS] [I_S] = [V_S]
```

**Symmetry of PEEC Blocks:**

| Block | Formula | Symmetric? | Reason |
|-------|---------|------------|--------|
| Z_LL | jω * L_ij | **Yes** | L_ij = L_ji (Neumann formula) |
| Z_SS | 1/(jω) * P_ij | **Yes** | P_ij = P_ji (potential coefficient) |
| Z_LS | jω * M_ij | **Yes** | Loop-Star coupling via mutual inductance |
| Z_SL | Z_LS^T | **Yes** | Reciprocity |

**Conclusion**: The PEEC Loop-Star system is **symmetric**.

### 1.2 MMM Part (Magnetic Material)

The MMM formulation uses the demagnetization tensor:

```
Z_MM * M = χ * H_ext
```

where `Z_MM[i][j] = δ_ij - χ_i * N_ij`

**Symmetry of N_ij (Demagnetization Tensor):**

```
N_ij[k][l] = (V_j / 4π) * (3*r_k*r_l - r²*δ_kl) / r^5
```

Since `r_ij = -r_ji` and the formula is quadratic in r:
- `N_ij^T = N_ji` (transpose equals reverse direction)
- Therefore `N_ij` is **symmetric** for any pair (i,j)

**Conclusion**: The MMM system is **symmetric** when χ is uniform.

### 1.3 Coupling Blocks (PEEC-MMM)

Original (non-symmetric) formulation:
- Z_LM: -jω * Φ_from_M (flux linkage)
- Z_ML: -χ/μ₀ * B_coil (H field from coil)

**Problem**: Z_LM ≠ Z_ML^T due to different units.

**Solution**: Variable scaling M' = sqrt(μ₀ * V) * M

```
Z_LM' = Z_LM / sqrt(μ₀ * V)
Z_ML' = Z_ML * sqrt(μ₀ * V)
```

By reciprocity theorem: **Z_LM' = Z_ML'^T**

**Verified numerically** (2026-01-11):
- Non-symmetric: Z = 4.292682e-03 + j1.926849e-03 Ohm
- Symmetric:     Z = 4.292682e-03 + j1.926849e-03 Ohm
- Difference:    6.5e-19 Ohm (machine precision)

## 2. Full Coupled System Structure

With symmetrization, the complete system becomes:

```
┌────────────────────────────────────────────────┐
│ [Z_LL   Z_LS   Z_LM'  0    ] [I_L ]   [V_L ]  │
│ [Z_LS^T Z_SS   0      Z_SE ] [I_S ] = [V_S ]  │
│ [Z_LM'^T 0     Z_MM   0    ] [M'  ]   [b'  ]  │
│ [0      Z_SE^T 0      Z_EE ] [P   ]   [D   ]  │
└────────────────────────────────────────────────┘
```

All blocks are symmetric or have symmetric transpose relationships.

## 3. CLN Model Reduction Strategy

### 3.1 Step 1: PEEC-only CLN (without magnetic materials)

**Goal**: Extract Cauer ladder circuit from Loop and Star parts separately.

```
PEEC System:
[Z_LL   Z_LS] [I_L]   [V_L]
[Z_LS^T Z_SS] [I_S] = [V_S]
```

**Loop part (Z_LL)**: Represents inductance
- Apply Lanczos algorithm to L matrix
- Extract: L_1 → R_1 → L_2 → R_2 → ... (RL ladder)

**Star part (Z_SS)**: Represents capacitance
- Apply Lanczos algorithm to P matrix
- Extract: C_1 → G_1 → C_2 → G_2 → ... (RC ladder)

**Coupling (Z_LS)**: Mutual inductance between Loop and Star
- Use ACA+ for low-rank approximation
- Represent as transformer coupling

### 3.2 Step 2: ACA+ Low-Rank Approximation

For large-scale problems, use HACApK (already available in Radia):

```cpp
// Existing HACApK interface
rad.SetHACApKParams(eps, leaf_size, eta)
rad.Solve(container, tol, maxiter, 2)  // Method 2 = HACApK
```

**ACA+ for L matrix**:
```
L ≈ U * V^T   (rank k << N)
```

This low-rank approximation enables efficient CLN extraction.

### 3.3 Step 3: Add Magnetic Material Coupling

Once PEEC CLN is established, add MMM coupling:

```
Full System CLN:
┌─────────────────────────────────────────┐
│  PEEC Ladder    ←→    MMM Ladder        │
│  (L-R chain)    Z_LM'  (Reluctance)     │
└─────────────────────────────────────────┘
```

**MMM as magnetic circuit**:
- Elements become reluctances: R_m = 1/(μ₀ * μ_r * A/l)
- N_ij coupling → magnetic transformer coupling

## 4. Implementation Plan

### Phase 1: Matrix Export API

Add APIs to export raw matrices for external CLN extraction:

```python
# Get PEEC matrices
L, R, P = rad.CndGetMatrices(conductor)

# Get MMM matrices
N = rad.GetDemagTensor(magnet)

# Get coupling matrices
Z_LM, Z_ML = rad.CplMagGetCoupling(solver)
```

### Phase 2: ACA+ for PEEC

Apply HACApK to partial inductance matrix:

```python
# Set ACA+ parameters
rad.SetHACApKParams(1e-4, 10, 2.0)

# Get low-rank factors
U, V, rank = rad.CndGetLowRankL(conductor)
```

### Phase 3: CLN Extraction (Lanczos)

Implement block Lanczos for symmetric matrix:

```python
def extract_cln(Z, order):
    """
    Extract Cauer ladder network from symmetric Z matrix.

    Returns: [(L1, R1), (L2, R2), ...] ladder elements
    """
    # Block Lanczos algorithm
    alpha, beta = lanczos(Z, order)

    # Convert to ladder elements
    return continued_fraction_to_ladder(alpha, beta)
```

### Phase 4: Circuit Output

Export to SPICE/Touchstone format:

```python
# SPICE subcircuit
rad.CplMagExportSPICE(solver, "model.cir", order=10)

# Touchstone S-parameters
rad.CplMagExportTouchstone(solver, "model.s1p", freqs)
```

## 5. Physical Interpretation

### 5.1 PEEC → Ladder Circuit

```
Port 1 ──[R1]──[L1]──[R2]──[L2]──...──[Rn]──[Ln]── Port 2
              |            |               |
             ===C1        ===C2           ===Cn
              |            |               |
             GND          GND             GND
```

### 5.2 MMM → Magnetic Circuit

```
Φ ──[Rm1]──[Rm2]──...──[Rmn]── Φ
     ↑       ↑           ↑
   NI₁     NI₂         NI_n  (MMF sources from coil)
```

### 5.3 Coupled System

```
Electrical Domain          Magnetic Domain
──────────────────         ─────────────────
    V                           Φ
    │                           │
   [Z_PEEC]  ←── Z_LM' ──→  [Z_MM]
    │                           │
    I                           M
```

## 6. References

1. G. Vecchi, "Loop-Star decomposition of basis functions", IEEE TAP, 1999
2. Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
3. K. Hollaus, "Effective Surface Impedance Method", 2024
4. A. Odabasioglu et al., "PRIMA: Passive Reduced-order Interconnect Macromodeling Algorithm", IEEE TCAD, 1998

## 7. Current Status (2026-01-11)

| Feature | Status | Notes |
|---------|--------|-------|
| Matrix symmetrization | **Done** | CplMagSetSymmetric() API |
| Symmetry verification | **Done** | test_symmetrization.py |
| PEEC Loop-Star structure | Defined | In rad_peec_mmm_coupled.h |
| ACA+ (HACApK) | Available | For MMM solver |
| CLN extraction | TODO | Need Lanczos implementation |
| SPICE export | TODO | Subcircuit format |

## Appendix A: Symmetrization Proof

**Reciprocity Relation:**

From Biot-Savart law and Maxwell's equations:

```
∫ H_coil · M dV = (1/μ₀) ∫ B_M · I dl
```

In matrix form:
```
Z_LM^T = μ₀ * V * Z_ML
```

**Scaling Transformation:**

Let M' = α * M where α = sqrt(μ₀ * V)

Then:
- Z_LM' = Z_LM / α
- Z_ML' = Z_ML * α = Z_ML * sqrt(μ₀ * V)

Check symmetry:
```
Z_LM'^T = (Z_LM / α)^T = Z_LM^T / α = (μ₀ * V * Z_ML) / sqrt(μ₀ * V)
        = sqrt(μ₀ * V) * Z_ML = Z_ML * α = Z_ML'
```

Therefore: **Z_LM'^T = Z_ML'** (symmetric!)

## Appendix B: Verification Code

```python
import radia as rad

# Create test problem
rad.UtiDelAll()
rad.FldUnits('m')

coil = rad.CndLoop([0,0,0], 0.05, [0,0,1], 'r', 2e-3, 2e-3, 5.8e7, 8, 36)
core = rad.ObjRecMag([0,0,0], [0.03,0.03,0.03], [0,0,0])
mat = rad.MatLin(1000)
rad.MatApl(core, mat)

# Non-symmetric solve
solver1 = rad.CplMagCreate(coil, core)
rad.CplMagSetSymmetric(solver1, 0)  # Non-symmetric
rad.CplMagSetFrequency(solver1, 1000)
rad.CplMagSetVoltage(solver1, 1.0, 0.0)
rad.CplMagSetMu(solver1, 1000, 0)
result1 = rad.CplMagSolve(solver1)

# Symmetric solve
solver2 = rad.CplMagCreate(coil, core)
rad.CplMagSetSymmetric(solver2, 1)  # Symmetric
rad.CplMagSetFrequency(solver2, 1000)
rad.CplMagSetVoltage(solver2, 1.0, 0.0)
rad.CplMagSetMu(solver2, 1000, 0)
result2 = rad.CplMagSolve(solver2)

# Compare
Z1, Z2 = result1['Z'], result2['Z']
print(f"Non-symmetric: Z = {Z1}")
print(f"Symmetric:     Z = {Z2}")
print(f"Difference:    {abs(Z1 - Z2)}")
# Output: Difference ~ 1e-18 (machine precision)
```
