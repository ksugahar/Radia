# VectorEddyCurrentFEMBEM Cross-Validation Analysis

Date: 2026-02-15

## Summary

Cross-validation of ShieldBEMSIBC vs VectorEddyCurrentFEMBEM revealed that
VectorFEMBEM loss is **frequency-independent** (~4.3 kW at all frequencies).
Root cause analysis identified **two independent bugs** that interact paradoxically.

## Test Setup

- Geometry: 20x20x10 mm aluminum block (sigma = 3.7e7 S/m, mu_r = 1)
- Mesh: maxh = 12 mm (18 volume elements, 94 H(curl) DOFs)
- Excitation: Uniform B_ext = [0, 0, 1] T
- Skin depth: 0.83 mm (10 kHz) to 0.08 mm (1 MHz)

## Bug 1: Missing curl-curl RHS term

**File**: `src/radia/ngbem_eddy.py`, `VectorEddyCurrentFEMBEM.solve()`, line 2138

The interior FEM equation for the scattered field A_s is:

```
curl(1/mu * curl A_s) + j*w*sigma*A_s = -curl(1/mu * curl A_inc) - j*w*sigma*A_inc
```

In weak form:

```
(1/mu * curl A_s, curl v) + jws*(A_s, v) = -(1/mu * curl A_inc, curl v) - jws*(A_inc, v)
```

So the RHS should be:

```
f_1 = -(a_curl + jws*a_mass) @ A_inc = -a_FEM @ A_inc
```

**Current code** (WRONG):
```python
rhs[:n1] = -1j * omega * sigma * self._a_mass_np @ A_inc_coeffs
```

**Missing**: `-self._a_curl_np @ A_inc_coeffs`

### Magnitude of the missing term

| Frequency | |w*s*M*A_inc| | |a_curl*A_inc| | curl/mass ratio |
|-----------|---------------|----------------|-----------------|
| 10 kHz    | 6.17e+05      | 1.40e+04       | 2.28%           |
| 100 kHz   | 6.17e+06      | 1.40e+04       | 0.23%           |
| 1 MHz     | 6.17e+07      | 1.40e+04       | 0.023%          |

The curl-curl term is small relative to the mass term (2.3% at 10 kHz, negligible at
higher frequencies). This is expected for high-conductivity materials where
`omega*sigma*mu >> 1/mu`.

## Bug 2 (Fundamental): Coarse mesh cannot resolve skin layer

Even with the correct RHS, the loss is **exactly zero** (P ~ 1e-26 W):

| Frequency | P_original (buggy RHS) | P_fixed (correct RHS) |
|-----------|------------------------|----------------------|
| 1 kHz     | 2.525e+03 W            | 1.169e-26 W          |
| 10 kHz    | 4.315e+03 W            | 1.040e-24 W          |
| 100 kHz   | 4.366e+03 W            | 6.581e-23 W          |
| 1 MHz     | 4.367e+03 W            | 1.064e-20 W          |

### Why P_fixed = 0

With the correct RHS `f_1 = -a_FEM @ A_inc`, the FEM equation becomes:

```
a_FEM * (A_s + A_inc) + B^T * j = 0
a_FEM * A_total + B^T * j = 0
```

On a coarse mesh (maxh = 12 mm >> delta = 0.08-0.83 mm), the H(curl) basis
cannot represent the exponential boundary layer `exp(-z/delta)`. The "best"
FEM solution is:

```
A_total = 0  (everywhere inside)
j = 0        (surface current)
```

This is physically the **perfect conductor limit** (sigma -> infinity).
The coarse mesh FEM cannot distinguish between sigma = 3.7e7 and sigma = infinity.

### Shielding diagnostic confirms this

| Frequency | |A_total|/|A_inc| | cos(A_scat, A_inc) |
|-----------|-------------------|--------------------|
| 1 kHz     | 3.05e-01          | -0.952             |
| 10 kHz    | 3.92e-02          | -0.999             |
| 100 kHz   | 3.94e-03          | -1.000             |
| 1 MHz     | 3.94e-04          | -1.000             |

A_scat -> -A_inc (perfect cancellation), A_total -> 0.

### Why P_original ~ const (the paradox)

Without the curl-curl RHS, the FEM equation is:

```
(a_curl + jws*M) * A_s + B^T * j = -jws * M * A_inc
```

At high omega, `jws*M` dominates, so `A_s ~ -A_inc`. But the curl-curl part
`a_curl * A_s ~ -a_curl * A_inc` creates an unbalanced residual:

```
A_total ~ a_curl * A_inc / (jws * M)
|A_total|^2 ~ |a_curl*A_inc|^2 / (w^2 * s^2 * ||M||^2)
```

Loss becomes:
```
P = 0.5 * w^2 * s * |A_total|^2 ~ 0.5 * |a_curl*A_inc|^2 / (s * ||M||^2)
```

This is **frequency-independent** (the w^2 in the loss formula exactly cancels
the 1/w^2 in |A_total|^2). The "nonzero loss" is entirely a numerical artifact
from the incomplete RHS.

## Root Cause Interaction

| Code state       | A_total            | Loss              | Correct? |
|------------------|--------------------|-------------------|----------|
| Current (buggy)  | Small, ~1/omega    | ~4.3 kW (const)   | NO (artifact from RHS imbalance) |
| Fixed (curl-curl)| Zero (to precision)| ~0 W              | NO (coarse mesh = perfect conductor) |
| Physical truth   | exp(-z/delta) layer| ~sqrt(omega)       | YES (needs delta-scale mesh) |

**Neither code state gives correct loss on a coarse mesh.** The volume FEM loss
formula `P = 0.5*w^2*s*|A_total|^2` fundamentally requires mesh elements finer
than the skin depth to resolve the boundary layer where the physical loss occurs.

## Fixes Applied

### Fix 1: Total field formulation (DONE)

Switched from scattered field (A_scat as unknown) to total field (A_total as
unknown) to avoid catastrophic cancellation:

**Old (scattered)**: A_scat + A_inc ~ 0 (both large, nearly cancel)
**New (total)**: A_total is direct unknown (no cancellation)

```python
# Old: Row 1 RHS = -jws*M*A_inc (buggy: missing curl-curl)
# New: Row 1 RHS = 0 (homogeneous interior)
# New: Row 2 RHS = +B @ A_inc (BEM drives system)
```

Result: A_total -> 0 (PEC limit on coarse mesh), j is frequency-independent
and nonzero (correct for PEC surface current).

### Fix 2: Surface current j cannot be used for SIBC loss (INVESTIGATED)

Investigation showed that j in the Johnson-Nedelec FEM-BEM coupling is an
auxiliary SLP representation density, NOT the physical surface current
K = n x H. The normalization differs by ~10^10 from the physical current.

**Why**: In VectorFEMBEM, the BEM equation relates j to the boundary trace of
the vector potential A (not to H or E). ShieldBEMSIBC uses the EFIE which
directly solves for the physical surface current K.

### Fix 3: Analytical SIBC loss estimate (DONE)

Added `compute_loss_sibc()` method that computes loss from the known incident
field H_inc on each boundary triangle:

```python
P = sum_faces 0.5 * Re(Zs) * |H_inc_tangential|^2 * area
```

This is the half-space SIBC approximation applied face-by-face.

### Cross-validation: compute_loss_sibc() vs ShieldBEMSIBC

| Frequency | P_analSIBC (W) | P_shield (W) | Ratio |
|-----------|----------------|--------------|-------|
| 1 kHz     | 2,617          | 678          | 3.87  |
| 10 kHz    | 8,274          | 9,005        | 0.92  |
| 100 kHz   | 26,165         | 19,604       | 1.33  |
| 1 MHz     | 82,741         | 50,796       | 1.63  |

The analytical SIBC gives the correct frequency scaling (sqrt(f)) and is
within a factor of 0.9-3.9x of ShieldBEMSIBC. The deviation comes from:
- Low freq (delta ~ thickness): half-space assumption breaks down
- High freq (delta << thickness): no edge/corner enhancement in flat approx

## Solver Selection Guide

| Regime | delta vs thickness | Recommended Solver | Notes |
|--------|--------------------|--------------------|-------|
| Thick skin | delta > thickness | VectorFEMBEM (fine mesh) | Need mesh resolution < delta |
| Moderate | delta ~ thickness | VectorFEMBEM (fine mesh) | Boundary layer meshing |
| Thin skin | delta << thickness | **ShieldBEMSIBC** | Mesh-independent, most accurate |
| Quick estimate | Any | VectorFEMBEM.compute_loss_sibc() | 1-4x of ShieldBEMSIBC |

**VectorFEMBEM is most useful** for:
1. Thick-skin problems where skin depth is resolvable by the mesh
2. Magnetic materials (mu_r >> 1) where ShieldBEMSIBC is less accurate
3. Computing A_total and H fields inside the conductor (not just loss)

**ShieldBEMSIBC is preferred** for:
1. Thin-skin shielding analysis (delta << thickness)
2. Loss computation on coarse meshes
3. Quick frequency sweeps (surface-only, no volume DOFs)

## Diagnostic Scripts

- `examples/ngbem_diagnostics/diagnose_vector_fembem.py` - Root cause diagnostic
- `examples/ngbem_diagnostics/validate_shield_vs_vector.py` - Cross-validation
- `examples/ngbem_diagnostics/test_sibc_loss.py` - SIBC loss method comparison

## References

- Johnson-Nedelec FEM-BEM coupling: C. Johnson and J.C. Nedelec, "On the
  coupling of boundary integral and finite element methods", Math. Comp. 35 (1980)
- Weggler stabilization: L. Weggler, "High order boundary element methods",
  PhD thesis, Saarland University (2011)
