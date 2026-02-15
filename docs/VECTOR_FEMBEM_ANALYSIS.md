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

## Required Fix

### Fix 1: Correct the RHS (minor, code correctness)

In `VectorEddyCurrentFEMBEM.solve()`, change line 2138 from:

```python
rhs[:n1] = -1j * self.omega * self.sigma * self._a_mass_np @ A_inc_coeffs
```

to:

```python
rhs[:n1] = -(self._a_curl_np
             + 1j * self.omega * self.sigma * self._a_mass_np
             ) @ A_inc_coeffs
```

### Fix 2: Alternative loss computation from surface current (major)

Instead of the volume loss formula, compute loss from the BEM surface current j
using the SIBC-like formula:

```python
P = 0.5 * Re(Zs) * integral_Gamma |j_tangential|^2 dS
```

This requires:
1. HDivSurface mass matrix `M_hdiv` (boundary L2 inner product)
2. Surface impedance `Zs = (1+j)/(sigma*delta)`

```python
def compute_loss_from_surface_current(self):
    """Loss from surface current (works on coarse mesh)."""
    Zs = (1 + 1j) / (self.sigma * self.delta)
    P = 0.5 * Zs.real * np.real(self._j_coeffs.conj() @ self._M_hdiv_np @ self._j_coeffs)
    return P
```

This would make VectorFEMBEM loss mesh-independent (like ShieldBEMSIBC) while
still using the 3-block FEM-BEM system to compute the surface current j correctly.

**However**, this requires that j from the FEM-BEM system correctly represents
the physical surface current. Current diagnostics show |j| ~ 1/omega
(decreasing), which is WRONG (should be ~constant). This is because the BEM
equation `S*[j;rho] = B*A_total` has A_total -> 0, so j -> 0.

### Fix 3: Reformulate to solve for A_total directly (fundamental)

Instead of the scattered field formulation (A_s = A_total - A_inc), reformulate
to solve for A_total directly with proper boundary conditions. This avoids the
catastrophic cancellation A_total = A_scat + A_inc ~ 0.

## Recommendations

1. **For thin-skin problems (delta << thickness)**: Use **ShieldBEMSIBC**.
   It handles skin depth analytically via SIBC, is mesh-independent, and
   is validated to work for mu_r >= 1.

2. **For moderate-skin problems (delta ~ thickness)**: VectorFEMBEM is
   appropriate only if the mesh resolves the skin depth. Use boundary
   layer meshing (fine elements near surface).

3. **Fix priority**: Fix 1 (curl-curl RHS) is a simple code correction.
   Fix 2 (surface current loss) requires rethinking the BEM coupling.
   Fix 3 (total field formulation) is the most robust but requires
   significant refactoring.

## Diagnostic Scripts

- `examples/ngbem_diagnostics/diagnose_vector_fembem.py` - Full diagnostic
- `examples/ngbem_diagnostics/validate_shield_vs_vector.py` - Cross-validation

## References

- Johnson-Nedelec FEM-BEM coupling: C. Johnson and J.C. Nedelec, "On the
  coupling of boundary integral and finite element methods", Math. Comp. 35 (1980)
- Weggler stabilization: L. Weggler, "High order boundary element methods",
  PhD thesis, Saarland University (2011)
