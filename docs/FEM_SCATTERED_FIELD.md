# FEM Scattered-Field Formulation for SIBC

## Overview

The scattered-field formulation decomposes the magnetic vector potential into
incident and scattered components:

```
A = A_inc + A_scat
```

- **A_inc** is the incident field, known analytically (uniform field problems)
  or computed from a free-space FEM solve (coil problems).
- **A_scat** is the scattered field, solved with a Robin boundary condition
  on the SIBC surface.

This decomposition avoids the systematic H_t extraction error that plagues
the total-field formulation when surface impedance boundary conditions are
applied.

## Formulation

For a conducting body with surface impedance Z_s in an incident field:

### Bilinear form (LHS)

```
int nu curl(A_scat) . curl(v) dx + (jw/Z_s) int A_scat . v ds(sibc) = RHS
```

where nu = 1/mu is the reluctivity and w = 2*pi*f is the angular frequency.

### Right-hand side

**Uniform field (sphere problem):** A_inc is an exact CoefficientFunction,
and volume terms cancel:

```
RHS = -(jw/Z_s) * <A_inc, v>_sibc - <n_cond x H_inc, v>_sibc
```

**Coil problem (two-step solve):** A_inc is obtained from a free-space FEM
solve (Step 1), then:

```
RHS = -robin * <A_inc, v>_sibc
```

The volume terms `int J . v dx - int nu B_inc . curl(v) dx` cancel because
A_inc satisfies the source PDE in the free-space domain.

## H_t Extraction

The tangential magnetic field intensity on the SIBC surface is extracted from
the total potential:

```
A_total = A_inc (CF) + A_scat (GridFunction)
H_t = |jw / Z_s| * A_t_rms
```

where A_t_rms is the RMS tangential component of A_total on the SIBC surface.

A_inc dominates the tangential component, so numerical error in A_scat does
not corrupt the H_t result. This is the key advantage of the scattered-field
decomposition.

## Why Total-Field Fails for H_t

The total-field formulation solves for A directly with a Robin boundary
condition. This approach gives correct inductance (volume integral) but
systematically wrong H_t values:

1. The Robin coefficient `|jw/Z_s|` is of order 10^9, while the curl-curl
   coefficient nu_0 is of order 10^6. This ratio mismatch causes the BND
   trace evaluation to return values dominated by the workpiece interior
   (where A approaches 0) rather than the exterior.

2. H_t is systematically -34% regardless of frequency or material parameters.

3. Inductance L (computed from volume integrals) is NOT affected by this
   issue -- only surface quantities like H_t are corrupted.

## Validation: Smythe Sphere

Conducting sphere R = 10 mm in uniform field B0, with Dowell Z_s.
H_ana is the analytical peak tangential H field (A/m) from Smythe.
Errors are relative to H_ana.

```
Mat      f[Hz]   R/d      H_ana     BEM    FEM-s    FEM-t
-----------------------------------------------------------------
copper      100   1.5     689.05  -40.2%   +0.1%  -34.3%
copper      500   3.4     842.34  -23.3%   +0.9%  -33.8%
copper     1000   4.8     878.48  -17.8%   +1.1%  -33.7%
copper     5000  10.7     930.18   -8.7%   +1.4%  -33.5%
copper     7000  12.7     936.92   -7.4%   +1.5%  -33.5%
copper    10000  15.1     942.96   -6.2%   +1.6%  -33.4%
copper    50000  33.8     960.33   -2.7%   +2.0%  -33.4%

 steel      100   2.8      37.87  -66.1%   -2.9%  -36.3%
 steel      500   6.3      81.34  -65.3%   -2.8%  -36.2%
 steel     1000   8.9     112.11  -64.7%   -2.7%  -36.1%
 steel     5000  19.9     225.39  -62.2%   -2.2%  -35.8%
 steel     7000  23.5     257.73  -61.3%   -2.1%  -35.7%
 steel    10000  28.1     295.32  -60.3%   -1.9%  -35.6%
 steel    50000  62.8     496.19  -52.7%   -1.0%  -35.0%
```

FEM-scattered achieves +/-3% accuracy across all conditions (R/delta from
1.5 to 62.8). The total-field formulation shows a consistent -34% bias.

BEM errors are large at low R/delta because the thin-skin assumption
underlying SIBC breaks down; this is expected and not a formulation error.

## Implementation

### Coil problem: `calc_fem_kelvin.py --formulation scattered`

Uses a two-step pardiso solver path:

1. **Step 1:** Solve the free-space problem (no Robin BC) to obtain A_inc
   as a GridFunction.
2. **Step 2:** Solve with Robin BC on the SIBC surface. The RHS is
   `-robin * <A_inc, v>_sibc`. This yields A_scat.
3. **Reconstruct:** `A_total = A_inc + A_scat`.

### Sphere problem: `sphere_uniform_field.py`

A_inc and H_inc are exact CoefficientFunctions (no discretization error):

```python
A_inc = CF((-B0/2 * y,  B0/2 * x,  0))
nxH_inc = CF((y * H0/R,  -x * H0/R,  0))
```

The single-step solve uses both terms in the RHS:

```
RHS = -(jw/Z_s) * <A_inc, v>_sibc - <n_cond x H_inc, v>_sibc
```

## EMMaterial Integration

All material properties (sigma, mu_r, BH curve) are managed by the
`EMMaterial` class in `src/radia/em_material.py`.

- **Presets:** steel (sigma = 2e6 S/m, mu_r = 100), copper (sigma = 5.8e7
  S/m, mu_r = 1), aluminum (sigma = 3.5e7 S/m, mu_r = 1).
- **Surface impedance:** `mat.dowell_Zs(frequency, R)` computes the Dowell
  surface impedance for a given frequency and local radius of curvature R.
- **ESIM solver:** `mat.create_esim_solver(frequency, half_thickness)` creates
  an Effective Surface Impedance Method solver for nonlinear or
  frequency-dependent materials.

## References

- Smythe, W.R., *Static and Dynamic Electricity*, 3rd ed., McGraw-Hill, 1968.
- Dowell, P.L., "Effects of eddy currents in transformer windings,"
  Proc. IEE, vol. 113, no. 8, pp. 1387-1394, 1966.
