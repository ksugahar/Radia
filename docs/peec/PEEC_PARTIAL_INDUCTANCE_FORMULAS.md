# PEEC Partial Inductance Formulas: Current State and Improvement Plan

**Date**: 2026-04-17
**Status**: Survey complete, implementation planned

---

## 1. Problem Statement

The current Radia PEEC implementation uses **filamentary approximation** for
partial inductance: each conductor segment is represented by its centerline,
and inductance is computed via the Neumann line integral.  This has two
known limitations:

1. **Self-inductance of short segments**: The Grover approximation
   `L = (mu_0/2pi) * l * [ln(2l/r) - 3/4]` has > 4% error at `h/r < 10`
   (Giussani et al. 2022, Fig. 5).  **Fixed 2026-04-17** by replacing with
   the exact Neumann formula (commit `3494fc16`).

2. **nwinc subdivision and circulating currents**: When `nwinc >= 3`,
   sub-filaments interact via filamentary mutual inductance, which does not
   account for the finite cross-section of the conductor.  At DC, this
   produces spurious circulating current patterns that reduce the apparent
   inductance by ~5% compared to Grover (which assumes uniform current).

   The root cause: filamentary Neumann mutual inductance between two
   parallel sub-filaments at distance `d` gives `M ~ ln(2l/d)`, but for
   finite-cross-section conductors the correct expression involves a 6-fold
   volume integral over both cross-sections (Ruehli 1972, eq. 6).

## 2. Current Implementation

| Quantity | Formula | Reference | Accuracy |
|----------|---------|-----------|----------|
| **Self L (circular)** | Exact Neumann: `l*asinh(l/r) - sqrt(l^2+r^2) + r + l/4` | Grover 1946, Aebischer 2014 [2] | Exact for all l/r |
| **Self L (rectangular)** | Rosa/Grover 6-term with atan corrections | Rosa 1908, Grover 1946 | Exact (6-fold closed-form) |
| **Mutual L (parallel)** | Filamentary Grover: `F(alpha,d) + F(beta,d) - ...` | Grover 1946 | Exact for **filaments**, approximate for finite cross-section |
| **Mutual L (general)** | 8-point Gauss-Legendre quadrature on Neumann integral | — | Numerical, ~0.1% for well-separated segments |
| **Mutual L (perpendicular)** | Returns 0 | Neumann integral | Exact |

### Verification (2026-04-17)

Circular torus (R = 30 mm, a = 3 mm):

| Method | L (nH) | vs Grover L_dc |
|--------|--------|----------------|
| Grover analytical | 99.23 | — |
| **3D FEM Kelvin (p=2, 1M DOF)** | **99.44** | **+0.21%** |
| PEEC nwinc=1 | 100.01 | +0.79% |
| PEEC nwinc=3 | 94.64 | -4.62% |
| PEEC nwinc=7 | 94.50 | -4.76% |

**FEM p=2 confirms Grover is accurate** for this geometry (R/a = 10).
PEEC nwinc=1 agrees well; nwinc >= 3 shows the circulating current artifact.

## 3. Improvement Plan (from Literature Survey)

### Phase 1: Exact Partial Mutual Inductance for Rectangular Cross-Sections

Replace the filamentary Neumann mutual inductance between sub-filaments
with the **exact 6-fold integral closed-form** from Ruehli 1972 and
Hoer & Love 1965.

For two parallel rectangular conductors of width `w`, height `h`,
length `l`, separated by distance `d`:

```
L_pp = (mu_0 / (4*pi*w1*h1*w2*h2)) * integral over both cross-sections
       of F(l, distance(p1, p2)) dp1 dp2
```

where `F(l, d) = l * asinh(l/d) - sqrt(l^2 + d^2) + d` is the
filamentary Neumann kernel.

The 6-fold integral reduces to a **closed-form sum of 36 terms**
involving `asinh`, `atan`, and algebraic functions (Hoer & Love 1965,
eq. 14-18; Piatek & Baron 2012, 2013).

**Impact**: With exact partial mutual L, `nwinc` subdivision becomes
physically correct at DC (no spurious circulating currents).  The
sub-filament mutual inductance correctly accounts for the averaging
over the finite cross-section area.

### Phase 2: GMD-Based Skin Effect Correction

Replace the DC internal inductance term (`l/4`) with a
frequency-dependent GMD correction (Aebischer 2017):

```
L(f) = L_ext(GMD_ext) + L_int(f)
```

where `GMD_ext` is the geometric mean distance of the cross-section
and `L_int(f)` decreases from `mu_0*l/(8*pi)` at DC to 0 at high
frequency as the skin effect concentrates current on the surface.

**Impact**: Correct frequency-dependent inductance without `nwinc`
subdivision.  `nwinc` remains useful for proximity effect between
adjacent conductors, but self-inductance frequency dependence is
handled analytically.

### Phase 3: Non-Parallel Mutual Inductance (Exact)

Replace 8-point Gauss quadrature for non-parallel filaments with
the analytical formula for skewed filaments (Grover 1946, Chapter 5;
Kalantarov-Zeitlin for circular filaments).

**Impact**: Removes numerical integration error for helical coils
where adjacent-turn segments are slightly non-parallel.

## 4. References

[1] A. E. Ruehli, "Inductance Calculations in a Complex Integrated
    Circuit Environment," IBM J. Res. Dev., vol. 16, no. 5,
    pp. 470-481, Sep. 1972. DOI: 10.1147/rd.165.0470

[2] H. A. Aebischer and B. Aebischer, "Improved Formulae for the
    Inductance of Straight Wires," Advanced Electromagnetics, vol. 3,
    no. 1, pp. 31-43, 2014. DOI: 10.7716/aem.v3i1.254

[3] H. A. Aebischer and B. Aebischer, "The GMD Method for Inductance
    Calculation Applied to Conductors with Skin Effect," Advanced
    Electromagnetics, vol. 6, no. 2, pp. 77-92, Sep. 2017.
    DOI: 10.7716/aem.v6i2.464

[4] C. Hoer and C. Love, "Exact Inductance Equations for Rectangular
    Conductors with Applications to More Complicated Geometries,"
    J. Res. Nat. Bur. Stand. C, vol. 69C, no. 2, pp. 127-137, 1965.

[5] Z. Piatek and B. Baron, "Exact Closed Form Formula for Self
    Inductance of Conductor of Rectangular Cross Section," Progress
    In Electromagnetics Research M, vol. 26, pp. 225-236, 2012.

[6] Z. Piatek, B. Baron, T. Szczegielniak, and A. Pasierbek, "Exact
    Closed Form Formula for Mutual Inductance of Conductors of
    Rectangular Cross Section," Przeglad Elektrotechniczny, vol. 89,
    no. 3a, pp. 61-64, 2013.

[7] L. M. Silveira, M. Kamon, and J. White, "Efficient Reduced-Order
    Modeling of Frequency-Dependent Coupling Inductances Associated
    with 3-D Interconnect Structures," IEEE Trans. Compon. Packag.
    Manuf. Technol. B, vol. 19, no. 2, pp. 283-288, May 1996.
    DOI: 10.1109/96.496030

[8] L. Giussani, L. Di Rienzo, M. Bechis, P. Cambareri, and
    C. de Falco, "Efficient PEEC Computation of Losses and Currents
    in Shields of Round Wires in Submarine Tripolar Cables," IEEE
    Trans. Magn., vol. 58, no. 9, Sep. 2022.
    DOI: 10.1109/TMAG.2022.3161840

[9] F. W. Grover, "Inductance Calculations: Working Formulas and
    Tables," Dover Publications, New York, 1946 (reprint of 1944
    original).

## 5. Local Paper Archive

All papers are stored in:
```
W:\03_文献・論文\00_電磁界解析\PEEC法\
```

| # | File | Ref |
|---|------|-----|
| 1 | `Inductance_Calculations_in_a_Complex_Integrated_Circuit_Environment.pdf` | [1] |
| 2 | `Improved_Formulae_for_the_Inductance_of_Straight_Wires_Aebischer_2014.pdf` | [2] |
| 3 | `The_GMD_Method_for_Inductance_Calculation_...Skin_Effect_Aebischer_2017.pdf` | [3] |
| 4 | `Exact_inductance_equations_for_rectangular_conductors_...Hoer_Love_1965.pdf` | [4] |
| 5 | `Exact_Closed_Form_Formula_for_Self_Inductance_...Piatek_Baron_2012.pdf` | [5] |
| 6 | `Efficient_reduced-order_modeling_...3-D_interconnect_structures.pdf` | [7] |
| 7 | `Efficient_PEEC_Computation_...Submarine_Tripolar_Cables.pdf` | [8] |
