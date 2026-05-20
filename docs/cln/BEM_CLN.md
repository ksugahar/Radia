# BEM-CLN: Per-Element Multipole CLN with Schur-F Termination

A multi-element extension of the single-conductor Schur-F augmented
Cauer Ladder Network (CLN) representation, suitable for assembling
fast eddy-current ROMs of small clusters of conductors using natural
open-boundary integral-equation coupling.

This document collects the polarizability-based assembly framework
that backs the IEEE Transactions on Magnetics submission "Universal
Cauer-SIBC Composition via Schur Complement" (Sugahara, Nagamine,
Hane, 2026), §V.G and §V.H.

## 1. Motivation

The single-conductor Schur-F augmentation
```
Y_R(s) = Y_CLN(s) + K_SIBC sqrt(s) / (s + d)
K_SIBC  = S sqrt(sigma / mu)
```
preserves the DC value, the canonical Cauer Taylor structure at
`s = 0`, and the non-rational SIBC asymptote `Y -> K_SIBC / sqrt(s)`
at `s -> infinity` for any single conductor.  The surface-area-only
scaling makes the augmentation parameter readable directly from CAD
metadata (Theorem 1 of the cited manuscript).

For systems with multiple conductors --- induction-heating workpieces
+ coils, paired transformer windings, accelerator-magnet yokes near
conductors --- the natural extension is BEM-style assembly: each
conductor is a CLN-Schur-F element, with inter-element coupling
delivered by integral-equation kernels.  This file documents that
extension.

## 2. Polarizability formulation (the right per-element variable)

The IGTE digest uses the cylinder admittance
`Y_cyl(s) = pi a^2 sigma * 2 I_1(gamma a) / (gamma a I_0(gamma a))`
with the DC value `Y_cyl(0) = pi a^2 sigma`.  For BEM-style coupling,
however, the natural per-element variable is the **polarizability**
```
alpha(s) = a^2 [1 - 2 I_1(gamma a) / (gamma a I_0(gamma a))]    (2D cylinder)
        = a^2 - Y_cyl(s) / (pi sigma).
```
This has the structurally correct limits **built in**:

| Limit       | Value      | Interpretation              |
|-------------|------------|------------------------------|
| `alpha(0)`        | `0`             | no eddy currents at DC      |
| `alpha(infinity)` | `a^2`           | full PEC diamagnetic shield |
| `|alpha(s)|`      | `<= a^2`        | bounded for all s           |

The bounded `|alpha| <= a^2` (or `<= V` in 3D) means the inter-element
coupling needs no phenomenological saturation factor --- a major
simplification over earlier prototype implementations.

The 3D analog for a cuboid of volume `V = a b c`:
```
alpha_3D(s) = V - Y_cuboid(s) / sigma
alpha_3D(0) = 0
alpha_3D(infinity) = V
```
where `Y_cuboid(s)` is the scalar-diffusion admittance of §V.B of the
manuscript.

The Schur-F-augmented per-element building block uses
`Y_cuboid_S(s) = Y_CLN_N(s) + K_SIBC^(cuboid) sqrt(s) / (s + d)` and
gives `alpha_S(s) = V - Y_cuboid_S(s) / sigma` automatically.

## 3. Element-to-element coupling (Green's function)

Each polarized element of strength `c_n = alpha(s) B^local_n`
produces a perturbation of the applied field at every other element
through the appropriate Green's function.

### 2D long-cylinder coupling (cross-section in x-y plane, B applied in y)
A 2D dipole `c cos(theta) / r` at the origin produces
`Delta B_y = c / D^2` at position `(D, 0)`.  The BEM-CLN system is
```
(I - alpha(s) G_2D) c = alpha(s) B_0 1
(G_2D)_ij = 1 / D_ij^2   for i != j, 0 on diagonal.
```

### 3D dipole coupling (cuboid pair, B applied in y, separation in x)
A 3D dipole `m_y` at origin produces `Delta B_y = -mu_0 m_y / (4 pi D^3)`
at `(D, 0, 0)`.  Combined with `c_n = alpha(s) B^local_n` and
`m = c` (same dimensions), the system reads
```
(I + alpha(s) G_3D) c = alpha(s) B_0 1
(G_3D)_ij = mu_0 / (4 pi D_ij^3)   for i != j, 0 on diagonal.
```
**Note the sign**: `+alpha G_3D` (opposite to 2D's `-alpha G_2D`),
reflecting the destructive `(3 cos^2(theta) - 1)` angular factor of
3D dipoles for moments perpendicular to the line connecting them.

### Stability
For both 2D and 3D, the assembly is well-posed provided
`alpha(s) lambda_max(G) < 1` (2D) or `alpha(s) lambda_max(G) > -1`
(3D).  For typical geometries with conductor spacing larger than
conductor size, both conditions hold automatically.

## 4. DOF accounting

Per element (one cylinder or cuboid):
```
DOF = N_Cauer + 1     (Cauer rungs + 1 Schur-F sqrt(s) block)
```
For N elements:
```
Total DOF = N (N_Cauer + 1)
Coupling matrix: N x N, dense (O(N^2) non-zeros) but small
```

Comparison with alternative methods:

| Method                  | DOF / element | Regime                |
|-------------------------|---------------|------------------------|
| MMM tetrahedron         | 3             | steady-state only      |
| MSC hexahedron          | 6             | steady-state only      |
| FEM-element CLN         | ~60 N_Cauer   | transient, bounded     |
| **BEM-CLN dipole+quad+Schur-F** | **8 (N_Cauer + 1)** | **transient, inhomogeneous + SIBC** |
| High-order curved tet (p=3) | ~60 + PML/Kelvin | transient, all regimes |

## 5. Verification status (Phase 1, 2, 2.5, 3 — May 2026)

| Phase | What was verified | Coupling | Status |
|-------|--------------------|----------|--------|
| Phase 1   | 2-cylinder 2D, phenomenological g(s,D) saturation | manual | Working, but Schur/Bessel ratio = 1.04 |
| Phase 2   | N=2,3,5,10 linear chain 2D, same phenomenological | manual | DOF scaling demonstrated |
| **Phase 2.5** | **Rigorous polarizability 2D**, no saturation | `1/D^2` (exact) | **Schur/Bessel = 0.998-0.999** |
| Phase 3 A | 3D quadrupole (Theorem 3) on 5x2x1 mm cuboid | n/a (single conductor) | Foster Mmax=99 r=0.88 |
| Phase 3 B | 3D 2-cuboid pair BEM-CLN | `mu_0/(4 pi D^3)` | Foster/Schur convergence verified |

Outstanding: Phase 3 B requires NGSolve A-formulation full-FEM cross-check
on the 2-cuboid geometry to verify the 3D rigorous coupling against the
ground-truth eddy-current solve.  Real induction-heating workpiece example
with N >= 5 elements is the next step.

## 6. Reference scripts

All Mathematica verification scripts live in
`examples/CLN/scripts/`:

| Script | Phase | What it verifies |
|--------|-------|------------------|
| `bem_cln_2cylinder.wls`           | Phase 1   | 2D 2-cylinder, phenomenological coupling (legacy) |
| `bem_cln_Ncylinder.wls`           | Phase 2   | 2D N-cylinder chain, phenomenological (legacy) |
| `bem_cln_2cylinder_rigorous.wls`  | Phase 2.5 | **Rigorous 2D polarizability, no saturation** |
| `bem_cln_2cuboid_3D.wls`          | Phase 3 B | **3D 2-cuboid BEM-CLN with 3D dipole coupling** |
| `cuboid3D_quadrupole.wls`         | (Th. 3 2D) | u_ext=(y^2-x^2)/2 on long cuboid |
| `cuboid3D_quadrupole_3D.wls`      | (Th. 3 3D) | u_ext=xy on 3D cuboid |

## 7. References

- Sugahara, Nagamine, Hane, "Universal Cauer-SIBC Composition via
  Schur Complement," IEEE Trans. Magn. (submitted), 2026.
- Takahashi, Hiruma, Fujiwara, Imamori, "Time-Domain Homogenization
  of Windings Using B-Input Cauer Ladder Network Method," IEEE Trans.
  Magn. 60(12), 2024.  IEEE Xplore document 10684727.
- Hiruma, Takahashi, Matsuo, "Homogenization Method Based on Cauer
  Ladder Network Representation of Unit Cell," IEEE Trans. Magn.
  60(12), 2024.  IEEE Xplore document 10736669.
- Hiruma, Igarashi, "Homogenization Method Based on Cauer Circuit
  via Unit Cell Approach," IEEE Trans. Magn. 56(2), 1-5, Feb. 2020.
