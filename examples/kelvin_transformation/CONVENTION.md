# Kelvin Transformation Convention (canonical)

This document declares the **single canonical convention** used by all
Kelvin-transformed examples and the centralized `radia.kelvin_source`
API in this repository.

## Authoritative reference

> **H. Nagamine, T. Yamaguchi, K. Sugahara**,
> "A Pullback-Based Formulation of Kelvin Transformation in
> Electromagnetic Field Analysis," CEFC 2026 (Thessaloniki), id 350.
>
> Derives the material transformation law for vector fields via
> **pullbacks of differential forms** and **bilinear energy functionals**
> (not a metric-based approach with ad-hoc scale factors). Validated
> numerically on a toroidal current loop (analytical exterior energy
> 3.333×10⁻⁸ J vs FEM on transformed domain 3.344×10⁻⁸ J, +0.33%).

See also:

> **K. Sugahara**, "Electromagnetic analysis of eddy current testing
> with Kelvin transformation," IEEE Trans. Magn. 58(9), 1–6, Sept. 2022
> (ref [3] in the Nagamine paper).

## Canonical rule

For **3D spherical Kelvin inversion** `k(r,θ,φ) = (R²/r, θ, φ)`
(conformal, orientation-reversing with `sgn(k) = −1`):

```
nu' = (r'/R)^2 * nu_0           [A-formulation, HCurl]
mu' = (R/r')^2 * mu_0           [Omega / H-formulation, H1 scalar]
```

where `r' = |r' − c_Kelvin|` is the distance in the **computational
(image) frame** from the Kelvin sphere center, and `R` is the Kelvin
sphere radius. The two factors are **inverses** of each other, and the
physical reciprocal relation `μ·ν = 1` holds pointwise — consistent
with Kelvin being a **physical coordinate transformation** of the same
material (not a numerical trick).

### Why these factors (Nagamine CEFC 2026 §II.A)

Pullback of orthonormal 1-form basis (eq. 5):
```
k*(e^{r'}) = −(R/r)² e^r      (same for e^θ, e^φ; conformal)
```

Hodge operator + inner product (eq. 7):
```
⋆(k*(B'))   = −(R/r)² · k*(⋆'(B'))          [B: 2-form]
g(k*(w'), k*(w')) = +(R/r)^4 · g'(w', w')    [w: 1-form]
```

Bilinear integrand (eq. 8):
```
ν ⟨dw, dA⟩ = ν ⟨dw', dA'⟩' · (R/r)^{2+2+4}
           = ν ⟨dw', dA'⟩' · (R/r)^8
```

Volume pullback: `dΩ = k*(−R⁶/r'⁶ dΩ')`, `sgn(k) = −1`.

Equating `W_m^E = ∫_{Ω_E} ν⟨dw,dA⟩ dΩ` with the transformed-domain
energy (eq. 9):

```
W_m^E = sgn(k) ∫_{Ω'} ν ⟨dw', dA'⟩' · (−r'²/R²) dΩ'
      = ∫_{Ω'} ν · (r'/R)^2 ⟨dw', dA'⟩' dΩ'
```

Comparing with `W_m' = ∫_{Ω'} ν' ⟨dw', dA'⟩' dΩ'` gives:

```
nu' = nu * (r'/R)^2
```

Reciprocal for mu: `mu' = mu / (r'/R)² = (R/r')² mu` (same as the
existing Ω-formulation convention).

### Axisymmetric (r, z) — same as 3D spherical

The axisymmetric case is the **same basic 3D spherical Kelvin** viewed
in the meridional (r, z) plane. `ρ' = sqrt(r² + (z - z_offset)²)` is
the 3D spherical distance from the Kelvin sphere center (on the axis
at `z = z_offset`). Use `kelvin_{nu,mu}_factor_axisym_cf` — the
formulas `(ρ'/R)²` and `(R/ρ')²` are identical to the 3D case.

### Cylindrical (2D) Kelvin — non-conformal

For cylindrical inversion `k(ρ,φ,z) = (R²/ρ, φ, z)` (Nagamine §II.B),
the transformation is **non-conformal** and yields an **anisotropic**
reluctivity tensor (eq. 12):

```
ν' = diag(1, 1, (ρ'/R)⁴) · ν
```

Only the axial (`z`) component picks up `(ρ'/R)⁴`; radial and
azimuthal components are **unchanged (factor = 1, identity)**.

- **2D H1/Ω-form** (scalar potential `φ(x,y)` with `grad(φ)` in (x,y)):
  bilinear form uses in-plane μ only → **Kelvin factor IS 1**
  (identity, not "no factor"). Use `kelvin_mu_factor_2d_Hx_Hy_cf()`.
- **2D A-form** (`A = A_z(x,y) ẑ`): bilinear form uses ν_z only →
  scalar factor `(ρ'/R)⁴`. Use `kelvin_nu_factor_2d_Az_cf(offset, R)`.

## Use in code

### FEM examples (NGSolve): use centralized helpers

```python
from radia.kelvin_source import (
    kelvin_material_factor_3d_cf,      # (r'/R)^2 for 3D spherical Kelvin
    kelvin_material_factor_axisym_cf,  # axisym (r,z) with Z-offset (3D sphere)
    kelvin_material_factor_2d_cf,      # 2D Cartesian (2D conformal)
    build_material_cf,                 # {material: value} CF builder
)

# A-formulation, axisym Z-offset (3D spherical):
nu_factor = kelvin_material_factor_axisym_cf(z_offset=0.15, R=0.06)
nu_cf = build_material_cf(
    mesh, nu0, nu_factor,
    overrides={"magnetic": nu0 / mu_r},
)

# Omega-formulation (3D spherical):
mu_factor_inv = kelvin_material_factor_3d_cf(center=(0, 0, 0), R=0.1)
# returns (r'/R)^2. For mu, use 1 / factor_inv = (R/r')^2.
```

### PEEC examples: use the A-evaluation helpers

```python
from radia.kelvin_source import (
    kelvin_map_3d,
    is_in_kelvin_exterior_domain,
    A_s_at_obs_with_kelvin,
)
```

## Historical note: 2026-04-15 / 2026-04-16 session

An earlier session (2026-04-15) performed an "empirical A/B test" in
`A-formulation/test_nu_convention.py` with two labels:

- `sugahara` (mislabeled): `ν = (R/ρ')² ν₀`
- `energy_invariant` (mislabeled): `ν = (ρ'/R)² ν₀`

and reported the first variant wins on a gapped-torus benchmark
(+0.90% vs +4.80% error). **Both labels are wrong**: the actual
Sugahara 2022 convention (per ref [3] in the Nagamine paper, and the
Nagamine derivation itself) is `ν = (ρ'/R)² ν₀`, which was labeled
`energy_invariant` in the test. The 2026-04-16 morning session
(mistakenly) flipped 5 A-formulation files to the wrong direction
based on these mislabeled labels; those fixes were reverted within
the same session upon reading the Nagamine paper.

The test_nu_convention.py label mismatch + the empirical result
(`(R/ρ')²` appearing to win on one coarse-mesh benchmark despite being
mathematically wrong) is filed as an open investigation — see
`examples/kelvin_transformation/docs/pullback_derivation_3D.md` §8.

## References

- **H. Nagamine, T. Yamaguchi, K. Sugahara** (CEFC 2026, id 350):
  pullback + bilinear energy functional derivation.
  Numerical validation: toroidal loop, ν' = (ρ'/R)² ν, FEM on Ω' gives
  exterior magnetic energy matching analytical dipole to +0.33%.
- **K. Sugahara** (IEEE TransMag 58(9), 2022): original reduced
  A-formulation result cited by Nagamine as ref [3].
- `docs/pullback_derivation_3D.md` — derivation from a different angle,
  cross-referenced with Nagamine.
