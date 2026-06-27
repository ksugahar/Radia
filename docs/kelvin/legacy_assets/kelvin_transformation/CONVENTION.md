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

### Axisymmetric (r, z) — same as 3D spherical, but mind the r-weight

The axisymmetric case is the **same basic 3D spherical Kelvin** viewed
in the meridional (r, z) plane. `ρ' = sqrt(r² + (z - z_offset)²)` is
the 3D spherical distance from the Kelvin sphere center (on the axis
at `z = z_offset`). Use `kelvin_{nu,mu}_factor_axisym_cf` — the
formulas `(ρ'/R)²` and `(R/ρ')²` are identical to the 3D case.

**NGSolve axisymmetric weight — easy to forget.** NGSolve has no
built-in axisymmetric mode; the `2πr` volume factor must be written
explicitly by the caller as `* r_coord * dx`:

```python
a += nu_cf * grad(u) * grad(v) * r_coord * dx   # correct
a += nu_cf * grad(u) * grad(v)           * dx   # WRONG: off by O(r)
```

The Kelvin helpers handle the **sphere-inversion pullback only**; they
do **not** absorb the r-weight. Dropping `* r_coord` is one of the
most common bugs in axisymmetric Kelvin setups.

### Cylindrical (2D) Kelvin — anisotropic tensor, factor depends on component

For cylindrical inversion `k(ρ,φ,z) = (R²/ρ, φ, z)` (Nagamine §II.B),
the transformation is **non-conformal** and yields an **anisotropic**
reluctivity / permeability tensor (eq. 12):

```
ν' = diag(1, 1, (ρ'/R)⁴) · ν              [reciprocal for μ']
μ' = diag(1, 1, (R/ρ')⁴) · μ
```

Only the axial (`z`) slot is modulated; in-plane slots are **identity**.
Which slot enters the bilinear form depends on **which B / H components
the problem uses**, not on whether it's "A-form" vs "Ω-form":

- **In-plane case (the common one)** — Kelvin factor = **1 (identity)**.
  Covers both
  * 2D Ω-form with scalar `φ(x,y)` (H in-plane), and
  * 2D A-form with `A = A_z(x,y) ẑ` (B in-plane, `B_z = 0`).
  
  Use `kelvin_factor_2d_inplane_cf()`. Deprecated aliases
  `kelvin_mu_factor_2d_Hx_Hy_cf()`, `kelvin_nu_factor_2d_Hx_Hy_cf()`, and
  `kelvin_nu_factor_2d_Az_cf(...)` all resolve to this (the `Az` name
  was misleading — the A_z scalar case still has in-plane B, factor = 1).

- **Axial case** — Kelvin factor = `(ρ'/R)⁴` (for `ν_zz`) or `(R/ρ')⁴`
  (for `μ_zz`). Covers the rare situation where B or H genuinely has a
  z-component, e.g. 2D A-form with in-plane `A = (A_x, A_y, 0)` giving
  `curl(A) = B_z ẑ`. Use `kelvin_{nu,mu}_factor_2d_axial_cf(offset, R)`.

## Use in code

### FEM examples (NGSolve): use centralized helpers

Most Radia-NGSolve problems are **3D**, so the isotropic 3D helpers
below cover the primary use case. Axisymmetric and 2D are reductions.

```python
from radia.kelvin_source import (
    # 3D (primary case, isotropic):
    kelvin_nu_factor_3d_cf,        # (rho'/R)^2 -- A-formulation, HCurl
    kelvin_mu_factor_3d_cf,        # (R/rho')^2 -- Omega / H-formulation, H1

    # Axisymmetric (r,z); 3D sphere viewed in the meridional plane:
    kelvin_nu_factor_axisym_cf,
    kelvin_mu_factor_axisym_cf,

    # 2D cylindrical (anisotropic tensor, see section above):
    kelvin_factor_2d_inplane_cf,        # = 1; common case
    kelvin_nu_factor_2d_axial_cf,       # (rho'/R)^4; rare (B_z ≠ 0)
    kelvin_mu_factor_2d_axial_cf,       # (R/rho')^4; rare (H_z ≠ 0)

    build_material_cf,                   # {material: value} CF builder
)

# 3D A-formulation:
nu_cf = build_material_cf(
    mesh, nu0,
    kelvin_nu_factor_3d_cf(center=(0, 0, 0), R=0.1),
    overrides={"magnetic": nu0 / mu_r},
)

# 3D Omega-formulation (reciprocal):
mu_cf = build_material_cf(
    mesh, mu0,
    kelvin_mu_factor_3d_cf(center=(0, 0, 0), R=0.1),
    overrides={"magnetic": mu_r * mu0},
)

# Axisymmetric A-formulation -- remember the r-weight!
from ngsolve import x as r_coord
nu_cf = build_material_cf(
    mesh, nu0,
    kelvin_nu_factor_axisym_cf(z_offset=0.15, R=0.06),
    overrides={"magnetic": nu0 / mu_r},
)
a += nu_cf * grad(u) * grad(v) * r_coord * dx   # * r_coord is mandatory
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
`docs/kelvin/KELVIN_TRANSFORMATION.md` §2 (and §9 references).

## References

- **H. Nagamine, T. Yamaguchi, K. Sugahara** (CEFC 2026, id 350):
  pullback + bilinear energy functional derivation.
  Numerical validation: toroidal loop, ν' = (ρ'/R)² ν, FEM on Ω' gives
  exterior magnetic energy matching analytical dipole to +0.33%.
- **K. Sugahara** (IEEE TransMag 58(9), 2022): original reduced
  A-formulation result cited by Nagamine as ref [3].
- `docs/pullback_derivation_3D.md` — derivation from a different angle,
  cross-referenced with Nagamine.
