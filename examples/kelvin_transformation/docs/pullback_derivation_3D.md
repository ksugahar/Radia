# 3D Kelvin Pullback of Differential Forms

This note derives the pullback formulas for differential forms (scalar
0-form, vector potential 1-form, magnetic flux 2-form) under the 3D
sphere-centered Kelvin map, and uses the formulas to check energy
invariance against the material-modulation conventions used in the
A-formulation FEM code.

## 1. Setup

Kelvin inversion centered at origin with radius `R`:

```
phi : r  ->  r' = (R^2 / |r|^2) r
```

Properties:
- Involutive: `phi(phi(r)) = r`.
- Maps the physical exterior `|r| > R` bijectively to the computational
  interior `|r'| < R`.
- The interface `|r| = R` is fixed pointwise.
- Image of physical infinity is `r' = 0` (the Kelvin "GND").

## 2. Jacobians

For the forward map `r' = (R^2/|r|^2) r`:

```
dr'^j / dr^i = R^2 * [delta^j_i / |r|^2 - 2 r^j r^i / |r|^4]
            = (R/|r|)^2 * (delta^j_i - 2 m^j m^i),     m = r/|r|
```

In computational variables (`|r| = R^2/rho'`, `m = n = r'/rho'`):

```
J^j_i := dr'^j / dr^i  =  (rho'/R)^2 * H^j_i
```

where `H = I - 2 n n^T` is a Householder reflection across the plane
perpendicular to the radial direction `n`.  Its determinant is `-1`
and `H^2 = I`.

By involution, the inverse Jacobian has the same form evaluated at the
other point:

```
(J^{-1})^j_i  =  dr^j / dr'^i  =  (R/rho')^2 * H^j_i
```

The volume Jacobian:

```
|det(J^{-1})| = (R/rho')^6     =>    dV_phys = (R/rho')^6 * dV_comp
```

## 3. 0-form (scalar)

A 0-form `f` pulls back by composition:

```
f_comp(r') = f_phys(r = phi(r'))
```

No prefactor; it is the same scalar value at the corresponding point.
For "harmonic continuation" conventions (classical Kelvin theorem),
authors sometimes redefine the scalar by an explicit `(R/rho')` factor;
that is a *choice of representation*, not a pullback.

## 4. 1-form (vector potential A)

A 1-form `omega = omega_i dr^i` pulls back by

```
(phi^* omega)_j (r') = (dr^i / dr'^j) * omega_i (r=phi(r'))
                   = (R/rho')^2 * H^i_j * omega_i(r_phys)
```

In matrix form (`H` symmetric):

```
A_comp(r') = (R / rho')^2 * H * A_phys(r_phys)
          = (R / rho')^2 * [ A_phys - 2 (A_phys . n) n ]
```

The factor decomposes into

- **(R/rho')^2 magnitude factor** (from the Jacobian determinant scaling
  per index), and
- **Householder reflection** that flips the radial component of A while
  preserving tangential components.

For tangential A (e.g. azimuthal A from a uniform `B_0 z_hat`
background), `H A = A` and only the scalar `(R/rho')^2` survives.

## 5. 2-form (magnetic flux B = curl A)

A 2-form `Omega = (1/2) Omega_{ij} dr^i ^ dr^j` pulls back by

```
(phi^* Omega)_{kl} = (dr^i/dr'^k)(dr^j/dr'^l) * Omega_{ij}
                  = (R/rho')^4 * H^i_k * H^j_l * Omega_{ij}
```

In 3D, `Omega` is dual to a vector via the Hodge star:
`Omega^k = (1/2) eps^{kij} Omega_{ij}`. Using the Levi-Civita identity

```
eps^{mkl} H^i_k H^j_l = det(H) * H^m_p * eps^{pij} = -H^m_p eps^{pij}
```

we obtain

```
B_comp = -(R/rho')^4 * H * B_phys
```

The pseudovector pickup of `det(H) = -1` is the geometric origin of the
sign flip familiar from the H-formulation rule
`H'_s = -(rho'/R)^2 H_s`. Note however the EXPONENT differs: A is
`(R/rho')^2` and B is `(R/rho')^4`.

## 6. Numerical verification (uniform B background)

Take `B_phys = B_0 z_hat`, so `A_phys = (B_0/2)(-y, x, 0)` (azimuthal,
tangential).

At `r' = (rho', 0, 0)`, `r_phys = (R^2/rho', 0, 0)` and `n = (1, 0, 0)`.
The Householder leaves `A_phys` invariant (it is tangential), so

```
A_comp(rho', 0, 0) = (R/rho')^2 * (B_0/2)(0, R^2/rho', 0)
                   = (B_0 / 2) * R^4 / rho'^3 * y_hat
```

Generalize to arbitrary r' (still tangential):

```
A_comp(r') = (B_0/2) * R^4 / rho'^4 * (-r'_y, r'_x, 0)
```

Compute curl_z directly:

```
A_comp_y = (B_0/2) R^4 r'_x / rho'^4
A_comp_x = -(B_0/2) R^4 r'_y / rho'^4

partial(A_comp_y)/partial(r'_x)  at  (rho', 0, 0)
   = (B_0/2) R^4 * [1/rho'^4 - 4 r'_x^2/rho'^6]
   = (B_0/2) R^4 * (1 - 4)/rho'^4 = -(3 B_0/2) R^4/rho'^4

-partial(A_comp_x)/partial(r'_y) at (rho', 0, 0)
   = (B_0/2) R^4 * 1/rho'^4

curl_z(A_comp) = -3/2 + 1/2 = -1, times (B_0 R^4 / rho'^4)
             = - (R/rho')^4 * B_0
```

This matches the 2-form pullback `B_comp = -(R/rho')^4 H B_phys`
(with `H = I` for tangential B at this evaluation point, the minus sign
coming from `det H = -1` in the Hodge identity). The 1-form pullback
applied to `A` with `(R/rho')^2` factor reproduces the 2-form pullback
in `B` with `(R/rho')^4` factor — a useful internal consistency check.

## 7. Energy invariance and the implied `nu_kelvin`

The physical magnetic energy

```
W = ∫_phys (1/2) nu_0 |B_phys|^2 dV_phys
```

should equal the computational integral with some `nu_kelvin`:

```
W = ∫_comp (1/2) nu_kelvin |B_comp|^2 dV_comp
```

Substituting `|B_comp|^2 = (R/rho')^8 |B_phys|^2` and
`dV_phys = (R/rho')^6 dV_comp`:

```
nu_0 * (R/rho')^6 * dV_comp = nu_kelvin * (R/rho')^8 * dV_comp
                             ----------------------------------
                                  on each computational cell
```

```
=> nu_kelvin = nu_0 * (rho'/R)^2
```

## 8. Canonical resolution: Nagamine CEFC 2026 (2026-04-16)

The derivation in sections 4-7 gives `nu_kelvin = (rho'/R)^2 * nu_0`.
This is the **canonical Nagamine / Sugahara 2022 convention** and is
the one to use in FEM code.

### 8.1. Authoritative reference

> **H. Nagamine, T. Yamaguchi, K. Sugahara**, "A Pullback-Based
> Formulation of Kelvin Transformation in Electromagnetic Field
> Analysis," CEFC 2026 (Thessaloniki), id 350.
>
> Derives the material transformation law for vector fields via
> pullbacks of differential forms and bilinear energy functionals (not
> a metric-based approach with ad-hoc scale factors). Confirms
> `nu' = (rho'/R)^2 * nu` for 3D spherical (conformal) Kelvin, citing
> Sugahara 2022 IEEE TransMag 58(9) as ref [3].

Key steps (Nagamine eqs. 5-9):

Spherical Kelvin map: `k(r, theta, phi) = (R^2/r, theta, phi)`,
orientation-reversing with `sgn(k) = -1`.

Pullback of orthonormal 1-form basis (conformal, same factor on each
angular component):
```
k*(e^{r'}) = -(R/r)^2 e^r     (same for e^theta, e^phi)
```

Hodge and inner product pullback (eq. 7):
```
star(k*(B')) = -(R/r)^2 * k*(star'(B'))           for 2-form B'
g(k*(w'), k*(w')) = +(R/r)^4 * g'(w', w')         for 1-forms
```

Bilinear integrand (eq. 8):
```
nu <dw, dA> = nu <dw', dA'>' * (R/r)^{2+2+4}
            = nu <dw', dA'>' * (R/r)^8
```

Volume: `dOmega = k*(-R^6/r'^6 dOmega')`. Substituting and using
`sgn(k) = -1`:

```
W_m^E = int_Omega_E nu <dw, dA> dOmega
      = sgn(k) int_Omega' nu <dw', dA'>' * (-r'^2/R^2) dOmega'
      = int_Omega' nu * (r'/R)^2 <dw', dA'>' dOmega'
```

Comparing with `W_m' = int_Omega' nu' <dw', dA'>' dOmega'`:

```
nu' = nu * (r'/R)^2        [CANONICAL: Nagamine eq. 9]
```

### 8.2. Numerical validation (Nagamine CEFC 2026 §III)

Toroidal current loop, major radius `a = 0.1 m`, wire radius
`b = 0.01 m`, magnetic moment `m = I * pi * a^2 = 1 A m^2`. Kelvin
sphere `R = 1 m`.

Exterior magnetic energy:

- Analytical (dipole approximation): `mu_0 m^2 / (12 pi R^3) = 3.333e-8 J`
- FEM on transformed domain Omega' using `nu' = (rho'/R)^2 * nu`: `3.344e-8 J`
- Error: `+0.33%`

This confirms that the correct FEM material coefficient in the
transformed domain is `nu' = (rho'/R)^2 * nu_0`, and that the
magnetic energy computed directly in Omega' coordinates reproduces
the physical exterior energy.

### 8.3. Cylindrical (2D non-conformal) Kelvin

Cylindrical inversion `k(rho, phi, z) = (R^2/rho, phi, z)` is NOT
conformal. The pullback yields an anisotropic reluctivity tensor
(Nagamine eq. 12):

```
nu' = diag(1, 1, (rho'/R)^4) * nu       (only axial component modulated)
```

Use cylindrical Kelvin examples directly with this tensor form, not
the scalar 3D factor.

### 8.4. Solution pullback vs material modulation

Two DIFFERENT factors arise in Kelvin, don't confuse them:

| Concept                | Factor                | Used in                            |
|------------------------|-----------------------|------------------------------------|
| Material `nu`          | `(rho'/R)^2 * nu_0`   | FEM bilinear form (Nagamine eq. 9) |
| Material `mu`          | `(R/rho')^2 * mu_0`   | Omega/H-form (reciprocal)          |
| Solution A (1-form)    | `(R/rho')^2 * H A`    | PEEC source pullback               |
| Solution B (2-form)    | `-(R/rho')^4 * H B`   | B-field pullback                   |

The material factor `(rho'/R)^2` for nu is the **bilinear energy
equivalence** result. The solution factors `(R/rho')^2` (for A) and
`(R/rho')^4` (for B) are **pullbacks of the field itself** and come
directly from the inverse Jacobian (see sections 4, 5 of this note).

### 8.5. Mislabel in `test_nu_convention.py` (historical note)

An earlier version of
`examples/kelvin_transformation/A-formulation/test_nu_convention.py`
defined two convention labels as follows:

| label (in old file)   | factor             |
|-----------------------|--------------------|
| `sugahara`            | `(R/rho')^2 * nu_0`|
| `energy_invariant`    | `(rho'/R)^2 * nu_0`|

The labels are **reversed** from the actual Sugahara 2022 result.
Per Nagamine's derivation (with Sugahara as co-author, citing ref
[3] = Sugahara 2022), the Sugahara convention is `(rho'/R)^2 * nu_0`,
which the old file labeled `energy_invariant`.

The 2026-04-15 session treated the mislabeled empirical result
(`(R/rho')^2` winning with +0.90% vs +4.80% on a gapped-torus
benchmark) as authoritative and "resolved" this note in favor of
`(R/rho')^2`. That resolution was incorrect. The present text (§8)
supersedes the earlier one.

### 8.6. Open investigation: empirical discrepancy in old test

The old test reported `(R/rho')^2 -> +0.90%` vs `(rho'/R)^2 -> +4.80%`
on a gapped-torus benchmark. Per Nagamine, `(rho'/R)^2` is correct,
so `(rho'/R)^2` should give the better result. The discrepancy is
likely due to FEM setup (GND / gauge regularization) and is logged
as a follow-up investigation (see
`memory/project_kelvin_e2c_deferred.md`). Re-running with proper GND
enforcement should restore agreement.

### 8.7. Centralized API

New code should import Nagamine-canonical factors from
`radia.kelvin_source` rather than re-deriving inline. See
`examples/kelvin_transformation/CONVENTION.md` for the factory
functions.

## 9. Validation script

See also `validate_radia_HB_kelvin.py` in the same directory for a
direct Radia-vs-FEM cross-check of the material factor + A 1-form
pullback.

## 9. Validation script

`examples/kelvin_transformation/A-formulation/validate_radia_HB_kelvin.py`
performs an A/B comparison on the same Kelvin mesh:

- Case A: full FEM (volume J, two-sphere Kelvin, canonical
  `nu_kelvin = (rho'/R)^2 nu_0` per Nagamine CEFC 2026).
- Case B: inject Radia ObjArcCur analytical A onto the same mesh, with
  the 1-form pullback in the Kelvin exterior, integrate
  `H . B = nu |curl A|^2`.

Empirical results (square 2a x 2a torus, R_coil=30mm, a_coil=3mm,
gap=5deg, R_K=60mm, offset=150mm, coarse mesh maxh=12mm), recorded
before the Nagamine resolution — kept as historical log; the
convention labels in the "inject" column refer to the factor applied
to A in the kelvin exterior, NOT the material factor:

| | L_inner [nH] | L_kext [nH] | Total |
|---|---|---|---|
| Case A FEM (volume J) | 78.8 | 2.6 | 81.4 |
| Case B (R/rho')^2 inject (Phase 1) | 78.8 | 73 | 152 |
| Case B (rho'/R)^2 inject (Phase 1) | 98.7 | 0.15 | 98.9 |
| Case B (rho'/R)^2 + Householder (Phase 2) | 98.7 | 0.62 | 99.4 |

Per Nagamine: the correct **A solution pullback** is `(R/rho')^2 H A`
(not `(rho'/R)^2`), and the correct **material nu** is
`(rho'/R)^2 nu_0`. The script in §9 was run in a mixed state where the
two factors were confused; results should be re-run with the corrected
`validate_radia_HB_kelvin.py` (factor fix committed 2026-04-16).

## 10. References

- **H. Nagamine, T. Yamaguchi, K. Sugahara**, "A Pullback-Based
  Formulation of Kelvin Transformation in Electromagnetic Field
  Analysis," CEFC 2026 (Thessaloniki), id 350. **Authoritative
  derivation**: pullback of orthonormal 1-form basis + bilinear energy
  functional gives `nu' = (rho'/R)^2 * nu` for 3D spherical (conformal)
  Kelvin; `nu' = diag(1,1,(rho'/R)^4) * nu` for 2D cylindrical
  (non-conformal). Validated numerically: toroidal loop analytical
  dipole energy 3.333e-8 J vs FEM 3.344e-8 J (+0.33%).
- **K. Sugahara**, "Electromagnetic analysis of eddy current testing
  with Kelvin transformation," IEEE Trans. Magn. 58(9), 1-6, Sept.
  2022 (cited as ref [3] in the Nagamine paper). Original A-formulation
  derivation of `nu' = (rho'/R)^2 * nu`.
- Sugahara, Nagamine, Kameari, "Kelvin Transformation for Open
  Boundary Problems in Reduced Potential Formulation"
  (digest, examples/kelvin_transformation/digest.pdf, 2026).
  Companion tutorial covering H-formulation rule
  `H'_s = -(rho'/R)^2 H_s`.
- A. Bossavit, *Computational Electromagnetism*, Academic Press, 1998
  (ref [4] in Nagamine) -- differential geometric framework for EM.
- Wong and Ciric, COMPEL 4(3), 1985; Freeman and Lowther, IEEE Trans.
  Magn. 24(6), 1988 (refs [1], [2] in Nagamine).
- Nabizadeh, Ramamoorthi, Chern, "Kelvin transformations for
  simulations on infinite domains", ACM Trans. Graphics 40(4), 2021
  (general k-form pullback framework).
