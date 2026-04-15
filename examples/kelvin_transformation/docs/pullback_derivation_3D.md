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

## 8. Open question: the code uses (R/rho')^2

The 3D HCurl A-formulation code in this repository
(`Coil_3D_A_HCurl_with_Kelvin.py`) uses

```python
kelvin_fac = R_K**2 / r_prime_sq          # = (R / rho')^2
nu_dict[m] = NU_0 * kelvin_fac            # nu_kelvin = (R/rho')^2 * nu_0
```

which is the **inverse** of the formula derived above. The empirical
result is nevertheless correct (gapped torus L = 89.44 nH vs analytical
torus 88.55 nH, +1.0% on a coarse mesh).

Possible reconciliations (any of which can be the right one):

1. The code's `nu_kelvin = (R/rho')^2 nu_0` is paired with a different A
   pullback than the pure 1-form formula above. If the FEM-solved
   `A_FEM_comp(r')` does NOT equal the 1-form pullback `(R/rho')^2 H
   A_phys`, but is instead something else (e.g. the literal
   `A_phys(r_phys)` value at the point — with ν compensating for the
   missing Jacobians), the bilinear form can still produce the right
   physical inductance.
2. The 2-form pullback derivation uses `det H = -1` for the Hodge
   identity. If the FEM convention treats B as a pure 2-form (no
   pseudovector dualisation), the sign and possibly an exponent change.
3. There is an additional source-side factor in the FEM that absorbs
   the `(rho'/R)^4` excess.

## 9. Validation script

`examples/kelvin_transformation/A-formulation/validate_radia_HB_kelvin.py`
performs an A/B comparison on the same Kelvin mesh:

- Case A: full FEM (volume J, Sugahara two-sphere convention,
  `nu_kelvin = (R/rho')^2 nu_0`).
- Case B: inject Radia ObjArcCur analytical A onto the same mesh, with
  the 1-form pullback in the Kelvin exterior, integrate
  `H . B = nu |curl A|^2`.

Empirical results (square 2a x 2a torus, R_coil=30mm, a_coil=3mm,
gap=5deg, R_K=60mm, offset=150mm, coarse mesh maxh=12mm):

| | L_inner [nH] | L_kext [nH] | Total |
|---|---|---|---|
| Case A FEM (volume J) | 78.8 | 2.6 | 81.4 |
| Case B (R/rho')^2 inject (Phase 1) | 78.8* | 73 | 152* (was 218 before fix) |
| Case B (rho'/R)^2 inject (Phase 1) | 98.7 | 0.15 | 98.9 |
| Case B (rho'/R)^2 + Householder (Phase 2) | 98.7 | 0.62 | 99.4 |

*The (R/rho')^2 row is from an earlier buggy implementation; the
Phase 2 pullback derived in section 4 of this note is the
mathematically correct 1-form formula.

The 25 % gap between Case A inner and Case B inner is unexplained by
mesh refinement alone; investigation continues. The (rho'/R)^2 form
is incorrect per the derivation above (it is the SCALAR-derived
gradient transformation factor, not the 1-form A pullback factor),
so the agreement Case B (rho'/R)^2 ~ Radia open-domain is partially
fortuitous and needs further analysis.

## 10. References

- Sugahara, Nagamine, Kameari, "Kelvin Transformation for Open
  Boundary Problems in Reduced Potential Formulation"
  (digest, examples/kelvin_transformation/digest.pdf, 2026).
  Provides the H-formulation rule
  `H'_s = -(rho'/R)^2 H_s` and the bilinear factor `(R/rho')^2`.
- Sugahara 2022, IEEE Trans. Magn. 58(9), "Electromagnetic analysis
  of eddy current testing with Kelvin transformation" -- A-formulation
  conformal-symmetry derivation (cited as ref [6] in the digest;
  contains the `(R/rho')^2` ν convention used in the code).
- Nabizadeh, Ramamoorthi, Chern, "Kelvin transformations for
  simulations on infinite domains", ACM Trans. Graphics 40(4), 2021.
