# Analytical Formulas — PDF cross-reference

This document maps each formula in [`src/radia/analytical_formulas/`](../src/radia/analytical_formulas/)
back to its origin in the IEE Japan review series

> 若尾真治, 五十嵐一, 藤原耕二, 野口聡, 松尾哲司, 亀有昭久,
> "Useful Formulas of Analytical Integration in Electromagnetic Field
> Computations (Part 1..9)", IEE Japan Joint Technical Meeting on
> Static Apparatus and Rotating Machinery (SA / RM), 2002–2007.

Numbering below is `Part N, eq M` exactly as printed in the original
research-meeting reports. The PDFs themselves are not redistributed with
the Radia repository (lab-internal copies).

## Scope

The package collects **closed-form** expressions only. It is the
reference layer that the rest of Radia (numerical solvers, BEM, FEM)
can be checked against. Each formula is one short Python module:
priority is correctness and traceability, not raw speed.

## Module index

### Group B — fields and global quantities

| Module | Symbols | PDF reference |
|--------|---------|---------------|
| [`ellipsoid`](../src/radia/analytical_formulas/ellipsoid.py)   | `demag_factor_prolate`, `demag_factor_oblate`, `demag_factor_rotational`, `ellipsoid_internal_field`, `ellipsoid_torque` | Part 5, eq 38–44 |
| [`ac_locus`](../src/radia/analytical_formulas/ac_locus.py)     | `ac_locus_axes`, `ac_locus_axes_batch` | Part 5, eq 29–37 |
| [`shielding`](../src/radia/analytical_formulas/shielding.py)   | `shielding_factor_cylinder`, `shielding_factor_sphere` | Part 1, eq 23–24 |
| [`rect_magnet_2d`](../src/radia/analytical_formulas/rect_magnet_2d.py) | `rect_magnet_2d_A`, `rect_magnet_2d_B` | Part 2, eq 2–3 |
| [`plate_eddy`](../src/radia/analytical_formulas/plate_eddy.py) | `plate_eddy_T`, `plate_eddy_J` | Part 1, eq 26–27 |

### Group C — coil geometries and numerical utilities

| Module | Symbols | PDF reference |
|--------|---------|---------------|
| [`solenoid_central`](../src/radia/analytical_formulas/solenoid_central.py)   | `fabri_F`, `solenoid_central_field`, `solenoid_axial_field` | Part 4 §4, eq 26–27 |
| [`three_phase_line`](../src/radia/analytical_formulas/three_phase_line.py)   | `vector_potential_z`, `field_xy`, `triangle_far_field_amplitude`, `planar_far_field_amplitude`, `helical_near_field_amplitude`, `helical_far_field_amplitude`, plus arrangement and balanced-current builders | Part 4 §5, Part 5 §3 |
| [`elliptic_integrals`](../src/radia/analytical_formulas/elliptic_integrals.py) | `K_hastings_2`, `E_hastings_2`, `K_hastings_4`, `E_hastings_4` | Part 3 §3, Tables 1–2 |
| [`gauss_legendre`](../src/radia/analytical_formulas/gauss_legendre.py)       | `gauss_legendre_nodes_weights`, `gauss_legendre_integrate`, `gauss_legendre_integrate_2d` | Part 3 §4, Table 3 |

### Group D — Part 6 / 8 / 9 extensions (eddy currents, AC, numerics, average B)

| Module | Symbols | PDF reference |
|--------|---------|---------------|
| [`plate_eddy`](../src/radia/analytical_formulas/plate_eddy.py) (extension)   | `plate_eddy_dissipation` (total Joule loss in the plate via numerical integration of the analytic J series) | Part 6 §3 |
| [`shielding`](../src/radia/analytical_formulas/shielding.py) (extensions)    | `shielding_factor_sphere_thin_ac`, `shielding_factor_cylinder_thin_ac` (AC, complex S), `spherical_shell_internal_field` (interior H, M, image dipole) | Part 6 §2, Part 8 §2 |
| [`conductor_impedance`](../src/radia/analytical_formulas/conductor_impedance.py) | `skin_depth`, `planar_surface_impedance` (Z_s = (1+j)/(σδ)), `cylinder_ac_impedance` (full Bessel solution), `cylinder_dc_resistance`, `cylinder_internal_inductance` | Part 6 §4–§5 |
| [`adaptive_quadrature`](../src/radia/analytical_formulas/adaptive_quadrature.py) | `patterson_nodes_weights` (n=0..3, 1/3/7/15 points), `adaptive_integrate` (node-reusing refinement) | Part 9 §2, Table 1 |
| [`cuboid_average_field`](../src/radia/analytical_formulas/cuboid_average_field.py) | `average_B_in_box` (closed-form C++ 64-corner sum; `method="numerical"` Gauss-Legendre kept for cross-checks), `average_demag_tensor` | Part 6 §7, eq 53–56 |
| [`induction_heating`](../src/radia/analytical_formulas/induction_heating.py) | `cylinder_axial_eddy_loss` (full Bessel form), `cylinder_axial_eddy_loss_small_ka` (Faraday eddy limit), `cylinder_axial_eddy_loss_thin_skin` (planar surface-impedance limit) | Smythe §11.07, Landau-Lifshitz §59, Jackson §5.18 |

Tests live in [`tests/analytical_formulas/`](../tests/analytical_formulas/);
runnable demonstrations in the notebook [`docs/analytical_formulas/analytical_formulas.ipynb`](analytical_formulas/analytical_formulas.ipynb).

## ellipsoid — rotational ellipsoid demag and torque

Geometry: semi-axes ``(a, a, c)`` with ``c`` the polar axis.
``c > a`` ⇒ prolate (cigar); ``c < a`` ⇒ oblate (disk); ``c = a`` ⇒
sphere (``N_x = N_y = N_z = 1/3``).

Polar demagnetization factor:

```
prolate  (eq 40)  N_z = (1 - e²)/e³ · (atanh e − e),   e = √(1 − (a/c)²)
oblate   (eq 39)  N_z = (1 + e²)/e³ · (e − atan e),    e = √((a/c)² − 1)
```

Equatorial factors are ``N_x = N_y = (1 − N_z) / 2`` and the trace is
1 (verified in
[`test_demag_factors_sum_to_unity`](../tests/analytical_formulas/test_ellipsoid.py)).

For aspect ratios within ``|c/a − 1| < 1e−4`` the implementation
switches to the Taylor series ``N_z = 1/3 − 2 ε /15 + 8 ε² /105 − ⋯``
with ``ε = (c/a)² − 1``, valid for both branches; this avoids the
``0/0`` form of the closed expressions at ``e → 0``.

Internal field (linear material, applied along principal axis ``i``,
eq 38):

```
H_i = H_0 / (1 + χ_r · N_i)
```

Torque (uniform field, axis tilted by ``α``, eq 44):

```
T_z = − μ_0 χ_r² (1 − 3 N_z)
        / [2 (1 + χ_r N_z) (2 + χ_r (1 − N_z))]
        · sin(2 α) H_0² V,    V = (4/3) π a² c.
```

Sign convention: ``T_z < 0`` is the restoring torque pulling ``α → 0``;
this is the prolate / soft-iron case. For oblate (``N_z > 1/3``) the
torque is repulsive (axis aligns perpendicular to the field).

## ac_locus — major / minor axis of an AC vector locus

Closed form (eq 37): for a complex 3-vector phasor ``B`` whose
time-domain version is ``b(t) = Re(B e^{j ω t})``,

```
s² = |B_x|² + |B_y|² + |B_z|²
z² = B_x² + B_y² + B_z²        (complex; the complex dot of B with itself)

|B|_max = √((s² + |z²|) / 2)
|B|_min = √((s² − |z²|) / 2)
```

The locus is a planar ellipse (3D mode) parameterised by time. Useful
to extract the peak field over a period for hysteresis-loss estimates
and the peak current for inductive heating uniformity constraints.

Tests cover (a) linear, circular and elliptical limit cases,
(b) invariance under physical-space rotations, (c) match against a
4 000-sample brute-force time sweep
([`test_axes_match_brute_force_time_sweep`](../tests/analytical_formulas/test_ac_locus.py)),
(d) batch / higher-rank shape handling.

## shielding — cylindrical / spherical magnetic shell

```
S_cyl = 4 μ_r / [(μ_r + 1)² − (a/b)² (μ_r − 1)²]                     (eq 23)
S_sph = 9 μ_r / [(μ_r + 2)(2 μ_r + 1) − 2 (a/b)³ (μ_r − 1)²]         (eq 24)
```

``S = |H_inside| / |H_applied|``; ``S → 1`` at ``μ_r → 1`` (no
shielding) and ``S → 0`` at ``μ_r → ∞`` for any finite wall thickness.
The cylinder formula assumes infinite length and a perpendicular
applied field; finite-cylinder caps are discussed in Rikitake [13] of
the PDF (not implemented).

## rect_magnet_2d — 2D uniformly-magnetised rectangular bar

Cross-section ``2 a × 2 b``, infinitely long along ``z``, magnetisation
``M = (M_x, M_y)`` in A/m (Radia convention; ``J = μ_0 M`` in Tesla).
With

```
u_i = x − x_i,  x_i ∈ {−a, +a},  i = 1, 2
v_j = y − y_j,  y_j ∈ {−b, +b},  j = 1, 2

G(u, v)      = u ln(u² + v²) + 2 v atan(u / v)
G_swap(u, v) = v ln(u² + v²) + 2 u atan(v / u)
```

```
A_z = (μ_0 / 4π) Σ_{i,j} (−1)^{i+j} [M_x G(u_i, v_j) − M_y G_swap(u_i, v_j)]
```

The PDF prints only the ``M_x`` term (eq 2); the ``M_y`` term is the
result of applying the same surface-magnetisation-current derivation to
the left and right faces.

```
B_x =  ∂A_z/∂y = (μ_0 / 2π) Σ (−1)^{i+j} [M_x atan(u_i/v_j) − ½ M_y ln(u_i² + v_j²)]
B_y = −∂A_z/∂x = −(μ_0 / 4π) Σ (−1)^{i+j} [M_x ln(u_i² + v_j²) − 2 M_y atan(v_j/u_i)]
```

Constants in ``∂(ln)/∂u`` and ``∂(atan)/∂v`` drop under the alternating
``(−1)^{i+j}`` sum.

Branch choice: ``np.arctan`` (principal branch ``[−π/2, π/2]``), not
``np.arctan2`` — the latter would introduce a ``2π·v`` jump at
``v = 0``.

## plate_eddy — eddy current in a thin rectangular plate

Slow-field limit of a plate of in-plane half-extents ``a × b`` and
thickness ``d`` in a uniform perpendicular ``B_z(t)``; valid when
``d / (μ σ Ḃ · min(a, b) / B) << 1`` (eq 25). Eddy current is purely
in-plane and is the curl of a single z-component vector potential
``T_z``:

```
T_z(x, y) = (σ Ḃ / 8) · {
    x² − a²
    + (32 a² / π³) Σ_{n=0}^∞ ((−1)^n / (2n+1)³)
          · cos(k_n x) cosh(k_n y) / cosh(k_n b)
}                                                   (eq 26)
```

with ``k_n = (2 n + 1) π / (2 a)``. Direct differentiation gives
(eq 27)

```
J_x =  (2 σ Ḃ a / π²) Σ ((−1)^n / (2n+1)²) cos(k_n x) sinh(k_n y) / cosh(k_n b)
J_y = −(σ Ḃ x / 4)
      + (2 σ Ḃ a / π²) Σ ((−1)^n / (2n+1)²) sin(k_n x) cosh(k_n y) / cosh(k_n b)
```

The closed-form ``− σ Ḃ x / 4`` term in ``J_y`` cancels exactly with
the cosh-cosh series at ``y = ± b`` so that ``J · n = 0`` on the y
edges (the residual is the Fourier triangle-wave reconstruction error
of ``x`` on ``[−a, a]``); on the x edges ``J_x`` vanishes term-by-term
because ``cos(k_n a) = 0``.

Truncation: terms decay as ``1/(2n+1)³`` (T) and ``1/(2n+1)²`` (J);
the ``cosh(k_n y) / cosh(k_n b)`` ratio decays exponentially in ``n``
for ``|y| < b``. The default ``n_terms = 200`` gives ``< 1e−10``
relative error away from the corner singularity. Use ``n_terms = 2000``
for sub-``1e−5`` accuracy at the boundary.

## solenoid_central — rectangular cross-section solenoid central / axial field

Geometry: rectangular cross-section solenoid with inner radius
``a_1``, outer radius ``a_2``, half-length ``b``, azimuthal current
density ``J``. With ``alpha = a_2 / a_1`` and ``beta = b / a_1``,

```
B_0 = mu_0 a_1 J F(alpha, beta)                                  (eq 26)
F(alpha, beta) = beta ln[(alpha + sqrt(alpha² + beta²))
                         / (1 + sqrt(1 + beta²))]                (eq 27)
```

The on-axis field at general ``z`` follows in closed form:

```
B_z(0, z) = (mu_0 J / 2) * [
    (b - z) ln((a_2 + sqrt(a_2² + (b - z)²))
              /(a_1 + sqrt(a_1² + (b - z)²)))
  + (b + z) ln((a_2 + sqrt(a_2² + (b + z)²))
              /(a_1 + sqrt(a_1² + (b + z)²)))
]
```

Limits: ``B_0 -> mu_0 J (a_2 - a_1)`` for an infinitely long
solenoid; ``B_0 ~ mu_0 a_1 J beta ln(alpha)`` for a thin pancake.

## three_phase_line — three-phase straight and helical lines

The phasor convention is the **peak** amplitude (eq 33 of Part 4):
``I_k(t) = Re(I_peak,k exp(j omega t))`` with
``I_peak,k = I_peak exp(- j 2 k pi / 3)``. For an RMS input
``I_rms`` pass ``I_peak = sqrt(2) I_rms``.

```
triangle (3 lines on circle a):       |B|(t) ~ 3 mu_0 a I_peak / (4 pi r²)   (eq 38-39)
planar   (3 lines spaced d):           |B|_max ~ sqrt(3) mu_0 d I_peak / (2 pi r²)   (eq 41-43)
helical  (a << p << r):               |B| ~ (3 pi mu_0 a I_peak / (4 sqrt(p r))) * exp(-2 pi r / p)   (eq 28)
```

Triangle far-field is circularly polarised
(``|B|_max = |B|_min``); planar is linearly polarised
(``|B|_min = 0``). The square-root-3 in the planar coefficient comes
from ``q* − q = j sqrt(3)`` of the three-phase set; this differs from
the OCR-extracted PDF text, but matches the rigorous Biot-Savart
calculation of [`field_xy`](../src/radia/analytical_formulas/three_phase_line.py)
to converging precision (verified in
[`test_three_phase_line.py`](../tests/analytical_formulas/test_three_phase_line.py)).

The hexagon arrangement that gives a 1/r**3 quadrupole asymptote
(Part 4, eq 44) requires a **specific** phase-on-vertex assignment
that is *not* satisfied by the simple
``hexagon_positions + balanced_six_phase_currents`` builder pair --
that combination gives a 1/r**2 dipole far-field (regression tested).
The 1/r**3 closed-form is therefore not exposed; users who need it
can construct the proper arrangement and integrate via
[`field_xy`](../src/radia/analytical_formulas/three_phase_line.py).

## elliptic_integrals — Hastings polynomial K(k), E(k)

```
K(k) ~ sum_{i=0}^{n} (a_i + b_i ln(1/x)) x**i,   x = 1 - k²
E(k) ~ sum_{i=0}^{n} (c_i + d_i ln(1/x)) x**i
```

Coefficients are listed in PDF Tables 1 (n = 2) and 2 (n = 4); they
are reproduced verbatim in
[`elliptic_integrals.py`](../src/radia/analytical_formulas/elliptic_integrals.py).
Accuracy:

|                   | claimed (PDF) | observed vs scipy |
|-------------------|---------------|-------------------|
| degree 2          | < 4e-5        | ~ 4e-5            |
| degree 4          | < 2e-8        | ~ 2e-8            |

`scipy.special.ellipk` / `ellipe` provide the same quantities to
machine precision via Carlson's reductions; the polynomial form is
preserved here as a documented reference and as a no-scipy fallback.

## gauss_legendre — quadrature nodes / weights up to n = 24

Standard ``n``-point rule on ``[-1, +1]``: roots of ``P_n`` as nodes,
``w_i = 2 (1 - x_i²) / (n P_{n-1}(x_i))²`` as weights, exact for
polynomials of degree ``2 n - 1``. Affine-mapped helpers
[`gauss_legendre_integrate`](../src/radia/analytical_formulas/gauss_legendre.py)
and `_integrate_2d` apply the rule on a general interval / rectangle.
The implementation delegates to ``numpy.polynomial.legendre.leggauss``
for the underlying computation, so all 24 entries of PDF Table 3 are
machine-precision rather than the printed-table truncation.

## plate_eddy (extended) — total Joule dissipation

```
P = (d / sigma) integral_{-a}^{+a} integral_{-b}^{+b} (J_x**2 + J_y**2) dx dy
```

with the analytic ``J_x``, ``J_y`` of Part 1 §6.1 evaluated by the
existing :func:`plate_eddy_J`. The PDF Part 6 §3 gives a compact
closed form

```
P = (4 sigma a**3 b d / 3) Bdot**2
    - (256 a**4 d / pi**5 sigma) sum_n tanh(lambda_n b) / lambda_n**5
```

but the OCR-extracted formula leaves the second-term sign and the
position of ``sigma`` ambiguous. We integrate the analytic ``J``
numerically by tensor-product Gauss-Legendre — exact up to
quadrature truncation, no ambiguity.

## shielding (extensions) — AC and shell-interior fields

Part 8 §2 gives the AC thin-shell shielding factors when the wall
thickness is much smaller than both the radius and the skin depth:

```
S_sph = 1 / (1 + j omega mu_0 sigma a Delta / 3)
S_cyl = 1 / (1 + j omega mu_0 sigma a Delta / 2)
```

The phase of ``S`` is the field's lag relative to the applied field.

Part 6 §2 gives the four uniform fields (interior, shell wall, image
moment) for a magnetic spherical shell:

```
H_hollow = 9 H0 / [9 + 2 (mu_r - 1)**2 (1 - (a/b)**3) / mu_r],
H_shell  = (3 (2 + mu_r) / (3 mu_r)) * H_hollow,
M_shell  = (mu_r - 1) * H_shell,
M_image  = -((mu_r - 1)(mu_r + 2) / (3 mu_r)) * H_hollow * V_outer.
```

## conductor_impedance — Part 6 §4 + §5

Skin depth, planar surface impedance, and the full Bessel solution
for the AC impedance of a solid cylindrical conductor. Limits:

* low frequency: ``Z = R_dc + j omega L_int``,
  ``L_int = mu_0 / (8 pi)`` per unit length.
* high frequency: ``Z = Z_s / (2 pi a)`` with the planar
  ``Z_s = (1 + j) / (sigma d_skin)``.

## adaptive_quadrature — Part 9 §2

Gauss-Patterson family of nested rules (n=0,1,2,3 give 1, 3, 7, 15
points) where the nodes of the lower-order rule are exactly contained
in the higher-order rules. Together with a dict-cache on
:func:`adaptive_integrate` this halves the function-evaluation cost
of refining quadrature accuracy. Polynomial exactness verified
against ``x^p`` for p up to 22 (n=3 rule). Higher-order tables
(n=4, 5 from PDF Table 1) were OCR-extracted but rejected pending
high-fidelity verification.

## cuboid_average_field — Part 6 §7 (closed-form C++ kernel, v4.22.0)

Spatial average of ``B`` over a target rectangular box from a uniform
``M`` source box. Useful for FEM-MMM coupling (cell-averaged ``B`` in
a rectangular FEM element), mesh-to-mesh transfer, and micromagnetics
on a regular cubic lattice (the original PDF use case).

The default path is a **closed-form 64-corner inclusion-exclusion
sum** of two 3-fold antiderivative primitives ``G1`` (diagonal demag)
and ``G2`` (off-diagonal demag), implemented in C++
(``src/core/rad_average_field.cpp``) with permutation rules that cover
all 9 entries of the average demag tensor. Performance is ~40 µs per
call (817× faster than the Gauss-Legendre numerical baseline).

Both ``G1`` and ``G2`` were derived in-house by direct sympy
integration of ``-∂²(1/R)/∂X²`` and ``-∂²(1/R)/∂X∂Y`` (the demag-tensor
kernels) over the displacement-box, after observing that the originally
cited parent reference Stafl 1967 §3.4 is a 2D rectangular-conductor
problem and not the 3D cuboid-magnetisation problem of Wakao Part 6 §7.
Verified at ``rel_err ≈ 1e-13`` against the Gauss-Legendre baseline
across diagonal-offset, side-by-side, distant, and self-source
geometries; the canonical self-cube demag tensor evaluates to ``-(1/3)·I``
exactly (trace = -1).

The closed-form path handles overlapping source-target boxes (including
the self-source case) via a ``mu_0 * M * V_overlap / V_T`` correction
term added to the H-field part returned by the corner sum.

**Numerical caveat.** The 64-corner sum suffers ULP cancellation when
``V_T ≪ V_S`` (the corner-primitive values are O(1) and the
alternating sum needs to evaluate to a much smaller number). For
true point-field evaluation, call ``CuboidMagnet.get_B(point)``
directly. For mismatched-scale box averages, ``average_B_in_box(...,
method="numerical", n_quad=8)`` falls back to tensor-product
Gauss-Legendre quadrature kept for diagnostic / cross-validation use.

References for the derivation: Wakao-Igarashi-Fujiwara-Kameari Part 6
§7 eq 53–56; Newell, Williams & Dunlop, J. Geophys. Res. 98 (1993)
9551–9555.

## induction_heating — canonical AC cylinder Joule loss (v4.23.0)

Closed-form time-averaged Joule loss per unit length of an infinite
conductive cylinder placed inside a long axial-current solenoid, the
canonical induction-furnace work-piece geometry. Closed-form via the
Kelvin functions (`scipy.special.{ber, bei, berp, beip}`):

```
P = (π / σ) H_0² (ka)
    · (ber(ka) ber'(ka) + bei(ka) bei'(ka))
    / (ber(ka)² + bei(ka)²),    k = sqrt(ω σ μ).
```

Two asymptotic forms validated against the full Bessel formula:

* `cylinder_axial_eddy_loss_small_ka` — thick-skin / Faraday-eddy
  limit: `P → π σ ω² μ² H_0² a^4 / 16`.  Matches the elementary
  derivation by Faraday's law plus Joule integral over a cylinder
  with uniform internal `B`.

* `cylinder_axial_eddy_loss_thin_skin` — high-frequency planar
  surface-impedance limit: `P → π a H_0² Re(Z_s)` with
  `Re(Z_s) = sqrt(ω μ / (2 σ))`.

Both limits agree with the full formula across the parameter sweep
in `tests/analytical_formulas/test_induction_heating.py` (tested
ka = 0.001 ... 100; small-ka rtol 1e-4, thin-skin rtol 0.02-0.10
depending on ka).

**Phase γ history.** v4.23.0's `induction_heating` module was the
*pivot* result of an aborted attempt to absorb Stafl 1967
"Electrodynamics of Electrical Machines" Chapter 4 §4.1-§4.4 closed-
forms. The PDF-OCR'd Stafl formulas systematically failed validation
against the canonical thick-skin and thin-skin limits (sphere ~1e10×
overestimate, transverse-cylinder *negative* losses, longitudinal-
cylinder ~800× off at low frequency).  The cylinder formula above
was therefore re-derived from first principles in-house and matches
the textbook references directly. Stafl §4.1 (transformer winding)
and §4.5 (filament near plate) remain unimplemented.

References: Smythe W. R., *Static and Dynamic Electricity*
(McGraw-Hill, 3rd ed. 1968), §11.07. Landau L. D., Lifshitz E. M.,
*Electrodynamics of Continuous Media* (Pergamon, 2nd ed. 1984), §59.
Jackson J. D., *Classical Electrodynamics* (Wiley, 3rd ed. 1999),
§5.18.

## Updating this document

When a new analytical formula is added to
[`src/radia/analytical_formulas/`](../src/radia/analytical_formulas/),
please:

1. add a row to the "Module index" table with the symbol names
   exported in `__init__.py`,
2. add a section that pins the PDF reference (Part N, eq M) and
   restates the formula in unambiguous notation, and
3. cross-link the corresponding test file under
   [`tests/analytical_formulas/`](../tests/analytical_formulas/) and the
   demonstration notebook
   [`docs/analytical_formulas/analytical_formulas.ipynb`](analytical_formulas/analytical_formulas.ipynb).
