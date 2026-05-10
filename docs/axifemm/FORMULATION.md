# axifemm — Formulation

A self-contained derivation of the Henrotte / Meeker axisymmetric finite-
element formulation that `radia.radia_axifemm` implements.  Read this
before reading the C++ source under `src/ext/axifemm/`.  Companion doc
[`AXIFEMM.md`](AXIFEMM.md) covers usage, API, and validation results.

## 1. Problem class

Axisymmetric magnetostatic + linear quasi-static eddy current.  Cylindrical
coordinates `(r, φ, z)`, all fields independent of `φ`.  Source current and
magnetisation are circumferential only:

$$
\mathbf{J} = J_\varphi(r, z)\,\hat\varphi,
\qquad
\mathbf{M} = M_\varphi(r, z)\,\hat\varphi.
$$

The induced eddy current is also circumferential under axial-symmetric
excitation, so the magnetic vector potential reduces to a single scalar
component:

$$
\mathbf{A}(r, z) = A_\varphi(r, z)\,\hat\varphi.
$$

## 2. Magnetic field from the vector potential

In cylindrical coordinates the curl of a φ-only vector field is

$$
\mathbf{B} = \nabla \times \mathbf{A}
= -\frac{\partial A_\varphi}{\partial z}\,\hat r
+ \frac{1}{r}\frac{\partial (r A_\varphi)}{\partial r}\,\hat z.
$$

Thus the radial component carries an explicit `∂/∂z` and the axial
component carries an `(1/r) ∂/∂r` of `r A_φ`.  The latter is the source
of the axis (`r = 0`) singularity that motivates the Henrotte basis
(see § 6).

It is convenient to introduce the **flux function**

$$
\psi(r, z) \equiv 2\pi r A_\varphi(r, z).
$$

Then `2π·∫B_z dr = ∂ψ/∂r` (a useful relationship for visualising flux
tubes), and the field becomes

$$
B_z = \frac{1}{2\pi r}\frac{\partial \psi}{\partial r},
\qquad
B_r = -\frac{1}{2\pi r}\frac{\partial \psi}{\partial z}.
$$

Both components carry a `1/r` factor when expressed via `ψ`.  The
Henrotte trick is to choose a basis so that the polynomial part of `ψ`
contains an explicit `r²` factor near the axis, exactly cancelling the
`1/r`.

## 3. Maxwell equations and the energy functional

Magnetostatic constitutive law inside a region of permeability `μ` and
magnetisation `M`:

$$
\mathbf{H} = \frac{1}{\mu}(\mathbf{B} - \mu_0 \mathbf{M}),
\qquad
\nabla \times \mathbf{H} = \mathbf{J}.
$$

Combining with `B = ∇ × A` gives the vector elliptic equation

$$
\nabla \times \left(\frac{1}{\mu}\nabla \times \mathbf{A}\right)
= \mathbf{J} + \nabla \times \left(\frac{\mu_0 \mathbf{M}}{\mu}\right).
$$

For axisymmetric `J = J_φ ĵ_φ` and `A = A_φ ĵ_φ` this reduces to a
scalar PDE on the meridian plane `Ω ⊂ {(r, z) : r ≥ 0}`:

$$
-\frac{\partial}{\partial r}\!\left[\frac{1}{\mu r}\frac{\partial(r A_\varphi)}{\partial r}\right]
-\frac{\partial}{\partial z}\!\left[\frac{1}{\mu}\frac{\partial A_\varphi}{\partial z}\right]
= J_\varphi - \text{div}\,(\mu_0 \mathbf{M}/\mu)\big|_\varphi.
$$

The corresponding energy functional (per radian of φ) on the meridian
plane is

$$
W[A_\varphi]
= \int_\Omega \left[\frac{1}{2\mu}|\mathbf{B}|^2 - J_\varphi A_\varphi\right]
  \cdot 2\pi r \, dr\, dz,
$$

with the Jacobian `2π r` from the φ-integration of the cylindrical
volume element.  Substituting `B_r, B_z` from § 2 and choosing the flux
function `ψ` as the unknown:

$$
W[\psi]
= \int_\Omega \frac{1}{8\pi^2 \mu r}
  \!\left[\!\left(\frac{\partial \psi}{\partial r}\right)^{\!2}
        + \!\left(\frac{\partial \psi}{\partial z}\right)^{\!2}\right]
  \cdot 2\pi r \, dr\, dz
- \int_\Omega \frac{J_\varphi \psi}{2\pi r} \cdot 2\pi r\, dr\, dz.
$$

The `r` in the Jacobian cancels one of the `1/r` factors, leaving a
weighted Dirichlet energy in `ψ`:

$$
W[\psi]
= \int_\Omega \frac{1}{4\pi \mu r}
  \!\left[\!\left(\frac{\partial \psi}{\partial r}\right)^{\!2}
        + \!\left(\frac{\partial \psi}{\partial z}\right)^{\!2}\right] dr\, dz
- \int_\Omega J_\varphi \psi \, dr\, dz.
$$

There remains one `1/r` factor multiplying the gradient term.  This is
the integrand singularity that any axisymmetric finite-element scheme
must cope with on axis-touching elements (`r → 0`).

## 4. Eddy currents add a σ-mass term

For a linear quasi-static eddy current with conductivity `σ(r, z)`,
Faraday's law gives `E = -∂A/∂t`, so the induced volume current is
`J_eddy = σ E = -σ ∂A_φ/∂t ĵ_φ`.  The corresponding term in the
weak form (multiplied by the test function and the `2π r` Jacobian) is

$$
\int_\Omega \sigma A_\varphi^{(\text{trial})} A_\varphi^{(\text{test})}
\cdot 2\pi r\, dr\, dz
= \int_\Omega \frac{\sigma}{2\pi r}\,\psi^{(\text{trial})} \psi^{(\text{test})}
   \cdot 2\pi r \, dr\, dz
= \int_\Omega \frac{\sigma}{2\pi r}\,\psi \psi'\, dr\, dz \cdot 2\pi
$$

Up to overall normalisation:

$$
M_\sigma[\psi, \psi'] = \frac{\sigma}{4\pi^2}
\int_\Omega \frac{\psi \psi'}{r^2}\, r\, dr\, dz
= \frac{\sigma}{4\pi^2}
\int_\Omega \frac{\psi \psi'}{r}\, dr\, dz.
$$

Or in the variable `s = r²` (see § 5), `dr = ds/(2r)`, so
`(1/r) dr dz = ds dz / (2 s)`:

$$
M_\sigma[\psi, \psi'] = \frac{\sigma}{8\pi^2}
\int_\Omega \frac{\psi \psi'}{s} \, ds\, dz.
$$

This is the σ-mass operator that `AxiHenrotteSigmaMassBFI` assembles.

## 5. Henrotte / Meeker change of variable: s = r²

Both the stiffness energy of § 3 and the σ-mass of § 4 contain a `1/r`
factor (or equivalently a `1/s` factor under `s = r²`).  Standard
finite-element shape functions polynomial in `(r, z)` with degree `p`
make the integrands

$$
\frac{(\partial_r \psi)^2}{r}\quad\text{and}\quad\frac{\psi^2}{r}
$$

singular as `r → 0`, which Gauss quadrature struggles to integrate
accurately.  In the limit, axis-touching elements lose effective
convergence rate from `p` to `p_eff = 1` even for `p ≥ 2`.

The Henrotte / Meeker observation [^Henrotte93][^MeekerFEMM]: choose
shape functions polynomial in `(s, z)` with `s = r²` instead.  Then

$$
\partial_r \psi = 2 r \, \partial_s \psi,
\qquad
(\partial_r \psi)^2 / r = 4 r \, (\partial_s \psi)^2 = 4 s^{1/2} \, (\partial_s \psi)^2,
$$

but the stiffness integrand transforms cleanly:

$$
\frac{(\partial_r \psi)^2}{r}\, dr\, dz
= 4 r (\partial_s \psi)^2 \, \frac{ds}{2r} \, dz
= 2 (\partial_s \psi)^2 \, ds \, dz,
$$

— the `r` and `1/r` cancel exactly, leaving a polynomial integrand in
`(s, z)` that integrates in **closed form** via Mathematica.

Similarly

$$
\frac{(\partial_z \psi)^2}{r}\, dr\, dz
= \frac{(\partial_z \psi)^2}{r} \cdot \frac{ds}{2r}\, dz
= \frac{(\partial_z \psi)^2}{2 s} \, ds\, dz,
$$

which keeps a `1/s` factor.  This second piece must be handled with
care on axis elements (`s_a = 0`) — see § 7.

## 6. Variational form on a single element

Let `Ω_e = [r_a, r_b] × [z_a, z_b]` be one axis-aligned rectangle, with
`s_a = r_a²`, `s_b = r_b²`.  Choose `ψ` to be polynomial in `(s, z)` of
total degree `p` (in NGSolve API: `H1Henrotte(mesh, order=p)`):

$$
\psi(s, z) = \sum_{m=1}^{N_p} c_m \, m(s, z),
$$

where the monomial set is

| `p` | `N_p` | Monomials `m(s, z)`                                  |
|----:|------:|------------------------------------------------------|
| 1   | 4     | `{1, s, z, s·z}`                                     |
| 2   | 9     | `{1, s, s², z, s·z, s²·z, z², s·z², s²·z²}`          |

The element stiffness in **monomial coefficients** is the symmetric
`N_p × N_p` matrix

$$
K_\phi^{(\text{mono})}_{m,n}
= \int_{s_a}^{s_b}\!\!\int_{z_a}^{z_b}
  \!\left[\frac{2}{\mu_z}\,\partial_s m \cdot \partial_s n
       + \frac{1}{2\mu_r s}\,\partial_z m \cdot \partial_z n\right]\!ds\, dz,
$$

where `μ_r, μ_z` are the (anisotropic) per-element permeabilities.  The
`1/s` in the second term is the residual `1/r` after the variable
change.  Both integrals are polynomial in `s` (for the first) and
admit a `log(s_b/s_a)` or `log(s_b)` closed form (for the second);
all entries are derived once symbolically and emitted as C++ code in
`q2_henrotte_generated.hpp` (auto-generated by
`packages/radia-axifemm/scripts/codegen_q2_henrotte.py`).

The σ-mass element matrix in monomial coefficients is

$$
M_{\sigma,\phi}^{(\text{mono})}_{m,n}
= \frac{\sigma}{4\pi^2}
  \int_{s_a}^{s_b}\!\!\int_{z_a}^{z_b}
    \frac{m(s,z)\, n(s,z)}{s}\, ds\, dz.
$$

Same `1/s` issue, same closed-form treatment.

### Conversion from monomial coefficients to nodal DOFs

The DOFs `V_j = A_φ` at each node form the API surface.  The monomial
coefficients are recovered from the DOFs via a per-element Vandermonde
inversion in the chosen monomial basis.  The conversion respects the
flux-function transformation `ψ = 2π r A_φ`:

$$
K_V = T \cdot V^{-T} \cdot K_\phi^{(\text{mono})} \cdot V^{-1} \cdot T,
\qquad
T = \mathrm{diag}(2\pi\, r_{\text{node}}),
$$

so that the assembled global system is in `A_φ`-DOFs at each mesh
vertex / edge midnode / face centre.  This conversion is done inside
`AxiHenrotteFE_Q2_AxisAligned::CalcShape` and the BFI integrators, not
by the user.

## 7. Axis-touching elements (`s_a = 0`)

When `r_a = 0`, the bottom edge of the element sits on the axis.  The
`1/s` factor in the second stiffness term and the σ-mass diverges as
`s → 0`.  The Henrotte basis avoids this by **dropping the
axis-incompatible monomials** that are non-zero at `s = 0`:

| `p` | Interior monomials | **Axis monomials** (axis-restricted basis)        |
|----:|--------------------|---------------------------------------------------|
| 1   | `{1, s, z, s·z}` (4) | `{s, s·z}` (2 — drops `1, z`)                     |
| 2   | 9 (above)          | `{s, s², s·z, s²·z, s·z², s²·z²}` (6 — drops the 3 monomials with no `s` factor) |

Geometrically, the dropped axis-side DOFs are forced to zero — they
correspond to `A_φ` not having an `s = const` mode at the axis (which
would imply `B_z` constant down to the axis with non-zero r component,
inconsistent with the assumption `A_φ ∝ r` near axis).  After
removing those DOFs the reduced `M_σ,φ` integrand has factors of `s` in
both `m` and `n` that cancel the `1/s`, restoring closed-form
integrability.

`AxiHenrotteFE_Q1_AxisAligned` and `AxiHenrotteFE_Q2_AxisAligned`
internally toggle between the interior basis (4 / 9 DOFs) and the
axis basis (2 / 6 DOFs) based on whether `s_a < ε`.  The dropped DOFs
are treated as a Dirichlet zero in NGSolve's free-DOF bookkeeping.

## 8. P1 triangle (FEMM `prob3big.cpp` direct port)

For unstructured triangle meshes — when axis-aligned quads are not
available — `axifemm` provides a P1-triangle fallback,
`AxiHenrotteFE_P1_Triangle` (3 DOFs/cell).  The basis is `{1, r², z}`
on each triangle, with shape functions linear in those variables (not
in `(r, z)`).  This reproduces FEMM's `StaticAxisymmetric()` exactly,
including the axis-touching cases (1- or 2-vertex on axis) that take
special-case stiffness formulas (see `axifemm_core.element_matrices`
in `tests/axifemm/_reference_python/` for the per-element reference
implementation).

The P1 triangle has lower per-DOF accuracy than Q2 quads
(`disk_convergence` example shows ~5–6 % gap to BEM-Foster reference
even at very fine meshes), so production use prefers Q2 quad whenever
a structured mesh is available.

## 9. Boundary conditions

Standard FEM boundary conditions in cylindrical coordinates:

- **Axis (`r = 0`)**: `A_φ = 0` is mathematically required (the
  vector potential has no φ-component on the symmetry axis).  This is
  imposed by Dirichlet on the axis edge label.
- **Far-field truncation (`r → ∞`)**: `A_φ → 0`, imposed by Dirichlet
  on the outer edges.  For applications needing a true open-domain
  treatment (no truncation), use Kelvin transformation (see
  `docs/kelvin/`) instead.
- **Symmetry plane (`z = 0` for an even-parity problem)**: `∂A_φ/∂z = 0`,
  Neumann (no boundary integral).

## 10. Summary of the assembled system

The complete weak form on the meridian plane Ω is

$$
\underbrace{K_V \cdot \mathbf{V}}_{\text{stiffness}}
+ \underbrace{M_\sigma \cdot \frac{d\mathbf{V}}{dt}}_{\text{eddy current}}
= \underbrace{\mathbf{F}}_{\text{source}},
$$

where `V` is the global vector of nodal `A_φ` DOFs.  In the API:

- `AxiHenrotteStiffnessBFI(mu_cf)` assembles `K_V` (per § 6, with
  `μ` from the CoefficientFunction).
- `AxiHenrotteSigmaMassBFI(sigma_cf)` assembles `M_σ` (per § 4 + 6,
  with `σ` from the CoefficientFunction).
- `H1Henrotte(mesh, order=p)` provides the appropriate FESpace
  (Q1 / Q2 quad, P1 triangle, with axis-element bookkeeping).
- The user assembles `K`, `M` with standard NGSolve `BilinearForm`
  and feeds them into any solver (direct, iterative, eigenvalue,
  Hiruma 3-term recurrence for Cauer ladder extraction).

## 10b. Heat-equation operator on the same FESpace (radia 4.31.0+)

The Henrotte basis is also the **mathematically natural function space
for axisymmetric scalar fields** like temperature `T(r, z)`.  The
reasoning is parity:

- Analytically extending an axisymmetric `T(r, z)` to `r < 0` gives
  `T(-r, z) = T(r, z)` — `T` is **an even function in r**.
- Even functions Taylor-expand only in even powers:
  `T = c_0 + c_2 r^2 + c_4 r^4 + ...`
- Standard P1 / P2 H1 basis on `(r, z)` includes `r`-linear and
  higher odd-`r` modes that **cannot occur in any axisymmetric
  distribution**.  The FE solver wastes DOFs trying to drive those
  spurious modes to zero through the `2 pi r` Jacobian weight.
- The Henrotte basis spans only even-`r` polynomials by construction
  (`{1, r², z}` for Q1, `{1, r², r⁴, z, r²z, r⁴z, ...}` for Q2),
  exactly matching the admissible function space.

The axisymmetric heat weak form (no eddy-current source for clarity):

$$
a_{\text{heat}}(T, v) = \int_\Omega k\,(\nabla T \cdot \nabla v)\,
                          \cdot 2\pi r\, dr\, dz
$$

After the change of variable `s = r²` (and the Jacobian
`2 pi r dr dz = pi ds dz` from § 5), this becomes

$$
a_{\text{heat}}(T, v) = \pi \int_\Omega k
   \big[\, 4 s\, \partial_s T\, \partial_s v
        + \partial_z T\, \partial_z v \,\big]\, ds\, dz.
$$

Note **no `1/s` factor** — heat is integrand-clean even on axis-
touching elements.  This is the architectural difference from the
magnetic stiffness (§ 6 has `1/(\mu_r s)` in the second term), and
it means **the full 9-monomial Q2 basis is admissible on axis-
touching elements** for heat (no axis basis reduction needed).

The transient heat capacity term is similarly clean:

$$
m_{\text{heat}}(T, v) = \pi \int_\Omega \rho c_p \, T \, v \, ds\, dz.
$$

Both element matrices are derived in closed form via sympy in
`packages/radia-axifemm/scripts/codegen_q_heat_henrotte.py` (or its
`C:/temp/codegen_axi_heat_henrotte.py` ad-hoc twin used during the
2026-05-10 derivation), emitted as
[`src/ext/axifemm/q_heat_henrotte_generated.hpp`](../../src/ext/axifemm/q_heat_henrotte_generated.hpp),
and assembled by the C++ integrators

- `radia.radia_axifemm.AxiHenrotteHeatStiffnessBFI(k_cf)` — `K_T`
- `radia.radia_axifemm.AxiHenrotteHeatMassBFI(rho_c_cf)` — `M_T`

DOF semantics are **nodal `T(vertex)`** directly — no
`T = diag(2 pi r_node)` flux-function transformation (T has no axis
boundary condition; standard nodal interpolation is correct).

Smoke test:
[`tests/axifemm/test_heat_henrotte_smoke.py`](../../tests/axifemm/test_heat_henrotte_smoke.py)
asserts that the assembled `K_T` is symmetric with one zero mode
(constant-T null space), `M_T` is SPD, and axis-touching elements
remain finite (no `1/s` divergence).

P1-triangle heat support is intentionally deferred to a follow-up;
production induction-heating workpiece thermal analyses use
structured quad meshes for accuracy anyway (Q2 quad converges much
faster than P1 triangle on the axis).

## 11. Cross-validation

The closed-form C++ implementation is cross-checked against three
independent paths:

1. **Pure-Python reference** — `tests/axifemm/_reference_python/`
   (axifemm_core for P1 triangle, axifemm_quad for Q1 quad,
   axifemm_quad_q2 for Q2 quad with Gauss 8×8 numerical quadrature).
   Test `tests/axifemm/test_python_reference_consistency.py` asserts
   the C++ stiffness eigenvalues match the Python prototype at
   machine precision for a single quad.

2. **Mathematica derivation** —
   `packages/radia-axifemm/scripts/validate_q2_codegen.py` runs the
   closed-form C++ vs the same monomial integrals re-derived in
   sympy at 50-digit precision.  Should match to ~1e-13 relative.

3. **Independent BEM-Foster reference** — the Cauer ladder time
   constants `τ_pair[k]` from a 1920-element BEM ring mesh
   (Mathematica, elliptic-integral kernel) are compared against
   the FE ladder via the Hiruma 3-term recurrence.  See
   `packages/radia-axifemm/tests/test_3way_cauer_cross_validation.py`
   and the table in [`AXIFEMM.md`](AXIFEMM.md#cauer-ladder-cross-validation-against-bem-nagamine-pipeline).

## 12. Why this matters for induction heating

Induction heating coils + axisymmetric workpieces (cylinders, disks,
tubes) generate **circumferential eddy currents** — exactly the case
where `A_φ` is the only nontrivial vector-potential component.  The
axifemm Henrotte / Meeker formulation:

- removes the standard P1 axis singularity in `B_r ∝ 1/r`,
- gives smooth per-element `B_z = const`, `B_r ∝ 1/r` — perfect for
  Cauer ladder time-constant extraction (the inductive transient of
  the workpiece is driven by these τ values),
- reduces the problem from 3-D to 2-D `(r, z)` without losing
  axisymmetric physics,
- runs at FEMM-grade accuracy per DOF, validated against FEMM `.mat`
  outputs to 0.1 % on the NMR benchmark
  (`examples/axifemm/nmr_validation/`).

This is why `radia.radia_axifemm` is a critical module of the radia
distribution and not a separate research package.

## References

[^Henrotte93]: F. Henrotte, B. Meys, A. Genon, W. Legros, "A new method
    for axisymmetric linear and nonlinear problems," *IEEE Transactions
    on Magnetics* **29**(2):1352–1355, March 1993.
    DOI [10.1109/20.250664](https://doi.org/10.1109/20.250664).

[^MeekerFEMM]: D. Meeker, *Finite Element Method Magnetics 4.2 — User's
    Manual + axisymmetric formulation notes*.  Reference C++ source:
    [`prob3big.cpp::StaticAxisymmetric()`](https://www.femm.info/wiki/Documentation).
