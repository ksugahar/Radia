"""
radia-core axifem knowledge base for the Radia + NGSolve MCP server.

`radia.axifem` in radia-core is an NGSolve add-on (registered FESpace
`axihenrotte`) that adds a Henrotte axisymmetric Q-element family which
NGSolve does not ship out of the box: polynomial basis in `s = r²` (not `r`)
representing the flux function `ψ = 2π r A_φ` (not `A_φ` itself). On the
Cu-disk eddy-current benchmark the `p=2` Q-element matches BEM-Foster to
0.27 % with ~6× fewer DOFs than NGSolve `H1 order=3`.

Read this when:
* Setting up an axisymmetric eddy-current / magnetostatic FEM problem with
  axis-touching elements.
* Comparing `axihenrotte p=k` to NGSolve `H1 order=k` and wondering why
  the latter loses 5–10 % accuracy near the symmetry axis.
* Cross-validating a Cauer-ladder / Hiruma 3-term workflow against BEM.
* Implementing a new FEM panel that needs the `1/r`-weighted axisymmetric
  weak form integrated in closed form rather than by Gauss quadrature.

The MCP server exposes this via axifemm_documentation(topic=...). Topics:
overview, api, basis_p1, basis_p2, vs_standard_h1, validation, kelvin,
file_layout, why_dropped_p3.
"""

AXIFEMM_OVERVIEW = """\
# radia-core axifem — Henrotte Axisymmetric Q-Element FE for NGSolve

## Status

| Order | DOFs / quad | Cu disk τ₁ vs BEM 224.31 µs (gap) |
|-------|-------------|------------------------------------|
| `p=1` Q-element | 4 (vertex only)            | 223.06 µs (0.55 %, fine mesh) |
| `p=2` Q-element | 9 (4 vertex + 4 edge + 1 face) | 223.69 µs (0.27 %, ~6× fewer DOFs than `p=1` for higher accuracy) |
| `p=3` Q-element | dropped — raw-monomial Vandermonde cond ≈ 10³⁰ exceeds double precision; would require shifted-Legendre basis to ship safely. Not on the roadmap. |

There is also a `p=1` triangle path (3 DOFs, a direct C++ axisymmetric
port) that matches a stored axisymmetric-magnet reference to 0.1 % on an
NMR-style benchmark.

## Key idea (why this is "an NGSolve feature NGSolve does not have")

NGSolve `H1(mesh, order=k)` represents `A_φ(r, z)` as a polynomial in `(r, z)`.
For the axisymmetric Maxwell weak form the volume element brings a factor of
`r dr dz` (toroidal Jacobian); the strong form contains `1/r` from
`B_z = (1/r) ∂(r A_φ)/∂r`. Result: the integrand mixes polynomials in `r`
with `1/r` weights, and Gauss quadrature loses 5–10 % accuracy near the
symmetry axis.

The Henrotte / Meeker trick:

* Substitute the **flux function** `ψ = 2π r A_φ`. Then
  `B_z = (1/(2π r)) ∂ψ/∂r = (1/π) ∂ψ/∂s` where `s = r²`.
  Pure polynomial in `(s, z)` — no `1/r` in the field formula.
* Choose the basis of `ψ` to be polynomial in `(s, z)`:
    p=1 (Q-element):  ψ ∈ span{1, s, z, s·z}             (4 DOFs)
    p=2 (Q-element):  ψ ∈ span{1, s, s², z, s·z, s²·z, z², s·z², s²·z²}  (9 DOFs)
* The `1/r`-weighted axisymmetric integrals become integrals of
  polynomials in `(s, z)` divided by `s` — closed-form integrable
  (Mathematica), with `log(s_b/s_a)` showing up for the truly
  log-singular components.
* After per-node `T = diag(2π r_node)` rescaling the assembled DOF
  vector stores `A_φ` directly, so user-facing eigenvalues
  (`τ_n = 1/λ_n`) come out in physical units with no extra factors.

NGSolve's standard FE library does not contain this basis (no FE family
polynomial in `s = r²`, no representation of `ψ` rather than `A_φ`).
`radia-core`'s `radia.axifem` module adds it via the public NGSolve add-on
pattern (custom `FESpace`, `DiffOp`, and `BilinearFormIntegrator`s), without
modifying NGSolve itself.
"""

AXIFEMM_API = """\
# Canonical API

```python
from ngsolve import Mesh, FESpace, BilinearForm, CoefficientFunction
from radia.axifem import AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI
import radia.axifem   # import once to register the FESpace

mesh = Mesh(...)                                    # axis-aligned quad mesh
fes  = FESpace("axihenrotte", mesh, order=2,
               dirichlet="axis|right|top|bot")      # p=2 Q-element

mu_cf    = CoefficientFunction(4 * 3.14159e-7)
sigma_cf = mesh.MaterialCF({"conductor": 5.8e7}, default=0.0)

a = BilinearForm(fes, symmetric=True)
a += AxiHenrotteStiffnessBFI(mu_cf)
a.Assemble()

m = BilinearForm(fes, symmetric=True)
m += AxiHenrotteSigmaMassBFI(sigma_cf)
m.Assemble()

# Solve K v = λ M v as you would with any NGSolve eigenproblem.
```

`FESpace("axihenrotte", mesh, order=k)` is the canonical entry point — it
follows the same pattern as `FESpace("h1ho", ...)`, `FESpace("hcurlho", ...)`,
etc. A back-compat helper `H1Henrotte(mesh, order=k, **flags)` is also
exported.

## Constraints

* **Mesh:** axis-aligned quadrilaterals only for `p=2` (the basis is
  polynomial in `s = r²` and `z` on `[s_a, s_b] × [z_a, z_b]`).  P1 triangle
  fallback exists at `p=1` for unstructured triangle meshes.
* **Dirichlet:** users MUST set `dirichlet="axis|..."` for any axis-touching
  geometry. The axis-restricted basis already produces zero rows / columns
  on the axis-side DOFs, so without an explicit Dirichlet flag the global
  matrix is singular.
* **Boundary 1D segments:** in `order=2` mode, a boundary segment exposes
  2 vertex DOFs and 1 edge midnode DOF, so `dirichlet="…"` on a boundary
  marks all three.
"""

AXIFEMM_BASIS_P1 = """\
# `p=1` Q-element (4 DOFs / quad)

Shape monomials in `(s = r², z)`:

```
{1, s, z, s·z}                    (interior)
{s, s·z}                          (axis-touching, sa = 0; the 2 axis-side
                                   nodes get zero shape functions)
```

DOFs: one per vertex. Lagrange-interpolatory at the 4 corners
`(s_a, z_a), (s_b, z_a), (s_b, z_b), (s_a, z_b)`.

The closed-form element matrices are derived in
`axifemm_quad.py:_element_matrices_quad_closed_form` and
`axifemm_quad.py:element_sigma_mass_quad`, validated against BEM v3 to
0.55 % on the Cu disk benchmark (τ₁ = 223.06 µs, ne = 15170, very-fine
mesh).
"""

AXIFEMM_BASIS_P2 = """\
# `p=2` Q-element (9 DOFs / quad)

Shape monomials in `(s = r², z)`:

```
interior (9): {1, s, s², z, s·z, s²·z, z², s·z², s²·z²}
axis     (6): {s, s², s·z, s²·z, s·z², s²·z²}
```

DOFs: 4 vertex + 4 edge midnode + 1 face center = 9 per quad.

## s-midpoint convention

Edge midnode positions use the **s-midpoint**, not the geometric
r-midpoint:

```
sm = (s_a + s_b) / 2            # arithmetic mean of s = r²
zm = (z_a + z_b) / 2
edge midnode at (sm, z_a)       # NOT at ((r_a + r_b)/2, z_a)
=> physical r_mid = √((r_a² + r_b²)/2)   ≠ (r_a + r_b)/2
```

Two adjacent elements with the same r-extent automatically agree on
edge-midnode placement and basis traces; no orientation flips are needed
for assembly continuity.

## Local DOF order (cyclic, with NGSolve edge permutation)

```
0..3:  4 corners (NGSolve QUAD vertex order: (sa,za),(sb,za),(sb,zb),(sa,zb))
4:     bottom edge midnode  (sm, za)    <- NGSolve edges[0]
5:     right  edge midnode  (sb, zm)    <- NGSolve edges[3]
6:     top    edge midnode  (sm, zb)    <- NGSolve edges[1]
7:     left   edge midnode  (sa, zm)    <- NGSolve edges[2]
8:     face center          (sm, zm)
```

NGSolve's QUAD local edge order is the *tensor-product* convention
(`bottom, top, left, right` = edges[0..3]); our Q-element local node
order is *cyclic CCW* (`bottom, right, top, left`). The FESpace's
`GetDofNrs` permutes accordingly. (Anyone refactoring the FESpace must
preserve this permutation, or assembly will silently produce wrong
edge-coupling and the eigenvalue spectrum collapses by a factor that
looks coincidentally like π — see project history.)

## Closed-form integration

```
K_phi[i,j]       = ∫ ∂m_i/∂s · ∂m_j/∂s / (π μ_Z) ds dz
                 + ∫ ∂m_i/∂z · ∂m_j/∂z / (4 π μ_R s) ds dz
M_sigma_phi[i,j] = σ / (4 π) · ∫ m_i · m_j / s ds dz
```

(Hessian-of-W convention. The 1/s integrand in M_sigma and in the W_r
component of K_phi produces `log(s_b/s_a)` for the constant component
and pure polynomials elsewhere; the axis case has all monomials carrying
an `s` factor so 1/s is integrable as polynomial.)

Per-element transformation to V-DOF (= A_φ at node):

```
K_V = T · V⁻ᵀ · K_phi · V⁻¹ · T,    T = diag(2 π r_node)
M_V = T · V⁻ᵀ · M_sigma_phi · V⁻¹ · T
```

with `V⁻¹` the (cached) per-element inverse Vandermonde of the
monomial basis at the 9 nodes (6 for axis case).
"""

AXIFEMM_VS_STANDARD_H1 = """\
# `axihenrotte p=2` vs NGSolve `H1 order=2` — the deeper differences

Beyond the obvious "the basis is in `r²`":

| Property | NGSolve `H1(mesh, order=2)` | `axihenrotte p=2` |
|---|---|---|
| Polynomial variable | `(r, z)` | `(s = r², z)` ⇒ degree 4 in `r` |
| Represented field | `A_φ` directly | `ψ = 2π r A_φ` (flux function); `A_φ` recovered via `/(2π r)` |
| Axis (`r = 0`) behaviour | `A_φ(0) ≠ 0` in general — must be Dirichlet'd explicitly. The `1/r` weight in axisymmetric integrals spoils Gauss-quadrature convergence near the axis (FEMM tradition reports 5–10 % loss). | `ψ` carries an `s = r²` factor in the axis-restricted basis, so `A_φ = ψ/(2π r) → 0` automatically as `r → 0`. No `1/r`-weighted Gauss quadrature; integrals are closed form. |
| Field formula | `B_z = (1/r) ∂(r A_φ)/∂r` (contains `1/r`) | `B_z = (1/π) ∂ψ/∂s` — pure derivative, no singularity |
| Element-matrix integration | NGSolve Gauss quadrature with `1/r` weight | Closed-form (Mathematica-derived per-monomial integrals, machine precision) |
| Effective convergence rate | Often `p_eff = 1` near the axis owing to weight singularity | Full `p = 2` in `(s, z)` |

Empirical comparison on the Cu disk eddy-current benchmark
(R = 10 mm, t = 2 mm, σ = 5.8 × 10⁷ S/m, BEM-Foster reference τ₁ = 224.31 µs):

| FE                                  | DOFs   | τ₁ (µs)   | Gap to BEM |
|-------------------------------------|--------|-----------|------------|
| NGSolve `H1 order=3` (v23 reference)| —      | 218.62    | 2.5 %      |
| `axihenrotte p=1` (Q1, fine mesh)   | 14904  | 223.06    | 0.55 %     |
| **`axihenrotte p=2` (Q2)**          | **9919** | **223.69** | **0.27 %** |

`axihenrotte p=2` outperforms NGSolve `H1 order=3` with fewer DOFs — that is
the value of the axisymmetric-specific basis.
"""

AXIFEMM_VALIDATION = """\
# Cross-validation references (per-element, machine precision)

* `axifemm/axifemm_quad.py` — validated Python prototype for `p=1`, gives
  τ₁ = 223.06 µs on the Cu disk (matches BEM v3 to 0.55 %).
* `axifemm/axifemm_quad_q2.py` — Python `p=2` Gauss-8×8 prototype; agrees
  with the Mathematica closed form to ~ 3.4 × 10⁻⁸ relative.
* `examples/axifemm/research/validate_q2_codegen.py` — runs both at the
  per-entry level after every `derive_quad_q2_henrotte.wls` re-run.
* `examples/axifemm/research/verification/test_q2_single_element.py` — end-to-end
  per-element BFI eigenvalue match (closed-form vs Gauss prototype, 1e-7
  tolerance).
* `examples/axifemm/research/verification/test_hiruma_disk_q2.py` — full disk Cu disk
  eddy-current Hiruma 3-term, expects τ₁ ≈ 223.7 µs.

## Cauer-ladder cross-validation against BEM (Phase 3-(3), Nagamine pipeline)

The Cauer ladder for the eddy-current problem follows Nagamine et al. 2026:

  in --R_0--+--R_2--+--R_4--+- ...
            |       |       |
           L_1     L_3     L_5  ...
            |       |       |
           gnd     gnd     gnd

R_{2k} (k = 0, 1, 2, ...) are the *series* resistors (even subscripts) and
L_{2k+1} are the *shunt* inductors (odd subscripts). Per-pair time constant:
    tau_pair[k] = L_{2k+1} / R_{2k}.

Two paths to the same ladder:

  (A) Nagamine BEM-Foster pipeline (independent reference):
      Mathematica bem_disk_axisym_cauer.wls computes 50 Foster eigenvalues
      and amplitudes on a 1920-element ring mesh, 20 alpha_n moments are
      derived, and Python disk_bem_cauer.py applies a 50-digit mpmath
      classical Cauer extraction. This is the mathematical equivalent of
      Nagamine's QD + equivalence-transform pipeline (Fig. 5 of his
      paper); we do NOT implement the verified-interval-arithmetic
      part, so our values are high-precision floats, not interval-
      rigorous bounds.

  (B) Differential-equation Henrotte FE + Hiruma 3-term Lanczos:
      C++ axihenrotte at order=1 / order=2, with Hiruma's 3-term
      recurrence reading off Cauer R, L coefficients directly:
          lambda_{2k+1} = w_{2k+1}^T K w_{2k+1} = 1 / R_{2k}    (conductance)
          lambda_{2k+2} = w_{2k+2}^T M w_{2k+2} =     L_{2k+1}  (inductance)
      so tau_pair[k] = lambda_{2k+1} * lambda_{2k+2}.

The Foster-amplitude normalisation differs between BEM and FE, so the
absolute R, L values differ by a common scale factor. The ratio
tau_pair[k] = L_{2k+1}/R_{2k} is normalisation-invariant and is the
comparison endpoint:

```
k   BEM Cauer    p=2 fine    p=1 very-fine    p=2/BEM gap    p=1/BEM gap
0   219.32 us    218.71      218.05          -0.28 %        -0.58 %
1    78.65       78.12        77.77          -0.68 %        -1.12 %
2    40.04       39.54        39.37          -1.24 %        -1.66 %
3    23.74       23.16        23.14          -2.46 %        -2.54 %
4    17.07       16.07        16.06          -5.86 %        -5.91 %
5    14.70       13.12        13.01         -10.77 %       -11.50 %
```

axihenrotte p=2 beats axihenrotte p=1 at every k (closer to the BEM
Cauer reference). The high-mode (k >= 4) divergence is the combined
effect of FE basis-order error at higher modes plus numerical
conditioning of the Cauer extraction at high stages (BEM itself starts
producing negative tau for k >= 6, addressed in Nagamine's verified-
interval pipeline which we have not implemented).

Reference:
  Nagamine, Yamaguchi, Sugahara, Hiruma, Mifune, Matsuo, "Verified
  Numerical Computations of the Cauer Network Representation of a Square
  Prism Conductor", manuscript 2026-05-04 (Japan J. Industrial Appl.
  Math. submission).

Test: `examples/axifemm/research/verification/test_3way_cauer_cross_validation.py`
Reference data (separate working tree, not in this repo):
  W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/
    bem_disk_axisym_cauer.wls     (Mathematica BEM + Foster amplitudes)
    disk_bem_cauer.py             (Python mpmath Cauer-I CFE)
    bem_disk_axisym_cauer.json
    bem_disk_axisym_cauer_python_results.json

## Hiruma 3-term ≠ Stoll χ-Foster Cauer-I (verified 2026-05-10)

The Hiruma 3-term Lanczos recurrence used here computes Cauer-I rungs of
  f_H(s) = bᵀ·(K - sM)⁻¹·b
which is the **Krylov-Padé impedance** representation. This is NOT the
same generating function as Sugahara/Kameari accumulation, which expands
  f_K(s) = uᵀ·M·(sM - K)⁻¹·M·u  with u = M⁻¹·b
giving the **χ-Foster (susceptibility)** Cauer-I — the one that matches
analytical Stoll Bessel for sphere and Stoll-equivalent BEM-Foster.

Sphere ground truth (Cu a=10 mm, B₀=1 T uniform; Mathematica Hankel-Padé
240 digits on Stoll spectrum):
  k=0  Stoll/Kameari τ = 694.142 μs  Hiruma 3-term τ = 728.85  (+5.0%)
  k=1  Stoll/Kameari τ = 154.604 μs  Hiruma 3-term τ = 171.51  (+10.9%)
  k=2  Stoll/Kameari τ =  64.075 μs  Hiruma 3-term τ =  71.68  (+11.9%)

Implication for axifemm: the τ₁ ≈ 223.7 μs cylinder reference (Hiruma) is
in the **Hiruma convention** and should not be mixed with BEM-Foster Cauer
values (which are in the Stoll/χ-Foster convention) without conversion.

To convert Hiruma matrices → Kameari/Stoll Cauer-I, swap the moment
formula:
  α_n^Kameari = uᵀ M (K⁻¹M)ⁿ u   (u = M⁻¹·b)
  α_n^Hiruma  = bᵀ (K⁻¹M)ⁿ K⁻¹ b
then run identical Hankel-Padé continued-fraction extraction.

POLICY (2026-05-10): Project repository (CauerLadderNetwork) has dropped
Hiruma in favour of Kameari accumulation as the canonical extractor. The
axifemm package retains its Hiruma routine for legacy verification, but
**new tests should target Kameari accumulation** so τ values are directly
comparable to Stoll analytical and BEM-Foster.

## Henrotte + Kelvin + CLN — workflow composition

The full eddy-current workflow on axisymmetric problems composes three pieces;
each is documented in its own knowledge file but the *combination* is the
intended canonical use:

  Henrotte axisym FE basis  +  z-offset Kelvin (or finite air box)  +  CLN extraction

* **Henrotte basis** (this file): polynomial in `s = r²` for `ψ = 2π r A_φ`.
  Produces `K`, `M` matrices in axisym geometry without `1/r`-Gauss errors
  near the axis.

* **Open-domain truncation: Kelvin works (Phase B3, 2026-05-12)**.
  The canonical 2D axisym z-offset Kelvin transformation now works
  end-to-end with `axihenrotte`. The historical blocker was that the
  `AxiHenrotteFESpace` did not expose enough `GetDofNrs(...)` overloads
  for `ngsolve.Periodic` to identify the `kelvin_int ↔ kelvin_ext` DOF
  pairs; that gap was closed by adding the `NodeId` /
  `GetVertexDofNrs` / `GetEdgeDofNrs` / `GetFaceDofNrs` overloads, after
  which the same recipe used for `H1` axisym Kelvin works for
  `H1Henrotte` (= `FESpace("axihenrotte", ...)`):

      fes = Periodic(H1Henrotte(mesh, order=2,
                                dirichlet="axis|axis_ext",
                                dirichlet_bbnd="GND"))
      mu_factor = kelvin_mu_factor_axisym_cf(z_offset=z_off, R=R_K)
      mu_cf     = build_material_cf(mesh, MU0, mu_factor,
                                    outer_keyword="kelvin")

  Verified on Cu sphere R=10 mm (Stoll Bessel ground truth):

  | mesh / curve            | τ₁ (µs) | gap to Stoll |
  |-------------------------|---------|--------------|
  | Stoll analytical        | 738.48  | —            |
  | axifemm p=2 + Kelvin    | 738.47  | -0.001 %     |
  | + Curve(2)              | 738.69  | +0.028 %     |

  This is the new canonical configuration for sphere/disk/cuboid axisym
  Cauer-ladder validation — see `kelvin` topic for full recipe and the
  documented `Periodic`-wrapping caveats.

  The finite air-box truncation (`A_φ = 0` on outer rectangle) remains
  valid and is still useful when you want to avoid building the Kelvin
  half-disc; convergence is just much slower (need R_air, Z_air ≈
  25×–50× the conductor extent for < 1 %).

* **CLN extraction** (`cln_3d.py` + `cln_notebooks/`): once you have
  `K, M, b` from the Henrotte assembly, choose a convention:
  - Hiruma 3-term Lanczos (Krylov-Padé impedance), or
  - Kameari accumulation (χ-Foster susceptibility, project canonical
    since 2026-05-10),
  and run the same Hankel-Padé continued-fraction Cauer extraction.

Cross-validation reference (2026-05-10, Cu disk R=10 mm, t=2 mm,
σ=5.8e7, B₀=1 T, leading Cauer rung τ_pair[0]):

| method                              | τ_pair[0] [μs] | gap to BEM |
|-------------------------------------|----------------|------------|
| BEM-Foster v3 (Nagamine pipeline)   | 219.32 (ref)   | —          |
| axihenrotte p=2 (Hiruma 3-term)     | 218.71         | -0.28 %    |
| axihenrotte p=1 very-fine (Hiruma)  | 218.05         | -0.58 %    |
| 3D HCurl (NGSolve + Kelvin)         | ≈ 218.7        | < 1 %      |
| Cylinder axisym VIM (144 cells)     | 211.85         | -3.4 %     |

The Cylinder VIM is run at a single coarse grid as a sanity check, not for
high accuracy; the axifemm `p=2` Q-element remains the recommended
production path for axisym Cu eddy-current Cauer extraction.

The same workflow on Cu sphere R=10 mm (Stoll Bessel ground truth, τ₁ = μ₀ σ R² / π² = 738.48 µs):

| method                              | τ_pair[0] [μs] | gap to Stoll |
|-------------------------------------|----------------|--------------|
| Stoll analytical (μ₀ σ R²/π²)       | 738.48 (ref)   | —            |
| **axifemm p=2 + z-offset Kelvin**   | **738.47**     | **-0.001 %** |
| axifemm p=2 + Kelvin + Curve(2)     | 738.69         | +0.028 %     |
| 3D HCurl (NGSolve + Kelvin)         | ≈ 694          | 0.027 % at L₁ (Stoll τ=694 convention) |
| Sphere axisym VIM (480 cells)       | 708.4          | +2.07 %      |

The axifemm + Kelvin result is the closest to Stoll across all available
axisym/3D methods on this benchmark (machine-precision agreement on the
leading rung). See `axifemm_documentation(topic="kelvin")` for the full
canonical recipe (Phase B3, commit 81f6415f).

## Hessian-of-W convention (load-bearing)

Both `K_phi` and `M_sigma_phi` are emitted in the **Hessian-of-W**
convention (matches `axifemm_quad.py:274`). The original
`derive_quad_q2_henrotte.wls` shipped with a coefficient-of-W convention
that was 2× too small for `K` and 2π× too small for `M`; this combined
multiplicatively to make the disk τ₁ off by ~ π and was diagnosed during
Phase A2 (commit 077e7b03). Anyone re-deriving the wls must keep the
Hessian convention, or every assembled quantity rescales by these
factors silently.
"""

AXIFEMM_FILE_LAYOUT = """\
# File layout

(the former standalone axifemm add-on was dissolved into radia-core on
2026-06-14: the C++ ships in the radia-core wheel as `radia.axifem`;
scripts/tests moved into the repo tree.)

```
src/ext/axifemm/                          # C++ source (built into the radia wheel)
  axi_henrotte_fe.{hpp,cpp}               # Q1, Q2 quad + P1 triangle FE classes
  axi_henrotte_fespace.{hpp,cpp}          # FESpace with order=1 / order=2 dispatch
  axi_henrotte_diffop.hpp                 # DifferentialOperators (value, gradient)
  axi_henrotte_integrators.{hpp,cpp}      # closed-form K and σ-mass BFI
  q2_henrotte_generated.hpp               # auto-generated, do not edit
  axifem.cpp                       # pybind11 entry -> radia.axifem
examples/axifemm/research/                # derivation + codegen scripts
  codegen_q2_henrotte.py                  # JSON → C++ codegen
  validate_q2_codegen.py                  # closed-form vs Gauss prototype
  q2_henrotte_test_values.json            # numerical reference values
examples/axifemm/research/verification/   # standalone __main__ verification scripts
  test_q2_single_element.py               # per-element BFI sanity check
  test_hiruma_disk_q1.py                  # disk Hiruma 3-term, p=1
  test_hiruma_disk_q2.py                  # disk Hiruma 3-term, p=2
  test_q2_assembly_diag.py                # 2-quad assembly diagnostic
tests/axifemm/                            # pytest golden tests (CI-collected)
  test_element_matrices.py, test_heat_*.py, test_python_reference_consistency.py
```

The Mathematica derivation lives upstream at
`W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/`
(`derive_quad_q2_henrotte.wls` and `quad_q2_henrotte_matrices.json`).

Comprehensive design + theory documentation: `docs/axifemm/AXIFEMM.md`.
"""

AXIFEMM_KELVIN = """\
# `axihenrotte` + z-offset Kelvin transformation (Phase B3, 2026-05-12)

The 2D axisymmetric z-offset Kelvin transformation works end-to-end with
the Henrotte basis as of the radia-core axifem Phase B3 work (commits 81f6415f /
6e963ab9). Same boilerplate as the canonical `H1` axisym Kelvin recipe
(`kelvin.a_formulation` topic) — substitute `H1Henrotte` for `H1` and
add a `Periodic(...)` wrap.

## Why this was blocked before Phase B3

`ngsolve.Periodic` traverses the `FESpace::GetDofNrs(NodeId, ...)` API
to identify the slave/master DOFs across the `kelvin_int ↔ kelvin_ext`
boundary. The original `AxiHenrotteFESpace` only implemented the
ElementId overload, so wrapping it in `Periodic` raised
`NotImplementedError` from the C++ side. Phase B3 added:

```
GetDofNrs(NodeId node, Array<DofId>& dnums)
GetVertexDofNrs(int vnr, Array<DofId>& dnums)
GetEdgeDofNrs(int enr, Array<DofId>& dnums)
GetFaceDofNrs(int fnr, Array<DofId>& dnums)
```

After this the same `Periodic` mechanism that works on stock `H1` works
on `H1Henrotte`, and the canonical Kelvin recipe carries over verbatim.

## Canonical recipe

```python
from ngsolve import Mesh, BilinearForm, Periodic, TaskManager
from netgen.occ import OCCGeometry, WorkPlane, Glue, MoveTo
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)
from radia.panels.add_kelvin import add_kelvin_2d_axisym
from radia.kelvin_source import (
    kelvin_mu_factor_axisym_cf, build_material_cf,
)

# 1. Build interior half-disc on x >= 0:
#    - conductor face (material "conductor"), bnd named "axis"
#      on x=0, anything else on the conductor boundary
#    - air_inner annulus (material "air_inner")
#    - outer arc edge of air_inner MUST be named "kelvin_int" before
#      passing to add_kelvin_2d_axisym
interior = build_my_interior_half_disk(...)   # user-defined

# 2. Append the Kelvin half-disc with z_offset.  Returns a Glue'd shape +
#    info dict with axis_labels = "axis|axis_ext" and the periodic pair
#    table already wired up.  GND vertex is at (0, z_offset).
shape, info = add_kelvin_2d_axisym(
    interior,
    R=R_K,                # Kelvin inversion radius
    z_offset=5 * R_K,     # see "z_offset = 5·R_K" rationale below
    maxh_kelvin=5e-3,
)
mesh = Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_kelvin))
# Optional: mesh.Curve(2)  -- see "Curve(2) trade-off" below

# 3. FESpace: H1Henrotte + Periodic.  Dirichlet on both axes + GND vertex.
fes = Periodic(H1Henrotte(
    mesh, order=2,
    dirichlet="axis|axis_ext",
    dirichlet_bbnd="GND",
))

# 4. Material coefficient: mu_0 inside, mu_0 * (rho'/R_K)^2 in the
#    "kelvin" region (NOTE the helper returns mu_factor, not nu_factor;
#    AxiHenrotteStiffnessBFI multiplies by mu, not nu).
mu_factor = kelvin_mu_factor_axisym_cf(z_offset=5 * R_K, R=R_K)
mu_cf = build_material_cf(mesh, MU0, mu_factor, outer_keyword="kelvin")
sigma_cf = mesh.MaterialCF({"conductor": SIGMA_CU}, default=0.0)

# 5. Assemble K and M.  check_unused=False is required because the
#    Periodic wrap introduces master DOFs that have no "owning" element.
a = BilinearForm(fes, symmetric=True, check_unused=False)
a += AxiHenrotteStiffnessBFI(mu_cf)
m_bf = BilinearForm(fes, symmetric=True, check_unused=False)
m_bf += AxiHenrotteSigmaMassBFI(sigma_cf)
with TaskManager():
    a.Assemble()
    m_bf.Assemble()

# 6. Extract the conductor subblock via Schur complement (air DOFs have
#    diag(M) = 0; they are eliminated to give S = K_cc - K_ca K_aa^-1 K_ac).
#    Solve generalised eigenproblem S v = lambda M_cc v, take tau = 1/lambda_min.
```

## Verified results (Cu sphere R = 10 mm, σ = 5.8×10⁷ S/m)

```
Stoll analytical                       τ₁ = 738.48 µs   (ref)
axihenrotte p=2 + Kelvin               τ₁ = 738.47 µs   -0.001 %
axihenrotte p=2 + Kelvin + Curve(2)    τ₁ = 738.69 µs   +0.028 %
```

Test: `examples/CLN/scripts/axifemm/test_p2_kelvin_sphere.py`.

## Critical gotchas

### `check_unused=False` on both BilinearForms
With `Periodic(fes)` wrapped around `H1Henrotte`, slave DOFs appear in
the FESpace ordering but are not owned by any element of the inner
space; without `check_unused=False` NGSolve raises "unused DOF" errors
during assembly. Both `a` (stiffness) and `m_bf` (mass) need the flag.

### `mu_factor` (not `nu_factor`) — convention difference vs `H1` recipe
The standard `H1` axisym Kelvin recipe in `kelvin.a_formulation` uses
`nu_kelvin = nu_0 * (rho'/R)^2` (reluctivity). `AxiHenrotteStiffnessBFI`
multiplies the integrand by `mu`, not `nu`, so the *correct* factor for
the Kelvin region is `mu_kelvin = mu_0 * (R/rho')^2` —
`kelvin_mu_factor_axisym_cf` returns exactly this (note `(R/rho')^2`,
NOT `(rho'/R)^2`). If you accidentally pass `kelvin_nu_factor_axisym_cf`
you get `(rho'/R)^2` and τ blows up by orders of magnitude near rho'=0.

### Element-centroid `mu` sampling — pick `z_offset = 5·R_K`
`AxiHenrotteStiffnessBFI` samples `mu` at the *element centroid*
(`axi_henrotte_integrators.cpp:710`), not per quadrature point. Since
`mu(rho')` diverges as `rho' → 0` near the GND vertex `(0, z_offset)`,
fine mesh near GND can over-sample large values. Putting the Kelvin
center at `z_offset = 5·R_K` keeps the GND vertex 4 Kelvin radii away
from the conductor and bounds the centroid sampling error well below
0.1 %. Smaller offsets (e.g. `2·R_K`) degrade results to the percent
level; larger offsets work but waste mesh on a region with vanishing
field. The Phase B3.2 follow-up would be a per-quadpt `mu` integrator
overload, but the centroid path already hits machine precision on the
sphere benchmark, so this is low priority.

### `Curve(2)` trade-off
NGSolve `mesh.Curve(2)` on the spherical conductor boundary improves
geometric accuracy 30× on coarse meshes (gap 0.867 % → 0.028 % on a
typical mesh), but on a *finer* mesh `Curve(2)` actually slightly
*increases* the gap (0.001 % straight → 0.028 % curved) because the
biquadratic Lagrange in (r, z) — not (s, z) — does not match the
Henrotte basis's even-power-of-r structure on curved edges. Use
`Curve(0)` (= straight edges, default) for sphere/disk benchmarks at
typical mesh resolution; reserve `Curve(2)` for coarse-mesh production
runs where the geometric error dominates the FE-basis error.

A proper curved Q2 element (basis space still in even powers of r, but
on a biquadratic isoparametric mapping) is Phase B5 work — currently a
Python prototype at
`examples/CLN/scripts/axifemm/axifemm_quad_q2_curved.py`, C++ port deferred.

### Geometry pre-requisites for `add_kelvin_2d_axisym`
Before calling the helper, the user-built interior must have:

  - all conductor edges on `x = 0` named exactly `"axis"`
  - all `air_inner` edges on `x = 0` named exactly `"axis"`
  - the outer arc of `air_inner` named exactly `"kelvin_int"`
  - no other edges named `"kelvin_int"` or starting with `"axis_"`

The helper then appends a mirror half-disc at `(0, z_offset)`, names its
arc `"kelvin_ext"`, identifies `kelvin_int ↔ kelvin_ext` as Periodic
pairs, names the new x=0 edges `"axis_ext"`, and registers the
`(0, z_offset)` vertex as `"GND"`. Mis-naming any of these breaks
either the Periodic identification (most common) or the Dirichlet trace
(silent — eigenvalues come out finite but wrong).

## What changed in radia-core axifem source

```
src/axi_henrotte_fespace.{hpp,cpp}
    + GetDofNrs(NodeId, Array<DofId>&)
    + GetVertexDofNrs / GetEdgeDofNrs / GetFaceDofNrs
src/axifem.cpp
    + pybind exports of the new overloads
```

Test: `examples/axifemm/research/verification/test_magnetized_sphere.py`
(the sphere verification; passes against Stoll to -0.001 %, replaces the
prior expected-failure xfail in Phase B2).
"""


AXIFEMM_WHY_DROPPED_P3 = """\
# Why `p=3` was attempted and dropped

`p=3` Q-element would have 16 DOFs / quad (4 corner + 8 edge + 4 face) and
in principle drive the Cu-disk gap below 0.1 %.

## Implementation status

The wls Mathematica derivation (256 + 256 + 144 + 144 = 800 entries) was
written, run, and codegen'd successfully (Phase A3, commit a84d53f1; later
reverted in 56f451fd). C++ FE class, BFI dispatch, and FESpace order=3
topology were all wired up cleanly — the implementation was structurally
correct.

## Why it was reverted

The 16×16 Vandermonde of the raw `{s^a z^b}` (a, b ∈ {0..3}) basis at
typical disk-mesh elements (s ≈ 10⁻⁷–10⁻⁶) has condition number ≈ 10³⁰,
which exceeds double precision (~10¹⁶). Gauss-Jordan inversion produces
noise that contaminates assembled K_V and M_V on every refinement: the
*coarse* Cu-disk mesh happened to give τ₁ = 223.67 µs (looks right
because the dominant mode is robust to noise), but mesh refinement made
results worse, eventually completely garbage at fine.

Empirically:

```
cond(V_4x4)  on p=1 mesh = 1e+2    ← Q-element p=1 fully OK
cond(V_9x9)  on p=2 mesh = 1e+8    ← p=2 borderline but persistent for typical meshes
cond(V_16x16) on p=3 mesh = 1e+30  ← p=3 broken before the BFI even runs
```

## What would un-block `p=3`

Switch the basis from raw monomials `{s^a z^b}` to a normalized /
orthogonal form. Three options, in increasing implementation cost:

1. **Normalised monomials** — substitute `ξ = (s − s_a)/(s_b − s_a)`,
   `η = (z − z_a)/(z_b − z_a)`, then redo the wls in `(ξ, η)`. Vandermonde
   condition drops to ~100. Same code structure, only the Mathematica
   integrals change shape (`1/(s_a + (s_b − s_a) ξ)` for the 1/s factor).
2. **Tensor-product 1D Lagrange** (recommended) — bypass Vandermonde
   inversion entirely. Compute four 4×4 Mathematica matrices `K_s, N_s,
   M_z, K_z`, then assemble `K_phi[(a,b),(a',b')] = K_s[a,a'] M_z[b,b'] /
   (πμ_Z) + N_s[a,a'] K_z[b,b'] / (4πμ_R)`. The 1D matrices come from
   integrating products of explicit 1D Lagrange polynomials (no inversion
   anywhere). Generalises trivially to `p=4..k`.
3. **Shifted Legendre on `[s_a, s_b] × [z_a, z_b]`** — orthogonal basis,
   most robust, and already standard in NGSolve's hierarchical FE
   internals; full re-derivation needed.

The pipeline (wls + codegen + FE class + BFI dispatch + FESpace topology)
is in place; only the *basis choice* needs replacing.

## Decision: `p=2` is the sweet spot

Given that `axihenrotte p=2` already beats NGSolve `H1 order=3` per DOF
on the canonical benchmark and reaches 0.27 % gap to BEM, `p=3` is not
on the roadmap. The relevant trade-off is "more accuracy per DOF" not
"more polynomial order"; for the latter, switching to BEM-Foster (already
shipped in `ngsolve_validation/bem_disk_axisym_*.wls`) is the appropriate
escape hatch.
"""


AXIFEMM_MAGNET = """\
# Permanent-magnet source term (axisymmetric magnetization edge-loop)

axihenrotte solves linear/eddy-current problems; a **permanent magnet** is
added purely as a `LinearForm` RHS — no FESpace/C++ change needed. This is the
standard axisymmetric magnetization edge-loop
(`H_c * (cos t * dr + sin t * dz)`) written as a continuous Galerkin form.

## Weak form

Magnetostatics with remanence: `curl(nu*B - nu*B_rem) = J`, so the magnet
contributes `RHS += int nu*B_rem . (curl v) dV`. In the A_phi convention
(DOFs store A_phi; `B_z = grad(u)[0] + u/r`, `B_r = -grad(u)[1]`), with
`nu*B_rem = H_c*(cos th, sin th)` and `th` measured from the r-axis
(`th = 90 deg` => axial magnetization):

    f += [ H_c*sin(th)*(r*grad(v)[0] + v) - H_c*cos(th)*r*grad(v)[1] ] dx_magnet

Stiffness (static, sigma=0):
    a += nu*(1/r)*(r*grad(u)[0]+u)*(r*grad(v)[0]+v)*dx + nu*r*grad(u)[1]*grad(v)[1]*dx

```python
th = math.radians(theta_deg)
reg = mesh.MaterialCF({"magnet": 1.0}, default=0.0)
f += reg * Hc * math.sin(th) * (x * grad(v)[0] + v) * dx
f += -reg * Hc * math.cos(th) * (x * grad(v)[1]) * dx
```

## Validation (uniformly-magnetized linear sphere, analytical)

Sphere radius a, rel. perm. mu_r, coercivity Hc, axial mag.:
    B_in = 2*mu0*mu_r*Hc/(mu_r+2)   (uniform inside)   ;  external = dipole.
For mu_r=2, Hc=3e5 -> B_in = 0.376991 T. The verified test
(`tests/test_magnetized_sphere.py`) gives, on H1Henrotte p=2:
    <B_z>_magnet (vol-avg)  = 0.376802 T   (-0.050 %)
    interior |B_r|          < 1.1e-4 T     (purely axial, as expected)
    external equator r=5a   = -1.497e-3 T  vs dipole -1.508e-3 (0.73 %)
Cross-checks: standard NGSolve H1 gives -0.010 % (same source term);
a stored axisymmetric Kelvin reference uses the same a/mu_r/Hc.

## Gotchas

* **Axis MUST be Dirichlet** (`dirichlet="axis|outer"`). A_phi=0 on r=0 is
  physical; without it the Henrotte {1,r²,z} basis leaves A_phi(0)!=0 and
  `B_z = grad(u)[0]+u/r` blows up at the axis. Name the FULL r=0 line
  (both magnet and air faces' `edges.Min(X)`) "axis".
* **On-axis B recovery** (r->0) via the `grad(gfu)` CoefficientFunction is
  unreliable for H1Henrotte (the `AxiHenrotteDiffOpGradient::Apply/ApplyTrans`
  fall back to the base-class stub — see the runtime "called base class
  apply" message). Sample interior fields OFF-axis (r>=0.3a) or use the
  volume-averaged metric; the SOLVE itself is correct on-axis.
* For an OPEN domain use the Kelvin reluctivity warp (see the "kelvin"
  topic): `nu_ext = (rho'/R)^2 * nu_0`. Truncation (Dirichlet far box)
  distorts the external dipole by O((a/R_far)^3)+mesh error.
"""

def get_axifemm_documentation(topic: str = "all") -> str:
    """
    Return radia-core axifem documentation for the requested topic.

    Topics:
      "all"             - All sections concatenated
      "overview"        - What it is and why NGSolve doesn't already have it
      "api"             - FESpace("axihenrotte", mesh, order=k) canonical usage
      "basis_p1"        - `p=1` Q-element (4 DOFs) basis details
      "basis_p2"        - `p=2` Q-element (9 DOFs) basis details + s-midpoint convention
      "vs_standard_h1"  - 6-property table comparing axihenrotte p=2 vs NGSolve H1 order=2
      "validation"      - Cross-validation references and Hessian-of-W convention
      "kelvin"          - Phase B3 z-offset Kelvin recipe (Periodic + H1Henrotte,
                          sphere -0.001 % vs Stoll, with gotchas: mu vs nu factor,
                          element-centroid mu sampling, Curve(2) trade-off)
      "magnet"          - Permanent-magnet source term (axisym magnetization edge-loop):
                          weak-form RHS, magnetized-sphere validation, axis-BC gotcha
      "file_layout"     - Where to find each piece (C++, Mathematica, tests)
      "why_dropped_p3"  - Why `p=3` was attempted, completed, and reverted
    """
    sections = {
        "overview":       AXIFEMM_OVERVIEW,
        "api":            AXIFEMM_API,
        "basis_p1":       AXIFEMM_BASIS_P1,
        "basis_p2":       AXIFEMM_BASIS_P2,
        "vs_standard_h1": AXIFEMM_VS_STANDARD_H1,
        "validation":     AXIFEMM_VALIDATION,
        "kelvin":         AXIFEMM_KELVIN,
        "magnet":         AXIFEMM_MAGNET,
        "file_layout":    AXIFEMM_FILE_LAYOUT,
        "why_dropped_p3": AXIFEMM_WHY_DROPPED_P3,
    }
    if topic == "all":
        return "\n\n".join(sections[k] for k in [
            "overview", "api", "basis_p1", "basis_p2",
            "vs_standard_h1", "validation", "kelvin", "magnet",
            "file_layout", "why_dropped_p3"
        ])
    if topic in sections:
        return sections[topic]
    available = ", ".join(sections.keys())
    raise ValueError(
        f"Unknown axifemm topic '{topic}'. Available: {available}, all"
    )
