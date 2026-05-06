"""
radia-axifemm knowledge base for the Radia + NGSolve MCP server.

`radia-axifemm` is an NGSolve add-on (registered FESpace `axihenrotte`) that
adds a FEMM/Henrotte axisymmetric Q-element family which NGSolve does not
ship out of the box: polynomial basis in `s = r²` (not `r`) representing the
flux function `ψ = 2π r A_φ` (not `A_φ` itself). On the Cu-disk eddy-current
benchmark the `p=2` Q-element matches BEM-Foster to 0.27 % with ~6× fewer
DOFs than NGSolve `H1 order=3`.

Read this when:
* Setting up an axisymmetric eddy-current / magnetostatic FEM problem with
  axis-touching elements.
* Comparing `axihenrotte p=k` to NGSolve `H1 order=k` and wondering why
  the latter loses 5–10 % accuracy near the symmetry axis.
* Cross-validating a Cauer-ladder / Hiruma 3-term workflow against BEM.
* Implementing a new FEM panel that needs the `1/r`-weighted axisymmetric
  weak form integrated in closed form rather than by Gauss quadrature.

The MCP server exposes this via axifemm_documentation(topic=...). Topics:
overview, api, basis_p1, basis_p2, vs_standard_h1, validation, file_layout,
why_dropped_p3.
"""

AXIFEMM_OVERVIEW = """\
# radia-axifemm — Henrotte Axisymmetric Q-Element FE for NGSolve

## Status

| Order | DOFs / quad | Cu disk τ₁ vs BEM 224.31 µs (gap) |
|-------|-------------|------------------------------------|
| `p=1` Q-element | 4 (vertex only)            | 223.06 µs (0.55 %, fine mesh) |
| `p=2` Q-element | 9 (4 vertex + 4 edge + 1 face) | 223.69 µs (0.27 %, ~6× fewer DOFs than `p=1` for higher accuracy) |
| `p=3` Q-element | dropped — raw-monomial Vandermonde cond ≈ 10³⁰ exceeds double precision; would require shifted-Legendre basis to ship safely. Not on the roadmap. |

There is also a `p=1` triangle path (3 DOFs, FEMM `prob3big.cpp` direct
port) that matches FEMM `.mat` outputs to 0.1 % on the FEMM NMR benchmark.

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
`radia-axifemm` adds it via the public NGSolve add-on pattern (custom
`FESpace`, `DiffOp`, and `BilinearFormIntegrator`s), without modifying
NGSolve itself.
"""

AXIFEMM_API = """\
# Canonical API

```python
from ngsolve import Mesh, FESpace, BilinearForm, CoefficientFunction
from radia_axifemm import AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI
import radia_axifemm   # import once to register the FESpace

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
* `packages/radia-axifemm/scripts/validate_q2_codegen.py` — runs both at the
  per-entry level after every `derive_quad_q2_henrotte.wls` re-run.
* `packages/radia-axifemm/tests/test_q2_single_element.py` — end-to-end
  per-element BFI eigenvalue match (closed-form vs Gauss prototype, 1e-7
  tolerance).
* `packages/radia-axifemm/tests/test_hiruma_disk_q2.py` — full disk Cu disk
  eddy-current Hiruma 3-term, expects τ₁ ≈ 223.7 µs.

## Cauer-I cross-validation against BEM (Phase 3-(3))

Beyond the leading τ₁ Foster comparison, the Cauer ladder per-stage time
constants τ_rung[n] = L_n × λ_{2n-1} have been cross-checked against an
**independent integral-equation BEM pipeline**: Mathematica
`bem_disk_axisym_cauer.wls` computes 50 Foster eigenvalues and amplitudes
on a 1920-element ring mesh, 20 α_n moments are derived, and a Python
script applies a 50-digit mpmath Cauer-I CFE (repeated Taylor inversion)
to extract τ_rung[n].

```
n   BEM Cauer    p=2 fine    p=1 very-fine    p=2/BEM gap    p=1/BEM gap
1   219.32 us    218.71      218.05          -0.28 %        -0.58 %
2    78.65       78.12        77.77          -0.68 %        -1.12 %
3    40.04       39.54        39.37          -1.24 %        -1.66 %
4    23.74       23.16        23.14          -2.46 %        -2.54 %
5    17.07       16.07        16.06          -5.86 %        -5.91 %
6    14.70       13.12        13.01         -10.77 %       -11.50 %
```

The comparison is between two FORMULATIONS, both implemented in C++:
  (A) integral-equation BEM (Mathematica + mpmath Cauer CFE) -- independent
      reference;
  (B) differential-equation Henrotte FE (axihenrotte order=1/2) +
      Hiruma 3-term recurrence -- the order=1 / order=2 entries are a
      convergence study WITHIN the same FE solver, not two independent
      methods.

Both formulations agree on the leading 3 Cauer rungs to ~1 % and on the
leading rung to 0.28 %. axihenrotte p=2 is closer to BEM than axihenrotte
p=1 at every stage.

Test: `packages/radia-axifemm/tests/test_3way_cauer_cross_validation.py`
Reference data (separate working tree, not in this repo):
  W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/
    bem_disk_axisym_cauer.wls     (Mathematica BEM + Foster amplitudes)
    disk_bem_cauer.py             (Python mpmath Cauer-I CFE)
    bem_disk_axisym_cauer.json
    bem_disk_axisym_cauer_python_results.json

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

```
packages/radia-axifemm/
  src/
    axi_henrotte_fe.{hpp,cpp}            # Q1, Q2 quad + P1 triangle FE classes
    axi_henrotte_fespace.{hpp,cpp}       # FESpace with order=1 / order=2 dispatch
    axi_henrotte_diffop.hpp              # DifferentialOperators (value, gradient)
    axi_henrotte_integrators.{hpp,cpp}   # closed-form K and σ-mass BFI
    q2_henrotte_generated.hpp            # auto-generated, do not edit
    radia_axifemm.cpp                    # pybind11 module entry
    __init__.py                          # Python re-exports
  scripts/
    codegen_q2_henrotte.py               # JSON → C++ codegen
    validate_q2_codegen.py               # closed-form vs Gauss prototype
    q2_henrotte_test_values.json         # numerical reference values
  tests/
    test_q2_single_element.py            # per-element BFI sanity check
    test_hiruma_disk_q1.py               # disk Hiruma 3-term, p=1
    test_hiruma_disk_q2.py               # disk Hiruma 3-term, p=2
    test_q2_assembly_diag.py             # 2-quad assembly diagnostic
```

The Mathematica derivation lives upstream at
`W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/`
(`derive_quad_q2_henrotte.wls` and `quad_q2_henrotte_matrices.json`).

Comprehensive design + theory documentation: `docs/axifemm/AXIFEMM.md`.
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


def get_axifemm_documentation(topic: str = "all") -> str:
    """
    Return radia-axifemm documentation for the requested topic.

    Topics:
      "all"             - All sections concatenated
      "overview"        - What it is and why NGSolve doesn't already have it
      "api"             - FESpace("axihenrotte", mesh, order=k) canonical usage
      "basis_p1"        - `p=1` Q-element (4 DOFs) basis details
      "basis_p2"        - `p=2` Q-element (9 DOFs) basis details + s-midpoint convention
      "vs_standard_h1"  - 6-property table comparing axihenrotte p=2 vs NGSolve H1 order=2
      "validation"      - Cross-validation references and Hessian-of-W convention
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
        "file_layout":    AXIFEMM_FILE_LAYOUT,
        "why_dropped_p3": AXIFEMM_WHY_DROPPED_P3,
    }
    if topic == "all":
        return "\n\n".join(sections[k] for k in [
            "overview", "api", "basis_p1", "basis_p2",
            "vs_standard_h1", "validation", "file_layout", "why_dropped_p3"
        ])
    if topic in sections:
        return sections[topic]
    available = ", ".join(sections.keys())
    raise ValueError(
        f"Unknown axifemm topic '{topic}'. Available: {available}, all"
    )
