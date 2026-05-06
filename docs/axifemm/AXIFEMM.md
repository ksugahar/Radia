# radia-axifemm — Henrotte Axisymmetric Q-Element FE for NGSolve

`radia-axifemm` is a NGSolve add-on package that adds a **Henrotte/Meeker
axisymmetric finite-element family** to NGSolve. It registers a FESpace named
`axihenrotte` and provides closed-form `BilinearFormIntegrator`s; the rest of
the workflow (mesh, dirichlet, eigsh, Hiruma 3-term) is plain NGSolve.

The basis lives on **axis-aligned quadrilateral meshes**; it is polynomial of
order `p` in the variable `s = r²` (not `r`), with the function being
represented being `ψ = 2π r A_φ` (the magnetic flux function), not `A_φ`
itself. After a per-node `T = diag(2π r)` transformation the assembled DOF
vector stores `A_φ` at each node directly. This is the key trick of the
Meeker / Henrotte axisymmetric formulation [^Henrotte93][^MeekerFEMM]; without
it, standard H1 elements lose 5–10 % accuracy near the symmetry axis owing to
the `1/r` weight in the axisymmetric integrals.

NGSolve does **not** ship this basis. `radia-axifemm` adds it without
modifying NGSolve itself, following the public `ngsolve-addon-template`
pattern (a custom `FESpace`, custom `DiffOps`, custom `BilinearFormIntegrator`s).

## Status

| Order | DOFs / quad | Status | Cu disk τ₁ vs BEM 224.31 µs |
|-------|-------------|--------|------------------------------|
| `p=1` (Q-element)   | 4 (vertex only)                       | shipping | 222.63 µs (gap **0.55 %**, fine mesh ne=15170) |
| `p=2` (Q-element)   | 9 (4 vertex + 4 edge mid + 1 face)    | shipping | 223.69 µs (gap **0.27 %**, ne=2530, ~ 6× fewer elements than `p=1` for better accuracy) |
| `p=3` (Q-element)   | —                                     | dropped  | raw-monomial Vandermonde cond ≈ 10³⁰ exceeds double precision; production-grade `p=3` requires switching to a shifted-Legendre basis on `[s_a, s_b] × [z_a, z_b]`. Not on the roadmap. |

There is also a `p=1` triangle path (3 DOFs/triangle, FEMM `prob3big.cpp`
direct port) for compatibility with general-shaped triangle meshes; this
matches FEMM `.mat` outputs to 0.1 % on the FEMM NMR benchmark.

## Quick start (canonical API)

```python
from ngsolve import Mesh, FESpace, BilinearForm, CoefficientFunction
from radia_axifemm import AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI
import radia_axifemm   # registers the FESpace

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

A convenience wrapper `H1Henrotte(mesh, order=k, **flags)` is also exported
for back-compat, but the canonical entry point is the standard NGSolve
`FESpace("axihenrotte", ...)` factory.

## How `axihenrotte p=2` differs from NGSolve `H1 order=2`

| Property | NGSolve `H1(mesh, order=2)` | `axihenrotte p=2` |
|---|---|---|
| Polynomial variable | `(r, z)` | **`(s = r², z)`** ⇒ degree 4 in `r` |
| Represented field | `A_φ` directly | **`ψ = 2π r A_φ`** (flux function); `A_φ` recovered via `/ (2π r)` |
| Axis (`r = 0`) behaviour | `A_φ(0) ≠ 0` in general — must Dirichlet explicitly. The `1/r` weight in the axisymmetric integrand spoils Gauss-quadrature convergence near the axis. | `ψ` carries an `s = r²` factor in the axis-restricted basis, so `A_φ = ψ/(2π r) → 0` automatically as `r → 0`. No `1/r`-weighted Gauss quadrature; integrals are closed form. |
| Field formula | `B_z = (1/r) ∂(r A_φ)/∂r` (contains `1/r`) | **`B_z = (1/π) ∂ψ/∂s`** — pure derivative, no singularity |
| Element-matrix integration | NGSolve Gauss quadrature with `1/r` weight | **Closed-form** (Mathematica-derived per-monomial integrals, machine precision) |
| Effective convergence rate | Often `p_eff = 1` near the axis owing to weight singularity | Full `p = 2` in `(s, z)` |

Empirical comparison on the Cu disk eddy-current benchmark
(R = 10 mm, t = 2 mm, σ = 5.8 × 10⁷ S/m, BEM-Foster reference τ₁ = 224.31 µs):

| FE                                  | DOFs   | τ₁ (µs)   | Gap to BEM |
|-------------------------------------|--------|-----------|------------|
| NGSolve `H1 order=3` (v23 reference)| —      | 218.62    | 2.5 %      |
| `axihenrotte p=1` (Q1, fine mesh)   | 14904  | 223.06    | 0.55 %     |
| **`axihenrotte p=2` (Q2)**          | **9919** | **223.69** | **0.27 %** |

`axihenrotte p=2` outperforms NGSolve `H1 order=3` with fewer DOFs; this is
the value of the axisymmetric-specific basis.

### 3-way Cauer-I cross-validation

Beyond the leading τ₁ comparison, the per-stage Cauer ladder time constants
`τ_rung[n] = L_n × λ_{2n-1}` (Hiruma's "stage time constant") have been
cross-checked against an **independent BEM-Foster-to-Cauer pipeline**
(Mathematica `bem_disk_axisym_cauer.wls` → 50-digit mpmath Cauer-I CFE on
the 20 leading α_n moments):

| n | BEM Cauer (µs) | Q2 fine (µs) | Q1 very-fine (µs) | Q2/BEM gap | Q1/BEM gap |
|---|---|---|---|---|---|
| 1 | 219.32 | 218.71 | 218.05 | **-0.28 %** | -0.58 % |
| 2 |  78.65 |  78.12 |  77.77 | **-0.68 %** | -1.12 % |
| 3 |  40.04 |  39.54 |  39.37 |  -1.24 %    | -1.66 % |
| 4 |  23.74 |  23.16 |  23.14 |  -2.46 %    | -2.54 % |
| 5 |  17.07 |  16.07 |  16.06 |  -5.86 %    | -5.91 % |
| 6 |  14.70 |  13.12 |  13.01 | -10.77 %    | -11.50 % |

Three **independent** numerical pipelines —
(i) Mathematica BEM with elliptic integrals + Foster-amplitude moments +
50-digit mpmath Cauer CFE,
(ii) Python Q1 axifemm prototype + Hiruma 3-term recurrence,
(iii) C++ `axihenrotte p=2` extension + Hiruma 3-term recurrence —
agree on the leading three Cauer rungs to ~ 1 % and on the leading rung to
0.28 %. `axihenrotte p=2` is closer to BEM than `axihenrotte p=1` at every
stage. This is **Phase 3-(3) cross-validation**, executed 2026-05-06; the
test lives at
[`tests/test_3way_cauer_cross_validation.py`](../../packages/radia-axifemm/tests/test_3way_cauer_cross_validation.py).

## Basis details

### `p=1` Q-element (4 DOFs)

Shape monomials in `(s, z)`:

```
{1, s, z, s·z}                    (interior)
{s, s·z}                          (axis-touching, sa = 0; the 2 axis-side
                                   nodes get zero shape functions)
```

DOFs: one per vertex. Lagrange-interpolatory at the four corners
`(s_a, z_a), (s_b, z_a), (s_b, z_b), (s_a, z_b)`.

### `p=2` Q-element (9 DOFs)

Shape monomials in `(s, z)`:

```
interior (9): {1, s, s², z, s·z, s²·z, z², s·z², s²·z²}
axis     (6): {s, s², s·z, s²·z, s·z², s²·z²}
```

DOFs: 4 vertex + 4 edge midnode + 1 face center = 9 per quad.

**Edge-midnode positions use the *s*-midpoint convention:** the bottom-edge
midnode of `[r_a, r_b] × [z_a, z_b]` lives at `s_m = (s_a + s_b)/2`, i.e.
physical `r = √((r_a² + r_b²)/2)` (quadratic mean of `r_a` and `r_b`), **not**
`(r_a + r_b)/2`. Two adjacent elements with the same `r`-extent automatically
agree on edge-midnode placement and basis traces; no orientation flips needed.

Internally the FE class stores a 9×9 inverse Vandermonde (or 6×6 padded for
axis elements) and the `BilinearFormIntegrator`s build per-element stiffness
and σ-mass matrices via

```
K_V = T · V⁻ᵀ · K_phi · V⁻¹ · T,    T = diag(2π r_node)
```

with `K_phi` provided in *closed form* by the Mathematica derivation
(`axifemm/derive_quad_q2_henrotte.wls` → `q2_henrotte_generated.hpp`).

### Convention: Hessian-of-W

Both `K_phi` and `M_sigma_phi` are emitted in the **Hessian-of-W** convention
(matching the validated `axifemm_quad.py` Q-element-`p=1` reference):

```
K_phi[i,j]       = ∫ ∂m_i/∂s · ∂m_j/∂s / (π μ_Z) ds dz
                 + ∫ ∂m_i/∂z · ∂m_j/∂z / (4 π μ_R s) ds dz
M_sigma_phi[i,j] = σ / (4 π) · ∫ m_i · m_j / s ds dz
```

The original `derive_quad_q2_henrotte.wls` shipped with a coefficient-of-W
convention that was 2× too small for `K` and 2π× too small for `M`; this was
diagnosed and corrected during Phase A2 (commit
[077e7b03](../../packages/radia-axifemm/)). Anyone re-deriving the wls
script must keep the Hessian convention.

## Cross-validation references (per-element, machine precision)

* `axifemm/axifemm_quad.py` — validated Python prototype for `p=1`, gives
  τ₁ = 223.06 µs on the Cu disk (matches BEM v3 to 0.55 %).
* `axifemm/axifemm_quad_q2.py` — Python `p=2` Gauss-8×8 prototype; agrees
  with the Mathematica closed form to ~ 3.4 × 10⁻⁸ relative.
* `packages/radia-axifemm/scripts/validate_q2_codegen.py` — runs both at the
  per-entry level after every `derive_quad_q2_henrotte.wls` re-run.

## File layout

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

## Why this is "an NGSolve feature NGSolve does not have"

NGSolve provides a vast standard FE library (`H1`, `HCurl`, `HDiv`, `L2`,
`VectorH1`, `Compound`, `HCurlDiv`, …) on triangles, quads, tets, hexes,
prisms, and pyramids of arbitrary order. None of these provides a basis
that is polynomial in `s = r²` rather than `r`, nor one that represents the
flux function `ψ = 2π r A_φ` rather than the field itself.

`axihenrotte` is what you need when:
* the problem is axisymmetric eddy-current or magnetostatic;
* axis touching elements occur in the mesh;
* you want closed-form per-element integration of the
  `1/r`-weighted axisymmetric integrand;
* you want the convergence rate that the FEMM / Henrotte tradition gets
  on this exact class of problem (≈ 2× more accurate than `H1 order=3`
  per DOF on the Cu-disk benchmark).

Within the NGSolve ecosystem it is a **drop-in FESpace** — once
`radia_axifemm` is imported, `FESpace("axihenrotte", mesh, order=k)`
works exactly like any other space.

[^Henrotte93]: F. Henrotte et al., "A new method for axisymmetric linear
    and nonlinear problems," *IEEE Transactions on Magnetics* 9(2):1352–1355,
    March 1993.
[^MeekerFEMM]: D. Meeker, FEMM 4.2 axisymmetric formulation notes
    (`prob3big.cpp` `StaticAxisymmetric`).
