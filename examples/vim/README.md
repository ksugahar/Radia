# FEEC VIM — the HDiv-type MMM/MSC

Research demonstration of the **FEEC (Finite Element Exterior Calculus) Volume Integral
Method** for the MMM/MSC magnetization problem. The goal: build an **HDiv-type** demag
operator whose magnetization "loops" are **field-null by construction**, complementing the
canonical **multipole-moment MMM** surface-charge backend with FEEC, curved/high-order geometry, and
symmetry-model capabilities.

> **Sibling research line in this directory:** a separate set of examples builds a **3D
> stream-function + cohomology COIL-DESIGN method** (Clebsch / vector-T potentials on the same
> FEEC de Rham complex, ACA+TSVD inverse, wire extraction, and an honest map of the
> non-convex / helicity frontier) — see
> [`README_streamfunction_cohomology.md`](README_streamfunction_cohomology.md). Same FEEC
> foundation, different application (coil design vs the demag operator above).

> **The question:** using NGSolve's H(div) basis inside the volume integral, are the
> magnetization loops (div-free circulations, `ker N`) FIELD-NULL — so no spurious loop
> modes arise, even on distorted elements at high μr?
>
> **The answer (this directory): YES, by construction** — and the resulting operator is
> symmetric, works on general/distorted/curved elements, and is high-order-capable.

## The two problems must be separated (see `memory/`)

| Problem | What it is | Does the FEEC basis fix it? |
|---|---|---|
| **A — loop-mode / de Rham** | On distorted hexes the constant-M ±1 loop is *not* field-null → carries spurious field → the 28%-wrong field at μr=1e5 ("the defect") | **YES — resolved by construction** (this directory) |
| **B — high-μr conditioning** | BiCGSTAB iters ~μr^1.5; the >15³ H-ILU factor wall. `A=(1/χ)I−N → −N` as μr→∞ | **NO** — near-null *spectrum*, formulation-independent (hex≈tet); H-ILU handles μr≤1e4, scaling is open |

## Why A is resolved (the criterion)

A magnetization `M` is **field-null** (a loop, in `ker N`) **iff it is charge-free**:
`ρ = −div M = 0` **and** `σ = M·n = 0` (both — div-free alone is not enough; a curl can
still carry surface charge). The FEEC loop `M = curl(interior H(curl))` (interior = zero
tangential trace) has **both**:

- `div(curl W) = 0` → no volume charge;
- `(curl W)·n = curl_S(W_tangential) = 0` → no surface charge (W interior ⇒ `W_t = 0`).

⇒ **charge-free ⇒ field-null everywhere, at any μr, on any element** — the contravariant
Piola map preserves *both* div and normal-trace (the de Rham commuting diagram), so
distortion cannot break it. This is the strong (everywhere) field-null vs the constant-M
basis's fragile (collocation-only) field-null that breaks under distortion. The operator
`N = BᵀGB` (B = charge map, G = Coulomb Gram) inherits `N·loop = 0` for **any** G.

## HDiv-type vs multipole-moment MMM — the clean win

**No fundamental inferiority.** The HDiv-type N is **symmetric**, works on **general**
elements (tet/hex/wedge), is **distorted-robust by construction**, and its accuracy is
comparable to multipole-moment MMM MSC. On top of parity, it has three accuracy-per-DOF advantages the
flat, hand-crafted multipole-moment MMM MSC **structurally cannot match**:

1. **de-Rham-exact on ANY mesh** — loops are field-null by construction (`4e-16`) on
   distorted hexes, vs the MSC retrofit's `6e-9`; no per-mesh element engineering. This also
   makes **symmetry models (1/2, 1/4, 1/8) automatic**: on a cut/reduced mesh the loops are
   just `ker(B)` (field-null `~4e-16`, count adapts 58→54→18→6 for sphere full→1/8), with **no
   cohomology-aware loop-star `installCycle`** — the "loop-removal is painful" problem of
   multipole-moment MMM MSC is eliminated (`test_hdiv_vim_symmetry_loops.py`). The symmetry *demag value*
   is the **image method** (`hdiv_demag_symmetry_image.py`): reflecting the reduced model's cap
   charge over the reduction planes (sign per z-mirror) reproduces the full demag from ~1/2, 1/4,
   1/8 the DOF (1/2 +0.08%, 1/4 +0.11%, 1/8 −0.32% vs full). So 1/4 & 1/8 are supported — loops
   automatic + demag via images, no hand-crafted loop-star.
2. **curved (isoparametric) geometry** — `mesh.Curve(p)`; the external field of a coarse
   sphere goes from flat `−10%` to `<0.3%` at the SAME ndof (`hdiv_demag_curved.py`). Flat
   `ObjHexahedron/ObjTetrahedron` cannot represent a curved boundary.
3. **polynomial high-order + curved → EXACT** — the surface demag Gram is the Laplace
   single-layer of `σ=M·n`, supplied high-order + curved + FMM by `ngsolve.bem`. On a coarse
   sphere the demag factor goes flat `+0.25%` (order-insensitive, faceting-floored) →
   **curved + order-2 `~1e-4%` (exact)** at fixed small ndof (`hdiv_demag_bem_singlelayer.py`).
   This is the accuracy-per-DOF win over flat lowest-order multipole-moment MMM MSC, on the demag factor
   directly — and it reuses NGSolve, no hand-rolled singular quadrature.

## What is built + validated — with REFERENCE HONESTY

What each accuracy number is measured against, precisely (this matters — Radia is a *trusted*
solver but NOT ground truth on a coarse mesh):

| Result | Script | Validated against | Status |
|---|---|---|---|
| Loops field-null on distorted hex | `ngsolve_loopfree_verify.py` | exact (charge-form field `3.7e-16`) | ✅ machine zero |
| Loop/star (Hodge) split | `hdiv_loop_star_split.py` | exact (`ker Q` charge-free `~1e-16`) | ✅ |
| Linear demag (sphere/cube → 1/3) | `hdiv_demag_solve` (`test_hdiv_vim_demag_solve.py`) | **ANALYTIC** 1/3 | ✅ `<0.5%` (C++ analytic charge Gram) |
| Nonlinear (damped Newton) | `solve_nonlinear_newton` (`test_hdiv_vim_tet_newton.py`) | **ANALYTIC** sphere fixed point | ✅ `<2%` deep-saturation (C++ charge Gram) |
| Nonlinear cross-check | `test_hdiv_vim_newton_vs_radia.py` | Radia MMM/MSC (`MatSatIsoTab`) | ✅ agree `<0.05%` (sphere) |
| Real BH table | `test_hdiv_vim_newton_table.py` | **ANALYTIC** uniform-sphere | ✅ `<1%` |
| Ellipsoid (D≠1/3) | `test_hdiv_vim_ellipsoid.py` | **ANALYTIC** prolate `N_z` | ✅ 2:1 `0.3%` |
| Volume Gram (analytic, div M≠0) | `test_hdiv_vim_volume_gram.py` | **ANALYTIC** (linear demag → 1/3) | ✅ the C++ analytic charge Gram (exact near AND far); required for non-uniform nonlinear |
| **Non-uniform nonlinear vs Radia** (cube + C-yoke) | `hdiv_cyoke_nonlinear.py`, `test_hdiv_vim_cyoke_nonlinear.py` | shipped **Radia** (both flat → valid) | ✅ volume-avg M_z agrees **<3%**, converges, mesh-stable. The C++ analytic charge Gram supplies the volume Gram required for div M≠0. |
| Scalable (C++ H-matrix + GMRES) | `test_hdiv_vim_newton_scalable.py` | **ANALYTIC** sphere fixed point + physical range | ✅ deep-saturation match, bounded iters |
| Distorted μr-independence | `test_hdiv_vim_solve.py` | iters bounded vs μr 10→1e4 | ✅ golden-locked |
| **Curved-mesh win** (elementary) | `hdiv_demag_curved.py` | **ANALYTIC** dipole / volume | ✅ external field flat `−10%` → Curve(3) `−0.26%` (~38× at same ndof) |
| **Curved + high-order demag** (production) | `hdiv_demag_bem_singlelayer.py` | **ANALYTIC** sphere 1/3 + spheroid + triaxial tensor | ✅ flat floored → curved + order-2 EXACT: sphere `~1e-4%`; prolate & oblate polar+transverse `<0.05%`; **triaxial ellipsoid** (a≠b≠c) all three distinct factors EXACT vs Osborn integral; sum rule `N_x+N_y+N_z=1` to `~1e-6`; Gram = `ngsolve.bem` single-layer |
| **Cubit curved-hex HDiv p-convergence** | `hdiv_cubit_hex_sphere_pconv.py` | **ANALYTIC** permeable sphere in dipole (`l=1`) and quadrupole (`l=2`) fields | ✅ Cubit `volume scheme sphere` 56-hex mesh, curved .vol order 1→3: dipole external-field error `23.5% → 0.21% → 0.13%`, quadrupole `35.6% → 0.59% → 0.19%`; on fixed order-3 geometry, HDiv p=0→3 gives quadrupole `32.1% → 0.63% → 0.36% → 0.19%`. This is projection + polynomial-charge field reconstruction, not the self-consistent high-order hex VIM solve. |
| **Cubit curved-hex harmonic HDiv-VIM solve** | `hdiv_cubit_hex_sphere_harmonic_vim.py` | **ANALYTIC** sphere harmonic demag coefficients `D_l=l/(2l+1)` | ✅ Self-consistent harmonic-subspace solve on the same Cubit 56-hex sphere. Order 3 gives dipole `D=0.333408` with external-field error `0.112%`, and quadrupole `D=0.399267` with external-field error `0.051%`. This validates the physically curl-free harmonic subspace; a general full-space high-order HDiv solve still needs a curl-free or mixed material constraint before using an unconstrained energy solve. |
| **Curved × nonlinear** — magnetization (honest) | `test_hdiv_vim_curved_nonlinear.py` | **ANALYTIC** spheroid M-H fixed point | ✅ curved nonlinear M exact (`<0.05%`); ⚠ but the curved win on the **magnetization** is MODEST (`~0.3%`) — the demag *ratio* cancels the volume faceting error. (Radia can't referee curved geometry — it facets; so curved nonlinear is validatable only on spheroids.) |
| **Curved × nonlinear** — field (the big win) | `hdiv_curved_nonlinear_field.py` | **ANALYTIC** dipole | ✅ external H field of a nonlinear soft-iron sphere: **flat `~+8.8%` at every point → Curve(3) `<0.4%` (~23×)**. The field inherits the ~9% volume error (dipole moment m=M·V); THIS is where curved × nonlinear pays off — the engineering deliverable (stray field), not M. |
| **Head-to-head vs shipped Radia** (B) | `compare_curved_vs_radia_field.py` | **ANALYTIC** dipole + shipped **Radia** | ✅ HDiv curved at the *coarsest* mesh (0.39%) beats shipped‑Radia‑**flat** at the *finest* (1.71%, 2042 tets); **~10–30× accuracy‑per‑resolution** at every h. Honest: accuracy‑per‑DOF (geometry‑driven); wall‑clock = the C++ lift (not done); Radia‑flat stands in for the also‑flat yano‑type. |

Multipole-moment MMM remains the large flat-cell moment baseline. It does not consume Cubit high-order curved
hex `.vol` nodes or NGSolve HDiv/Piola shape functions, so the curved-hex p-convergence and
harmonic HDiv-VIM checks above are HDiv-only evidence.

**Which quantity discriminates the curved win, and why it matters:** with the *crude
sub-point* Gram (`hdiv_demag_curved.py`) the demag FACTOR does NOT cleanly discriminate — its
~2% quadrature *bias* masks the ~0.25% geometry signal — so for that elementary method use the
EXTERNAL FIELD (geometry-only error). This is a limitation of the crude Gram, **not** of the
demag factor: with the **proper** Gram (the `ngsolve.bem` Laplace single-layer) the demag
factor discriminates cleanly **and p-converges** — flat floors at `+0.25%`, curved + order-2
is exact (`hdiv_demag_bem_singlelayer.py`). (An earlier note here wrongly attributed the
non-discrimination to "near-isotropy"; the verify-first single-layer result corrected it.)

## Golden tests

`tests/feec/` (full feec suite **85 passing**) — Newton, Newton-vs-Radia, Wilton Gram, volume Gram,
scalable, ellipsoid, BH table, distorted robustness, curved win (elementary), curved + high-order demag
exact on sphere + full spheroid tensor + general triaxial ellipsoid (sum rule) via the `ngsolve.bem`
single-layer, curved × nonlinear (honest modest magnetization win + the ~23× field win), and the
head-to-head accuracy-per-resolution win vs the shipped Radia solver (`test_curved_vs_radia_field.py`).

## Benchmark result management

Benchmark JSON files are part of the example evidence and are managed in this directory, next to the
script that generated them. Do not leave the only copy under a remote scratch directory such as
`C:/temp`: run large jobs locally on mdx for speed and memory isolation, then copy the completed JSON
back into `examples/vim/` before publishing.

Recommended pattern for mdx:

```powershell
# On mdx, run from a local scratch directory.
python C:/temp/radia_mdx_vim/hdiv_matvec_scaling.py

# Copy the resulting JSON back to this directory on the LAB checkout.
# The canonical tracked result is examples/vim/hdiv_matvec_scaling.json.
```

For multipole-moment MMM large-scale evidence, the dense nonlinear baseline
`bench_multipole_moment_scaling.py` intentionally stops at small/medium DOF. The 165k-DOF path is the
HACApK method-2 storage/solve benchmark:

```powershell
python bench_moment_storage_scaling.py --cases 60x46x10 --methods 2 --timeout 7200
```

`60x46x10` gives `60*46*10*6 = 165600` multipole-moment MMM DOF. Dense LU (`method 0`) is intentionally not
run at that size because its dense matrix alone would exceed the mdx memory budget.
The mdx result published on 2026-06-24 is
`results_moment_storage_scaling_165600_mdx_20260624.json`: 165600 DOF, 13389.5 MB peak memory,
830.75 s solve time, 1 BiCGSTAB iteration, and converged.

The 3-way solver comparison (method 0 dense LU vs method 1 dense-K BiCGSTAB vs method 2
HACApK-BiCGSTAB, the collocation-MMMM coarse tier of 2026-07-02) is
`bench_moment_solvers.py`: per-case subprocess isolation (Benchmark Policy memory accuracy),
compact-cube AND elongated-bar geometries (near-field vs ACA-compressible), cold + warm solve
timing (the warm solve exposes the method-2 cross-solve K cache = the optimization-inner-loop
per-iteration cost), and external-B cross-method correctness columns.

```powershell
python bench_moment_solvers.py --sweep        # writes results_moment_solvers.json
```

## Detailed home

The narrative + decisions live in the radia-mcp **`hdiv_vim`** MCP knowledge
(`overview` / `status` / `nonlinear`); `memory/` holds the Problem-A/B investigation record.

## Honest open items (productionization alongside multipole-moment MMM)

The **surface** Gram, curved + high-order + FMM-scalable, is now SOLVED by reusing the
`ngsolve.bem` Laplace single-layer (`hdiv_demag_bem_singlelayer.py`) — no hand-rolled singular
quadrature needed. The remaining work:

1. **Full operator on the single-layer** — wire `N = BᵀVB` (`B` = HDiv charge map, `V` =
   `ngsolve.bem` single-layer) and run a self-consistent linear solve on a curved/high-order
   body, validated vs analytic. The demag-factor proof (`<σ,Vσ>/V_vol` → 1/3) is done; the
   operator + solve is the next assembly.
2. **Volume charge for non-uniform M** — `ngsolve.bem` is boundary-only, so the nonlinear case
   (`div M ≠ 0` → volume charge `ρ`) still needs the Newtonian volume potential (`phi_tet`,
   already built) on curved/high-order cells; combine with the single-layer for the full
   curved+high-order nonlinear operator.
3. **C++ maturity** — the production charge Gram (single-layer surface + `phi_tet` volume) and
   the Newton loop in C++ behind a Radia API. The big lift that turns the validated prototype
   into a shipped production backend.
