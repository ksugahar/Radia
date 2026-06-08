# FEEC VIM — the HDiv-type MMM/MSC (replacing the yano-type distortion elements)

Research demonstration of the **FEEC (Finite Element Exterior Calculus) Volume Integral
Method** for the MMM/MSC magnetization problem. The goal: build an **HDiv-type** demag
operator whose magnetization "loops" are **field-null by construction**, and — if it matches
or beats the **yano-type** element-engineering (Yano's hand-crafted elements that suppress
the loop-star / `A_ls` component on distorted hexes, preserved in the private ELF repo) —
retire the yano-type from public Radia.

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

## HDiv-type vs yano-type — the clean win

**No fundamental inferiority.** The HDiv-type N is **symmetric**, works on **general**
elements (tet/hex/wedge), is **distorted-robust by construction**, and its accuracy is
comparable to yano-type. On top of parity, it has three accuracy-per-DOF advantages the
flat, hand-crafted yano-type **structurally cannot match**:

1. **de-Rham-exact on ANY mesh** — loops are field-null by construction (`4e-16`) on
   distorted hexes, vs the MSC retrofit's `6e-9`; no per-mesh element engineering.
2. **curved (isoparametric) geometry** — `mesh.Curve(p)`; the external field of a coarse
   sphere goes from flat `−10%` to `<0.3%` at the SAME ndof (`hdiv_demag_curved.py`). Flat
   `ObjHexahedron/ObjTetrahedron` cannot represent a curved boundary.
3. **polynomial high-order** — the FEEC construction is order-agnostic (loops stay
   field-null at all orders); the p-convergence speed win is the major continuation.

## What is built + validated — with REFERENCE HONESTY

What each accuracy number is measured against, precisely (this matters — Radia is a *trusted*
solver but NOT ground truth on a coarse mesh):

| Result | Script | Validated against | Status |
|---|---|---|---|
| Loops field-null on distorted hex | `ngsolve_loopfree_verify.py` | exact (charge-form field `3.7e-16`) | ✅ machine zero |
| Loop/star (Hodge) split | `hdiv_loop_star_split.py` | exact (`ker Q` charge-free `~1e-16`) | ✅ |
| Linear demag (sphere/cube → 1/3) | `hdiv_demag_tet.py` | **ANALYTIC** 1/3 | ✅ `<0.15%` (Wilton surface Gram) |
| Nonlinear (damped Newton) | `hdiv_demag_tet_nonlinear.py` | **ANALYTIC** sphere fixed point | ✅ `<0.05%` deep-saturation |
| Nonlinear cross-check | `test_hdiv_vim_newton_vs_radia.py` | Radia MMM/MSC (`MatSatIsoTab`) | ✅ agree `<0.05%` (sphere) |
| Real BH table | `test_hdiv_vim_newton_table.py` | **ANALYTIC** uniform-sphere | ✅ `<0.2%` |
| Ellipsoid (D≠1/3) | `test_hdiv_vim_ellipsoid.py` | **ANALYTIC** prolate `N_z` | ✅ 2:1 `0.3%` |
| Volume Gram (`phi_tet`) | `test_hdiv_vim_volume_gram.py` | Radia (cube nonlinear) | ⚠ 13%→6.2% **agreement** (no analytic truth; cross-method difference, not a verified error) |
| Scalable (C++ H-matrix + GMRES) | `test_hdiv_vim_newton_scalable.py` | dense reference | ✅ machine precision |
| Distorted μr-independence | `test_hdiv_vim_solve.py` | iters bounded vs μr 10→1e4 | ✅ golden-locked |
| **Curved-mesh win** | `hdiv_demag_curved.py` | **ANALYTIC** dipole / volume | ✅ external field flat `−10%` → Curve(3) `−0.26%` (~38× at same ndof) |

The **demag FACTOR does NOT discriminate the curved win** (a coarse inscribed polyhedron is
already near-isotropic → `D_z ≈ 1/3` regardless of faceting); the win lives in
geometry-derived field quantities (external field/flux/force) — see `hdiv_demag_curved.py`'s
docstring and `test_hdiv_vim_curved.py`'s `test_demag_factor_does_not_discriminate`.

## Golden tests

`tests/feec/test_hdiv_vim_*.py` (full feec suite **65 passing**) — Newton, Newton-vs-Radia,
Wilton Gram, volume Gram, scalable, ellipsoid, BH table, distorted robustness, curved win.

## Detailed home

The narrative + decisions live in the radia-mcp **`hdiv_vim`** MCP knowledge
(`overview` / `status` / `nonlinear`); `memory/` holds the Problem-A/B investigation record.

## Honest open items (the productionization to actually retire yano-type)

1. **Curved demag OPERATOR** — `build_demag`/`wilton_surface_block`/`phi_tet` still use FLAT
   vertices; the accurate curved Gram reuses `src/radia/bem/sibc_hacapk.py::_ss_block_curved_trafo`
   (Galerkin single-layer, proper singular self/near on curved elements). An **integration,
   not research** — the curved-sampling primitive (`_trafo_sample`) and the geometry-win
   proof (vs analytic truth) are already done in `hdiv_demag_curved.py`.
2. **High-order Gram** — Graglia polynomial-density weakly-singular quadrature; the
   p-convergence accuracy-per-DOF win over lowest-order yano-type.
3. **C++ maturity** — Wilton + `phi_tet` + curved in the C++ HACApK charge Gram, the full
   Newton loop in C++, and a Radia API. This is the big lift that turns the validated
   prototype into the shipped replacement.
