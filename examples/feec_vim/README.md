# FEEC VIM — does NGSolve's basis solve the loop-mode problem?

Research demonstration of the **FEEC (Finite Element Exterior Calculus) Volume Integral
Method** for the MMM/MSC magnetization problem, focused on **one question**:

> Using NGSolve's H(curl)/H(div) basis inside the volume integral, are the magnetization
> "loops" (div-free circulations, `ker N`) FIELD-NULL — so no spurious loop modes arise,
> even on distorted elements at high μr?

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
basis's fragile (collocation-only) field-null that breaks under distortion.

## Files

- `ngsolve_loopfree_verify.py` — builds a genuinely distorted hex mesh in **real NGSolve
  6.2.2604**, forms `M = curl(interior H(curl))`, and confirms charge-free → field-null:
  - `||div M|| / ||M|| = 1.05e-13`
  - `||M·n||_bnd / ||M|| = 0` (exact)
  - charge-form demag field at external points `= 3.7e-16` (machine zero) → **field-null**
  - (contrast) uniform `M=(0,0,1)`: charge `σ=0.58`, field `1.8e-2` — a "star", not a loop
  - `ngsolve_loopfree_verify.json` — the recorded numbers
- The symbolic twin: `radia_mcp/mathematica/basis_functions/vim_loopfree.wls` (+ the full
  FEEC suite: `h1`, `hcurl`, `hdiv`, `simplex_ho`, `derham`, `cohomology`, `vim_field`).

## Status / honest next steps (before any C++)

Problem A is **settled at the formulation level** (symbolic + real NGSolve). Problem A is
*also* already handled in production Radia by the local-null-vector `installCycle` fix — so
the FEEC VIM's marginal value is **high-order accuracy + a principled (not patched) loop
space**, not fixing an open production defect.

Remaining formulation work (all in Python/NGSolve, before C++):

1. **loop/star Hodge split + star-block conditioning** — confirm the charge-free loop
   space is exactly `curl(interior H(curl)) ⊕ cohomology H¹` (dimension count) and that the
   star (charge-carrying) block is well-posed once loops are removed.
2. **VIM RHS / collocation correspondence** — the evaluation points must match the
   high-order basis (the "矢野 element" lesson).
3. **a small end-to-end VIM solve** on a distorted multi-element patch vs a trusted
   reference — the genuine remaining build, whose hard part is **singular near-field
   quadrature** of the 1/r kernel for the high-order basis.

Problem B (high-μr conditioning / scaling) is orthogonal and not addressed by the basis.
