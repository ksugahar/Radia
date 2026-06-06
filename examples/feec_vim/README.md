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

## Plan: build the HDiv-type MMM/MSC; if it works, retire the yano-type (2026-06-06)

The **HDiv-type MMM/MSC** uses NGSolve's H(div) basis for the magnetization; its loops are
field-null **by construction** (de Rham), replacing the **yano-type** element-engineering
(Yano's hand-crafted elements that suppress the A_ls / loop-star component on distorted
hexes). Plan: build the HDiv-type; **if it works**, retire the yano-type from public Radia
and **preserve it in the ELF repo**. (Lowest-order HDiv-type == the shipped Radia MSC +
`installCycle` de-Rham loops, which already works; the new build is the **high-order** VIM.)

Progress:

1. ✅ **loop/star (Hodge) split on the real HDiv space** — `hdiv_loop_star_split.py`:
   the charge map `Q : M ↦ (−div M, M·n|_bnd)` splits HDiv into **star(charge-carrying) ⊕
   loop(charge-free)** (distorted 2×2×2, HDiv order 1: 240 = 159 star + 81 loop; `rank Q =
   160 charge-dofs − 1` Gauss charge-conservation). Loops from `ker Q` are charge-free
   (`div`, `M·n` ~1e-16) and **field-null** (`|H|`~3e-17) — the **operator-level**
   loop-mode-free statement: loops are in `ker(N_demag)` exactly, no element engineering.
   The 159-dim **star** space is where the VIM is solved. (The loop space is defined by the
   charge map directly; a naive `curl(HCurl_p)` count does not apply — that lands in HDiv
   order *p−1*, a different space.)
2. **next — assemble the demag operator on the star space + a small validated solve**:
   the genuine remaining build, whose hard part is the **weakly-singular 1/r charge-Galerkin
   self/near term** (Duffy / Graglia surface + volume moments). Validate against a known
   demag (e.g. sphere factor 1/3) and against Radia MSC at lowest order.
3. **RHS / collocation correspondence** — evaluation points matched to the high-order basis
   (the "矢野 element" lesson).

Problem B (high-μr conditioning / scaling) is orthogonal and not addressed by the basis;
H-ILU handles μr≤1e4 (shipped).
