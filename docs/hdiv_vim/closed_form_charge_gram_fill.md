# Closed-form charge-Gram fill for the HDiv-VIM hex/wedge near block

**Design + decision record (2026-07-05).** The HDiv-VIM RT1 hex/wedge charge Gram
`G_ab = INT_A INT_B c_a(x) c_b(y)/|x-y|` (the demag operator `N = B^T G B`; `c` = div of the RT1 basis, a Q1
trilinear charge) is currently filled per-entry by a 6-sub-tet Duffy-graded quadrature (`QuadBlockHex`). This
note records the decision to replace that per-entry fill with a **semi-closed form** — a closed box Newtonian
potential (inner) plus a **symmetric** tensor Gauss outer — which (1) removes the ~1e-5 reflection-symmetry
defect exactly and (2) lowers the per-entry cost. The closed forms are derived and verified here; the C++ port
into `src/core/rad_hacapk_hdiv.cpp` is gated (see *Implementation gate*).

## Decision

- Adopt the **semi-closed charge-Gram entry** for **disjoint** hex/wedge pairs (near + far blocks). The
  coincident/self block keeps the existing exact `QuadBlockHex` path (it is already exact to 3.66e-15).
- This is direction **"A" (closed-form per-entry fill)**, chosen over "B" (structured-hex BTTB/FFT, rejected as
  too narrow — it only helps regular lattices).
- The reflection-symmetry defect's strict xfail
  (`validation_test/feec/test_hdiv_radfld_contract.py`) stays as the documented known-limitation until the port
  lands; the 1e-5 field-parity error is accepted in the interim (two orders below the ~1e-3 engineering
  tolerance; the demag factor and the parity-averaged `M_avg` are already correct).

## Why: LAB build profile (the bottleneck)

Coarse LAB profile (`build_time` from `H.stats`, structured hex block):

| n | charge DOF | C++ build_time | ms/DOF | blocks (dense/lowrank) |
|---|-----------:|---------------:|-------:|------------------------|
| 6 | 2592 | 17.6 s | 6.8 | 762 / 370 |
| 8 | 5632 | 40.2 s | 7.1 | 1600 / 1152 |
| 10 | 10400 | 109.8 s | 10.6 | 3435 / 2242 |

- **The Gram build is 90 % of the full HDiv-VIM hex solve** (n=8: 59.2 s build of 66.1 s; the mass-Riesz CG
  iterate is ~10 %). Build is the lever.
- **The hex/wedge Gram has no near/far quadrature split.** `far_quad` / `ho_far_factor` are TET-only;
  hex/wedge route to `_build_charge_gram_hex` -> `QuadBlockHex`, which applies the full 6-sub-tet
  (`glout_n=6`, `glin_n=5`) Duffy fill to **every** entry, dense and ACA-sampled alike. So the closed form
  replaces the *entire* per-entry cost, and it can *add* the near/far adaptivity the hex path lacks today.

## Why it works: the semi-closed Gram (recipe)

For box A (x-charge `c_a(x) = prod_{k in Sa} x_k`) and box B (y-charge `c_b(y) = prod_{k in Sb} y_k`), with
`Sa, Sb` subsets of `{1,2,3}` spanning the Q1 = 8-dim trilinear charge space `{1,x,y,z,xy,xz,yz,xyz}`:

```
G_ab   = sum over SYMMETRIC tensor-Gauss nodes x of A:  w_x * c_a(x) * Phi_b(x)
Phi_b(x) = INT_B c_b(y)/|x-y| dy
         = sum_{T subset Sb} (-1)^|T| * (prod_{k in Sb\T} x_k) * P_T(x, B)
P_T(x,B) = - sum over the 8 corners of B: (+/- sign) * F_T( x - corner )     # inclusion-exclusion, u = x - y
```

`F_T` is the 3-fold antiderivative with `d^3 F_T / du dv dw = (prod_{k in T} u_k)/|u|` (listed below). The inner
box potential is thus **closed** (no singular quadrature). The **symmetric** outer rule (nodes symmetric about
A's centre) commutes with the box reflections, so a symmetry-forbidden (transverse) entry is **machine-zero by
construction** — this is the exact reflection-symmetry fix, with no 6-sub-tet split, no Duffy, no 16x.

## Per-entry cost estimate (LAB)

The flat hex `QuadBlockHex` already uses a **closed** inner (`PhiAtHO_Analytic`, the analytic moment) times a
Gauss-Duffy outer; the difference is therefore the **outer** rule: `QuadBlockHex` = 6 sub-tets x `glout_n=6`^3
≈ **1296 Duffy-graded outer points/host** (same for near and far). The semi-closed symmetric tensor outer
converges as (relative error vs a high-order reference):

| pair | 1e-6 reached at | outer points |
|------|-----------------|-------------:|
| TOUCH const-const (hardest near) | n=6 | 216 |
| TOUCH const-const (1e-7) | n=8 | 512 |
| TOUCH xyz-xyz | n=6 | 216 |
| SEPARATED const-const (far) | n=4 | 64 |

So the semi-closed uses **~4x fewer outer points on near/touching pairs (343 vs 1296)** and **~20x fewer on
far pairs (64 vs 1296)** — the far win because the hex path has no far reduction today. The inner cost per outer
point is comparable (both closed).

**Empirical wall-clock (LAB, fair numpy-vs-numpy at matched 1e-6, same closed inner in both).** Timing the
symmetric tensor outer (A) against a *lean* Kuhn 6-sub-tet outer (B, order tuned to 1e-6 — already leaner than
the C++'s fixed 1296-pt Duffy outer), per entry:

| pair / charge | A symmetric outer | B 6-sub-tet outer | A faster by |
|---------------|-------------------|-------------------|------------:|
| TOUCH const-const | n=7, 343 pts | q=5, 750 pts | 1.45x |
| TOUCH xyz-xyz | n=5, 125 pts | q=5, 750 pts | 2.59x |
| SEP const-const | n=4, 64 pts | q=4, 384 pts | 1.75x |
| SEP xyz-xyz | n=4, 64 pts | q=4, 384 pts | 1.89x |

The closed inner was cross-checked against `scipy.integrate.tplquad` (0 / 6e-16). A is **1.45–2.6x faster per
entry than even a lean sub-tet outer**; against the *actual* C++ QuadBlockHex (fixed 1296 outer pts, no far
reduction) A uses ~4x (near) to ~20x (far) fewer outer points, so the real-C++ advantage is larger than the
measured numpy ratio. Absolute seconds are numpy; the RATIO is language-neutral (same inner, same primitives).
The delivered C++ wall-clock factor still needs the C++ port measured on **mdx** (LAB timing is contended), but
the algorithm-level win is now empirically confirmed, not only counted.

## The verified closed antiderivatives (C99; `r = sqrt(u*u+v*v+w*w)`)

The 7 moment forms are verified `d^3 F_T/du dv dw == (prod_{k in T} u_k)/r` to 1e-16 AND their box potentials
corner-sum-match NIntegrate (the sufficient FTC condition + the actual corner-sum both hold). The const is the
standard prism potential — see *Self-verify lesson*.

```c
// r = sqrt(u*u + v*v + w*w)
// T = {}   const  -> the standard box (prism) Newtonian potential INT_B 1/|x-y| dy
//                    (MacMillan 1930 / Nagy 2000 / Waldvogel 1979). Use the existing analytic const box
//                    potential, or Integrate[1/Sqrt[u^2+v^2+w^2], w,u,v]. DO NOT hand-transcribe a compact
//                    asinh/atan2 one-liner without validating its corner-sum (see Self-verify lesson).
F_u   = (-21*u*u*w - 2*w*w*w + 12*v*w*r + 12*u*u*u*(atan(w/u) - atan((v*w)/(u*r)))
         + 9*v*(u*u+v*v)*atanh(w/r) + 6*w*(3*u*u+w*w)*log(v+r) - 3*v*(-3*u*u+v*v)*log(w+r))/36;
F_v   = (4*u*w*r - 4*v*v*v*atan((u*w)/(v*r)) + 2*w*(3*v*v+w*w)*atanh(u/r)
         + 3*u*(u*u+v*v)*atanh(w/r) - u*(u*u-3*v*v)*log(w+r))/12;
F_w   = (4*u*v*r - 4*w*w*w*atan((u*v)/(w*r)) + 3*v*(v*v+w*w)*atanh(u/r)
         + 2*u*(u*u+3*w*w)*atanh(v/r) - v*(v*v-3*w*w)*log(u+r))/12;
F_uv  = (w*r*(5*(u*u+v*v) + 2*w*w) + 3*(u*u+v*v)*(u*u+v*v)*atanh(w/r))/24;
F_uw  = (-3*w*w*w*w + 8*v*v*v*r + 20*v*w*w*r + u*u*(-18*w*w + 20*v*r) + 12*(u*u+w*w)*(u*u+w*w)*log(v+r))/96;
F_vw  = (2*u*r*(2*u*u + 5*(v*v+w*w)) + 6*(-u*u*u*u + (v*v+w*w)*(v*v+w*w))*atanh(u/r)
         + 3*u*u*u*u*(log(u+r) - log(-u+r)))/48;
F_uvw = r*r*r*r*r/15;
```

Port notes: the single-arg `atan`/`atanh`/`log` are plain C99. Removable near-face singularities (a shared-face
pair drives an offset coordinate -> 0, e.g. `w->0` in `atan(uv/(w r))`) need a limit/epsilon guard — the
potential is finite there. Keep the outer Gauss nodes **symmetric** about A's centre (mandatory for the exact
transverse cancellation). Use a **low** outer order (~2-3/dim) for well-separated pairs, higher (~6-8) near
touching — this is the near/far adaptivity the hex path lacks today.

## Scope + the self block

Disjoint pairs only. The coincident/self block (A=B) is **out of scope** for the closed corner-sum: the box
potential antiderivative picks up an interior branch offset (measured `closed - numeric = 3*pi/2` at the box
centre — the classic Waldvogel interior branch), so the naive corner-sum is wrong there. It is left to the
existing `QuadBlockHex` self term, which is already exact (3.66e-15). A fully-closed self block would need the
octant-dependent branch corrections; it is not needed.

## Validation oracle

`packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/hdiv_charge_gram.wls` is the verified reference.
`wolframscript -file hdiv_charge_gram.wls` -> ALL PASS (full Q1 basis, disjoint blocks). The C++ port must
reproduce, on disjoint hex pairs: the box-box entries to ~1e-9 (e.g. unit-cube stack const-const = 0.9808850)
and the transverse exact-symmetry entries to machine-zero. Cross-check the C++ Gram against
`chargeGramSemiclosed[Sa, Sb, A, B, n]` before wiring it into the solve; then the existing
`validation_test/feec/test_hdiv_*` goldens gate the result.

## Implementation gate (codex + mdx)

The C++ port into `src/core/rad_hacapk_hdiv.cpp` (`QuadBlockHex` / `_build_charge_gram_hex`) is **codex's**
(source-edit ownership). Before it lands: measure the real per-entry factor on **mdx** (LAB timing is
contended). The operation-count estimate above (~2.5x near / ~20x far fewer outer points + the free symmetry
fix) is the go signal; the mdx wall-clock is the confirmation. This is a constant-factor build speedup — it does
not change the HACApK O(N log N) / ~N^1.23 structure.

## Self-verify lesson

The first compact hand-written `asinh`/`atan2` "Waldvogel" const `F0` tried here was **wrong**: its corner-sum
gave 1.070 vs the true 0.804 (its `d^3` was 0.42 off; the single-arg Nagy form failed identically). It was
caught only by checking the **quantity actually used** (the corner-sum / box potential), not just `d^3`. Verify
a closed form by the value it produces in situ, not only by a derivative identity — and prefer the established
prism-potential form for the const rather than a re-derived one-liner.
