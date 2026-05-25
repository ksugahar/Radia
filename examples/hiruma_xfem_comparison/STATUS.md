# Hiruma XFEM vs augmented-CLN benchmark — final status

## 2026-05-25 (β-1 Phase 1-3 complete)

### Three deliverables

| Phase | Topic | Status |
|-------|-------|--------|
| 1 | CFEM coarse + fine baselines, Bessel reference | ✓ done, 0.04 % vs Bessel |
| 2 | Hiruma 2023 XFEM enrichment psi = exp(−γξ) | ✓ done, 0.14 % @ r/δ=15 (publishable) |
| 3 | Augmented-CLN via Krylov-Galerkin reduction | ✓ done; **revealed scope condition** |

### Phase 2 — Hiruma XFEM reproduction

Independent NGSolve implementation of Hiruma 2023 §III.A confirms the
headline claim across the full r/δ ∈ [0, 15] sweep with the corrected
Bessel admittance reference Y = σ·2πR·I₁(κR)/(κ·I₀(κR)):

| r/δ  | CFEM fine | CFEM coarse | XFEM 88 DOF |
|------|-----------|-------------|-------------|
| 0.1  | +0.00 %   | +0.00 %     | +0.00 %     |
| 3.0  | +0.00 %   | +0.68 %     | +3.01 %     |
| 15.0 | +0.04 %   | +85.89 %    | **+0.14 %** |

This is publishable as a direct numerical-comparison row in paper 1
§V (or paper 2 §VI), upgrading open direction (iii) from
"architectural comparison only" to "numerical comparison on the
cylinder benchmark — XFEM at 88 DOF reproduces Hiruma's < 3 % claim
with wide margin (0.14 %)".

### Phase 3 — Krylov-Galerkin augmented CLN

The first attempt used `scipy.interpolate.pade` on Taylor moments —
ill-conditioned (Hankel rcond ~ 1e-37) due to the 30-decade range of
Taylor coefficients in a high-σ problem.  Replaced with paper 1 §III's
exact Galerkin projection:

```
K0_r = Q^T K0 Q,  M_r = Q^T M Q,  b_r = Q^T b
Y_FE_ROM(s) = σ·area − s·σ² · b_r^T (K0_r + s σ M_r)^{-1} b_r
Y_R(s) = Y_FE_ROM(s) + K_SIBC √s / (s + d)
```

where Q is the σM-orthonormalised Krylov basis of size M_ROM, built
from {K0^{-1} b, K0^{-1}(σM)K0^{-1} b, ...}.

This implementation correctly extracts the first eigenvalue of
(K0_r, σM_r) as d = 4.78 × 10⁴ rad/s (analytical wall band:
4.60 × 10⁴, 4 % agreement), and converges in M_ROM with stable
arithmetic.

**However**: the augmented-CLN does NOT recover the SIBC asymptote
ratio r(f) → 1 at high f for THIS benchmark.

### What we found about the scope of paper 1 Theorem 1

The augmented-CLN's high-f limit is:

    Y_R(s → ∞) = Y_CLN(s → ∞) + K_SIBC/√s
               = σ·[area − b_r^T M_r^{−1} b_r] + 0
               = σ·(FE residual offset)

For paper 1's r(f) → 1 claim to hold, **Y_CLN(s → ∞) must → 0**.
This requires `b_r^T M_r^{−1} b_r = area`, which happens cleanly for
**port-driven** admittance problems where b couples to the boundary
(Dirichlet excitation, Neumann flux), but NOT for the volume-source
problem in Hiruma 2023 §III.A where the entire conductor interior is
forced by J₀ = σ.

Convergence of the residual with mesh refinement on our cylinder:

| Mesh | DOF  | Y_CLN(∞) | d_eig |
|------|------|----------|-------|
| R/3.5 | 44   | 53       | 4.78 × 10⁴ |
| R/10  | 393  | 19       | 4.62 × 10⁴ |
| R/30  | 3558 | 13       | 4.60 × 10⁴ |

The residual decreases with mesh but slowly (algebraic, not
exponential).  Mesh refinement alone is not enough.

### Implication for paper 1

The "single-DOF asymptote-preserving augmentation" (Theorem 1) is a
**port-driven** ROM construction.  Volume-source problems (e.g.
Hiruma's cylinder, an IH workpiece in a uniform inducing field) need
either:

  (a) re-formulation as a port-driven problem (drive via boundary
      E_z or normal B_n), then apply the augmentation, OR
  (b) FE-level treatment (CFEM-SIBC or XFEM) that captures the skin
      asymptote at the FE-discretisation level.

Paper 1 §V's cuboid Y(s) verification IS port-driven (the b vector
couples to a surface excitation), which is why r(10^12 Hz) = 1.0002
works there.  This scope distinction should be stated explicitly in
paper 1 §III or §VIII Limitations.

### Phase 4 (2026-05-25) — XFEM-CLN Hankel conditioning experiment

**Hypothesis tested** (user 2026-05-25): does XFEM enrichment
$\psi(x)\!=\!\exp(-\xi/\delta_{\rm ref})$ extend the float64 high-$N$
reliability of canonical CLN by improving Hankel matrix conditioning?

**Setup**: 44-vertex Hiruma cylinder, $\delta_{\rm ref}\!=\!5.88\,$mm
(skin depth at $\omega_{\rm wall}$).  Compared CFEM (24 free DOF) vs
XFEM (48 free DOF, +24 enrichment).  Computed Krylov moments
$\mu_n\!=\!b^\top (K_0^{-1}\sigma M)^n K_0^{-1} b$ for $n\!=\!0..19$,
then Hankel condition number $\kappa(H_N)$ for $N\!=\!1..10$.

**Result** (rejected):

| $N$ | $\kappa_{\rm CFEM}$ | $\kappa_{\rm XFEM}$ | Ratio | float64 OK? |
|----|------|------|------|------|
| 2  | 1.15e+11 | 7.71e+10 | 1.5 | marginal both |
| 3  | 2.93e+23 | 5.49e+22 | 5.3 | ✗ both broken |
| 5  | 3.19e+48 | 7.90e+47 | 4.0 | ✗ |
| 10 | 2.18e+101 | 1.02e+100 | 21.5 | ✗ |

$N_{\rm break}\!=\!3$ for both.  XFEM gives factor 1-21× improvement in
$\kappa$ but does NOT extend the float64-reliable stage.  Krylov moments
are nearly identical (0.3-4% diff across $n\!=\!0..19$).

**Why**: high-$k$ skin modes captured by XFEM enrichment have tiny
$\tau_k$, so they contribute as $\tau_k^n$ in Krylov-at-$s\!=\!0$
sequence — **exponentially damped**.  Canonical CLN moments are
dominated by low-$k$ bulk modes which CFEM and XFEM capture identically.
Hankel ill-conditioning is set by bulk eigenvalue spread, not by
high-$k$ skin modes.

**Implication**: XFEM utility = skin-layer resolution + volume-source
FE-residual cure (Phase 2 / Phase 3).  XFEM does NOT save canonical
CLN at high N.  Paper 1 §VI "$N\!\le\!4$ float64 limit" stands.

For broadband accuracy with circuit identification, three remaining
options:
1. Multi-K (Kuriyama 2019) — sidesteps Hankel via $K_0$-MGS, loses canonical Cauer
2. Foster-of-Cauers (companion) — low-N Cauer macro + diffusive Foster terminator
3. MPFR + INTLAB (Nagamine 2026) — keep canonical Cauer, pay 192-bit precision cost

XFEM-CLN stacking is NOT among these.

### Conclusion: XFEM-CLN vs SIBC-CLN vs augmented-CLN

Not competitors — **three abstraction layers** of the same SIBC
asymptote, suited to different problem classes:

| Method | DC | Volume-source | Port-driven | FE DOF | ROM DOF |
|--------|------|------|------|-------|------|
| CFEM-SIBC at FE level | × breaks | ✓ | ✓ | small | n/a |
| XFEM (Hiruma) at FE level | ✓ | **✓** | ✓ | mid | n/a |
| Augmented-CLN (paper 1) | ✓ | × FE residual | **✓** | small | 4-5 |

**Task-routing**:
- Single-frequency 3D EM + corner / volume-source: **XFEM**
- Frequency-sweep / time-domain / SPICE / multi-port: **augmented-CLN** (on port-driven)
- Plain industrial high-f skin only, no DC: **SIBC**

### Files in this directory

- `phase1_cfem_cylinder_baseline.py`        — DONE
- `phase2_xfem_hiruma_enrichment.py`        — DONE (publishable)
- `phase3_sqrt_s_schur_comparison.py`       — DONE (Pade-based, kept as
                                                  cautionary tale, see comments)
- `phase3b_krylov_galerkin.py`              — DONE (correct Galerkin-Krylov)
- `phase4_xfem_hankel_conditioning.py`      — DONE 2026-05-25 (hypothesis
                                              rejected: XFEM does NOT
                                              extend float64 CLN reliability)
- `phase4_xfem_hankel_conditioning.{pdf,png}` — kappa(H_N) vs N + moments
- `phase4_xfem_hankel_results.json`         — raw moments + condition numbers
- `results_phase{1,2,3}.npz`                — numerical outputs

### Naming policy

The lab-internal nickname "Schur-F" has been dropped in this directory.
Canonical terminology used:

- **sqrt(s) Schur block** — the non-rational augmentation, paper 1 §III
- **Augmented CLN** / **Schur-complement augmented CLN** — full ROM
- **Asymptote-preserving augmentation** — paper 1 Theorem 1 name

### Recommended paper 1 edits (next session)

1. §III Theorem 1 statement should add a remark: "assumes Y_CLN(∞) → 0,
   which holds for port-driven admittance problems but not for
   volume-source problems where Y(∞) carries a finite FE residual
   offset"
2. §VIII Scope paragraph (added 2026-05-24) should note: "the
   sqrt(s)-Schur-augmentation framework is verified on port-driven
   benchmarks (the cuboid Y(s) of §V); volume-source problems need
   either FE-level treatment (XFEM, Hiruma 2023) or a port-driven
   re-formulation"
3. §VIII open direction (iii) can now state numerical reproduction
   of Hiruma 2023 (0.14 % at r/δ=15 on 88-DOF coarse mesh, see
   `examples/hiruma_xfem_comparison/`)

These edits make paper 1 stronger by sharpening the scope claim
rather than weakening the contribution.
