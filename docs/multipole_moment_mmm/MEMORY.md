# Multipole-Moment MMM Memory

Date policy: 2026_06_26

This file records short "do not repeat this mistake" notes for the
multipole-moment MMM implementation.  Temporary JSON files are scratch logs;
the durable decisions and key numbers must be written here.

## Terminology: Retired Yano-MSC vs Live Collocation MMMM

Decision, 2026_06_28:

- "Yano-MSC" names the retired historical MSC implementation path.  Do not use
  that name for the live solver.
- The live surface-charge soft-iron solver is the canonical collocation MMMM
  path: multipole-moment MMM with MSC face charges, HACApK matrix-vector
  acceleration, and two-sided loop/co-loop deflation for loop-heavy cases.
  *(Superseded 2026_07_02: deflation / loop-free were REMOVED -- see "Coarse
  Tier: HACApK-BiCGSTAB, Loop-Free Abandoned" below.)*
- `demag_backend="yano"` must not remain as a compatibility alias.  Keeping it
  conflates the retired Yano-MSC path with live collocation MMMM and is a bug
  source.
- The explicit backend name for the live surface-charge path is
  `collocation_mmmm`.
- Galerkin MMMM is retired as a production branch because it duplicates the
  symmetric `B^T G B` direction already covered by HDiv-VIM while lacking
  nonlinear and image-symmetry coverage.

## Coarse Tier: HACApK-BiCGSTAB, Loop-Free Abandoned

Decision + measurements, 2026_07_02:

- Collocation MMMM **gives up loop-free** (Sugahara): the internal M is
  field-correct but loop-polluted, acceptable for the COARSE / optimization
  tier.  Accurate + hysteresis work uses the loop-free HDiv-VIM (the PRIMARY
  soft-iron method).  The two-sided loop/co-loop deflation and the moment
  H-LU were REMOVED and must not come back (no-pivot H-LU is wrong AND slow
  at compression < 1; `memory/collocation_loopfree_abandoned.md`).
- `rad.Solve(..., method=2)` on PURE-HEX moment routes to **HACApK-BiCGSTAB,
  matvec-only**: the chi-free geometry K is a `RadHACApKMomentSystem`
  H-matrix, CACHED CROSS-SOLVE on `radTApplication` together with the chi-free
  localL/diagK blocks (validity = interaction ptr + hacapk eps/leaf/eta +
  bit-exact centroid compare).  tet/wedge/mixed method-2 -> dense moment LU.
- mdx 3-way benchmark (2026_07_02, retired
  `results_moment_solvers_mdx_20260702.json`, inventoried in
  `docs/hdiv_vim/vim_examples_retirement_results.json`, 27/27 cases,
  mu_r = 200, eps = 1e-4): COLD solve -- cube 24.6k DOF: dense-K BiCGSTAB
  33.7 s vs H-matrix 4.4 s (7.7x); C-yoke 14.9k DOF: 13.9 s vs 2.5 s (5.4x);
  H-matrix alone reaches 48k DOF in 9.3 s / 2.5 GB.  WARM solve (the
  optimization-inner-loop per-iteration cost, cross-solve cache hit): cube
  24.6k DOF 4.0 s -> 0.062 s (65x); near-FLAT 0.03-0.11 s across all sizes.
  Memory: C-yoke 14.9k DOF 3.8 GB (dense K) vs 0.6 GB (H).  Correctness:
  method-2 external B matches the reference method to <= 5e-5 (the ACA eps);
  iterations 51-55 on the C-yoke, N-independent (block Jacobi holds at
  engineering mu_r).
- Verification recipe for "is the H-matrix engaged?":
  `rad.GetSolveStats()['linear_iterations'] > 0` (LU = 0) and a `hacapk_eps`
  sweep on an ELONGATED geometry must change the field (a compact cube is
  near-field dominated -> mostly dense blocks -> eps-insensitive).
- Kernel default flip (2026_07_02, same day): the ANALYTIC closed-form moment
  kernel (`moment_analytic_kernel`, d2efb88d) is now the DEFAULT -- exact
  (removes the Gauss error) AND 1.5x faster on the dominant H-matrix build
  (mdx knob matrix: ctype 28k DOF 2.97 -> 1.95 s, cube 24.6k 2.76 -> 1.84 s).
  Gauss stays selectable via `moment_analytic_kernel=False` for cross-checks.
  Knob-matrix negatives: `hacapk_eps=1e-3` shifts the field ~3% (keep 1e-4);
  `hacapk_leaf` 16/64 is a wash vs 32.

Speed review round 3, 2026_07_03 (verdict: sufficiently optimized -- CLOSED,
no code change):

- WARM re-solve runs **0 BiCGSTAB iterations** (`iterations_solve2 = 0` in all
  12 mdx method-2 cases): the flat magnetization array persists on the CACHED
  interaction across `rad.Solve` calls (ResetM does not clear it), so the
  Krylov start vector is the previous solution and the initial-residual check
  passes immediately.  Warm cost = one verification H-matvec + O(N) overhead
  (LAB 16.5k DOF: 0.030 s; ObjBckg python-callback RHS measured 0.1-0.4 ms).
- The mdx bench used tol 1e-10; at the PRODUCTION default `bicgstab_tol=1e-4`
  cold iterations drop ~50 -> ~14, so a cold solve is ~72% H-build / ~27%
  Krylov (LAB 16.5k DOF: 2.8 s = 2.03 + 0.76).  Cold is BUILD-bounded; the
  per-entry cost is already the analytic closed form and the sampled-entry
  count is HACApK's (as-is policy) -- no material lever left.
- Rejected micro-levers: skip Gauss face-sample fill when analytic (compute
  ~ms; fixed-size embedded arrays free no memory; conditional cache content
  would re-create the `cross-solve-cache-config-flag-key-and-lifecycle`
  hazard); GMRES default (stays opt-in); ObjBckg callback batching
  (negligible); dense-K cross-solve cache for methods 0/1 (redundant with
  method 2).  High-mu_r coarse-space preconditioning stays deprioritized (see
  "Engineering Benchmark Range" below -- the 3-mode coarse correction already
  failed it empirically).
- Timing-noise lesson (LAB): a probe row showing 14 iters SLOWER than 54 iters
  was CPU contention (load jumped 5% -> 57% mid-run), not code -- re-run
  timing anomalies on a quiet machine before believing them.

## Anderson+Picard Default (Hysteresis Stabilization)

Decision + measurements, 2026_07_03 (Sugahara: "Anderson+Picardで行きます"):

- `moment_anderson_depth` DEFAULTS to 1; the safeguarded Anderson(1) mixing is
  gated on ANY moment Picard solve (method 0 LU / 1 dense-K / 2 H-matrix) via
  the new `ctx.last_solve_was_moment` flag (previously method-2-only, default
  off).  Opt out: `rad.SolverConfig(moment_anderson_depth=0)`.
- Why: plain B-input Picard DIVERGES on a strongly-coupled 4x4x4 hex hysteresis
  block at the first steep DESCENDING-branch step (`relax_param=0.3` does NOT
  rescue -- mixing direction, not step size, is what matters).  Anderson(1)
  completes the full loop: method 2 median 4.5 iters/step (max 22), method 0
  median 10 (max 86), final Mz m0-vs-m2 rel ~6e-5.  Weakly-coupled goldens
  (single hex, 5-tet cube) cannot see this failure mode.
- Safeguarded acceptance (keep accelerated iterate only if the residual drops)
  means linear / well-behaved solves are unchanged: 74 regression tests green
  on the default flip.
- Routing fact locked the same day: on all-moment hysteresis bodies BOTH
  `b_input_newton=True` and `b_input_hantila=True` route to this same moment
  B-input Picard(+Anderson) (verified bit-identical); the dense 3-DOF
  Newton/Hantila survive only for genuine 3-DOF dipoles (none exist for soft
  iron post-unification).
- Golden: `validation_test/hysteresis/test_binput_moment.py::
  test_E_coupled_block_loop_default_anderson` (locks the default value, the
  full-loop completion on m0+m2, hysteretic branch separation, and m0==m2).

## Engineering Benchmark Range

Decision, 2026_06_26:

- Do not let `mu_r >= 10000` linear high-permeability cases dominate
  engineering conclusions.  They are useful as numerical conditioning stress
  tests, but they are not the main production target for real soft-iron
  devices.
- For performance claims, prioritize nonlinear BH-curve cube/C-yoke runs and
  linear reference cases around ordinary engineering permeability
  (`mu_r ~ 100` to `5000`, with `mu_r ~ 1000` as the common lock/reference
  scale).
- If a preconditioner only helps the `mu_r >= 10000` stress test but slows or
  does not improve the nonlinear/`mu_r ~ 1000` cases, remove it or keep it as
  a separately justified new research branch.  Do not present it as the
  default path.

Two-stage smoke record, 2026_06_26:

| Stage | Case | Method 2 block-Jacobi inner iters | Removed 3-mode coarse correction inner iters | Result |
|---|---:|---:|---:|---|
| A, engineering | linear `mu_r=1000`, 3x3x2 block | 19 | 28 | worse |
| A, engineering | nonlinear BH, 3x4x2 block | 532 | 568 | worse |
| B, stress | linear `mu_r=10000`, 3x3x2 block | 31 | 36 | worse |
| B, stress | linear `mu_r=100000`, 3x3x2 block | 59 | 69 | worse |

Conclusion:

- The 3-mode global dipole coarse correction did not reduce BiCGSTAB
  iterations in either Stage A or Stage B.
- It was deleted from the implementation and public configuration surface.
- The raw scratch JSON was `C:\temp\radia_moment_two_stage_eval_2026_06_26.json`;
  this markdown table is the durable record.

## Method 1 Pure-Hex BiCGSTAB Must Stay Matrix-Free

Decision, 2026_06_26:

- Pure 6-DOF hex multipole-moment method 1 uses `MomentSystemBlock6x6` in split
  form `A(chi)x = Lx + chi*Kx` and stores only the RHS, current iterate, work
  vectors, and 6x6 element-block Jacobi inverses.
- Do not rebuild the dense `SystemMatrix` for pure-hex BiCGSTAB just because
  method 0 still needs dense LU.  That repeats the old comparison-path shortcut
  and hides the real method-1 memory behavior.
- Wedge/pyramid and mixed hex-wedge method 1 also stay matrix-free through
  `MomentSystemBlockAny`; do not bring back a dense variable-DOF method-1 path.
- There is no scalar or identity preconditioner substitute: if an element block
  inverse cannot be built, fail loud and fix the block/preconditioner issue.
- Smoke record: `C:\temp\radia_moment_method1_matrixfree_smoke_2026_06_26.json`
  compares sequential method 0 vs pure-hex method 1 at 24 and 108 DOF;
  external-B relative differences were about `6.3e-13` and `3.4e-12`.
- Follow-up smoke record: `C:\temp\radia_moment_accel_followup_2026_06_26.json`
  compares mixed hex+wedge method 1 against method 0; external-B relative
  difference was about `1.7e-11`.

## Three-Mode Coarse Correction Was Removed

Decision, 2026_06_26:

- The former `SolverConfig(moment_two_level_precond=True)` 3-mode global dipole
  correction was deleted from the implementation and public configuration
  surface.
- The small 108-DOF linear smoke in
  `C:\temp\radia_moment_accel_followup_2026_06_26.json` preserved the solution
  (`relB ~ 1.0e-11`) but did not reduce inner iterations (`14 -> 14`) and
  increased solve time from coarse-space overhead.
- The two-stage smoke summarized above was worse in every small case:
  `19->28`, `532->568`, `31->36`, and `59->69` iterations.
- Passing `moment_two_level_precond` now fails loud.
- Do not reintroduce this 3-mode correction as a casual performance knob.
- A future global/hierarchical preconditioner must be a new design and must
  first pass Stage A engineering cases before any Stage B stress-test gain is
  considered useful.

## Inexact BiCGSTAB Was Removed

Observation, re-confirmed on LAB during the 2026_06_26 acceleration pass:

- Inexact BiCGSTAB can reduce accumulated inner iterations, but the benefit is
  limited compared with nonlinear outer-loop and preconditioner improvements.
- On the 144-DOF nonlinear 3x4x2 lock case, exact method-2 BiCGSTAB used 446
  accumulated inner iterations and 60 outer Picard iterations.
- The tightened inexact schedule reduced this to 344 inner iterations but still
  60 outer iterations, with external-B difference about `4.2e-12`.
- A too-loose early schedule (`1e-4` to `1e-3` class when `bicgstab_tol=1e-9`)
  caused the nonlinear method-2 lock test to miss convergence. Do not repeat
  that schedule.
- Therefore, the implementation was removed from `SolverConfig` and the
  method-2 solve path.  Do not reintroduce it as a default or a casual
  performance knob.

Practical rule:

- Treat inexact BiCGSTAB as a failed optimization branch unless a larger
  benchmark proves that outer convergence and final fields are unchanged.
- Prefer work on nonlinear acceleration, GMRES comparison, and redesigned
  global/hierarchical preconditioning when seeking larger gains.

Relevant measurement:

- `C:\temp\radia_moment_krylov_accel_smoke_2026_06_26.json`
